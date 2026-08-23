/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifdef _WIN32
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif


// Qt includes

// Slicer includes
#include "qSlicerDynamicPETModuleWidget.h"
#include "ui_qSlicerDynamicPETModuleWidget.h"

#include <QApplication>
#include <QEventLoop>
#include <QStandardItemModel>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QDoubleSpinBox>
#include <QTimer>

#include <vtkWeakPointer.h>

#include <cmath>
#include <limits>
#include <algorithm>
#include <array>

#ifdef _WIN32
#include <sstream>
#include <iomanip>

// strptime replacement for Windows
inline char* strptime(
    const char* s,
    const char* f,
    struct tm* tm)
{
    std::istringstream input(s);
    input >> std::get_time(tm, f);

    if (input.fail())
    {
        return nullptr;
    }

    return (char*)(s + input.tellg());
}
#endif

namespace
{

enum class IFCurveDomain
{
  WholeBlood = 0,
  TotalPlasma = 1,
  ParentPlasma = 2
};

enum class PBIFTemplateDomain
{
  MatchPatient = 0,
  WholeBlood = 1,
  TotalPlasma = 2
};

enum class ActivityUnit
{
  BqPerMl = 0,
  KBqPerMl = 1,
  SUVbw = 2
};

struct InputFunctionResult
{
  // Plasma curve before optional parent-fraction correction.
  // If plasmaIsParent is true, this curve is already parent plasma.
  std::vector<double> framePlasma;
  std::vector<double> nativePlasmaTimesSec;
  std::vector<double> nativePlasmaValues;

  // Total whole blood used by the vascular TCM term.
  std::vector<double> frameWholeBlood;
  std::vector<double> nativeWholeBloodTimesSec;
  std::vector<double> nativeWholeBloodValues;

  // Parent fraction remains independent until method-specific sampling.
  std::vector<double> parentFractionTimesSec;
  std::vector<double> parentFractionValues;

  // Ready-to-use PET-frame plasma input for MTGA.
  std::vector<double> frameModelPlasma;

  // Whether the original patient IF observation at each PET frame is retained.
  // A false value does not create a hole in the continuous IF: the missing IF
  // observation is reconstructed across neighboring retained samples for
  // integration, but MTGA excludes that frame from regression.
  std::vector<bool> frameKeep;

  // Number of complete PET frames for which the final input function is
  // supported. Terminally removed IDIF observations shorten this support
  // instead of invalidating the whole input function.
  size_t supportFrameCount{0};

  bool plasmaIsParent{false};
  bool applyParentFraction{false};
  bool hasWholeBlood{false};

  bool pbifApplied{false};
  double pbifScale{1.0};
  IFCurveDomain pbifCalibrationDomain{IFCurveDomain::WholeBlood};
  std::vector<double> pbifPatientCalibrationTimesSec;
  std::vector<double> pbifPatientCalibrationValues;
  std::vector<double> pbifScaledValues;

  bool sourceProcessingApplied{false};
  QString sourceProcessingLabel;
  std::vector<double> processedSourcePreviewTimesSec;
  std::vector<double> processedSourcePreviewValues;
  FengParameters fengParameters;
};

struct PreviewCurve
{
  QString name;
  std::vector<double> times;
  std::vector<double> values;
  bool pointsOnly{false};
};

static std::vector<double> lowessPredict(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& evalX,
    double span,
    int robustIterations = 2)
{
    const size_t n = x.size();
    if (n == 0 || n != y.size())
    {
        return {};
    }
    if (n < 3)
    {
        std::vector<double> out;
        out.reserve(evalX.size());
        for (double t : evalX)
        {
            const auto it = std::lower_bound(x.begin(), x.end(), t);
            if (it == x.begin()) out.push_back(y.front());
            else if (it == x.end()) out.push_back(y.back());
            else
            {
                const size_t r = static_cast<size_t>(std::distance(x.begin(), it));
                const size_t l = r - 1;
                const double f = (t - x[l]) / (x[r] - x[l]);
                out.push_back(y[l] + f * (y[r] - y[l]));
            }
        }
        return out;
    }

    span = std::max(0.05, std::min(1.0, span));
    const size_t neighborhood = std::min(
        n,
        std::max<size_t>(3, static_cast<size_t>(std::ceil(span * static_cast<double>(n)))));

    auto localPrediction = [&](double x0, const std::vector<double>& robustWeights)
    {
        std::vector<double> distances(n);
        for (size_t i = 0; i < n; ++i)
        {
            distances[i] = std::abs(x[i] - x0);
        }
        std::vector<double> sortedDistances = distances;
        std::nth_element(
            sortedDistances.begin(),
            sortedDistances.begin() + static_cast<std::ptrdiff_t>(neighborhood - 1),
            sortedDistances.end());
        double bandwidth = sortedDistances[neighborhood - 1];
        if (bandwidth <= 1e-12)
        {
            bandwidth = *std::max_element(distances.begin(), distances.end());
        }
        if (bandwidth <= 1e-12)
        {
            return y.front();
        }

        double sw = 0.0;
        double swx = 0.0;
        double swy = 0.0;
        double swxx = 0.0;
        double swxy = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            const double u = distances[i] / bandwidth;
            if (u >= 1.0)
            {
                continue;
            }
            const double oneMinusCube = 1.0 - u * u * u;
            double w = oneMinusCube * oneMinusCube * oneMinusCube;
            if (i < robustWeights.size())
            {
                w *= robustWeights[i];
            }
            const double dx = x[i] - x0;
            sw += w;
            swx += w * dx;
            swy += w * y[i];
            swxx += w * dx * dx;
            swxy += w * dx * y[i];
        }

        if (sw <= 1e-16)
        {
            return 0.0;
        }

        const double denom = sw * swxx - swx * swx;
        if (std::abs(denom) <= 1e-16)
        {
            return swy / sw;
        }

        // Centering x around x0 makes the fitted intercept the prediction at x0.
        return (swy * swxx - swx * swxy) / denom;
    };

    std::vector<double> robustWeights(n, 1.0);
    std::vector<double> fittedAtData(n, 0.0);

    for (int iteration = 0; iteration <= robustIterations; ++iteration)
    {
        for (size_t i = 0; i < n; ++i)
        {
            fittedAtData[i] = localPrediction(x[i], robustWeights);
        }

        if (iteration == robustIterations)
        {
            break;
        }

        std::vector<double> absResiduals(n, 0.0);
        for (size_t i = 0; i < n; ++i)
        {
            absResiduals[i] = std::abs(y[i] - fittedAtData[i]);
        }
        std::vector<double> tmp = absResiduals;
        std::nth_element(tmp.begin(), tmp.begin() + static_cast<std::ptrdiff_t>(n / 2), tmp.end());
        const double medianAbsResidual = tmp[n / 2];
        const double scale = 6.0 * medianAbsResidual;
        if (scale <= 1e-16)
        {
            break;
        }

        for (size_t i = 0; i < n; ++i)
        {
            const double u = absResiduals[i] / scale;
            if (u >= 1.0)
            {
                robustWeights[i] = 0.0;
            }
            else
            {
                const double oneMinusSquare = 1.0 - u * u;
                robustWeights[i] = oneMinusSquare * oneMinusSquare;
            }
        }
    }

    std::vector<double> out;
    out.reserve(evalX.size());
    for (double x0 : evalX)
    {
        out.push_back(std::max(0.0, localPrediction(x0, robustWeights)));
    }
    return out;
}

bool segmentNameLessCaseInsensitive(
    const std::string& a,
    const std::string& b)
{
  const QString qa = QString::fromStdString(a);
  const QString qb = QString::fromStdString(b);

  const int result =
    QString::compare(qa, qb, Qt::CaseInsensitive);

  if (result != 0)
    return result < 0;

  return QString::compare(
           qa, qb, Qt::CaseSensitive) < 0;
}


double statisticDispersionSigma(
    const VoxelStatistics& stats,
    const std::string& statistic)
{
  if (statistic == "Mean")
  {
    return stats.stddev;
  }

  if (statistic == "Median")
  {
    // For Gaussian data: IQR = 1.34898 * sigma.
    return stats.iqr / 1.3489795003921634;
  }

  if (statistic == "Peak")
  {
    return stats.peakStddev;
  }

  return std::numeric_limits<double>::quiet_NaN();
}

double inverseVarianceWeightFromSigma(double sigma)
{
  if (!std::isfinite(sigma) || sigma <= 1e-12)
  {
    // Spatial dispersion is a weighting proxy, not a known measurement
    // variance.  Do not turn a zero/undefined proxy into an infinite weight.
    return 1.0;
  }

  return 1.0 / (sigma * sigma);
}

void normalizePositiveWeights(std::vector<double>& weights)
{
  double sum = 0.0;
  size_t count = 0;

  for (const double weight : weights)
  {
    if (std::isfinite(weight) && weight > 0.0)
    {
      sum += weight;
      ++count;
    }
  }

  if (count == 0 || sum <= 0.0)
  {
    return;
  }

  const double mean = sum / static_cast<double>(count);
  for (double& weight : weights)
  {
    if (std::isfinite(weight) && weight > 0.0)
    {
      weight /= mean;
    }
  }
}

QString formatTCMBoundStatus(unsigned int flags)
{
  QStringList hits;

  auto append =
      [&](unsigned int lowerFlag,
          unsigned int upperFlag,
          const QString& name)
      {
        if (flags & lowerFlag)
        {
          hits << name + ":lower";
        }
        if (flags & upperFlag)
        {
          hits << name + ":upper";
        }
      };

  append(TCM_BOUND_K1_LOWER, TCM_BOUND_K1_UPPER, "K1");
  append(TCM_BOUND_K2_LOWER, TCM_BOUND_K2_UPPER, "k2");
  append(TCM_BOUND_K3_LOWER, TCM_BOUND_K3_UPPER, "k3");
  append(TCM_BOUND_K4_LOWER, TCM_BOUND_K4_UPPER, "k4");
  append(TCM_BOUND_VB_LOWER, TCM_BOUND_VB_UPPER, "vb");
  append(TCM_BOUND_TD_LOWER, TCM_BOUND_TD_UPPER, "td");
  append(TCM_BOUND_KA_LOWER, TCM_BOUND_KA_UPPER, "ka");
  append(TCM_BOUND_FA_LOWER, TCM_BOUND_FA_UPPER, "fA");

  return hits.join(", ");
}

std::vector<std::string> sortedSegmentIDs(
    const std::map<std::string, std::string>& segmentNames)
{
  std::vector<std::string> ids;
  ids.reserve(segmentNames.size());

  for (const auto& [segmentID, displayName] : segmentNames)
  {
    ids.push_back(segmentID);
  }

  std::sort(
    ids.begin(),
    ids.end(),
    [&](const std::string& a, const std::string& b)
    {
      return segmentNameLessCaseInsensitive(
        segmentNames.at(a),
        segmentNames.at(b));
    });

  return ids;
}

void padFittedCurveToFullFrames(
    double*& fittedCurve,
    size_t fittedFrameCount,
    size_t fullFrameCount)
{
  if (!fittedCurve || fittedFrameCount >= fullFrameCount)
  {
    return;
  }

  double* padded = new double[fullFrameCount];
  const double nan = std::numeric_limits<double>::quiet_NaN();

  for (size_t i = 0; i < fullFrameCount; ++i)
  {
    padded[i] = i < fittedFrameCount ? fittedCurve[i] : nan;
  }

  delete[] fittedCurve;
  fittedCurve = padded;
}

}


//-----------------------------------------------------------------------------
class qSlicerDynamicPETModuleWidgetPrivate: public Ui_qSlicerDynamicPETModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerDynamicPETModuleWidget);

protected:
  qSlicerDynamicPETModuleWidget* const q_ptr;
  void setDoubleField(QLineEdit* le, double lo, double hi, int decimals);
  void setIntField(QLineEdit* le, int lo, int hi);
  void previewParametricVoxelSelection();
  void updateBodySupportUI();
  void updateLiverParameterUI();
  void previewInputFunction();
  void previewPBIF();
  void previewParentFraction();
  void previewCompanionWholeBlood();
  void resetInputFunctionData(
    bool clearExternalData);
  std::string selectedIFInterpolation() const;
  IFCurveDomain selectedIFCurveDomain() const;
  PBIFTemplateDomain selectedPBIFTemplateDomain() const;
  ActivityUnit selectedDisplayActivityUnit() const;
  ActivityUnit selectedExternalIFActivityUnit() const;
  ActivityUnit selectedCompanionWholeBloodActivityUnit() const;
  ActivityUnit petStoredActivityUnit() const;
  QString activityUnitLabel(ActivityUnit unit) const;
  bool getCommonSUVbwFactor(
      double& factor,
      QString* errorMessage = nullptr) const;
  bool convertActivityValue(
      double value,
      ActivityUnit from,
      ActivityUnit to,
      double& converted,
      QString* errorMessage = nullptr) const;
  bool convertActivityVector(
      const std::vector<double>& values,
      ActivityUnit from,
      ActivityUnit to,
      std::vector<double>& converted,
      QString* errorMessage = nullptr) const;
  void updateQuantitativeUnitUI();
  void removeInputFunctionPreview();
  void removePreviewGroup(const std::string& groupName);
  void showCurvePreview(
      const std::string& groupName,
      const QString& title,
      const QString& yAxisTitle,
      const std::vector<PreviewCurve>& curves);
  std::vector<unsigned char>
      MTGAOptimizedSelection;

  std::vector<double>
      MTGAOptimizedKiValues;

  std::vector<double>
      MTGAOptimizedDVValues;

  std::string MTGAOptimizedKiNodeID;
  std::string MTGAOptimizedDVNodeID;
  std::string MTGAOptimizedRGBNodeID;

  std::vector<std::string> TCMOptimizedNodeIDs;
  std::string TCMOptimizedModelSelectionNodeID;

public:
  std::vector<std::string> segmentDisplayOrder;
  qSlicerDynamicPETModuleWidgetPrivate(qSlicerDynamicPETModuleWidget& object);
  ~qSlicerDynamicPETModuleWidgetPrivate()=default;
  void init();
  void populatePatientComboBox();
  void populateStudyComboBox(vtkIdType patientID);
  void populateNodeComboBox(QComboBox* comboBox, vtkIdType parentItemID, const char * requiredNodeType, const std :: string requiredModality);
  void populateSegmentCheckboxes(vtkIdType SegItemID);
  void populatePlotSegmentCheckboxes();
  void populateIF();
  void updateInputFunctionStatus();
  void updateROIModelingAvailability();
  void updateParametricImagingAvailability();
  std::string selectedIFStatistic();
  bool buildCurrentSegmentInputFunction(
      std::vector<double>& values,
      QString* errorMessage = nullptr,
      std::vector<bool>* keepMask = nullptr);

  void updateInputFunctionUI();

  void invalidateInputFunctionResults();

  bool loadExternalInputFunctionCSV(
      const QString& filePath,
      QString* errorMessage = nullptr);
  bool loadCompanionWholeBloodCSV(
      const QString& filePath,
      QString* errorMessage = nullptr);
  bool loadPBIFCSV(
      const QString& filePath,
      QString* errorMessage = nullptr);
  bool loadParentFractionCSV(
      const QString& filePath,
      QString* errorMessage = nullptr);
  bool loadTwoColumnCurveCSV(
      const QString& filePath,
      const QString& secondColumnName,
      double minimumValue,
      double maximumValue,
      double insertedZeroTimeValue,
      std::vector<double>& timesSec,
      std::vector<double>& values,
      bool& zeroAnchorAdded,
      QString* errorMessage = nullptr) const;

  bool buildCurrentInputFunction(
      InputFunctionResult& result,
      bool requireWholeBlood,
      QString* errorMessage = nullptr);

  bool buildCurrentInputFunctionWeights(
      std::vector<double>& weights,
      bool weighted,
      QString* errorMessage = nullptr,
      bool excludeRemovedInputFrames = true);

  bool hasValidInputFunction(
      QString* errorMessage = nullptr,
      bool requireWholeBlood = false);

  bool ensureParametricPETData(
      QString* errorMessage = nullptr);

  double interpolateInputFunction(
      const std::vector<double>& times,
      const std::vector<double>& values,
      double targetTime,
      const std::string& interpolationType) const;

  double averageInputFunctionOverInterval(
      const std::vector<double>& times,
      const std::vector<double>& values,
      double startTime,
      double endTime,
      const std::string& interpolationType) const;
  double integrateInputFunctionOverInterval(
      const std::vector<double>& times,
      const std::vector<double>& values,
      double startTime,
      double endTime,
      const std::string& interpolationType) const;
  double integrateFrameAverageCurveOverInterval(
      const std::vector<double>& frameValues,
      double startTime,
      double endTime);
  double evaluateFrameCurve(
      const std::vector<double>& frameValues,
      double targetTime,
      const std::string& interpolationType);
  double averagePlasmaTimesParentFractionOverInterval(
      const std::vector<double>& plasmaTimes,
      const std::vector<double>& plasmaValues,
      bool plasmaIsFrameCurve,
      const std::string& plasmaInterpolation,
      const std::vector<double>& parentTimes,
      const std::vector<double>& parentValues,
      double startTime,
      double endTime);
  bool pbrAtTime(
      double timeSec,
      double& pbr,
      QString* errorMessage = nullptr) const;

  QString externalIFPath;
  std::vector<double> externalIFTimesSec;
  std::vector<double> externalIFConcentrations;
  bool externalIFZeroAnchorAdded{false};

  QString externalWholeBloodPath;
  std::vector<double> externalWholeBloodTimesSec;
  std::vector<double> externalWholeBloodConcentrations;
  bool externalWholeBloodZeroAnchorAdded{false};

  QString pbifPath;
  std::vector<double> pbifTimesSec;
  std::vector<double> pbifTemplateValues;
  bool pbifZeroAnchorAdded{false};

  QString parentFractionPath;
  std::vector<double> parentFractionTimesSec;
  std::vector<double> parentFractionValues;
  bool parentFractionZeroAnchorAdded{false};
  bool suvbwFactorValidated{false};

  bool inputFunctionCacheValid{false};
  InputFunctionResult cachedInputFunction;
  void populateVOI(std :: string ifID);
  void populateVOIMTGA(std :: string ifID);
  void populateResultsVOI();
  void populateResultsVOIMTGA();
  void populateResultsTable(std :: string segmentID);
  void populateResultsMTGATable(std :: string segmentID);
  void populateModelsTCM(std :: string segmentID);
  void populateModelsMTGA(std :: string segmentID);
  void populateTimeBarMTGA();
  void populateTimeBarMTGAImg();
  void populateModelCombo(QComboBox* comboToFill,
                          const std::string& otherSelectedModel,
                          const std::string& currentSelectedModel,
                          const std::string& segmentID);
  void populateModelComboTCM(QComboBox* comboToFill,
                             const std::string& otherSelectedModel,
                             const std::string& currentSelectedModel,
                             const std::string& segmentID);
  void setPostTACEnabled(bool enabled);
  void updateMTGAOutputUI();
  void updateTCMOutputUI();
  void updateMTGAOptimizationUI();
  enum class MTGAOptimizedClass : unsigned char
  {
    Excluded = 0,
    Patlak = 1,
    Reversible = 2
  };
  enum class BodySupportSource
  {
      Auto = 0,
      CT = 1,
      PET = 2,
      Union = 3,
      Intersection = 4
  };
  void generateMTGAOptimizedResult();

  void refreshMTGAOptimizedRGB();

  vtkMRMLScalarVolumeNode*
  createMTGAOptimizedScalarVolume(
      const std::vector<double>& values,
      const QString& name,
      vtkMRMLScalarVolumeNode* refPETNode,
      vtkMRMLSubjectHierarchyNode* shNode,
      vtkIdType refPetID);

  void removeMTGAOptimizedSceneNodes();

  void populateTCMOptimizationModels();
  void updateTCMOptimizationUI();

  void generateTCMOptimizedResult();
  void removeTCMOptimizedSceneNodes();

  void outputMTGAParametricResult(
      const std::string& modelID,
      vtkSlicerDynamicPETLogic* logic,
      vtkMRMLScalarVolumeNode* refPETNode,
      vtkMRMLSubjectHierarchyNode* shNode,
      vtkIdType refPetID);

  void outputTCMParametricResult(
      const std::string& modelID,
      vtkSlicerDynamicPETLogic* logic,
      vtkMRMLScalarVolumeNode* refPETNode,
      vtkMRMLSubjectHierarchyNode* shNode,
      vtkIdType refPetID);

  bool exportParametricMapDICOM(
      vtkMRMLScalarVolumeNode* refPETNode,
      const std::vector<double>& values,
      const std::string& method,
      const std::string& modelID,
      const std::string& field,
      const QString& outputDirectory,
      int seriesNumber,
      const QString& unitCode,
      const QString& unitMeaning,
      const QString& derivationDetails = QString());

  std::map<std::string, QString> MTGAImgFitSignatures;
  std::map<std::string, QString> TCMImgFitSignatures;

  bool parametricFitRunning{false};

  // Common parametric-imaging voxel selection.
  std::vector<unsigned char> parametricVoxelMask;
  std::vector<int> parametricFitVoxelIndices;

  vtkSmartPointer<vtkOrientedImageData>
      parametricBodySupportImage;

  vtkSmartPointer<vtkOrientedImageData>
      parametricFinalFitMaskImage;

  vtkIdType parametricVoxelSelectionCTID{
      vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID};

  vtkIdType parametricVoxelSelectionPETID{
      vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID};

  bool ensureParametricVoxelSelection();
  void invalidateParametricVoxelSelection();

  void resetParametricImagingSelections();
  void setPETItemID(vtkIdType newPetID);

};

//-----------------------------------------------------------------------------
// qSlicerDynamicPETModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerDynamicPETModuleWidgetPrivate::qSlicerDynamicPETModuleWidgetPrivate(qSlicerDynamicPETModuleWidget& object): q_ptr(&object)
{
  Q_Q(qSlicerDynamicPETModuleWidget);
}


void qSlicerDynamicPETModuleWidgetPrivate::setDoubleField(QLineEdit* le, double lo, double hi, int decimals)
{

  QObject::connect(le, &QLineEdit::editingFinished, le, [le, lo, hi, decimals]() {
    const QLocale loc = le->locale(); // respect UI locale (comma/dot)
    bool ok = false;
    double x = loc.toDouble(le->text(), &ok);
    if (!ok) {
      le->setText(loc.toString(lo, 'f', decimals));
      return;
    }
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    le->setText(loc.toString(x, 'f', decimals));
  });

  return ;
}

void qSlicerDynamicPETModuleWidgetPrivate::setIntField(QLineEdit* le, int lo, int hi)
{
  QObject::connect(le, &QLineEdit::editingFinished, le, [le, lo, hi]() {
    bool ok = false;
    int val = le->text().toInt(&ok);
    if (!ok) {
      le->setText(QString::number(lo));
      return;
    }
    if (val < lo) val = lo;
    if (val > hi) val = hi;
    le->setText(QString::number(val));
  });

  return ;
}

void qSlicerDynamicPETModuleWidgetPrivate::
previewParametricVoxelSelection()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->ensureParametricVoxelSelection())
    {
        QMessageBox::warning(
            q,
            QObject::tr("Patient support"),
            QObject::tr(
                "Could not create the parametric fitting support mask."));

        return;
    }

    vtkMRMLScene* scene =
        q->mrmlScene();

    if (!scene)
    {
        return;
    }

    vtkMRMLSubjectHierarchyNode* shNode =
        vtkMRMLSubjectHierarchyNode::
            GetSubjectHierarchyNode(scene);

    if (!shNode)
    {
        return;
    }

    vtkMRMLScalarVolumeNode* petNode =
        vtkMRMLScalarVolumeNode::SafeDownCast(
            shNode->GetItemDataNode(
                q->petID));

    if (!petNode)
    {
        return;
    }

    vtkSlicerDynamicPETLogic* logic =
        vtkSlicerDynamicPETLogic::SafeDownCast(
            q->logic());

    if (!logic)
    {
        return;
    }

    logic->CreateOrUpdateBodySupportPreview(
        this->parametricFinalFitMaskImage,
        petNode);
}

void qSlicerDynamicPETModuleWidgetPrivate::
updateBodySupportUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const bool hasCT =
        q->ctID !=
        vtkMRMLSubjectHierarchyNode::
            INVALID_ITEM_ID;

    const bool hasPET =
        q->petID !=
        vtkMRMLSubjectHierarchyNode::
            INVALID_ITEM_ID;

    const bool supportEnabled =
        this->BodySupportEnabledCheckBoxImg->
            isChecked();

    int sourceMode =
        this->BodySupportSourceImg->
            currentIndex();

    // If CT disappeared while a CT-dependent mode was selected,
    // safely return to Auto.
    if (!hasCT &&
        (sourceMode == 1 ||
         sourceMode == 3 ||
         sourceMode == 4))
    {
        this->BodySupportSourceImg->
            blockSignals(true);

        this->BodySupportSourceImg->
            setCurrentIndex(0);

        this->BodySupportSourceImg->
            blockSignals(false);

        sourceMode = 0;

        this->invalidateParametricVoxelSelection();
    }

    // Enable/disable individual source choices.
    QStandardItemModel* model =
        qobject_cast<QStandardItemModel*>(
            this->BodySupportSourceImg->
                model());

    if (model)
    {
        if (model->item(0))
            model->item(0)->setEnabled(hasPET); // Auto

        if (model->item(1))
            model->item(1)->setEnabled(hasPET && hasCT); // CT

        if (model->item(2))
            model->item(2)->setEnabled(hasPET); // PET

        if (model->item(3))
            model->item(3)->setEnabled(hasPET && hasCT); // Union

        if (model->item(4))
            model->item(4)->setEnabled(hasPET && hasCT); // Intersection
    }

    const bool sourceUsesCT =
        supportEnabled &&
        hasCT &&
        (
            sourceMode == 0 || // Auto -> CT if available
            sourceMode == 1 || // CT
            sourceMode == 3 || // Union
            sourceMode == 4    // Intersection
        );

    const bool sourceUsesPET =
        supportEnabled &&
        hasPET &&
        (
            sourceMode == 2 || // PET
            sourceMode == 3 || // Union
            sourceMode == 4 || // Intersection
            (sourceMode == 0 && !hasCT) // Auto without CT
        );

    this->BodySupportSourceImg->
        setEnabled(
            supportEnabled &&
            hasPET);

    this->BodySupportCTThresholdLabelImg->
        setEnabled(sourceUsesCT);

    this->BodySupportCTThresholdImg->
        setEnabled(sourceUsesCT);

    this->BodySupportMarginLabelImg->
        setEnabled(sourceUsesCT);

    this->BodySupportMarginImg->
        setEnabled(sourceUsesCT);

    this->BodySupportPETCompositeLabelImg->
        setEnabled(sourceUsesPET);

    this->BodySupportPETCompositeImg->
        setEnabled(sourceUsesPET);

    // Shared CT/PET operation.
    this->BodySupportFillHolesCheckBoxImg->
        setEnabled(
            supportEnabled &&
            hasPET);

    this->BodySupportAdvancedCollapsibleButtonImg->
        setEnabled(
            supportEnabled &&
            hasPET);

    this->BodySupportPreviewButtonImg->
        setEnabled(hasPET);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
removePreviewGroup(
    const std::string& groupName)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    vtkMRMLScene* scene = q->mrmlScene();
    if (!scene)
    {
        return;
    }

    std::vector<std::string> nodeIDs;

    for (int i = 0;
         i < scene->GetNumberOfNodes();
         ++i)
    {
        vtkMRMLNode* node =
            scene->GetNthNode(i);

        if (!node || !node->GetID())
        {
            continue;
        }

        const char* group =
            node->GetAttribute(
                "SlicerDynamicPET.IFPreviewGroup");

        if (group && groupName == group)
        {
            nodeIDs.push_back(node->GetID());
        }
    }

    vtkCollection* viewNodes =
        scene->GetNodesByClass(
            "vtkMRMLPlotViewNode");

    if (viewNodes)
    {
        for (int i = 0;
             i < viewNodes->GetNumberOfItems();
             ++i)
        {
            vtkMRMLPlotViewNode* viewNode =
                vtkMRMLPlotViewNode::SafeDownCast(
                    viewNodes->GetItemAsObject(i));

            if (!viewNode ||
                !viewNode->GetPlotChartNodeID())
            {
                continue;
            }

            const std::string activeChartID =
                viewNode->GetPlotChartNodeID();

            if (std::find(
                    nodeIDs.begin(),
                    nodeIDs.end(),
                    activeChartID) !=
                nodeIDs.end())
            {
                viewNode->SetPlotChartNodeID(nullptr);
            }
        }

        viewNodes->Delete();
    }

    for (const std::string& id : nodeIDs)
    {
        vtkMRMLNode* node =
            scene->GetNodeByID(id);

        if (node)
        {
            scene->RemoveNode(node);
        }
    }
}

void
qSlicerDynamicPETModuleWidgetPrivate::
showCurvePreview(
    const std::string& groupName,
    const QString& title,
    const QString& yAxisTitle,
    const std::vector<PreviewCurve>& curves)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    vtkMRMLScene* scene = q->mrmlScene();
    if (!scene || curves.empty())
    {
        return;
    }

    this->removePreviewGroup(groupName);

    double maximumTimeSec = 0.0;
    for (const PreviewCurve& curve : curves)
    {
        for (double timeSec : curve.times)
        {
            if (std::isfinite(timeSec))
            {
                maximumTimeSec = std::max(maximumTimeSec, timeSec);
            }
        }
    }

    const bool useMinutes = maximumTimeSec >= 300.0;
    const std::string timeColumnName =
        useMinutes ? "Time (min)" : "Time (s)";
    const double timeScale = useMinutes ? (1.0 / 60.0) : 1.0;

    vtkNew<vtkMRMLPlotChartNode> chartNode;
    chartNode->SetName(
        ("DynamicPET." + groupName + ".Chart").c_str());
    chartNode->SetTitle(title.toStdString().c_str());
    chartNode->SetXAxisTitle(timeColumnName.c_str());
    chartNode->SetYAxisTitle(
        yAxisTitle.toStdString().c_str());
    chartNode->SetAttribute(
        "SlicerDynamicPET.IFPreviewGroup",
        groupName.c_str());

    scene->AddNode(chartNode);

    int curveIndex = 0;

    for (const PreviewCurve& curve : curves)
    {
        if (curve.times.empty() ||
            curve.times.size() != curve.values.size())
        {
            continue;
        }

        vtkNew<vtkMRMLTableNode> tableNode;
        tableNode->SetName(
            ("DynamicPET." +
             groupName +
             ".Table." +
             std::to_string(curveIndex)).c_str());
        tableNode->SetAttribute(
            "SlicerDynamicPET.IFPreviewGroup",
            groupName.c_str());
        scene->AddNode(tableNode);

        vtkNew<vtkDoubleArray> timeArray;
        timeArray->SetName(timeColumnName.c_str());

        vtkNew<vtkDoubleArray> valueArray;
        valueArray->SetName("Value");

        for (size_t i = 0;
             i < curve.times.size();
             ++i)
        {
            timeArray->InsertNextValue(
                curve.times[i] * timeScale);
            valueArray->InsertNextValue(
                curve.values[i]);
        }

        tableNode->AddColumn(timeArray);
        tableNode->AddColumn(valueArray);

        vtkNew<vtkMRMLPlotSeriesNode> seriesNode;
        seriesNode->SetName(
            curve.name.toStdString().c_str());
        seriesNode->SetAttribute(
            "SlicerDynamicPET.IFPreviewGroup",
            groupName.c_str());
        seriesNode->SetPlotType(
            vtkMRMLPlotSeriesNode::PlotTypeScatter);
        seriesNode->SetAndObserveTableNodeID(
            tableNode->GetID());
        seriesNode->SetXColumnName(timeColumnName.c_str());
        seriesNode->SetYColumnName("Value");

        if (curve.pointsOnly)
        {
            seriesNode->SetLineStyle(
                vtkMRMLPlotSeriesNode::LineStyleNone);
        }
        else
        {
            seriesNode->SetMarkerStyle(
                vtkMRMLPlotSeriesNode::MarkerStyleNone);
        }

        scene->AddNode(seriesNode);
        seriesNode->SetUniqueColor();

        chartNode->AddAndObservePlotSeriesNodeID(
            seriesNode->GetID());

        ++curveIndex;
    }

    vtkMRMLLayoutNode* layoutNode =
        vtkMRMLLayoutNode::SafeDownCast(
            scene->GetFirstNodeByClass(
                "vtkMRMLLayoutNode"));

    if (layoutNode)
    {
        layoutNode->SetViewArrangement(
            vtkMRMLLayoutNode::
                SlicerLayoutConventionalPlotView);
    }

    vtkMRMLPlotViewNode* plotViewNode =
        vtkMRMLPlotViewNode::SafeDownCast(
            scene->GetFirstNodeByClass(
                "vtkMRMLPlotViewNode"));

    if (plotViewNode)
    {
        plotViewNode->SetPlotChartNodeID(
            chartNode->GetID());
    }
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateLiverParameterUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const bool liverSelected =
        std::find(
            q->modelsID.begin(),
            q->modelsID.end(),
            std::string("Liver DBIF")) !=
        q->modelsID.end();

    const std::array<QWidget*, 8> liverWidgets =
    {
        this->liverKaLabel,
        this->liverKaInitial,
        this->liverKaLower,
        this->liverKaUpper,
        this->liverFaLabel,
        this->liverFaInitial,
        this->liverFaLower,
        this->liverFaUpper
    };

    for (QWidget* widget : liverWidgets)
    {
        if (widget)
        {
            widget->setVisible(
                liverSelected);
        }
    }
}


void
qSlicerDynamicPETModuleWidgetPrivate::
previewInputFunction()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    InputFunctionResult result;
    QString error;

    if (!this->buildCurrentInputFunction(
            result,
            false,
            &error))
    {
        QMessageBox::warning(
            q,
            QObject::tr("Input Function"),
            error);
        return;
    }

    std::vector<double> sourceTimes;
    std::vector<double> sourceValues;

    const bool segmentSource =
        this->IFSourceSelector->currentIndex() == 0;

    if (segmentSource)
    {
        std::vector<double> allSourceValues;
        std::vector<bool> sourceKeep;

        if (!this->buildCurrentSegmentInputFunction(
                allSourceValues,
                &error,
                &sourceKeep))
        {
            QMessageBox::warning(
                q,
                QObject::tr("Input Function"),
                error);
            return;
        }

        sourceTimes.reserve(allSourceValues.size());
        sourceValues.reserve(allSourceValues.size());

        for (size_t i = 0; i < allSourceValues.size(); ++i)
        {
            if (i < sourceKeep.size() && !sourceKeep[i])
            {
                continue;
            }

            sourceTimes.push_back(
                q->timePoints[i]);
            sourceValues.push_back(allSourceValues[i]);
        }
    }
    else
    {
        sourceTimes = this->externalIFTimesSec;
        sourceValues = this->externalIFConcentrations;

        if (this->externalIFZeroAnchorAdded &&
            !sourceValues.empty())
        {
            sourceValues[0] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    const ActivityUnit displayUnit =
        this->selectedDisplayActivityUnit();

    const ActivityUnit sourceUnit =
        segmentSource
        ? this->petStoredActivityUnit()
        : this->selectedExternalIFActivityUnit();

    std::vector<double> sourceDisplayValues;
    if (!this->convertActivityVector(
            sourceValues,
            sourceUnit,
            displayUnit,
            sourceDisplayValues,
            &error))
    {
        QMessageBox::warning(
            q,
            QObject::tr("Input Function"),
            error);
        return;
    }

    auto toDisplay =
        [&](const std::vector<double>& nativeValues,
            std::vector<double>& displayValues) -> bool
        {
            return this->convertActivityVector(
                nativeValues,
                this->petStoredActivityUnit(),
                displayUnit,
                displayValues,
                &error);
        };

    const size_t previewFrameCount =
        std::min(
            result.supportFrameCount > 0
                ? result.supportFrameCount
                : q->durations.size(),
            q->durations.size());

    std::vector<double> frameTimes;
    frameTimes.reserve(previewFrameCount);
    for (size_t i = 0;
         i < previewFrameCount;
         ++i)
    {
        frameTimes.push_back(
            q->timePoints[i]);
    }

    auto clippedFrameValues =
        [&](const std::vector<double>& values)
        {
            const size_t n =
                std::min(previewFrameCount, values.size());
            return std::vector<double>(
                values.begin(),
                values.begin() +
                    static_cast<std::ptrdiff_t>(n));
        };

    std::vector<PreviewCurve> curves;
    curves.push_back(
        {QObject::tr("Original source"),
         sourceTimes,
         sourceDisplayValues,
         true});

    if (result.sourceProcessingApplied &&
        !result.processedSourcePreviewTimesSec.empty() &&
        result.processedSourcePreviewTimesSec.size() ==
            result.processedSourcePreviewValues.size())
    {
        std::vector<double> processedDisplayValues;
        if (!toDisplay(
                result.processedSourcePreviewValues,
                processedDisplayValues))
        {
            QMessageBox::warning(q, QObject::tr("Input Function"), error);
            return;
        }
        curves.push_back(
            {result.sourceProcessingLabel,
             result.processedSourcePreviewTimesSec,
             processedDisplayValues,
             false});
    }

    if (result.hasWholeBlood)
    {
        std::vector<double> displayValues;
        const bool useNative =
            !result.nativeWholeBloodTimesSec.empty();
        const std::vector<double> frameValues =
            clippedFrameValues(result.frameWholeBlood);
        const std::vector<double>& nativeValues =
            useNative
            ? result.nativeWholeBloodValues
            : frameValues;

        if (!toDisplay(nativeValues, displayValues))
        {
            QMessageBox::warning(q, QObject::tr("Input Function"), error);
            return;
        }

        curves.push_back(
            {QObject::tr("Total whole blood"),
             useNative
                 ? result.nativeWholeBloodTimesSec
                 : frameTimes,
             displayValues,
             false});
    }

    if (!result.plasmaIsParent)
    {
        std::vector<double> displayValues;
        const bool useNative =
            !result.nativePlasmaTimesSec.empty();
        const std::vector<double> frameValues =
            clippedFrameValues(result.framePlasma);
        const std::vector<double>& nativeValues =
            useNative
            ? result.nativePlasmaValues
            : frameValues;

        if (!toDisplay(nativeValues, displayValues))
        {
            QMessageBox::warning(q, QObject::tr("Input Function"), error);
            return;
        }

        curves.push_back(
            {QObject::tr("Total plasma"),
             useNative
                 ? result.nativePlasmaTimesSec
                 : frameTimes,
             displayValues,
             false});
    }

    if (result.plasmaIsParent ||
        result.applyParentFraction)
    {
        std::vector<double> displayValues;
        const std::vector<double> modelFrameValues =
            clippedFrameValues(result.frameModelPlasma);
        if (!toDisplay(modelFrameValues, displayValues))
        {
            QMessageBox::warning(q, QObject::tr("Input Function"), error);
            return;
        }

        curves.push_back(
            {QObject::tr("Parent plasma / model input (frame average)"),
             frameTimes,
             displayValues,
             false});
    }

    double pbr0 = 0.0;
    QString pbrText;
    if (this->pbrAtTime(0.0, pbr0, nullptr))
    {
        pbrText = QObject::tr(" | PBR(0)=%1").arg(pbr0, 0, 'g', 6);
    }

    this->showCurvePreview(
        "InputFunctionPreview",
        QObject::tr("Input Function - %1%2")
            .arg(this->activityUnitLabel(displayUnit))
            .arg(pbrText),
        this->activityUnitLabel(displayUnit),
        curves);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
previewPBIF()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (this->pbifTimesSec.size() < 2 ||
        this->pbifTemplateValues.size() !=
            this->pbifTimesSec.size())
    {
        QMessageBox::warning(
            q,
            QObject::tr("PBIF calibration"),
            QObject::tr(
                "No valid PBIF template CSV has been loaded."));
        return;
    }

    std::vector<double> rawPBIF =
        this->pbifTemplateValues;

    if (this->pbifZeroAnchorAdded &&
        !rawPBIF.empty())
    {
        rawPBIF[0] =
            std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<PreviewCurve> curves;
    curves.push_back(
        {QObject::tr("Unscaled PBIF samples"),
         this->pbifTimesSec,
         rawPBIF,
         true});

    const double previewStepSec =
        std::max(0.1, this->timeStepEdit->text().toDouble());
    const double previewEndSec = this->pbifTimesSec.back();

    std::vector<double> densePBIFTimes;
    std::vector<double> densePBIFValues;
    for (double t = 0.0; t <= previewEndSec + 1e-9; t += previewStepSec)
    {
        densePBIFTimes.push_back(t);
        densePBIFValues.push_back(
            this->interpolateInputFunction(
                this->pbifTimesSec,
                this->pbifTemplateValues,
                t,
                this->selectedIFInterpolation()));
    }

    curves.push_back(
        {QObject::tr("PBIF interpolation used by model"),
         densePBIFTimes,
         densePBIFValues,
         false});

    QString title =
        QObject::tr("PBIF template");
    QString yAxisTitle =
        QObject::tr("Template value (input units)");

    if (this->PBIFOptionCheckBox->isChecked())
    {
        InputFunctionResult result;
        QString error;

        if (!this->buildCurrentInputFunction(
                result,
                false,
                &error))
        {
            QMessageBox::warning(
                q,
                QObject::tr("PBIF calibration"),
                error);
            return;
        }

        if (!result.pbifApplied)
        {
            QMessageBox::warning(
                q,
                QObject::tr("PBIF calibration"),
                QObject::tr(
                    "No valid PBIF calibration is available."));
            return;
        }

        // The unscaled PBIF may be normalized or expressed in arbitrary
        // template units. Do not overlay it with the patient curve on the
        // same quantitative y axis once calibration is active.
        curves.clear();

        const ActivityUnit displayUnit =
            this->selectedDisplayActivityUnit();

        std::vector<double> patientDisplayValues;
        if (!this->convertActivityVector(
                result.pbifPatientCalibrationValues,
                this->petStoredActivityUnit(),
                displayUnit,
                patientDisplayValues,
                &error))
        {
            QMessageBox::warning(q, QObject::tr("PBIF calibration"), error);
            return;
        }

        curves.push_back(
            {QObject::tr("Patient calibration curve"),
             result.pbifPatientCalibrationTimesSec,
             patientDisplayValues,
             false});

        std::vector<double> denseScaledPBIFNative;
        denseScaledPBIFNative.reserve(densePBIFValues.size());
        for (double value : densePBIFValues)
        {
            denseScaledPBIFNative.push_back(value * result.pbifScale);
        }

        std::vector<double> denseScaledPBIF;
        if (!this->convertActivityVector(
                denseScaledPBIFNative,
                this->petStoredActivityUnit(),
                displayUnit,
                denseScaledPBIF,
                &error))
        {
            QMessageBox::warning(q, QObject::tr("PBIF calibration"), error);
            return;
        }

        curves.push_back(
            {QObject::tr("Scaled PBIF used by model"),
             densePBIFTimes,
             denseScaledPBIF,
             false});

        const QString calibrationDomainName =
            result.pbifCalibrationDomain == IFCurveDomain::WholeBlood
            ? QObject::tr("Whole blood")
            : QObject::tr("Total plasma");

        title =
            QObject::tr("PBIF calibration - %1, scale %2, AUC %3 to %4 s")
                .arg(calibrationDomainName)
                .arg(result.pbifScale, 0, 'g', 8)
                .arg(this->PBIFCalibrationStartSpinBox->value())
                .arg(this->PBIFCalibrationEndSpinBox->value());

        yAxisTitle =
            this->activityUnitLabel(this->selectedDisplayActivityUnit());
    }

    this->showCurvePreview(
        "PBIFPreview",
        title,
        yAxisTitle,
        curves);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
previewParentFraction()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (this->parentFractionTimesSec.size() < 2 ||
        this->parentFractionValues.size() !=
            this->parentFractionTimesSec.size())
    {
        QMessageBox::warning(
            q,
            QObject::tr("Parent fraction"),
            QObject::tr(
                "No valid parent-fraction CSV has been loaded."));
        return;
    }

    std::vector<double> rawValues =
        this->parentFractionValues;

    if (this->parentFractionZeroAnchorAdded &&
        !rawValues.empty())
    {
        rawValues[0] =
            std::numeric_limits<double>::quiet_NaN();
    }

    const double endTime =
        this->parentFractionTimesSec.back();
    const int n = 300;

    std::vector<double> denseTimes;
    std::vector<double> denseValues;
    denseTimes.reserve(n);
    denseValues.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        const double t =
            endTime *
            static_cast<double>(i) /
            static_cast<double>(n - 1);

        denseTimes.push_back(t);
        denseValues.push_back(
            this->interpolateInputFunction(
                this->parentFractionTimesSec,
                this->parentFractionValues,
                t,
                "linear"));
    }

    this->showCurvePreview(
        "ParentFractionPreview",
        QObject::tr("Parent Fraction"),
        QObject::tr("Parent fraction"),
        {
          {QObject::tr("Measured parent fraction"),
           this->parentFractionTimesSec,
           rawValues,
           true},
          {QObject::tr("Linear interpolation"),
           denseTimes,
           denseValues,
           false}
        });
}

void
qSlicerDynamicPETModuleWidgetPrivate::
previewCompanionWholeBlood()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (this->externalWholeBloodTimesSec.size() < 2 ||
        this->externalWholeBloodConcentrations.size() !=
            this->externalWholeBloodTimesSec.size())
    {
        QMessageBox::warning(
            q,
            QObject::tr("Whole blood"),
            QObject::tr(
                "No valid companion whole-blood CSV has been loaded."));
        return;
    }

    std::vector<double> rawValues =
        this->externalWholeBloodConcentrations;

    if (this->externalWholeBloodZeroAnchorAdded &&
        !rawValues.empty())
    {
        rawValues[0] =
            std::numeric_limits<double>::quiet_NaN();
    }

    QString error;
    std::vector<double> rawDisplay;
    std::vector<double> curveDisplay;
    const ActivityUnit displayUnit = this->selectedDisplayActivityUnit();

    if (!this->convertActivityVector(
            rawValues,
            this->selectedCompanionWholeBloodActivityUnit(),
            displayUnit,
            rawDisplay,
            &error) ||
        !this->convertActivityVector(
            this->externalWholeBloodConcentrations,
            this->selectedCompanionWholeBloodActivityUnit(),
            displayUnit,
            curveDisplay,
            &error))
    {
        QMessageBox::warning(q, QObject::tr("Whole blood"), error);
        return;
    }

    this->showCurvePreview(
        "CompanionWholeBloodPreview",
        QObject::tr("Companion Whole Blood"),
        this->activityUnitLabel(displayUnit),
        {
          {QObject::tr("Whole-blood samples"),
           this->externalWholeBloodTimesSec,
           rawDisplay,
           true},
          {QObject::tr("Whole blood"),
           this->externalWholeBloodTimesSec,
           curveDisplay,
           false}
        });
}

void
qSlicerDynamicPETModuleWidgetPrivate::
resetInputFunctionData(
    bool clearExternalData)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    q->IFID.clear();

    if (clearExternalData)
    {
        this->externalIFPath.clear();
        this->externalIFTimesSec.clear();
        this->externalIFConcentrations.clear();
        this->externalIFZeroAnchorAdded = false;
        this->IFCSVPathEdit->clear();

        this->externalWholeBloodPath.clear();
        this->externalWholeBloodTimesSec.clear();
        this->externalWholeBloodConcentrations.clear();
        this->externalWholeBloodZeroAnchorAdded = false;
        this->IFWholeBloodPathEdit->clear();

        this->pbifPath.clear();
        this->pbifTimesSec.clear();
        this->pbifTemplateValues.clear();
        this->pbifZeroAnchorAdded = false;
        this->PBIFPathEdit->clear();
        this->PBIFOptionCheckBox->blockSignals(true);
        this->PBIFOptionCheckBox->setChecked(false);
        this->PBIFOptionCheckBox->blockSignals(false);

        this->parentFractionPath.clear();
        this->parentFractionTimesSec.clear();
        this->parentFractionValues.clear();
        this->parentFractionZeroAnchorAdded = false;
        this->ParentFractionPathEdit->clear();
        this->MetaboliteCorrectionCheckBox->blockSignals(true);
        this->MetaboliteCorrectionCheckBox->setChecked(false);
        this->MetaboliteCorrectionCheckBox->blockSignals(false);
    }

    this->removeInputFunctionPreview();
    this->invalidateInputFunctionResults();
    this->updateInputFunctionUI();
    this->updateParametricImagingAvailability();
}

void
qSlicerDynamicPETModuleWidgetPrivate::
removeInputFunctionPreview()
{
    this->removePreviewGroup(
        "InputFunctionPreview");
    this->removePreviewGroup(
        "PBIFPreview");
    this->removePreviewGroup(
        "ParentFractionPreview");
    this->removePreviewGroup(
        "CompanionWholeBloodPreview");
}

std::string
qSlicerDynamicPETModuleWidgetPrivate::
selectedIFInterpolation() const
{
    return
        this->IFInterpolationSelector->
            currentIndex() == 1
        ? "const"
        : "linear";
}

IFCurveDomain
qSlicerDynamicPETModuleWidgetPrivate::
selectedIFCurveDomain() const
{
    // Image-derived vascular TACs are whole blood by definition
    // for this pipeline. The external selector is ignored here.
    if (this->IFSourceSelector->currentIndex() == 0)
    {
        return IFCurveDomain::WholeBlood;
    }

    const int index =
        this->IFCurveTypeSelector->currentIndex();

    if (index == 1)
    {
        return IFCurveDomain::TotalPlasma;
    }

    if (index == 2)
    {
        return IFCurveDomain::ParentPlasma;
    }

    return IFCurveDomain::WholeBlood;
}

PBIFTemplateDomain
qSlicerDynamicPETModuleWidgetPrivate::
selectedPBIFTemplateDomain() const
{
    const int index =
        this->PBIFDomainSelector->currentIndex();

    if (index == 1)
    {
        return PBIFTemplateDomain::WholeBlood;
    }

    if (index == 2)
    {
        return PBIFTemplateDomain::TotalPlasma;
    }

    return PBIFTemplateDomain::MatchPatient;
}

ActivityUnit
qSlicerDynamicPETModuleWidgetPrivate::
selectedDisplayActivityUnit() const
{
    return this->QuantitativeDisplayUnitSelector->currentIndex() == 1
        ? ActivityUnit::BqPerMl
        : ActivityUnit::SUVbw;
}

ActivityUnit
qSlicerDynamicPETModuleWidgetPrivate::
selectedExternalIFActivityUnit() const
{
    const int index = this->IFCSVUnitSelector->currentIndex();
    if (index == 1)
    {
        return ActivityUnit::KBqPerMl;
    }
    if (index == 2)
    {
        return ActivityUnit::SUVbw;
    }
    return ActivityUnit::BqPerMl;
}

ActivityUnit
qSlicerDynamicPETModuleWidgetPrivate::
selectedCompanionWholeBloodActivityUnit() const
{
    const int index = this->IFWholeBloodUnitSelector->currentIndex();
    if (index == 1)
    {
        return ActivityUnit::KBqPerMl;
    }
    if (index == 2)
    {
        return ActivityUnit::SUVbw;
    }
    return ActivityUnit::BqPerMl;
}

ActivityUnit
qSlicerDynamicPETModuleWidgetPrivate::
petStoredActivityUnit() const
{
    Q_Q(const qSlicerDynamicPETModuleWidget);
    return q->dPETvalueType == "SUVbw"
        ? ActivityUnit::SUVbw
        : ActivityUnit::BqPerMl;
}

QString
qSlicerDynamicPETModuleWidgetPrivate::
activityUnitLabel(ActivityUnit unit) const
{
    if (unit == ActivityUnit::SUVbw)
    {
        return QObject::tr("SUVbw");
    }
    if (unit == ActivityUnit::KBqPerMl)
    {
        return QObject::tr("kBq/mL");
    }
    return QObject::tr("Bq/mL");
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
getCommonSUVbwFactor(
    double& factor,
    QString* errorMessage) const
{
    Q_Q(const qSlicerDynamicPETModuleWidget);

    factor = 0.0;

    if (!this->suvbwFactorValidated || q->suvFactors.empty())
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "A validated SUVbw conversion factor is not available for the selected PET sequence.");
        }
        return false;
    }

    const double reference = q->suvFactors.front();
    if (!std::isfinite(reference) || reference <= 0.0)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "The stored SUVbw conversion factor is invalid.");
        }
        return false;
    }

    for (double current : q->suvFactors)
    {
        if (!std::isfinite(current) || current <= 0.0)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The SUVbw conversion factor is missing or invalid in one or more PET frames.");
            }
            return false;
        }

        const double tolerance =
            1e-9 * std::max(1.0, std::abs(reference));
        if (std::abs(current - reference) > tolerance)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "Frame-dependent SUVbw factors are not supported by the current input-function unit conversion.");
            }
            return false;
        }
    }

    factor = reference;
    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
convertActivityValue(
    double value,
    ActivityUnit from,
    ActivityUnit to,
    double& converted,
    QString* errorMessage) const
{
    if (!std::isfinite(value))
    {
        converted = value;
        return true;
    }

    if (from == to)
    {
        converted = value;
        return true;
    }

    double bqPerMl = value;

    if (from == ActivityUnit::KBqPerMl)
    {
        bqPerMl = value * 1000.0;
    }
    else if (from == ActivityUnit::SUVbw)
    {
        double factor = 0.0;
        if (!this->getCommonSUVbwFactor(factor, errorMessage))
        {
            return false;
        }
        bqPerMl = value / factor;
    }

    if (to == ActivityUnit::BqPerMl)
    {
        converted = bqPerMl;
        return true;
    }

    if (to == ActivityUnit::KBqPerMl)
    {
        converted = bqPerMl / 1000.0;
        return true;
    }

    double factor = 0.0;
    if (!this->getCommonSUVbwFactor(factor, errorMessage))
    {
        return false;
    }

    converted = bqPerMl * factor;
    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
convertActivityVector(
    const std::vector<double>& values,
    ActivityUnit from,
    ActivityUnit to,
    std::vector<double>& converted,
    QString* errorMessage) const
{
    converted.clear();
    converted.reserve(values.size());

    for (double value : values)
    {
        double output = 0.0;
        if (!this->convertActivityValue(
                value,
                from,
                to,
                output,
                errorMessage))
        {
            converted.clear();
            return false;
        }
        converted.push_back(output);
    }

    return true;
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateQuantitativeUnitUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const bool hasPET =
        q->petID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID &&
        q->sequencePETNode != nullptr;

    this->QuantitativeDisplayUnitLabel->setEnabled(hasPET);
    this->QuantitativeDisplayUnitSelector->setEnabled(hasPET);

    if (!hasPET)
    {
        return;
    }

    const ActivityUnit nativeUnit = this->petStoredActivityUnit();

    this->QuantitativeDisplayUnitSelector->blockSignals(true);
    this->QuantitativeDisplayUnitSelector->setCurrentIndex(
        nativeUnit == ActivityUnit::SUVbw ? 0 : 1);
    this->QuantitativeDisplayUnitSelector->blockSignals(false);

    QStandardItemModel* model = qobject_cast<QStandardItemModel*>(
        this->QuantitativeDisplayUnitSelector->model());

    if (model)
    {
        QStandardItem* suvItem = model->item(0);
        QStandardItem* bqItem = model->item(1);

        if (suvItem)
        {
            suvItem->setEnabled(
                nativeUnit == ActivityUnit::SUVbw || this->suvbwFactorValidated);
        }
        if (bqItem)
        {
            bqItem->setEnabled(
                nativeUnit == ActivityUnit::BqPerMl || this->suvbwFactorValidated);
        }
    }
}

//-----------------------------------------------------------------------------
void qSlicerDynamicPETModuleWidgetPrivate::init()
{
  Q_Q(qSlicerDynamicPETModuleWidget);
  this->setupUi(q);

  this->suvbwFactorValidated = false;
  this->updateQuantitativeUnitUI();

  this->IFSourceSelector->setCurrentIndex(0);

  QStandardItemModel* ifSourceModel =
      qobject_cast<QStandardItemModel*>(
          this->IFSourceSelector->model());

  if (ifSourceModel)
  {
    // Segment-derived IDIF supported now.
    if (ifSourceModel->item(0))
    {
      ifSourceModel->item(0)->setEnabled(true);
    }

    // Next implementation step.
    if (ifSourceModel->item(1))
    {
      ifSourceModel->item(1)->setEnabled(true);
    }

    QObject::connect(
        this->IFSourceSelector,
        QOverload<int>::of(
            &QComboBox::currentIndexChanged),
        q,
        [this, q](int source)
        {
            if (source == 0)
            {
                const int index =
                    this->IFSelector->
                        currentIndex();

                q->IFID =
                    index >= 0
                    ? this->IFSelector->
                        itemData(index)
                        .toString()
                        .toStdString()
                    : std::string();
            }
            else
            {
                q->IFID.clear();
            }

            const std::string excludedVOI =
                source == 0
                ? q->IFID
                : std::string();

            this->populateVOI(
                excludedVOI);

            this->populateVOIMTGA(
                excludedVOI);

            this->invalidateInputFunctionResults();
            this->updateInputFunctionUI();
        });

    QObject::connect(
        this->IFCSVBrowseButton,
        &QPushButton::clicked,
        q,
        [this, q]()
        {
            const QString path =
                QFileDialog::
                    getOpenFileName(
                        q,
                        QObject::tr(
                            "Select input-function CSV"),
                        QString(),
                        QObject::tr(
                            "CSV files (*.csv)"));

            if (path.isEmpty())
            {
                return;
            }

            QString error;

            if (!this->
                loadExternalInputFunctionCSV(
                    path,
                    &error))
            {
                QMessageBox::warning(
                    q,
                    QObject::tr(
                        "Input Function"),
                    error);

                return;
            }

            if (this->externalIFZeroAnchorAdded)
            {
                QMessageBox::warning(
                    q,
                    QObject::tr("Input Function"),
                    QObject::tr(
                        "The first input-function sample occurs after 0 seconds.\n\n"
                        "SlicerDynamicPET added an assumed concentration of 0 at "
                        "t = 0 s and will interpolate between this anchor and the "
                        "first measured sample.\n\n"
                        "This assumption can influence the early input-function area."));
            }

            this->IFCSVPathEdit->
                setText(path);

            this->invalidateInputFunctionResults();
            this->updateInputFunctionUI();
        });


    QObject::connect(
        this->IFCurveTypeSelector,
        QOverload<int>::of(
            &QComboBox::currentIndexChanged),
        q,
        [this](int)
        {
            this->invalidateInputFunctionResults();
            this->updateInputFunctionUI();
        });

    auto onPBRChanged =
        [this]()
        {
            this->invalidateInputFunctionResults();
            this->updateInputFunctionUI();
        };

    QObject::connect(
        this->pbrp1Edit,
        &QLineEdit::editingFinished,
        q,
        onPBRChanged);

    QObject::connect(
        this->pbrp2Edit,
        &QLineEdit::editingFinished,
        q,
        onPBRChanged);

    QObject::connect(
        this->pbrp3Edit,
        &QLineEdit::editingFinished,
        q,
        onPBRChanged);

  }

  // Future parent-fraction models are visible but disabled.
  if (QStandardItemModel* model =
          qobject_cast<QStandardItemModel*>(
              this->ParentFractionProcessingSelector->model()))
  {
      for (int i = 1; i < model->rowCount(); ++i)
      {
          if (model->item(i))
          {
              model->item(i)->setEnabled(false);
          }
      }
  }

  this->IFCurveTypeSelector->setCurrentIndex(0);
  this->IFSourceProcessingSelector->setCurrentIndex(0);
  this->ParentFractionProcessingSelector->setCurrentIndex(0);
  this->PBIFDomainSelector->setCurrentIndex(0);

  QObject::connect(
      this->IFWholeBloodBrowseButton,
      &QPushButton::clicked,
      q,
      [this, q]()
      {
          const QString path =
              QFileDialog::getOpenFileName(
                  q,
                  QObject::tr(
                      "Select whole-blood CSV"),
                  QString(),
                  QObject::tr(
                      "CSV files (*.csv)"));

          if (path.isEmpty())
          {
              return;
          }

          QString error;
          if (!this->loadCompanionWholeBloodCSV(
                  path,
                  &error))
          {
              QMessageBox::warning(
                  q,
                  QObject::tr("Whole blood"),
                  error);
              return;
          }

          this->IFWholeBloodPathEdit->setText(path);

          if (this->externalWholeBloodZeroAnchorAdded)
          {
              QMessageBox::warning(
                  q,
                  QObject::tr("Whole blood"),
                  QObject::tr(
                      "The first whole-blood sample occurs after 0 seconds.\n\n"
                      "An assumed point (0 s, 0) was added internally. "
                      "The original CSV file was not modified."));
          }

          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
      });

  QObject::connect(
      this->IFWholeBloodPreviewButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
          this->previewCompanionWholeBlood();
      });

  QObject::connect(
      this->PBIFBrowseButton,
      &QPushButton::clicked,
      q,
      [this, q]()
      {
          const QString path =
              QFileDialog::getOpenFileName(
                  q,
                  QObject::tr(
                      "Select PBIF template CSV"),
                  QString(),
                  QObject::tr(
                      "CSV files (*.csv)"));

          if (path.isEmpty())
          {
              return;
          }

          QString error;
          if (!this->loadPBIFCSV(
                  path,
                  &error))
          {
              QMessageBox::warning(
                  q,
                  QObject::tr("PBIF calibration"),
                  error);
              return;
          }

          this->PBIFPathEdit->setText(path);

          if (this->pbifZeroAnchorAdded)
          {
              QMessageBox::warning(
                  q,
                  QObject::tr("PBIF calibration"),
                  QObject::tr(
                      "The PBIF template begins after 0 seconds.\n\n"
                      "An assumed point (0 s, 0) was added internally. "
                      "The original CSV file was not modified."));
          }

          this->PBIFCalibrationStartSpinBox->setValue(0.0);

          double endTime =
              this->pbifTimesSec.empty()
              ? 0.0
              : this->pbifTimesSec.back();

          if (!q->timePoints.empty())
          {
              endTime = std::min(
                  endTime,
                  q->timePoints.back());
          }

          this->PBIFCalibrationEndSpinBox->setValue(
              endTime);

          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
      });

  QObject::connect(
      this->PBIFPreviewButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
          this->previewPBIF();
      });

  QObject::connect(
      this->ParentFractionBrowseButton,
      &QPushButton::clicked,
      q,
      [this, q]()
      {
          const QString path =
              QFileDialog::getOpenFileName(
                  q,
                  QObject::tr(
                      "Select parent-fraction CSV"),
                  QString(),
                  QObject::tr(
                      "CSV files (*.csv)"));

          if (path.isEmpty())
          {
              return;
          }

          QString error;
          if (!this->loadParentFractionCSV(
                  path,
                  &error))
          {
              QMessageBox::warning(
                  q,
                  QObject::tr("Parent fraction"),
                  error);
              return;
          }

          this->ParentFractionPathEdit->setText(path);

          if (this->parentFractionZeroAnchorAdded)
          {
              QMessageBox::warning(
                  q,
                  QObject::tr("Parent fraction"),
                  QObject::tr(
                      "The first parent-fraction measurement occurs after 0 seconds.\n\n"
                      "An assumed point (0 s, 1.0) was added internally. "
                      "The original CSV file was not modified."));
          }

          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
      });

  QObject::connect(
      this->ParentFractionPreviewButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
          this->previewParentFraction();
      });

  auto invalidateIFPipeline =
      [this]()
      {
          // Clear the cached/precomputed IF before rebuilding status. This is
          // important for nonlinear Feng processing, which must be recomputed
          // exactly once after a pipeline setting changes rather than showing
          // a stale cached source during the immediate UI refresh.
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
          this->updateROIModelingAvailability();
          this->updateParametricImagingAvailability();
      };

  QObject::connect(
      this->IFCSVUnitSelector,
      QOverload<int>::of(&QComboBox::currentIndexChanged),
      q,
      [invalidateIFPipeline](int)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->IFWholeBloodUnitSelector,
      QOverload<int>::of(&QComboBox::currentIndexChanged),
      q,
      [invalidateIFPipeline](int)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->QuantitativeDisplayUnitSelector,
      QOverload<int>::of(&QComboBox::currentIndexChanged),
      q,
      [this, q](int)
      {
          // Display conversion never changes fitting data/results.
          this->removeInputFunctionPreview();
          this->updateInputFunctionStatus();
          if (!q->segmentTACs.empty())
          {
              q->onPlotbutton();
          }
      });

  QObject::connect(
      this->PBIFOptionCheckBox,
      &QCheckBox::toggled,
      q,
      [invalidateIFPipeline](bool)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->PBIFDomainSelector,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [invalidateIFPipeline](int)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->PBIFCalibrationStartSpinBox,
      QOverload<double>::of(
          &QDoubleSpinBox::valueChanged),
      q,
      [invalidateIFPipeline](double)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->PBIFCalibrationEndSpinBox,
      QOverload<double>::of(
          &QDoubleSpinBox::valueChanged),
      q,
      [invalidateIFPipeline](double)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->MetaboliteCorrectionCheckBox,
      &QCheckBox::toggled,
      q,
      [invalidateIFPipeline](bool)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->IFSourceProcessingSelector,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [invalidateIFPipeline](int)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->IFLowessSpanSpinBox,
      QOverload<double>::of(
          &QDoubleSpinBox::valueChanged),
      q,
      [invalidateIFPipeline](double)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->ParentFractionProcessingSelector,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [invalidateIFPipeline](int)
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->IFInterpolationSelector,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this](int)
      {
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
      });

  QObject::connect(
      this->IFPreviewButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
          this->previewInputFunction();
      });


  this->StuSelector->setEnabled(false);
  this->CTSelector->setEnabled(false);
  this->PETSelector->setEnabled(false);
  this->SegSelector->setEnabled(false);
  this->segmentSelectAll->setEnabled(false);
  this->saveExcelButton->setEnabled(false);

  this->setPostTACEnabled(false);
  this->updateParametricImagingAvailability();

  this->VOIsegmentSelectAll->setEnabled(false);
  this->VOIMTGAsegmentSelectAll->setEnabled(false);
  this->saveTCMfittedExcelButton->setEnabled(false);
  this->saveMTGAfittedExcelButton->setEnabled(false);
  this->saveTCMExcelButton->setEnabled(false);
  this->saveMTGAExcelButton->setEnabled(false);

  this->TCMResultsTable->setSortingEnabled(true);
  this->MTGAResultsTable->setSortingEnabled(true);

  // Make connections
  QObject::connect( this->PatSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onPatChanged(int)));
  QObject::connect( this->StuSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onStuChanged(int)));
  QObject::connect( this->CTSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onCTChanged(int)) );
  QObject::connect( this->PETSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onPETChanged(int)) );
  QObject::connect( this->SegSelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onSegChanged(int)));
  QObject::connect( this->TACbutton, SIGNAL(clicked(bool)),
    q, SLOT(onTACbutton()));
  QObject::connect( this->segmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onSelectAllbutton()));
  QObject::connect( this->direxcel, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelPathChanged(const QString&)));
  QObject::connect( this->fileexcel, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelPathChanged(const QString&)));
  QObject::connect( this->saveExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveExcelbutton()));
  QObject::connect( this->plotButton, SIGNAL(clicked(bool)),
    q, SLOT(onPlotbutton()));
  QObject::connect( this->plotTCMButton, SIGNAL(clicked(bool)),
    q, SLOT(onPlotTCMbutton()));
  QObject::connect( this->plotMTGAButton, SIGNAL(clicked(bool)),
    q, SLOT(onPlotMTGAbutton()));
  // QObject::connect( this->PlotErrorCheckbox, SIGNAL(toggled(bool)),
  //   q, SLOT(onPlotbutton()));
  QObject::connect(
      this->IFSelector,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      &qSlicerDynamicPETModuleWidget::
          onIFSelectionChanged);
  QObject::connect(
      this->IFStatSelector,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this, q](int)
      {
        this->updateInputFunctionStatus();

        // ROI results depend on the chosen IF statistic.
        q->clearFITdata();
        q->clearFITMTGAdata();

        // Parametric fits also depend on it.
        this->MTGAImgFitSignatures.clear();
        this->TCMImgFitSignatures.clear();

        q->MTGAImgOutcomes.clear();
        q->TCMImgOutcomes.clear();

        q->enableFITbutton();
        q->enableFITMTGAbutton();
        q->enableFITMTGAImgbutton();
        q->enableFITTCMImgbutton();
      });
  QObject::connect( this->VOIsegmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onVOISelectAllbutton()));
  QObject::connect( this->VOIMTGAsegmentSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onVOIMTGASelectAllbutton()));
  QObject::connect( this->ModelsSelectAll, SIGNAL(clicked(bool)),
    q, SLOT(onModelsAllbutton()));
  QObject::connect( this->ModelsMTGASelectAll, SIGNAL(clicked(bool)),
      q, SLOT(onModelsMTGAAllbutton()));
  QObject::connect( this->ModelsTCMSelectAll, SIGNAL(clicked(bool)),
      q, SLOT(onModelsTCMSelectAllbutton()));
  QObject::connect( this->ModelsSelectAllMTGAImg, SIGNAL(clicked(bool)),
      q, SLOT(onModelsSelectAllMTGAImgbutton()));
  QObject::connect( this->ModelsSelectAllTCMImg, SIGNAL(clicked(bool)),
      q, SLOT(onModelsSelectAllTCMImgbutton()));
  QObject::connect( this->FITbutton, SIGNAL(clicked(bool)),
    q, SLOT(onFITbutton()));
  QObject::connect( this->FITbuttonTCMImg, SIGNAL(clicked(bool)),
    q, SLOT(onFITTCMImgbutton()));
  QObject::connect( this->RESETbutton, SIGNAL(clicked(bool)),
    q, SLOT(onResetbutton()));
  QObject::connect( this->FITMTGAbutton, SIGNAL(clicked(bool)),
    q, SLOT(onFITMTGAbutton()));
  QObject::connect( this->FITbuttonMTGAImg, SIGNAL(clicked(bool)),
    q, SLOT(onFITMTGAImgbutton()));
  QObject::connect( this->VOISelector, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onVOISelectionChanged(int)));
  QObject::connect( this->VOISelectorMTGA, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onVOIMTGASelectionChanged(int)));
  QObject::connect( this->direxceltcm, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelTCMPathChanged(const QString&)));
  QObject::connect( this->direxcelmtga, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelMTGAPathChanged(const QString&)));
  QObject::connect( this->fileexceltcm, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelTCMPathChanged(const QString&)));
  QObject::connect( this->fileexcelmtga, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelMTGAPathChanged(const QString&)));
  QObject::connect( this->direxceltcmfitted, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelTCMfittedPathChanged(const QString&)));
  QObject::connect( this->direxcelmtgafitted, SIGNAL(currentPathChanged(const QString&)),
    q, SLOT(onExcelMTGAfittedPathChanged(const QString&)));
  QObject::connect( this->fileexceltcmfitted, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelTCMfittedPathChanged(const QString&)));
  QObject::connect( this->fileexcelmtgafitted, SIGNAL(textChanged(const QString&)),
    q, SLOT(onExcelMTGAfittedPathChanged(const QString&)));
  QObject::connect( this->saveTCMExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveTCMExcelbutton()));
  QObject::connect( this->saveMTGAExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveMTGAExcelbutton()));
  QObject::connect( this->saveTCMfittedExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveTCMfittedExcelbutton()));
  QObject::connect( this->saveMTGAfittedExcelButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveMTGAfittedExcelbutton()));
  QObject::connect(this->olsFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onOLSclicked()));
  QObject::connect(this->olsFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onOLSImgclicked()));
  QObject::connect(this->weightedFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onWLSclicked()));
  QObject::connect(this->weightedFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onWLSImgclicked()));
  QObject::connect(this->robustFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onRLSclicked()));
  QObject::connect(this->robustFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onRLSImgclicked()));
  QObject::connect(this->standardFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onStdFitclicked()));
  QObject::connect(this->standardFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onStdFitImgclicked()));
  QObject::connect(this->weightFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onWFitclicked()));
  QObject::connect(this->weightFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onWFitImgclicked()));
  QObject::connect(this->MTGAModel1, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onMTGAModelBox(int)));
  QObject::connect(this->MTGAModel2, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onMTGAModelBox(int)));
  QObject::connect(this->TCMModel1, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onTCMModelBox(int)));
  QObject::connect(this->TCMModel2, SIGNAL(currentIndexChanged(int)),
    q, SLOT(onTCMModelBox(int)));


  // MTGA controls
  this->setDoubleField(this->framingNormEdit, 0.01, 3600.0, 2);
  this->setDoubleField(this->huberTuneEdit,  1e-3, 10.0,   6);
  this->setDoubleField(this->tolEdit,        1e-12, 1e-1, 12);
  this->setIntField   (this->maxIterEdit,    1,     100000);

  // MTGA imaging controls
  this->setDoubleField(this->framingNormEditImg, 0.01, 3600.0, 2);
  this->setDoubleField(this->huberTuneEditImg,  1e-3, 10.0,   6);
  this->setDoubleField(this->tolEditImg,        1e-12, 1e-1, 12);
  this->setIntField   (this->maxIterEditImg,    1,     100000);
  #ifdef HAVE_OPENMP
  this->setIntField(this->numThreadsMTGA, 1, omp_get_max_threads());
  this->numThreadsMTGA->setText(QString::number(omp_get_max_threads()));
  #else
  this->setIntField(this->numThreadsMTGA, 1, 1);
  #endif

  // TCM params
  this->setDoubleField(this->k1Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k1Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k1Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->k2Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k2Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k2Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->k3Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k3Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k3Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->k4Initial, 0.0, 10.0, 6);
  this->setDoubleField(this->k4Lower,   0.0, 10.0, 6);
  this->setDoubleField(this->k4Upper,   0.0, 10.0, 6);

  this->setDoubleField(this->vbInitial, 0.0, 1.0, 6);
  this->setDoubleField(this->vbLower,   0.0, 1.0, 6);
  this->setDoubleField(this->vbUpper,   0.0, 1.0, 6);

  this->setDoubleField(this->tdInitial, -10.0, 600., 3);
  this->setDoubleField(this->tdLower,   -10.0, 600., 3);
  this->setDoubleField(this->tdUpper,   -10.0, 600., 3);

  // Liver DBIF-only parameters.
  this->setDoubleField(this->liverKaInitial, 0.0, 10.0, 6);
  this->setDoubleField(this->liverKaLower,   0.0, 10.0, 6);
  this->setDoubleField(this->liverKaUpper,   0.0, 10.0, 6);

  this->setDoubleField(this->liverFaInitial, 0.0, 1.0, 6);
  this->setDoubleField(this->liverFaLower,   0.0, 1.0, 6);
  this->setDoubleField(this->liverFaUpper,   0.0, 1.0, 6);

  this->setDoubleField(this->decayConstEdit, 1e-6, 10.0, 10);
  this->setDoubleField(this->timeStepEdit,   0.001, 60.0, 6);

  this->setDoubleField(this->pbrp1Edit, -10.0, 10.0, 6);
  this->setDoubleField(this->pbrp2Edit,   0.0, 10.0, 6);
  this->setDoubleField(this->pbrp3Edit,   0.0, 10.0, 6);

  this->setIntField(this->maxIterTCMEdit, 1, 100000);

  // TCM imaging params
  this->setDoubleField(this->k1InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k1LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k1UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->k2InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k2LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k2UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->k3InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k3LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k3UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->k4InitialImg, 0.0, 10.0, 6);
  this->setDoubleField(this->k4LowerImg,   0.0, 10.0, 6);
  this->setDoubleField(this->k4UpperImg,   0.0, 10.0, 6);

  this->setDoubleField(this->vbInitialImg, 0.0, 1.0, 6);
  this->setDoubleField(this->vbLowerImg,   0.0, 1.0, 6);
  this->setDoubleField(this->vbUpperImg,   0.0, 1.0, 6);

  this->setDoubleField(this->tdInitialImg, -10.0, 600.0, 3);
  this->setDoubleField(this->tdLowerImg,   -10.0, 600.0, 3);
  this->setDoubleField(this->tdUpperImg,   -10.0, 600.0, 3);

  this->setDoubleField(this->decayConstEditImg, 1e-6, 10.0, 10);

  this->setIntField(this->maxIterTCMEditImg, 1, 100000);
  #ifdef HAVE_OPENMP
  this->setIntField(this->numThreadsTCM, 1, omp_get_max_threads());
  this->numThreadsTCM->setText(QString::number(omp_get_max_threads()));
  #else
  this->setIntField(this->numThreadsTCM, 1, 1);
  #endif


  this->TACCollapsibleButton->setCollapsed(true);
  this->TCMCollapsibleButton->setCollapsed(true);
  this->MTGACollapsibleButton->setCollapsed(true);
  this->MTGAStatTestButton->setCollapsed(true);
  this->MTGAStatTestButton->setEnabled(false);
  this->TCMStatTestButton->setCollapsed(true);
  this->TCMStatTestButton->setEnabled(false);
  for (const QString& name : q->checkboxNames)
  {
    QCheckBox* cb = new QCheckBox(name, this->PlotStatsCheckContents);
    this->PlotStatsCheckLayout->addWidget(cb);
  }

  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  mainContext.evalScript(R"PYTHON(
try:
    import pandas as pd
except ImportError:
    import slicer
    slicer.util.pip_install("pandas")
    import pandas as pd

try:
    import xlsxwriter
except ImportError:
    import slicer
    slicer.util.pip_install("xlsxwriter")

import importlib
import importlib.metadata
import sys
import slicer

DPE_HIGHDICOM_REQUIRED_VERSION = "0.28.1"


def DPE_get_highdicom():
    try:
        installed_version = (
            importlib.metadata.version("highdicom")
        )
    except importlib.metadata.PackageNotFoundError:
        installed_version = None

    if installed_version != DPE_HIGHDICOM_REQUIRED_VERSION:

        highdicom_already_loaded = any(
            name == "highdicom"
            or name.startswith("highdicom.")
            for name in sys.modules
        )

        slicer.util.pip_install(
            "--upgrade --force-reinstall --no-deps "
            "highdicom=="
            + DPE_HIGHDICOM_REQUIRED_VERSION
        )

        importlib.invalidate_caches()

        # If another highdicom version was already imported,
        # replacing files on disk does not safely replace the
        # Python classes already resident in this Slicer process.
        if highdicom_already_loaded:
            raise RuntimeError(
                "SlicerDynamicPET installed highdicom "
                + DPE_HIGHDICOM_REQUIRED_VERSION
                + ", but another highdicom version was already "
                  "loaded in this Slicer session.\n\n"
                  "Please restart Slicer once."
            )

    import highdicom as hd

    if hd.__version__ != DPE_HIGHDICOM_REQUIRED_VERSION:
        raise RuntimeError(
            "SlicerDynamicPET requires highdicom "
            + DPE_HIGHDICOM_REQUIRED_VERSION
            + ", but Python loaded version "
            + str(hd.__version__)
            + " from:\n"
            + str(hd.__file__)
        )

    return hd


hd = DPE_get_highdicom()

try:
    import numpy as np
except ImportError:
    import slicer
    slicer.util.pip_install("numpy")

def DPE_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Time(s)", "FrameStart_s", "FrameMid_s", "FrameEnd_s", "Duration_s", "Mean", "Median", "StDev","IQR","Min", "Max", "Q1", "Q3", "Peak", "PeakStDev", "PeakVoxelCount", "VoxelCount","Volume(mm3)","Volume(cm3)"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_saveTCM_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Model", "K1", "k2", "k3", "k4", "ka", "fA", "vb", "td", "Ki", "DV", "AIC", "BIC", "MASE", "chi^2_nu", "BoundHits"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_saveMTGA_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[list[str]]] - sheet name to 2D table
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Model", "Ki", "DV", "Intercept", "R2", "AIC", "MASE"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_generic_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[dict]] - sheet name to list of row dicts
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            if not data:
                continue  # skip empty sheets

            # Create DataFrame from list of dicts — columns inferred automatically
            df = pd.DataFrame(data)

            # Optional: ensure "Time(s)" is first column if present
            if "Time(s)" in df.columns:
                taccols = [x for x in df.columns if "TAC" in x]
                if len(taccols)>0:
                  cols = ["Time(s)"] + taccols + [c for c in df.columns if not np.isin(c, ["Time(s)"]+taccols)]
                else:
                  cols = ["Time(s)"] + [c for c in df.columns if c != "Time(s)"]
                df = df[cols]

            df.to_excel(writer, sheet_name=sheet, index=False)

def DPE_genericMTGA_save_multisheet_excel(filepath, sheet_data_dict):
    """
    filepath: str - full path to xlsx
    sheet_data_dict: dict[str, list[dict]] - sheet name to list of row dicts
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            if not data:
                continue  # skip empty sheets

            # Create DataFrame from list of dicts — columns inferred automatically
            df = pd.DataFrame(data)

            # Optional: ensure "Time(s)" is first column if present
            cols = []
            if "Patlak_x" in df.columns:
              cols += ["Patlak_x", "Patlak_y", "Patlak_fitted"]
            if "Logan_x" in df.columns:
              cols += ["Logan_x", "Logan_y", "Logan_fitted"]
            if "RE_x" in df.columns:
              cols += ["RE_x", "RE_y", "RE_fitted"]
            if len(cols)>0:
              df = df[cols]

            df.to_excel(writer, sheet_name=sheet, index=False)

# --------------------------------------------------------------------------
# DICOM PM spatial-source cache.
#
# Only one PET source is retained. A different geometry UID set
# automatically replaces the previous cache.
# --------------------------------------------------------------------------

_DPE_PMAP_SOURCE_CACHE = {
    "key": None,
    "source_images": None,
}


def DPE_clear_parametric_map_source_cache():
    _DPE_PMAP_SOURCE_CACHE["key"] = None
    _DPE_PMAP_SOURCE_CACHE["source_images"] = None

def DPE_export_parametric_map(
    volume_node_id,
    geometry_instance_uids,
    all_instance_uids,
    output_path,
    series_description,
    series_number,
    quantity_code,
    quantity_meaning,
    method_code,
    method_meaning,
    unit_code,
    unit_meaning,
    derivation_details
):
    import os
    import numpy as np
    import slicer
    import vtk
    hd = DPE_get_highdicom()

    try:
        # ------------------------------------------------------------
        # 1. Retrieve temporary Slicer parametric volume
        # ------------------------------------------------------------
        volume_node = slicer.mrmlScene.GetNodeByID(
            str(volume_node_id)
        )

        if volume_node is None:
            return {
                "ok": False,
                "error":
                    "Temporary parametric volume node was not found."
            }

        # ------------------------------------------------------------
        # 2. Separate spatial construction sources from provenance.
        # ------------------------------------------------------------

        geometry_uid_list = str(
            geometry_instance_uids
        ).split()

        all_uid_list = str(
            all_instance_uids
        ).split()

        if not geometry_uid_list:
            return {
                "ok": False,
                "error":
                    "Source PET geometry UID list is empty."
            }

        if not all_uid_list:
            return {
                "ok": False,
                "error":
                    "Source PET provenance UID list is empty."
            }


        # ------------------------------------------------------------
        # 3. Read only the DICOM objects needed to define spatial
        # geometry.
        #
        # Classic PET:
        #   one temporal frame -> one complete slice stack.
        #
        # Enhanced PET:
        #   one temporal frame -> one multiframe DICOM object.
        #
        # These objects are cached and reused by every parameter map
        # exported from this PET.
        # ------------------------------------------------------------

        source_cache_key = tuple(
            geometry_uid_list
        )

        if (
            _DPE_PMAP_SOURCE_CACHE["key"]
                == source_cache_key
            and
            _DPE_PMAP_SOURCE_CACHE["source_images"]
                is not None
        ):
            source_images = (
                _DPE_PMAP_SOURCE_CACHE[
                    "source_images"
                ]
            )

        else:

            source_paths = []
            seen_paths = set()

            for uid in geometry_uid_list:

                path = (
                    slicer.dicomDatabase
                    .fileForInstance(uid)
                )

                if (
                    path
                    and os.path.isfile(path)
                    and path not in seen_paths
                ):
                    seen_paths.add(path)
                    source_paths.append(path)

            if not source_paths:
                return {
                    "ok": False,
                    "error":
                        "Could not resolve the PET spatial "
                        "reference DICOM instances from the "
                        "Slicer DICOM database."
                }

            # Metadata only: source pixel values are not required
            # to construct the derived parametric volume.
            import pydicom

            source_images = [
                pydicom.dcmread(
                    path,
                    stop_before_pixels=True
                )
                for path in source_paths
            ]

            _DPE_PMAP_SOURCE_CACHE["key"] = (
                source_cache_key
            )

            _DPE_PMAP_SOURCE_CACHE[
                "source_images"
            ] = source_images


        # ------------------------------------------------------------
        # Spatial source type.
        # ------------------------------------------------------------

        source_is_multiframe = [
            int(
                getattr(
                    source,
                    "NumberOfFrames",
                    1
                )
            ) > 1
            for source in source_images
        ]

        has_multiframe_sources = any(
            source_is_multiframe
        )

        has_singleframe_sources = any(
            not value
            for value in source_is_multiframe
        )

        if (
            has_multiframe_sources
            and has_singleframe_sources
        ):
            return {
                "ok": False,
                "error":
                    "PET spatial reference contains a mixture "
                    "of single-frame and multiframe DICOM images."
            }

        # ------------------------------------------------------------
        # Validate spatial source datasets.
        # ------------------------------------------------------------

        if (
            has_multiframe_sources
            and len(source_images) != 1
        ):
            return {
                "ok": False,
                "error":
                    "Enhanced PET spatial reference must contain "
                    "exactly one multiframe DICOM instance."
            }


        # ------------------------------------------------------------
        # Normalize mandatory Type-2 patient/study attributes.
        #
        # Type 2 attributes must exist, but may legitimately have
        # an empty value when unknown.
        # ------------------------------------------------------------

        required_type2_attributes = {
            "PatientName": "",
            "PatientID": "",
            "PatientBirthDate": "",
            "PatientSex": "",
            "StudyDate": "",
            "StudyTime": "",
            "ReferringPhysicianName": "",
            "StudyID": "",
            "AccessionNumber": "",
        }

        for source in source_images:
            for attribute_name, empty_value in \
                    required_type2_attributes.items():

                if not hasattr(source, attribute_name):
                    setattr(
                        source,
                        attribute_name,
                        empty_value
                    )


        first_source = source_images[0]


        required_type1_attributes = [
            "StudyInstanceUID",
            "SeriesInstanceUID",
            "SOPInstanceUID",
            "SOPClassUID",
        ]

        for attribute_name in required_type1_attributes:
            if (
                not hasattr(first_source, attribute_name)
                or not str(
                    getattr(
                        first_source,
                        attribute_name
                    )
                ).strip()
            ):
                return {
                    "ok": False,
                    "error":
                        "Source PET is missing mandatory DICOM "
                        "attribute "
                        + attribute_name
                        + "."
                }


        if not hasattr(
                first_source,
                "FrameOfReferenceUID"
        ):
            return {
                "ok": False,
                "error":
                    "Source PET does not contain "
                    "FrameOfReferenceUID."
            }


        frame_of_reference_uid = str(
            first_source.FrameOfReferenceUID
        )


        # All geometry source images must share the same
        # patient coordinate system.
        for source in source_images:

            if (
                hasattr(source, "FrameOfReferenceUID")
                and
                str(source.FrameOfReferenceUID)
                    != frame_of_reference_uid
            ):
                return {
                    "ok": False,
                    "error":
                        "PET spatial reference contains more than "
                        "one FrameOfReferenceUID."
                }


        constructor_source_images = source_images
        # ------------------------------------------------------------
        # 4. Get parametric values from Slicer
        #
        # Slicer NumPy ordering:
        #   [K, J, I] == [slice, row, column]
        # ------------------------------------------------------------
        pixel_array = (
            slicer.util.arrayFromVolume(volume_node)
            .copy()
            .astype(np.float32)
        )

        if pixel_array.ndim != 3:
            return {
                "ok": False,
                "error":
                    "Parametric map is not a 3D scalar volume."
            }

        # ------------------------------------------------------------
        # 5. Construct KJI -> RAS affine from Slicer's IJK -> RAS
        # ------------------------------------------------------------
        ijk_to_ras_vtk = vtk.vtkMatrix4x4()

        volume_node.GetIJKToRASMatrix(
            ijk_to_ras_vtk
        )

        ijk_to_ras = np.array(
            [
                [
                    ijk_to_ras_vtk.GetElement(r, c)
                    for c in range(4)
                ]
                for r in range(4)
            ],
            dtype=np.float64
        )

        # highdicom's Volume array axes are:
        #
        #   axis 0 = slice  = K
        #   axis 1 = row    = J
        #   axis 2 = column = I
        #
        # Slicer's matrix columns are I, J, K.
        kji_to_ras = np.eye(
            4,
            dtype=np.float64
        )

        kji_to_ras[:3, 0] = ijk_to_ras[:3, 2]
        kji_to_ras[:3, 1] = ijk_to_ras[:3, 1]
        kji_to_ras[:3, 2] = ijk_to_ras[:3, 0]
        kji_to_ras[:3, 3] = ijk_to_ras[:3, 3]

        # highdicom accepts the source affine convention explicitly
        # and converts RAS -> DICOM LPS internally.
        parametric_volume = hd.Volume(
            array=pixel_array,
            affine=kji_to_ras,
            coordinate_system="PATIENT",
            frame_of_reference_uid=
                frame_of_reference_uid,
            from_reference_convention="RAS"
        )

        # ------------------------------------------------------------
        # 6. Real-world quantity definition
        # ------------------------------------------------------------
        quantity = hd.sr.CodedConcept(
            value=str(quantity_code),
            scheme_designator="99SDPET",
            meaning=str(quantity_meaning)
        )

        unit = hd.sr.CodedConcept(
            value=str(unit_code),
            scheme_designator="UCUM",
            meaning=str(unit_meaning)
        )

        finite_values = pixel_array[
            np.isfinite(pixel_array)
        ]

        if finite_values.size == 0:
            return {
                "ok": False,
                "error":
                    "Parametric map contains no finite values."
            }

        value_min = float(
            finite_values.min()
        )

        value_max = float(
            finite_values.max()
        )

        # Avoid a degenerate mapping range.
        if value_max <= value_min:
            value_max = value_min + 1.0e-12

        mapping = hd.pm.RealWorldValueMapping(
            lut_label=str(quantity_code)[:16],
            lut_explanation=
                str(quantity_meaning)[:64],
            value_range=(
                value_min,
                value_max
            ),
            quantity_definition=quantity,
            unit=unit
        )

        # ------------------------------------------------------------
        # 7. Display window
        # ------------------------------------------------------------
        window_width = max(
            value_max - value_min,
            1.0e-12
        )

        window_center = (
            value_min + value_max
        ) / 2.0

        # ------------------------------------------------------------
        # 8. Construct standards-based DICOM PM
        #
        # highdicom uses Volume geometry to match the PM frames
        # against source DICOM frames/images.
        # ------------------------------------------------------------
        pm = hd.pm.ParametricMap(
            source_images=constructor_source_images,

            pixel_array=parametric_volume,

            series_instance_uid=hd.UID(),

            series_number=int(series_number),

            sop_instance_uid=hd.UID(),

            instance_number=1,

            manufacturer="SlicerDynamicPET",

            manufacturer_model_name=
                "SlicerDynamicPET",

            software_versions="development",

            device_serial_number=
                "SlicerDynamicPET",

            contains_recognizable_visual_features=False,

            real_world_value_mappings=[
                mapping
            ],

            voi_lut_transformations=[
                hd.VOILUTTransformation(
                    window_center=
                        window_center,
                    window_width=
                        window_width
                )
            ],

            series_description=
                str(series_description)[:64]
        )

        # ------------------------------------------------------------
        # Enhanced dynamic PET:
        #
        # highdicom required one multiframe instance for geometric
        # PM construction, but the kinetic result was derived from
        # ALL temporal Enhanced PET instances.
        #
        # Record all of those source SOP instances at image level.
        # ------------------------------------------------------------

        # ------------------------------------------------------------
        # Record ALL temporal PET source SOP instances as provenance.
        #
        # These DICOM files do not need to be opened. Their SOP
        # Instance UIDs are already retained on the Slicer sequence.
        # ------------------------------------------------------------

        from pydicom.dataset import Dataset
        from pydicom.sequence import Sequence

        source_sop_class_uid = str(
            first_source.SOPClassUID
        )

        source_references = []

        for source_uid in all_uid_list:

            reference = Dataset()

            reference.ReferencedSOPClassUID = (
                source_sop_class_uid
            )

            reference.ReferencedSOPInstanceUID = (
                str(source_uid)
            )

            source_references.append(
                reference
            )

        pm.SourceImageSequence = Sequence(
            source_references
        )

        # Keep model provenance human-readable.
        derivation_parts = [
            str(method_meaning),
            "methodCode=" + str(method_code),
            "quantity=" + str(quantity_meaning),
        ]

        details = str(derivation_details).strip()

        if details:
            derivation_parts.append(details)

        pm.DerivationDescription = (
            "; ".join(derivation_parts)
        )[:1024]

        # ------------------------------------------------------------
        # 9. Write PM file
        # ------------------------------------------------------------
        output_path = os.path.abspath(
            str(output_path)
        )

        os.makedirs(
            os.path.dirname(output_path),
            exist_ok=True
        )

        if os.path.isfile(output_path):
            os.remove(output_path)

        pm.save_as(
            output_path,
            enforce_file_format=True
        )

        if not os.path.isfile(output_path):
            return {
                "ok": False,
                "error":
                    "Parametric Map construction completed "
                    "but no DICOM file was written."
            }

        return {
            "ok": True,

            "path":
                output_path,

            "geometry_source_count":
                len(source_images),

            "provenance_source_count":
                len(all_uid_list),

            "source_sop_class":
                str(first_source.SOPClassUID),

            "pm_frames":
                int(pm.NumberOfFrames),

            "study_uid":
                str(pm.StudyInstanceUID),

            "frame_of_reference_uid":
                str(pm.FrameOfReferenceUID)
        }

    except Exception as exc:
        import traceback

        return {
            "ok": False,
            "error":
                str(exc)
                + "\n\n"
                + traceback.format_exc()
        }
)PYTHON");

  for (const QString& name : q->ModelsNamesMTGA)
  {
    QCheckBox* cb = new QCheckBox(name, this->ModelsMTGACheckContents);
    this->ModelsMTGACheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsMTGAChanged()));
    QCheckBox* cb2 = new QCheckBox(name, this->ModelsMTGACheckContents);
    this->ModelsCheckLayoutMTGAImg->addWidget(cb2);
    QObject::connect(cb2, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsMTGAImgChanged()));
  }

  for (const QString& name : q->ModelsNamesTCM)
  {
    QCheckBox* cb = new QCheckBox(name, this->ModelsCheckContents);
    this->ModelsCheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsChanged()));
    QCheckBox* cb2 = new QCheckBox(name, this->ModelsCheckContents);
    this->ModelsCheckLayoutTCMImg->addWidget(cb2);
    QObject::connect(cb2, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsTCMImgChanged()));
  }

  // ROI-only optimization-derived liver dual-blood-input model.
  // Keep it out of the voxelwise model list.
  QCheckBox* liverDBIFCheckBox =
      new QCheckBox(
          "Liver DBIF",
          this->ModelsCheckContents);

  liverDBIFCheckBox->setToolTip(
      QObject::tr(
          "ROI-only optimization-derived dual-blood-input liver model. "
          "Uses the selected arterial/aortic whole-blood input and "
          "estimates portal-vein dispersion (ka) and arterial fraction (fA). "
          "The published 7-parameter form is used; input delay is fixed at 0."));

  this->ModelsCheckLayout->
      addWidget(
          liverDBIFCheckBox);

  QObject::connect(
      liverDBIFCheckBox,
      SIGNAL(stateChanged(int)),
      q,
      SLOT(onModelsChanged()));

  this->updateLiverParameterUI();

  for (int i = 0;
       i < q->StatsNames.size();
       ++i)
  {
    this->StatSelector->addItem(
        q->StatsNames[i],
        q->StatsNames[i]);

    this->StatSelectorMTGA->addItem(
        q->StatsNames[i],
        q->StatsNames[i]);

    this->IFStatSelector->addItem(
        q->StatsNames[i],
        q->StatsNames[i]);
  }

  this->StatSelector->setCurrentIndex(0);
  this->StatSelectorMTGA->setCurrentIndex(0);
  this->IFStatSelector->setCurrentIndex(0);

  // --------------------------------------------------------------------------
  // Parametric imaging output controls
  // --------------------------------------------------------------------------

  QObject::connect(
      this->MTGASaveDICOMCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool)
      {
        this->updateMTGAOutputUI();
      });

  QObject::connect(
      this->MTGAShowInSlicerCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool)
      {
        this->updateMTGAOutputUI();
      });

  QObject::connect(
      this->MTGADICOMDirectoryImg,
      &ctkPathLineEdit::currentPathChanged,
      q,
      [this](const QString&)
      {
        this->updateMTGAOutputUI();
      });


  QObject::connect(
      this->TCMSaveDICOMCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool)
      {
        this->updateTCMOutputUI();
      });

  QObject::connect(
      this->TCMShowInSlicerCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool)
      {
        this->updateTCMOutputUI();
      });

  QObject::connect(
      this->TCMDICOMDirectoryImg,
      &ctkPathLineEdit::currentPathChanged,
      q,
      [this](const QString&)
      {
        this->updateTCMOutputUI();
      });


  // --------------------------------------------------------------------------
  // MTGA model-selection controls
  // --------------------------------------------------------------------------

  QObject::connect(
      this->MTGAUseVuongCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool checked)
      {
        this->MTGAVuongAlphaSpinBoxImg
            ->setEnabled(checked);

        this->updateMTGAOptimizationUI();
      });

  QObject::connect(
      this->MTGAReversibleModelComboImg,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this](int)
      {
        this->updateMTGAOptimizationUI();
      });

  QObject::connect(
      this->MTGASelectionCriterionComboImg,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this](int)
      {
        this->updateMTGAOptimizationUI();
      });

  QObject::connect(
      this->GenerateMTGAOptimizedImgButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
        this->generateMTGAOptimizedResult();
      });

  QObject::connect(
      this->RefreshMTGARGBButtonImg,
      &QPushButton::clicked,
      q,
      [this]()
      {
        this->refreshMTGAOptimizedRGB();
      });


  // --------------------------------------------------------------------------
  // TCM model-selection controls
  // --------------------------------------------------------------------------

  QObject::connect(
      this->TCMUseStatTestsCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool checked)
      {
        this->TCMStatAlphaSpinBoxImg
            ->setEnabled(checked);

        this->updateTCMOptimizationUI();
      });

  QObject::connect(
      this->TCMOptimizationModelsSelectAllImg,
      &QPushButton::clicked,
      q,
      [this]()
      {
        int numberOfModels = 0;
        bool allChecked = true;

        for (int i = 0;
             i < this->TCMOptimizationModelsCheckLayoutImg->count();
             ++i)
        {
          QLayoutItem* item =
              this->TCMOptimizationModelsCheckLayoutImg->itemAt(i);

          QCheckBox* cb =
              qobject_cast<QCheckBox*>(item->widget());

          if (!cb)
          {
            continue;
          }

          ++numberOfModels;

          if (!cb->isChecked())
          {
            allChecked = false;
          }
        }

        const bool newState =
            !(numberOfModels > 0 && allChecked);

        for (int i = 0;
             i < this->TCMOptimizationModelsCheckLayoutImg->count();
             ++i)
        {
          QLayoutItem* item =
              this->TCMOptimizationModelsCheckLayoutImg->itemAt(i);

          QCheckBox* cb =
              qobject_cast<QCheckBox*>(item->widget());

          if (!cb)
          {
            continue;
          }

          cb->blockSignals(true);
          cb->setChecked(newState);
          cb->blockSignals(false);
        }

        this->updateTCMOptimizationUI();
      });

  QObject::connect(
      this->TCMSelectionCriterionComboImg,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this](int)
      {
        this->updateTCMOptimizationUI();
      });

  QObject::connect(
      this->GenerateTCMOptimizedImgButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
        this->generateTCMOptimizedResult();
      });

  const std::vector<QCheckBox*> tcmOptimizationParameterBoxes =
  {
    this->TCMOptK1CheckBoxImg,
    this->TCMOptk2CheckBoxImg,
    this->TCMOptk3CheckBoxImg,
    this->TCMOptk4CheckBoxImg,
    this->TCMOptvbCheckBoxImg,
    this->TCMOpttdCheckBoxImg,
    this->TCMOptKiCheckBoxImg,
    this->TCMOptDVCheckBoxImg
  };

  for (QCheckBox* cb : tcmOptimizationParameterBoxes)
  {
    QObject::connect(
        cb,
        &QCheckBox::toggled,
        q,
        [this](bool)
        {
          this->updateTCMOptimizationUI();
        });
  }

  this->updateMTGAOutputUI();
  this->updateTCMOutputUI();

  this->updateMTGAOptimizationUI();
  this->populateTCMOptimizationModels();

  QObject::connect(
      this->BodySupportPreviewButtonImg,
      &QPushButton::clicked,
      q,
      [this]()
      {
          this->previewParametricVoxelSelection();
      });

  QObject::connect(
      this->BodySupportCTThresholdImg,
      QOverload<double>::of(
          &QDoubleSpinBox::valueChanged),
      q,
      [this](double)
      {
          this->invalidateParametricVoxelSelection();
      });

  QObject::connect(
      this->BodySupportMarginImg,
      QOverload<double>::of(
          &QDoubleSpinBox::valueChanged),
      q,
      [this](double)
      {
          this->invalidateParametricVoxelSelection();
      });

  QObject::connect(
      this->BodySupportFillHolesCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool)
      {
          this->invalidateParametricVoxelSelection();
      });
  QObject::connect(
      this->BodySupportSourceImg,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this](int)
      {
          this->invalidateParametricVoxelSelection();
          this->updateBodySupportUI();
      });

  QObject::connect(
      this->BodySupportEnabledCheckBoxImg,
      &QCheckBox::toggled,
      q,
      [this](bool)
      {
          this->invalidateParametricVoxelSelection();
          this->updateBodySupportUI();
      });

  QObject::connect(
      this->BodySupportPETCompositeImg,
      QOverload<int>::of(
          &QComboBox::currentIndexChanged),
      q,
      [this](int)
      {
          this->invalidateParametricVoxelSelection();
      });

  // Initialize enabled/disabled states once setupUi() and all
  // connections are complete.
  this->updateBodySupportUI();

  this->updateInputFunctionUI();
}

void qSlicerDynamicPETModuleWidgetPrivate::populatePatientComboBox() {
  Q_Q(qSlicerDynamicPETModuleWidget);

  vtkIdType currentSelectedID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = this->PatSelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->PatSelector->itemData(currentIndex).value<vtkIdType>();
  }

  this->PatSelector->blockSignals(true);  // Optional: prevent signal emission
  this->PatSelector -> clear();
  this->PatSelector->addItem(QString::fromStdString("None"), QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));


  vtkMRMLScene* scene = q->mrmlScene();

  if (!scene) {
    q->patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    this->populateStudyComboBox(q->patID);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  vtkIdType rootID = shNode->GetSceneItemID();
  if (rootID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID){
    q->patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    this->populateStudyComboBox(q->patID);
    return;
  }

  std::function<void(vtkIdType)> visit;
  std::vector<vtkIdType> patients;
  visit = [&](vtkIdType itemID)
  {
    if (itemID != rootID && shNode->HasItemAttribute(itemID, "Level"))
    {
      std :: string level = shNode->GetItemAttribute(itemID, "Level");
      if (level == "Patient")
      {
        patients.push_back(itemID);
      }
    }

    std::vector<vtkIdType> children;
    shNode->GetItemChildren(itemID, children);
    for (vtkIdType childID : children)
    {
      visit(childID);
    }
  };

  visit(rootID);
  int restoredIndex = 0;
  for (vtkIdType id : patients)
  {
    std::string name = shNode->GetItemName(id);
    this->PatSelector->addItem(QString::fromStdString(name), QVariant::fromValue(id));

    if (id == currentSelectedID)
    {
      restoredIndex = this->PatSelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->PatSelector->setCurrentIndex(restoredIndex);
  }

  this->PatSelector->blockSignals(false);  // Re-enable signals

  vtkIdType passonID = restoredIndex>0 ? currentSelectedID : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  q->patID = passonID;
  this->populateStudyComboBox(q->patID);

  return;
}


void qSlicerDynamicPETModuleWidgetPrivate::populateStudyComboBox(vtkIdType patientID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);
  // Save current selection by study ID
  vtkIdType currentSelectedStudyID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = this->StuSelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedStudyID = this->StuSelector->itemData(currentIndex).value<vtkIdType>();
  }

  this->StuSelector->blockSignals(true);
  this->StuSelector->clear();
  // Add a "None" option first
  this->StuSelector->addItem(QString::fromStdString("None"), QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));
  if (patientID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    this->StuSelector->setEnabled(false);
    // "None" selected — ignore or reset state
    q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    this->populateNodeComboBox(this->CTSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               q->stuID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }
  this->StuSelector->setEnabled(true);


  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    this->StuSelector->setEnabled(false);
    q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    // "None" selected — ignore or reset state
    this->populateNodeComboBox(this->CTSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               q->stuID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    this->StuSelector->setEnabled(false);
    q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    // "None" selected — ignore or reset state
    this->populateNodeComboBox(this->CTSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "CT"
                              );
    this->populateNodeComboBox(this->PETSelector,
                               q->stuID,
                               "vtkMRMLScalarVolumeNode",
                               "PT"
                              );
    this->populateNodeComboBox(this->SegSelector,
                               q->stuID,
                               "vtkMRMLSegmentationNode",
                               ""
                              );
    return;
  }

  // Retrieve the children of the given patient
  std::vector<vtkIdType> children;
  shNode->GetItemChildren(patientID, children);

  // Index to restore, default to the "None" option, which is at index 0
  int restoredIndex = 0;
  for (vtkIdType childID : children)
  {
    if (shNode->HasItemAttribute(childID, "Level"))
    {
      std::string level = shNode->GetItemAttribute(childID, "Level");
      if (level == "Study")
      {
        std::string name = shNode->GetItemName(childID);
        this->StuSelector->addItem(QString::fromStdString(name), QVariant::fromValue(childID));

        // Check if this study was previously selected
        if (childID == currentSelectedStudyID)
        {
          restoredIndex = this->StuSelector->count() - 1;
        }
      }
    }
  }

  // Restore the previous selection (if still present)
  this->StuSelector->setCurrentIndex(restoredIndex);
  this->StuSelector->blockSignals(false);

  vtkIdType passonID = restoredIndex>0 ? currentSelectedStudyID : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  q->stuID = passonID;
  this->populateNodeComboBox(this->CTSelector,
                             q->stuID,
                             "vtkMRMLScalarVolumeNode",
                             "CT"
                            );
  this->populateNodeComboBox(this->PETSelector,
                             q->stuID,
                             "vtkMRMLScalarVolumeNode",
                             "PT"
                            );
  this->populateNodeComboBox(this->SegSelector,
                             q->stuID,
                             "vtkMRMLSegmentationNode",
                             ""
                            );
  return;
}

void qSlicerDynamicPETModuleWidgetPrivate::populateNodeComboBox(
  QComboBox* comboBox,
  vtkIdType parentItemID,
  const char * requiredNodeType,
  const std :: string requiredModality = ""  // Optional: empty string disables filtering
)
{
  Q_Q(qSlicerDynamicPETModuleWidget);
  // Save current selection
  vtkIdType currentSelectedItemID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  int currentIndex = comboBox->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedItemID = comboBox->itemData(currentIndex).value<vtkIdType>();
  }

  comboBox->blockSignals(true);
  comboBox->clear();
  comboBox->addItem("None", QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));

  if (parentItemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    comboBox->setEnabled(false);
    comboBox->blockSignals(false);
    if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
    } else {
      if (requiredModality=="CT")
        q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      if (requiredModality=="PT")
        this->setPETItemID(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID);
      q->enableTACbutton();
    }
    return;
  }

  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    comboBox->setEnabled(false);
    comboBox->blockSignals(false);
    if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
    } else {
      if (requiredModality=="CT")
        q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      if (requiredModality=="PT")
        this->setPETItemID(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID);
      q->enableTACbutton();
    }
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    comboBox->setEnabled(false);
    comboBox->blockSignals(false);
    if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
    } else {
      if (requiredModality=="CT")
        q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      if (requiredModality=="PT")
        this->setPETItemID(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID);
      q->enableTACbutton();
    }
    return;
  }

  if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
    if (q->petID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
      comboBox->setEnabled(false);
      comboBox->blockSignals(false);
      q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
      this->populateSegmentCheckboxes(q->segID);
      return;
    }
  }
  comboBox->setEnabled(true);
  int restoredIndex = 0;

  std::function<void(vtkIdType)> collectItems;
  collectItems = [&](vtkIdType itemID)
  {
    vtkMRMLNode* dataNode =
        shNode->GetItemDataNode(itemID);

    if (dataNode &&
        dataNode->IsA(requiredNodeType))
    {
        const char* internalAttribute =
            dataNode->GetAttribute(
                "SlicerDynamicPET.InternalNode");

        const bool isInternalDynamicPETNode =
            internalAttribute &&
            std::string(internalAttribute) == "1";

        if (!isInternalDynamicPETNode)
        {
            bool hasatt =
                shNode->HasItemAttribute(
                    itemID,
                    "DICOM.Modality");

            std::string modalityAttr =
                hasatt
                ? shNode->GetItemAttribute(
                      itemID,
                      "DICOM.Modality")
                : "";

            bool modalityMatches =
                requiredModality.empty() ||
                requiredModality ==
                    modalityAttr;

            if (requiredModality == "PT")
            {
                std::vector<std::string>
                    pt_attributes =
                        dataNode->
                            GetAttributeNames();

                auto hasseq =
                    std::find(
                        pt_attributes.begin(),
                        pt_attributes.end(),
                        "Sequences.BaseName");

                modalityMatches =
                    modalityMatches &&
                    hasseq !=
                        pt_attributes.end();
            }

            if (modalityMatches)
            {
                std::string name =
                    shNode->GetItemName(
                        itemID);

                comboBox->addItem(
                    QString::fromStdString(name),
                    QVariant::fromValue(itemID));

                if (itemID ==
                    currentSelectedItemID)
                {
                    restoredIndex =
                        comboBox->count() - 1;
                }
            }
        }
    }

    std::vector<vtkIdType> children;
    shNode->GetItemChildren(itemID, children);
    for (vtkIdType childID : children)
    {
      collectItems(childID);
    }
  };

  collectItems(parentItemID);
  comboBox->setCurrentIndex(restoredIndex);
  comboBox->blockSignals(false);

  vtkIdType passonID = restoredIndex>0 ? currentSelectedItemID : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  if (std::string(requiredNodeType)=="vtkMRMLSegmentationNode") {
    q->segID = passonID;
    this->populateSegmentCheckboxes(q->segID);
  } else {
    if (requiredModality=="CT")
      q->ctID = passonID;
    if (requiredModality=="PT")
      this->setPETItemID(passonID);
    q->enableTACbutton();
  }
}


void qSlicerDynamicPETModuleWidgetPrivate::populateSegmentCheckboxes(vtkIdType SegItemID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->segmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->segmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Step 2: Clear existing checkboxes
  this->SegmentCheckContents->blockSignals(true);
  QLayoutItem* item;
  while ((item = this->segmentCheckLayout->takeAt(0)) != nullptr)
  {
    delete item->widget();
    delete item;
  }

  if (SegItemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // this->segmentCheckLayout->setEnabled(false);
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkMRMLScene* scene = q->mrmlScene();
  if (!scene)
  {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode)
  {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(SegItemID));
  // Step 3: Repopulate based on new segmentation node
  if (!segNode) {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  if (q->sequencePETNode == nullptr || q->segSequenceNode == nullptr) {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  vtkSegmentation* segmentation = segNode->GetSegmentation();
  if (!segmentation) {
    this->SegmentCheckContents->blockSignals(false);
    q->segmentIDs.clear();
    q->enableTACbutton();
    this->segmentSelectAll->setEnabled(false);
    return;
  }

  q->segmentIDs.clear();
  std :: vector<std :: string> segmentIDs = segmentation->GetSegmentIDs();

  std::sort(
    segmentIDs.begin(),
    segmentIDs.end(),
    [&](const std::string& a, const std::string& b)
    {
      vtkSegment* segmentA = segmentation->GetSegment(a);
      vtkSegment* segmentB = segmentation->GetSegment(b);

      if (!segmentA)
        return segmentB != nullptr;

      if (!segmentB)
        return false;

      return segmentNameLessCaseInsensitive(
        segmentA->GetName(),
        segmentB->GetName());
    });
  this->segmentDisplayOrder = segmentIDs;
  for (vtkIdType i = 0; i < static_cast<vtkIdType>(segmentIDs.size()); ++i)
  {
    std::string segmentID = segmentIDs[i];
    vtkSegment* vtksegment = segmentation->GetSegment(segmentID);
    std::string segmentName = vtksegment->GetName();

    QCheckBox* checkbox = new QCheckBox(QString::fromStdString(segmentName));
    checkbox->setProperty("SegmentID", QString::fromStdString(segmentID));

    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    checkbox->setChecked(wasSelected);
    this->segmentCheckLayout->addWidget(checkbox);
    QObject::connect(checkbox, SIGNAL(stateChanged(int)),
                 q, SLOT(onSegmentsChanged()));
    if (wasSelected) {
      q->segmentIDs.push_back(QString::fromStdString(segmentID));
    }
  }
  q->enableTACbutton();

  this->segmentCheckLayout->addStretch();
  this->segmentSelectAll->setEnabled(true);
  this->SegmentCheckContents->blockSignals(false);
}


void qSlicerDynamicPETModuleWidgetPrivate::populatePlotSegmentCheckboxes()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslyPlotSelectedIDs;
  for (int i = 0; i < this->PlotsegmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->PlotsegmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslyPlotSelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Step 2: Clear existing checkboxes
  this->PlotSegmentCheckContents->blockSignals(true);
  QLayoutItem* item;
  while ((item = this->PlotsegmentCheckLayout->takeAt(0)) != nullptr)
  {
    delete item->widget();
    delete item;
  }

  if (q->segmentTACsnames.empty() || q->segmentTACs.empty()) {
    this->TACCollapsibleButton->setEnabled(false);
    this->SegmentCheckContents->blockSignals(false);
    return;
  }
  this->TACCollapsibleButton->setEnabled(true);

  for (const std::string& segmentID : this->segmentDisplayOrder)
  {
    auto nameIt = q->segmentTACsnames.find(segmentID);

    // TAC may not have been computed for every segment
    if (nameIt == q->segmentTACsnames.end())
      continue;

    const std::string& segmentName = nameIt->second;

    QCheckBox* checkbox = new QCheckBox(QString::fromStdString(segmentName));
    checkbox->setProperty("SegmentID", QString::fromStdString(segmentID));

    bool wasSelected = previouslyPlotSelectedIDs.contains(QString::fromStdString(segmentID));
    checkbox->setChecked(wasSelected);
    this->PlotsegmentCheckLayout->addWidget(checkbox);
    // QObject::connect(checkbox, SIGNAL(stateChanged(int)),
    //              q, SLOT(onPlotSegmentsChanged()));
  }

  this->PlotsegmentCheckLayout->addStretch();
  this->PlotSegmentCheckContents->blockSignals(false);
}

void qSlicerDynamicPETModuleWidgetPrivate::populateTimeBarMTGA()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (q->numberOfTimepoints <= 0 ||
      q->timePoints.empty() ||
      q->durations.size() != q->timePoints.size())
  {
    return;
  }

  const int lastFrame = q->numberOfTimepoints;

  this->timeOffsetSlider->setMinimum(1);
  this->timeOffsetSlider->setMaximum(lastFrame);
  if (this->timeOffsetSlider->value() < 1 ||
      this->timeOffsetSlider->value() > lastFrame)
  {
    this->timeOffsetSlider->setValue(1);
  }

  this->timeEndSlider->setMinimum(1);
  this->timeEndSlider->setMaximum(lastFrame);
  if (this->timeEndSlider->value() < 1 ||
      this->timeEndSlider->value() > lastFrame)
  {
    this->timeEndSlider->setValue(lastFrame);
  }

  this->TCMEndSlider->setMinimum(1);
  this->TCMEndSlider->setMaximum(lastFrame);
  if (this->TCMEndSlider->value() < 1 ||
      this->TCMEndSlider->value() > lastFrame)
  {
    this->TCMEndSlider->setValue(lastFrame);
  }

  frameEdit->setReadOnly(true);
  timeSecEdit->setReadOnly(true);
  timeMinEdit->setReadOnly(true);

  auto updateMTGAEndInfo = [this, q]()
  {
    const int index = this->timeEndSlider->value();
    if (index < 1 || index > static_cast<int>(q->timePoints.size()))
      return;
    const double endSec = q->timePoints[index - 1];
    this->timeEndInfoEdit->setText(
        QObject::tr("Frame %1 | end %2 s (%3 min)")
            .arg(index)
            .arg(endSec, 0, 'f', 2)
            .arg(endSec / 60.0, 0, 'f', 2));
  };

  auto updateTCMEndInfo = [this, q]()
  {
    const int index = this->TCMEndSlider->value();
    if (index < 1 || index > static_cast<int>(q->timePoints.size()))
      return;
    const double endSec = q->timePoints[index - 1];
    this->TCMEndInfoEdit->setText(
        QObject::tr("Frame %1 | complete through %2 s (%3 min)")
            .arg(index)
            .arg(endSec, 0, 'f', 2)
            .arg(endSec / 60.0, 0, 'f', 2));
  };

  QObject::disconnect(
      this->timeOffsetSlider,
      SIGNAL(valueChanged(int)),
      q,
      SLOT(onSliderChanged(int)));
  QObject::connect(
      this->timeOffsetSlider,
      SIGNAL(valueChanged(int)),
      q,
      SLOT(onSliderChanged(int)));

  QObject::disconnect(this->timeEndSlider, nullptr, q, nullptr);
  QObject::connect(
      this->timeEndSlider,
      &QSlider::valueChanged,
      q,
      [this, q, updateMTGAEndInfo](int index)
      {
        if (index < this->timeOffsetSlider->value())
        {
          this->timeOffsetSlider->setValue(index);
        }
        updateMTGAEndInfo();
        q->enableFITMTGAbutton();
      });

  QObject::disconnect(this->TCMEndSlider, nullptr, q, nullptr);
  QObject::connect(
      this->TCMEndSlider,
      &QSlider::valueChanged,
      q,
      [q, updateTCMEndInfo](int)
      {
        updateTCMEndInfo();
        q->enableFITbutton();
      });

  q->onSliderChanged(this->timeOffsetSlider->value());
  updateMTGAEndInfo();
  updateTCMEndInfo();
}

void qSlicerDynamicPETModuleWidgetPrivate::populateTimeBarMTGAImg()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (q->numberOfTimepoints <= 0 ||
      q->timePoints.empty() ||
      q->durations.size() != q->timePoints.size())
  {
    return;
  }

  const int lastFrame = q->numberOfTimepoints;

  this->timeOffsetSliderImg->setMinimum(1);
  this->timeOffsetSliderImg->setMaximum(lastFrame);
  if (this->timeOffsetSliderImg->value() < 1 ||
      this->timeOffsetSliderImg->value() > lastFrame)
  {
    this->timeOffsetSliderImg->setValue(1);
  }

  this->timeEndSliderImg->setMinimum(1);
  this->timeEndSliderImg->setMaximum(lastFrame);
  if (this->timeEndSliderImg->value() < 1 ||
      this->timeEndSliderImg->value() > lastFrame)
  {
    this->timeEndSliderImg->setValue(lastFrame);
  }

  this->TCMEndSliderImg->setMinimum(1);
  this->TCMEndSliderImg->setMaximum(lastFrame);
  if (this->TCMEndSliderImg->value() < 1 ||
      this->TCMEndSliderImg->value() > lastFrame)
  {
    this->TCMEndSliderImg->setValue(lastFrame);
  }

  this->frameEditImg->setReadOnly(true);
  this->timeSecEditImg->setReadOnly(true);
  this->timeMinEditImg->setReadOnly(true);

  auto updateMTGAEndInfo = [this, q]()
  {
    const int index = this->timeEndSliderImg->value();
    if (index < 1 || index > static_cast<int>(q->timePoints.size()))
      return;
    const double endSec = q->timePoints[index - 1];
    this->timeEndInfoEditImg->setText(
        QObject::tr("Frame %1 | end %2 s (%3 min)")
            .arg(index)
            .arg(endSec, 0, 'f', 2)
            .arg(endSec / 60.0, 0, 'f', 2));
  };

  auto updateTCMEndInfo = [this, q]()
  {
    const int index = this->TCMEndSliderImg->value();
    if (index < 1 || index > static_cast<int>(q->timePoints.size()))
      return;
    const double endSec = q->timePoints[index - 1];
    this->TCMEndInfoEditImg->setText(
        QObject::tr("Frame %1 | complete through %2 s (%3 min)")
            .arg(index)
            .arg(endSec, 0, 'f', 2)
            .arg(endSec / 60.0, 0, 'f', 2));
  };

  QObject::disconnect(
      this->timeOffsetSliderImg,
      SIGNAL(valueChanged(int)),
      q,
      SLOT(onSliderImgChanged(int)));
  QObject::connect(
      this->timeOffsetSliderImg,
      SIGNAL(valueChanged(int)),
      q,
      SLOT(onSliderImgChanged(int)));

  QObject::disconnect(this->timeEndSliderImg, nullptr, q, nullptr);
  QObject::connect(
      this->timeEndSliderImg,
      &QSlider::valueChanged,
      q,
      [this, q, updateMTGAEndInfo](int index)
      {
        if (index < this->timeOffsetSliderImg->value())
        {
          this->timeOffsetSliderImg->setValue(index);
        }
        updateMTGAEndInfo();
        q->enableFITMTGAImgbutton();
      });

  QObject::disconnect(this->TCMEndSliderImg, nullptr, q, nullptr);
  QObject::connect(
      this->TCMEndSliderImg,
      &QSlider::valueChanged,
      q,
      [q, updateTCMEndInfo](int)
      {
        updateTCMEndInfo();
        q->enableFITTCMImgbutton();
      });

  q->onSliderImgChanged(this->timeOffsetSliderImg->value());
  updateMTGAEndInfo();
  updateTCMEndInfo();
}

void qSlicerDynamicPETModuleWidgetPrivate::populateIF()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  const std::string previousID = q->IFID;

  this->IFSelector->blockSignals(true);
  this->IFSelector->clear();

  this->IFSelector->addItem(
      QObject::tr("None"),
      QString());

  int restoredIndex = 0;

  if (!q->segmentTACsnames.empty() &&
      !q->segmentTACs.empty())
  {
    for (const std::string& segmentID :
         this->segmentDisplayOrder)
    {
      auto it =
          q->segmentTACsnames.find(segmentID);

      if (it ==
          q->segmentTACsnames.end())
      {
        continue;
      }

      const std::string& displayName =
          it->second;

      this->IFSelector->addItem(
          QString::fromStdString(displayName),
          QString::fromStdString(segmentID));

      if (segmentID == previousID)
      {
        restoredIndex =
            this->IFSelector->count() - 1;
      }
    }
  }

  this->IFSelector->setCurrentIndex(
      restoredIndex);

  this->IFSelector->blockSignals(false);

  if (restoredIndex > 0)
  {
    q->IFID =
        this->IFSelector
            ->itemData(restoredIndex)
            .toString()
            .toStdString();
  }
  else
  {
    q->IFID.clear();
  }

  // Both ROI-modeling branches consume the same
  // canonical input function.
  this->populateVOI(q->IFID);
  this->populateVOIMTGA(q->IFID);

  this->updateInputFunctionUI();
}


std::string
qSlicerDynamicPETModuleWidgetPrivate::
selectedIFStatistic()
{
  const int index =
      this->IFStatSelector->currentIndex();

  if (index < 0)
  {
    return std::string();
  }

  return
      this->IFStatSelector
          ->itemData(index)
          .toString()
          .toStdString();
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
buildCurrentSegmentInputFunction(
    std::vector<double>& values,
    QString* errorMessage,
    std::vector<bool>* keepMask)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  values.clear();
  if (keepMask)
  {
    keepMask->clear();
  }

  if (q->IFID.empty())
  {
    if (errorMessage)
    {
      *errorMessage =
          QObject::tr(
              "No input-function segment is selected.");
    }

    return false;
  }

  const auto it =
      q->segmentTACs.find(q->IFID);

  if (it == q->segmentTACs.end())
  {
    if (errorMessage)
    {
      *errorMessage =
          QObject::tr(
              "The selected input-function TAC is not available.");
    }

    return false;
  }

  const std::string statistic =
      this->selectedIFStatistic();

  if (statistic.empty())
  {
    if (errorMessage)
    {
      *errorMessage =
          QObject::tr(
              "No input-function statistic is selected.");
    }

    return false;
  }

  const std::vector<VoxelStatistics>& stats =
      it->second;

  if (stats.size() !=
      q->timePoints.size())
  {
    if (errorMessage)
    {
      *errorMessage =
          QObject::tr(
              "Input-function frame count does not match the dynamic PET.");
    }

    return false;
  }

  values.reserve(stats.size());

  for (const VoxelStatistics& frameStats :
       stats)
  {
    double value = 0.0;

    if (statistic == "Mean")
    {
      value = frameStats.mean;
    }
    else if (statistic == "Median")
    {
      value = frameStats.median;
    }
    else if (statistic == "Peak")
    {
      value = frameStats.peak;
    }
    else
    {
      if (errorMessage)
      {
        *errorMessage =
            QObject::tr(
                "Unknown input-function statistic.");
      }

      values.clear();
      return false;
    }

    values.push_back(value);
    if (keepMask)
    {
      keepMask->push_back(frameStats.keep);
    }
  }

  return true;
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateInputFunctionUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const int source =
        this->IFSourceSelector->currentIndex();

    const bool segmentSource = source == 0;
    const bool csvSource = source == 1;

    const bool primarySourceSelected =
        segmentSource
        ? !q->IFID.empty()
        : (this->externalIFTimesSec.size() >= 2 &&
           this->externalIFTimesSec.size() ==
               this->externalIFConcentrations.size());

    // Keep source selection available, but do not expose an active
    // processing pipeline when there is no actual IF selected.
    this->IFAdvancedCollapsibleButton->setEnabled(
        primarySourceSelected);

    if (!primarySourceSelected)
    {
        this->IFAdvancedCollapsibleButton->setCollapsed(true);
    }

    const bool lowessSelected =
        this->IFSourceProcessingSelector->currentIndex() == 1;
    this->IFLowessSpanLabel->setVisible(lowessSelected);
    this->IFLowessSpanSpinBox->setVisible(lowessSelected);

    this->IFLabel->setVisible(segmentSource);
    this->IFSelector->setVisible(segmentSource);
    this->IFStatLabel->setVisible(segmentSource);
    this->IFStatSelector->setVisible(segmentSource);

    this->IFCSVFileLabel->setVisible(csvSource);
    this->IFCSVPathEdit->setVisible(csvSource);
    this->IFCSVUnitSelector->setVisible(csvSource);
    this->IFCSVBrowseButton->setVisible(csvSource);

    // Segment IDIFs are always whole blood. Biological-domain choice
    // is only meaningful for an external curve.
    this->IFCurveTypeLabel->setVisible(csvSource);
    this->IFCurveTypeSelector->setVisible(csvSource);

    const IFCurveDomain domain =
        this->selectedIFCurveDomain();

    const bool parentPlasmaSource =
        csvSource &&
        domain == IFCurveDomain::ParentPlasma;

    this->IFWholeBloodFileLabel->
        setVisible(parentPlasmaSource);
    this->IFWholeBloodPathEdit->
        setVisible(parentPlasmaSource);
    this->IFWholeBloodUnitSelector->
        setVisible(parentPlasmaSource);
    this->IFWholeBloodBrowseButton->
        setVisible(parentPlasmaSource);
    this->IFWholeBloodPreviewButton->
        setVisible(parentPlasmaSource);

    // Parent plasma is already metabolite corrected and cannot be used
    // to derive total whole blood with PBR.
    const bool pbrRelevant =
        !parentPlasmaSource;

    this->pbrpLabel->setEnabled(pbrRelevant);
    this->pbrp1Edit->setEnabled(pbrRelevant);
    this->pbrp2Edit->setEnabled(pbrRelevant);
    this->pbrp3Edit->setEnabled(pbrRelevant);

    if (parentPlasmaSource)
    {
        this->PBIFOptionCheckBox->blockSignals(true);
        this->PBIFOptionCheckBox->setChecked(false);
        this->PBIFOptionCheckBox->blockSignals(false);

        this->MetaboliteCorrectionCheckBox->blockSignals(true);
        this->MetaboliteCorrectionCheckBox->setChecked(false);
        this->MetaboliteCorrectionCheckBox->blockSignals(false);
    }

    this->PBIFOptionCheckBox->
        setEnabled(!parentPlasmaSource);

    const bool pbifEnabled =
        !parentPlasmaSource &&
        this->PBIFOptionCheckBox->isChecked();

    this->PBIFFileLabel->setEnabled(pbifEnabled);
    this->PBIFPathEdit->setEnabled(pbifEnabled);
    this->PBIFBrowseButton->setEnabled(pbifEnabled);
    this->PBIFPreviewButton->setEnabled(
        this->pbifTimesSec.size() >= 2);
    this->PBIFDomainLabel->setEnabled(pbifEnabled);
    this->PBIFDomainSelector->setEnabled(pbifEnabled);
    this->PBIFCalibrationIntervalLabel->setEnabled(pbifEnabled);
    this->PBIFCalibrationStartSpinBox->setEnabled(pbifEnabled);
    this->PBIFCalibrationToLabel->setEnabled(pbifEnabled);
    this->PBIFCalibrationEndSpinBox->setEnabled(pbifEnabled);

    this->MetaboliteCorrectionCheckBox->
        setEnabled(!parentPlasmaSource);

    const bool parentCorrectionEnabled =
        !parentPlasmaSource &&
        this->MetaboliteCorrectionCheckBox->isChecked();

    this->ParentFractionFileLabel->
        setEnabled(parentCorrectionEnabled);
    this->ParentFractionPathEdit->
        setEnabled(parentCorrectionEnabled);
    this->ParentFractionBrowseButton->
        setEnabled(parentCorrectionEnabled);
    this->ParentFractionPreviewButton->
        setEnabled(
            this->parentFractionTimesSec.size() >= 2);
    this->ParentFractionProcessingLabel->
        setEnabled(parentCorrectionEnabled);
    this->ParentFractionProcessingSelector->
        setEnabled(parentCorrectionEnabled);

    this->IFWholeBloodPreviewButton->setEnabled(
        parentPlasmaSource &&
        this->externalWholeBloodTimesSec.size() >= 2);

    QString error;
    InputFunctionResult result;

    const bool validForPlasmaModel =
        this->buildCurrentInputFunction(
            result,
            false,
            &error);

    this->IFPreviewButton->setEnabled(
        validForPlasmaModel);

    if (validForPlasmaModel)
    {
        QString status;

        if (segmentSource)
        {
            QString name =
                QString::fromStdString(q->IFID);

            const auto it =
                q->segmentTACsnames.find(q->IFID);

            if (it != q->segmentTACsnames.end())
            {
                name =
                    QString::fromStdString(
                        it->second);
            }

            status =
                QObject::tr(
                    "Valid IDIF: %1 | %2 | Whole blood")
                    .arg(name)
                    .arg(
                        QString::fromStdString(
                            this->selectedIFStatistic()));
        }
        else
        {
            QString domainName =
                QObject::tr("Whole blood");

            if (domain == IFCurveDomain::TotalPlasma)
            {
                domainName =
                    QObject::tr("Total plasma");
            }
            else if (domain == IFCurveDomain::ParentPlasma)
            {
                domainName =
                    QObject::tr("Parent plasma");
            }

            status =
                QObject::tr(
                    "Valid external input function: %1 samples | %2")
                    .arg(this->externalIFTimesSec.size())
                    .arg(domainName);

            if (this->externalIFZeroAnchorAdded)
            {
                status +=
                    QObject::tr(
                        " | assumed (0 s, 0) anchor");
            }
        }

        status += QObject::tr(" | model scale=%1")
            .arg(this->activityUnitLabel(this->petStoredActivityUnit()));

        if (result.sourceProcessingApplied)
        {
            status += QObject::tr(" | processing=%1")
                .arg(result.sourceProcessingLabel);
            if (this->IFSourceProcessingSelector->currentIndex() == 2)
            {
                status += QObject::tr(" | Feng tau=%1 min, lambdas=%2/%3/%4 min^-1")
                    .arg(result.fengParameters.tau, 0, 'g', 4)
                    .arg(result.fengParameters.lambda1, 0, 'g', 4)
                    .arg(result.fengParameters.lambda2, 0, 'g', 4)
                    .arg(result.fengParameters.lambda3, 0, 'g', 4);
            }
        }

        if (this->selectedDisplayActivityUnit() != this->petStoredActivityUnit())
        {
            status += QObject::tr(" | display=%1")
                .arg(this->activityUnitLabel(this->selectedDisplayActivityUnit()));
        }

        if (!parentPlasmaSource)
        {
            double pbr0 = 0.0;
            if (this->pbrAtTime(0.0, pbr0, nullptr))
            {
                status += QObject::tr(" | PBR(0)=%1")
                    .arg(pbr0, 0, 'g', 6);
            }
        }

        if (result.pbifApplied)
        {
            status +=
                QObject::tr(
                    " | PBIF scale=%1")
                    .arg(
                        result.pbifScale,
                        0,
                        'g',
                        6);
        }

        if (result.applyParentFraction)
        {
            status +=
                QObject::tr(
                    " | parent-fraction corrected");
        }

        if (parentPlasmaSource &&
            !result.hasWholeBlood)
        {
            status +=
                QObject::tr(
                    " | TCM requires companion whole blood");
        }

        this->IFStatusLabel->setText(status);
    }
    else
    {
        this->IFStatusLabel->setText(error);
    }

    // External CSV and processed IF stages do not yet carry
    // propagated uncertainty for voxelwise WLS.
    if (csvSource ||
        pbifEnabled ||
        parentCorrectionEnabled)
    {
        if (this->weightedFitCheckBoxImg->isChecked())
        {
            this->weightedFitCheckBoxImg->setChecked(false);
        }

        if (this->weightFitCheckBoxImg->isChecked())
        {
            this->weightFitCheckBoxImg->setChecked(false);
        }

        this->weightedFitCheckBoxImg->setEnabled(false);
        this->weightFitCheckBoxImg->setEnabled(false);
    }
    else
    {
        this->weightedFitCheckBoxImg->setEnabled(true);
        this->weightFitCheckBoxImg->setEnabled(true);
    }
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateInputFunctionStatus()
{
    this->updateInputFunctionUI();
    this->updateROIModelingAvailability();
    this->updateParametricImagingAvailability();
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateROIModelingAvailability()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const int roiIndex =
        this->PlotsTabWidget->indexOf(
            this->ROIModelingWidget);

    if (roiIndex < 0)
    {
        return;
    }

    const bool tacReady =
        !q->segmentTACs.empty() &&
        !q->segmentTACsnames.empty();

    bool inputReady = false;

    if (tacReady)
    {
        QString ignoredError;
        inputReady = this->hasValidInputFunction(
            &ignoredError,
            false);
    }

    this->PlotsTabWidget->setTabEnabled(
        roiIndex,
        tacReady && inputReady);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateParametricImagingAvailability()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const int imagingIndex =
        this->PlotsTabWidget->indexOf(this->ImagingWidget);

    if (imagingIndex < 0)
    {
        return;
    }

    bool enabled = false;

    if (q->petID !=
            vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID &&
        q->sequencePETNode != nullptr)
    {
        QString ignoredError;
        enabled = this->hasValidInputFunction(
            &ignoredError,
            false);
    }

    this->PlotsTabWidget->setTabEnabled(
        imagingIndex,
        enabled);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
invalidateInputFunctionResults()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    // Any pipeline edit makes the currently displayed QC curves stale.
    this->inputFunctionCacheValid = false;
    this->cachedInputFunction = InputFunctionResult{};
    this->removeInputFunctionPreview();

    q->clearFITdata();
    q->clearFITMTGAdata();

    this->MTGAImgFitSignatures.clear();
    this->TCMImgFitSignatures.clear();

    q->MTGAImgOutcomes.clear();
    q->TCMImgOutcomes.clear();

    q->enableFITbutton();
    q->enableFITMTGAbutton();
    q->enableFITMTGAImgbutton();
    q->enableFITTCMImgbutton();
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
loadTwoColumnCurveCSV(
    const QString& filePath,
    const QString& secondColumnName,
    double minimumValue,
    double maximumValue,
    double insertedZeroTimeValue,
    std::vector<double>& timesSec,
    std::vector<double>& values,
    bool& zeroAnchorAdded,
    QString* errorMessage) const
{
    QFile file(filePath);

    if (!file.open(
            QIODevice::ReadOnly |
            QIODevice::Text))
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Could not open the selected CSV file.");
        }
        return false;
    }

    QTextStream stream(&file);

    if (stream.atEnd())
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "The CSV file is empty.");
        }
        return false;
    }

    QString header =
        stream.readLine().trimmed();

    if (!header.isEmpty() &&
        header.at(0).unicode() == 0xFEFF)
    {
        header.remove(0, 1);
    }

    const QStringList headerColumns =
        header.split(
            ',',
            Qt::KeepEmptyParts);

    if (headerColumns.size() != 2 ||
        headerColumns[0].trimmed() != "time_s" ||
        headerColumns[1].trimmed() != secondColumnName)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Invalid CSV header.\n\nExpected exactly:\ntime_s,%1")
                    .arg(secondColumnName);
        }
        return false;
    }

    std::vector<double> parsedTimes;
    std::vector<double> parsedValues;

    int lineNumber = 1;

    while (!stream.atEnd())
    {
        ++lineNumber;

        const QString line =
            stream.readLine().trimmed();

        if (line.isEmpty())
        {
            continue;
        }

        const QStringList columns =
            line.split(
                ',',
                Qt::KeepEmptyParts);

        if (columns.size() != 2)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Invalid CSV row at line %1. "
                        "Exactly two values are required.")
                        .arg(lineNumber);
            }
            return false;
        }

        bool timeOK = false;
        bool valueOK = false;

        const double time =
            columns[0].trimmed().toDouble(&timeOK);

        const double value =
            columns[1].trimmed().toDouble(&valueOK);

        if (!timeOK ||
            !valueOK ||
            !std::isfinite(time) ||
            !std::isfinite(value))
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Non-numeric or non-finite value "
                        "at CSV line %1.")
                        .arg(lineNumber);
            }
            return false;
        }

        if (time < 0.0)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Negative time at CSV line %1.")
                        .arg(lineNumber);
            }
            return false;
        }

        if (value < minimumValue ||
            value > maximumValue)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Value outside the allowed range "
                        "at CSV line %1.")
                        .arg(lineNumber);
            }
            return false;
        }

        if (!parsedTimes.empty() &&
            time <= parsedTimes.back())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "CSV times must be strictly increasing. "
                        "Problem at line %1.")
                        .arg(lineNumber);
            }
            return false;
        }

        parsedTimes.push_back(time);
        parsedValues.push_back(value);
    }

    if (parsedTimes.size() < 2)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "The CSV must contain at least two samples.");
        }
        return false;
    }

    bool addedAnchor = false;

    if (parsedTimes.front() > 1e-6)
    {
        parsedTimes.insert(
            parsedTimes.begin(),
            0.0);
        parsedValues.insert(
            parsedValues.begin(),
            insertedZeroTimeValue);
        addedAnchor = true;
    }
    else
    {
        parsedTimes.front() = 0.0;
    }

    timesSec = std::move(parsedTimes);
    values = std::move(parsedValues);
    zeroAnchorAdded = addedAnchor;

    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
loadExternalInputFunctionCSV(
    const QString& filePath,
    QString* errorMessage)
{
    std::vector<double> times;
    std::vector<double> values;
    bool zeroAnchorAdded = false;

    if (!this->loadTwoColumnCurveCSV(
            filePath,
            "concentration",
            0.0,
            std::numeric_limits<double>::infinity(),
            0.0,
            times,
            values,
            zeroAnchorAdded,
            errorMessage))
    {
        return false;
    }

    this->externalIFPath = filePath;
    this->externalIFTimesSec = std::move(times);
    this->externalIFConcentrations = std::move(values);
    this->externalIFZeroAnchorAdded = zeroAnchorAdded;

    // A companion whole-blood curve is specifically paired with an
    // already-parent-plasma source. A newly loaded primary curve must
    // never inherit that auxiliary curve silently.
    this->externalWholeBloodPath.clear();
    this->externalWholeBloodTimesSec.clear();
    this->externalWholeBloodConcentrations.clear();
    this->externalWholeBloodZeroAnchorAdded = false;
    this->IFWholeBloodPathEdit->clear();

    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
loadCompanionWholeBloodCSV(
    const QString& filePath,
    QString* errorMessage)
{
    std::vector<double> times;
    std::vector<double> values;
    bool zeroAnchorAdded = false;

    if (!this->loadTwoColumnCurveCSV(
            filePath,
            "concentration",
            0.0,
            std::numeric_limits<double>::infinity(),
            0.0,
            times,
            values,
            zeroAnchorAdded,
            errorMessage))
    {
        return false;
    }

    this->externalWholeBloodPath = filePath;
    this->externalWholeBloodTimesSec = std::move(times);
    this->externalWholeBloodConcentrations = std::move(values);
    this->externalWholeBloodZeroAnchorAdded = zeroAnchorAdded;

    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
loadPBIFCSV(
    const QString& filePath,
    QString* errorMessage)
{
    std::vector<double> times;
    std::vector<double> values;
    bool zeroAnchorAdded = false;

    if (!this->loadTwoColumnCurveCSV(
            filePath,
            "template_value",
            0.0,
            std::numeric_limits<double>::infinity(),
            0.0,
            times,
            values,
            zeroAnchorAdded,
            errorMessage))
    {
        return false;
    }

    this->pbifPath = filePath;
    this->pbifTimesSec = std::move(times);
    this->pbifTemplateValues = std::move(values);
    this->pbifZeroAnchorAdded = zeroAnchorAdded;

    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
loadParentFractionCSV(
    const QString& filePath,
    QString* errorMessage)
{
    std::vector<double> times;
    std::vector<double> values;
    bool zeroAnchorAdded = false;

    if (!this->loadTwoColumnCurveCSV(
            filePath,
            "parent_fraction",
            0.0,
            1.0,
            1.0,
            times,
            values,
            zeroAnchorAdded,
            errorMessage))
    {
        return false;
    }

    this->parentFractionPath = filePath;
    this->parentFractionTimesSec = std::move(times);
    this->parentFractionValues = std::move(values);
    this->parentFractionZeroAnchorAdded = zeroAnchorAdded;

    return true;
}

double
qSlicerDynamicPETModuleWidgetPrivate::
interpolateInputFunction(
    const std::vector<double>& times,
    const std::vector<double>& values,
    double targetTime,
    const std::string& interpolationType) const
{
    if (times.size() < 2 ||
        times.size() != values.size() ||
        targetTime < times.front() ||
        targetTime > times.back())
    {
        return
            std::numeric_limits<double>::
                quiet_NaN();
    }

    if (targetTime == times.back())
    {
        return values.back();
    }

    const auto upper =
        std::upper_bound(
            times.begin(),
            times.end(),
            targetTime);

    if (upper == times.begin())
    {
        return values.front();
    }

    const size_t right =
        static_cast<size_t>(
            std::distance(
                times.begin(),
                upper));

    const size_t left =
        right - 1;

    if (interpolationType == "const")
    {
        return values[left];
    }

    const double t1 = times[left];
    const double t2 = times[right];
    const double y1 = values[left];
    const double y2 = values[right];

    return y1 +
        (y2 - y1) *
        (targetTime - t1) /
        (t2 - t1);
}

double
qSlicerDynamicPETModuleWidgetPrivate::
integrateInputFunctionOverInterval(
    const std::vector<double>& times,
    const std::vector<double>& values,
    double startTime,
    double endTime,
    const std::string& interpolationType) const
{
    if (times.size() < 2 ||
        times.size() != values.size() ||
        endTime <= startTime ||
        startTime < times.front() ||
        endTime > times.back())
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> knots;
    knots.push_back(startTime);

    for (double t : times)
    {
        if (t > startTime &&
            t < endTime)
        {
            knots.push_back(t);
        }
    }

    knots.push_back(endTime);

    double integral = 0.0;

    for (size_t i = 0;
         i + 1 < knots.size();
         ++i)
    {
        const double a = knots[i];
        const double b = knots[i + 1];

        if (interpolationType == "const")
        {
            const double y =
                this->interpolateInputFunction(
                    times,
                    values,
                    0.5 * (a + b),
                    interpolationType);

            if (!std::isfinite(y))
            {
                return
                    std::numeric_limits<double>::quiet_NaN();
            }

            integral += y * (b - a);
        }
        else
        {
            const double ya =
                this->interpolateInputFunction(
                    times,
                    values,
                    a,
                    interpolationType);

            const double yb =
                this->interpolateInputFunction(
                    times,
                    values,
                    b,
                    interpolationType);

            if (!std::isfinite(ya) ||
                !std::isfinite(yb))
            {
                return
                    std::numeric_limits<double>::quiet_NaN();
            }

            integral +=
                0.5 * (ya + yb) * (b - a);
        }
    }

    return integral;
}

double
qSlicerDynamicPETModuleWidgetPrivate::
averageInputFunctionOverInterval(
    const std::vector<double>& times,
    const std::vector<double>& values,
    double startTime,
    double endTime,
    const std::string& interpolationType) const
{
    const double integral =
        this->integrateInputFunctionOverInterval(
            times,
            values,
            startTime,
            endTime,
            interpolationType);

    if (!std::isfinite(integral) ||
        endTime <= startTime)
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    return integral / (endTime - startTime);
}

double
qSlicerDynamicPETModuleWidgetPrivate::
integrateFrameAverageCurveOverInterval(
    const std::vector<double>& frameValues,
    double startTime,
    double endTime)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (frameValues.size() != q->durations.size() ||
        q->durations.empty() ||
        endTime <= startTime ||
        startTime < 0.0 ||
        q->timePoints.empty() ||
        endTime > q->timePoints.back())
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    double integral = 0.0;
    double frameStart = 0.0;

    for (size_t i = 0;
         i < q->durations.size();
         ++i)
    {
        const double frameEnd =
            frameStart + q->durations[i];

        const double overlapStart =
            std::max(startTime, frameStart);

        const double overlapEnd =
            std::min(endTime, frameEnd);

        if (overlapEnd > overlapStart)
        {
            integral +=
                frameValues[i] *
                (overlapEnd - overlapStart);
        }

        frameStart = frameEnd;
    }

    return integral;
}

double
qSlicerDynamicPETModuleWidgetPrivate::
evaluateFrameCurve(
    const std::vector<double>& frameValues,
    double targetTime,
    const std::string& interpolationType)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (frameValues.size() != q->durations.size() ||
        frameValues.empty() ||
        q->timePoints.size() != frameValues.size() ||
        targetTime < 0.0 ||
        targetTime > q->timePoints.back())
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    if (interpolationType == "const")
    {
        double frameStart = 0.0;

        for (size_t i = 0;
             i < q->durations.size();
             ++i)
        {
            const double frameEnd =
                frameStart + q->durations[i];

            if (targetTime < frameEnd ||
                i + 1 == q->durations.size())
            {
                return frameValues[i];
            }

            frameStart = frameEnd;
        }
    }

    if (frameValues.size() == 1)
    {
        return frameValues.front();
    }

    std::vector<double> midpoints;
    midpoints.reserve(frameValues.size());

    for (size_t i = 0;
         i < frameValues.size();
         ++i)
    {
        midpoints.push_back(
            q->timePoints[i] -
            0.5 * q->durations[i]);
    }

    size_t left = 0;
    size_t right = 1;

    if (targetTime <= midpoints.front())
    {
        left = 0;
        right = 1;
    }
    else if (targetTime >= midpoints.back())
    {
        left = frameValues.size() - 2;
        right = frameValues.size() - 1;
    }
    else
    {
        const auto upper =
            std::upper_bound(
                midpoints.begin(),
                midpoints.end(),
                targetTime);

        right = static_cast<size_t>(
            std::distance(
                midpoints.begin(),
                upper));
        left = right - 1;
    }

    const double t1 = midpoints[left];
    const double t2 = midpoints[right];
    const double y1 = frameValues[left];
    const double y2 = frameValues[right];

    const double value =
        y1 +
        (y2 - y1) *
        (targetTime - t1) /
        (t2 - t1);

    // Match the existing KMAP fine-sampling behavior.
    return std::max(0.0, value);
}

double
qSlicerDynamicPETModuleWidgetPrivate::
averagePlasmaTimesParentFractionOverInterval(
    const std::vector<double>& plasmaTimes,
    const std::vector<double>& plasmaValues,
    bool plasmaIsFrameCurve,
    const std::string& plasmaInterpolation,
    const std::vector<double>& parentTimes,
    const std::vector<double>& parentValues,
    double startTime,
    double endTime)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (endTime <= startTime ||
        parentTimes.size() < 2 ||
        parentTimes.size() != parentValues.size() ||
        startTime < parentTimes.front() ||
        endTime > parentTimes.back())
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    if (plasmaIsFrameCurve)
    {
        if (plasmaValues.size() != q->durations.size())
        {
            return
                std::numeric_limits<double>::quiet_NaN();
        }
    }
    else if (plasmaTimes.size() < 2 ||
             plasmaTimes.size() != plasmaValues.size() ||
             startTime < plasmaTimes.front() ||
             endTime > plasmaTimes.back())
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> knots;
    knots.push_back(startTime);

    for (double t : parentTimes)
    {
        if (t > startTime &&
            t < endTime)
        {
            knots.push_back(t);
        }
    }

    if (plasmaIsFrameCurve)
    {
        if (plasmaInterpolation == "const")
        {
            for (double t : q->timePoints)
            {
                if (t > startTime &&
                    t < endTime)
                {
                    knots.push_back(t);
                }
            }
        }
        else
        {
            for (size_t i = 0;
                 i < q->durations.size();
                 ++i)
            {
                const double midpoint =
                    q->timePoints[i] -
                    0.5 * q->durations[i];

                if (midpoint > startTime &&
                    midpoint < endTime)
                {
                    knots.push_back(midpoint);
                }
            }
        }
    }
    else
    {
        for (double t : plasmaTimes)
        {
            if (t > startTime &&
                t < endTime)
            {
                knots.push_back(t);
            }
        }
    }

    knots.push_back(endTime);

    std::sort(knots.begin(), knots.end());
    knots.erase(
        std::unique(
            knots.begin(),
            knots.end()),
        knots.end());

    auto plasmaAt =
        [&](double t)
        {
            if (plasmaIsFrameCurve)
            {
                return this->evaluateFrameCurve(
                    plasmaValues,
                    t,
                    plasmaInterpolation);
            }

            return this->interpolateInputFunction(
                plasmaTimes,
                plasmaValues,
                t,
                plasmaInterpolation);
        };

    auto parentAt =
        [&](double t)
        {
            return this->interpolateInputFunction(
                parentTimes,
                parentValues,
                t,
                "linear");
        };

    constexpr double inverseSqrt3 =
        0.57735026918962576451;

    double integral = 0.0;

    for (size_t i = 0;
         i + 1 < knots.size();
         ++i)
    {
        const double a = knots[i];
        const double b = knots[i + 1];
        const double midpoint = 0.5 * (a + b);
        const double halfWidth = 0.5 * (b - a);

        const double x1 =
            midpoint -
            inverseSqrt3 * halfWidth;

        const double x2 =
            midpoint +
            inverseSqrt3 * halfWidth;

        const double p1 = plasmaAt(x1);
        const double p2 = plasmaAt(x2);
        const double f1 = parentAt(x1);
        const double f2 = parentAt(x2);

        if (!std::isfinite(p1) ||
            !std::isfinite(p2) ||
            !std::isfinite(f1) ||
            !std::isfinite(f2))
        {
            return
                std::numeric_limits<double>::quiet_NaN();
        }

        integral +=
            halfWidth *
            (p1 * f1 + p2 * f2);
    }

    return integral / (endTime - startTime);
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
pbrAtTime(
    double timeSec,
    double& pbr,
    QString* errorMessage) const
{
    const double p1 =
        this->pbrp1Edit->text().toDouble();
    const double p2 =
        this->pbrp2Edit->text().toDouble();
    const double p3 =
        this->pbrp3Edit->text().toDouble();

    pbr =
        p1 *
        std::exp(-p2 * timeSec / 60.0) +
        p3;

    if (!std::isfinite(pbr) ||
        pbr <= 1e-12)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Plasma/whole-blood ratio is non-positive "
                    "at t=%1 s.")
                    .arg(timeSec);
        }
        return false;
    }

    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
buildCurrentInputFunction(
    InputFunctionResult& result,
    bool requireWholeBlood,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    result = InputFunctionResult{};

    if (this->inputFunctionCacheValid)
    {
        result = this->cachedInputFunction;
        if (requireWholeBlood && !result.hasWholeBlood)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "TCM requires total whole blood for its vascular term.");
            }
            return false;
        }
        return true;
    }

    if (q->durations.empty() ||
        q->timePoints.empty() ||
        q->durations.size() != q->timePoints.size())
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Select a valid dynamic PET dataset first.");
        }
        return false;
    }

    const int source =
        this->IFSourceSelector->currentIndex();

    if (source != 0 && source != 1)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "The selected input-function source "
                    "is not supported.");
        }
        return false;
    }

    const double petEndTimeSec =
        q->timePoints.back();

    const std::string inputInterpolation =
        this->selectedIFInterpolation();

    const IFCurveDomain sourceDomain =
        this->selectedIFCurveDomain();

    auto fillFrameAverages =
        [&](const std::vector<double>& times,
            const std::vector<double>& values,
            const std::string& interpolation,
            std::vector<double>& output) -> bool
        {
            output.clear();
            output.reserve(q->durations.size());

            double frameStart = 0.0;

            for (double duration : q->durations)
            {
                const double frameEnd =
                    frameStart + duration;

                const double value =
                    this->averageInputFunctionOverInterval(
                        times,
                        values,
                        frameStart,
                        frameEnd,
                        interpolation);

                if (!std::isfinite(value))
                {
                    return false;
                }

                output.push_back(value);
                frameStart = frameEnd;
            }

            return true;
        };

    // ------------------------------------------------------------
    // 1. Patient-specific source -> explicit whole blood / plasma.
    // ------------------------------------------------------------

    if (source == 0)
    {
        // Segment-derived IDIF is always total whole blood.
        // Keep the observation mask separate from the continuous IF. If an
        // IDIF point is removed in the TAC plot, reconstruct the IF through
        // that gap for integration, while retaining frameKeep=false so MTGA
        // can exclude that regression observation.
        std::vector<double> observedWholeBlood;
        if (!this->buildCurrentSegmentInputFunction(
                observedWholeBlood,
                errorMessage,
                &result.frameKeep))
        {
            return false;
        }

        if (observedWholeBlood.size() != q->durations.size() ||
            result.frameKeep.size() != observedWholeBlood.size())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "IDIF frame count does not match PET timing.");
            }
            return false;
        }

        result.frameWholeBlood = observedWholeBlood;

        const bool hasRemovedInputPoint =
            std::find(
                result.frameKeep.begin(),
                result.frameKeep.end(),
                false) != result.frameKeep.end();

        // A missing first IDIF observation would require unknown pre-injection
        // history. Keep that boundary mandatory. Missing observations at the
        // end, however, simply shorten the usable acquisition.
        if (result.frameKeep.empty() ||
            !result.frameKeep.front())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "The first IDIF observation must be retained. "
                        "Input history before the first available PET frame "
                        "is not extrapolated.");
            }
            return false;
        }

        auto lastKeepIt =
            std::find(
                result.frameKeep.rbegin(),
                result.frameKeep.rend(),
                true);

        if (lastKeepIt == result.frameKeep.rend())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "No retained IDIF observations are available.");
            }
            return false;
        }

        const size_t lastKeptIndex =
            result.frameKeep.size() - 1 -
            static_cast<size_t>(
                std::distance(
                    result.frameKeep.rbegin(),
                    lastKeepIt));

        result.supportFrameCount = lastKeptIndex + 1;

        const size_t retainedCount =
            static_cast<size_t>(
                std::count(
                    result.frameKeep.begin(),
                    result.frameKeep.begin() +
                        static_cast<std::ptrdiff_t>(
                            result.supportFrameCount),
                    true));

        if (retainedCount < 2)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "At least two IDIF observations must remain.");
            }
            return false;
        }

        if (hasRemovedInputPoint)
        {
            std::vector<double> retainedTimes;
            std::vector<double> retainedValues;
            retainedTimes.reserve(result.supportFrameCount);
            retainedValues.reserve(result.supportFrameCount);

            for (size_t i = 0;
                 i < result.supportFrameCount;
                 ++i)
            {
                if (!result.frameKeep[i])
                {
                    continue;
                }

                retainedTimes.push_back(
                    q->timePoints[i] -
                    0.5 * q->durations[i]);
                retainedValues.push_back(
                    observedWholeBlood[i]);
            }

            // Reconstruct only internal missing observations. Trailing removed
            // points are outside support and therefore are never used by a fit.
            for (size_t i = 0;
                 i < result.supportFrameCount;
                 ++i)
            {
                if (result.frameKeep[i])
                {
                    continue;
                }

                const double midpoint =
                    q->timePoints[i] -
                    0.5 * q->durations[i];

                result.frameWholeBlood[i] =
                    this->interpolateInputFunction(
                        retainedTimes,
                        retainedValues,
                        midpoint,
                        inputInterpolation);
            }
        }

        const int sourceProcessing =
            this->IFSourceProcessingSelector->currentIndex();

        if (sourceProcessing == 1)
        {
            std::vector<double> retainedTimes;
            std::vector<double> retainedValues;
            std::vector<double> evaluationTimes;
            retainedTimes.reserve(result.supportFrameCount);
            retainedValues.reserve(result.supportFrameCount);
            evaluationTimes.reserve(result.supportFrameCount);

            for (size_t i = 0; i < result.supportFrameCount; ++i)
            {
                const double midpoint =
                    q->timePoints[i] - 0.5 * q->durations[i];
                evaluationTimes.push_back(midpoint);
                if (result.frameKeep[i])
                {
                    retainedTimes.push_back(midpoint);
                    retainedValues.push_back(observedWholeBlood[i]);
                }
            }

            const std::vector<double> smoothed =
                lowessPredict(
                    retainedTimes,
                    retainedValues,
                    evaluationTimes,
                    this->IFLowessSpanSpinBox->value(),
                    2);

            if (smoothed.size() != result.supportFrameCount)
            {
                if (errorMessage)
                {
                    *errorMessage =
                        QObject::tr("LOWESS input-function smoothing failed.");
                }
                return false;
            }

            for (size_t i = 0; i < result.supportFrameCount; ++i)
            {
                result.frameWholeBlood[i] = smoothed[i];
                result.processedSourcePreviewTimesSec.push_back(q->timePoints[i]);
                result.processedSourcePreviewValues.push_back(smoothed[i]);
            }
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel =
                QObject::tr("LOWESS (span %1)")
                    .arg(this->IFLowessSpanSpinBox->value(), 0, 'f', 2);
        }
        else if (sourceProcessing == 2)
        {
            vtkSlicerDynamicPETLogic* logic =
                vtkSlicerDynamicPETLogic::SafeDownCast(q->logic());
            if (!logic)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr("DynamicPET logic is unavailable.");
                }
                return false;
            }

            std::vector<double> retainedTimes;
            std::vector<double> retainedValues;
            std::vector<double> retainedFrameStarts;
            std::vector<double> retainedFrameEnds;
            for (size_t i = 0; i < result.supportFrameCount; ++i)
            {
                if (!result.frameKeep[i])
                {
                    continue;
                }
                retainedTimes.push_back(
                    q->timePoints[i] - 0.5 * q->durations[i]);
                retainedValues.push_back(observedWholeBlood[i]);
                retainedFrameStarts.push_back(q->timePoints[i] - q->durations[i]);
                retainedFrameEnds.push_back(q->timePoints[i]);
            }

            std::vector<double> fittedObservations;
            std::string fengError;
            if (!logic->FitFengInputFunction(
                    retainedTimes,
                    retainedValues,
                    &retainedFrameStarts,
                    &retainedFrameEnds,
                    true,
                    result.fengParameters,
                    fittedObservations,
                    &fengError))
            {
                if (errorMessage)
                {
                    *errorMessage =
                        QObject::tr("Feng input-function fit failed: %1")
                            .arg(QString::fromStdString(fengError));
                }
                return false;
            }

            for (size_t i = 0; i < result.supportFrameCount; ++i)
            {
                result.frameWholeBlood[i] =
                    logic->AverageFengInputFunction(
                        q->timePoints[i] - q->durations[i],
                        q->timePoints[i],
                        result.fengParameters);
            }

            const double supportEndSec =
                q->timePoints[result.supportFrameCount - 1];
            // Dense representation of the analytic Feng curve is independent
            // of the user-selected TCM integration step. The TCM backend may
            // subsequently sample this continuous representation more finely.
            const double denseStepSec = 0.25;

            for (double t = 0.0; t < supportEndSec; t += denseStepSec)
            {
                result.nativeWholeBloodTimesSec.push_back(t);
                result.nativeWholeBloodValues.push_back(
                    logic->EvaluateFengInputFunction(t, result.fengParameters));
            }
            if (result.nativeWholeBloodTimesSec.empty() ||
                result.nativeWholeBloodTimesSec.back() < supportEndSec)
            {
                result.nativeWholeBloodTimesSec.push_back(supportEndSec);
                result.nativeWholeBloodValues.push_back(
                    logic->EvaluateFengInputFunction(
                        supportEndSec,
                        result.fengParameters));
            }

            result.processedSourcePreviewTimesSec =
                result.nativeWholeBloodTimesSec;
            result.processedSourcePreviewValues =
                result.nativeWholeBloodValues;
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel = QObject::tr("Feng model");
        }

        result.framePlasma.resize(
            result.frameWholeBlood.size());

        for (size_t i = 0;
             i < result.frameWholeBlood.size();
             ++i)
        {
            const double midpoint =
                q->timePoints[i] -
                0.5 * q->durations[i];

            double pbr = 0.0;
            if (!this->pbrAtTime(
                    midpoint,
                    pbr,
                    errorMessage))
            {
                return false;
            }

            result.framePlasma[i] =
                result.frameWholeBlood[i] * pbr;
        }

        if (!result.nativeWholeBloodTimesSec.empty() &&
            result.nativeWholeBloodTimesSec.size() ==
                result.nativeWholeBloodValues.size())
        {
            result.nativePlasmaTimesSec =
                result.nativeWholeBloodTimesSec;
            result.nativePlasmaValues.resize(
                result.nativeWholeBloodValues.size());
            for (size_t i = 0; i < result.nativeWholeBloodValues.size(); ++i)
            {
                double pbr = 0.0;
                if (!this->pbrAtTime(
                        result.nativeWholeBloodTimesSec[i],
                        pbr,
                        errorMessage))
                {
                    return false;
                }
                result.nativePlasmaValues[i] =
                    result.nativeWholeBloodValues[i] * pbr;
            }
        }

        result.hasWholeBlood = true;
        result.plasmaIsParent = false;

        // None/LOWESS retain the established frame-aware fine-sampling path.
        // Feng supplies its analytic curve densely so TCM can evaluate the
        // fitted source without reconstructing the bolus from frame means.
    }
    else
    {
        if (this->externalIFTimesSec.size() < 2 ||
            this->externalIFTimesSec.size() !=
                this->externalIFConcentrations.size())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "No valid external input-function CSV "
                        "has been loaded.");
            }
            return false;
        }

        if (this->externalIFTimesSec.back() + 1e-6 <
            petEndTimeSec)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Input-function coverage is insufficient.\n\n"
                        "PET end time: %1 s\n"
                        "Input-function end time: %2 s\n\n"
                        "Extrapolation is not performed.")
                        .arg(petEndTimeSec)
                        .arg(this->externalIFTimesSec.back());
            }
            return false;
        }

        // External activity curves are converted to the PET sequence native
        // quantitative scale before any biological-domain processing. This is
        // the key guard that prevents mixing image SUVbw with Bq/mL blood data.
        std::vector<double> externalPatientValues;
        if (!this->convertActivityVector(
                this->externalIFConcentrations,
                this->selectedExternalIFActivityUnit(),
                this->petStoredActivityUnit(),
                externalPatientValues,
                errorMessage))
        {
            return false;
        }

        std::vector<double> processedPatientTimes =
            this->externalIFTimesSec;
        std::vector<double> processedPatientValues =
            externalPatientValues;

        const int sourceProcessing =
            this->IFSourceProcessingSelector->currentIndex();

        if (sourceProcessing == 1)
        {
            processedPatientValues =
                lowessPredict(
                    processedPatientTimes,
                    processedPatientValues,
                    processedPatientTimes,
                    this->IFLowessSpanSpinBox->value(),
                    2);
            if (processedPatientValues.size() != processedPatientTimes.size())
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr("LOWESS input-function smoothing failed.");
                }
                return false;
            }
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel =
                QObject::tr("LOWESS (span %1)")
                    .arg(this->IFLowessSpanSpinBox->value(), 0, 'f', 2);
        }
        else if (sourceProcessing == 2)
        {
            if (sourceDomain == IFCurveDomain::ParentPlasma)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "Feng source modeling is intended for whole-blood or total-plasma input functions, not an already metabolite-corrected parent-plasma curve. Use None or LOWESS for this source.");
                }
                return false;
            }

            vtkSlicerDynamicPETLogic* logic =
                vtkSlicerDynamicPETLogic::SafeDownCast(q->logic());
            if (!logic)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr("DynamicPET logic is unavailable.");
                }
                return false;
            }

            std::vector<double> fittedObservations;
            std::string fengError;
            if (!logic->FitFengInputFunction(
                    processedPatientTimes,
                    processedPatientValues,
                    nullptr,
                    nullptr,
                    false,
                    result.fengParameters,
                    fittedObservations,
                    &fengError))
            {
                if (errorMessage)
                {
                    *errorMessage =
                        QObject::tr("Feng input-function fit failed: %1")
                            .arg(QString::fromStdString(fengError));
                }
                return false;
            }

            processedPatientTimes.clear();
            processedPatientValues.clear();
            // Dense representation of the analytic Feng curve is independent
            // of the user-selected TCM integration step. The TCM backend may
            // subsequently sample this continuous representation more finely.
            const double denseStepSec = 0.25;
            for (double t = 0.0; t < petEndTimeSec; t += denseStepSec)
            {
                processedPatientTimes.push_back(t);
                processedPatientValues.push_back(
                    logic->EvaluateFengInputFunction(t, result.fengParameters));
            }
            processedPatientTimes.push_back(petEndTimeSec);
            processedPatientValues.push_back(
                logic->EvaluateFengInputFunction(
                    petEndTimeSec,
                    result.fengParameters));
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel = QObject::tr("Feng model");
        }

        if (result.sourceProcessingApplied)
        {
            result.processedSourcePreviewTimesSec = processedPatientTimes;
            result.processedSourcePreviewValues = processedPatientValues;
        }

        if (sourceDomain == IFCurveDomain::WholeBlood)
        {
            result.nativeWholeBloodTimesSec =
                processedPatientTimes;
            result.nativeWholeBloodValues =
                processedPatientValues;

            result.nativePlasmaTimesSec =
                processedPatientTimes;
            result.nativePlasmaValues.resize(
                processedPatientValues.size());

            for (size_t i = 0;
                 i < processedPatientValues.size();
                 ++i)
            {
                double pbr = 0.0;
                if (!this->pbrAtTime(
                        processedPatientTimes[i],
                        pbr,
                        errorMessage))
                {
                    return false;
                }

                result.nativePlasmaValues[i] =
                    processedPatientValues[i] * pbr;
            }

            result.hasWholeBlood = true;
            result.plasmaIsParent = false;
        }
        else if (sourceDomain == IFCurveDomain::TotalPlasma)
        {
            result.nativePlasmaTimesSec =
                processedPatientTimes;
            result.nativePlasmaValues =
                processedPatientValues;

            result.nativeWholeBloodTimesSec =
                processedPatientTimes;
            result.nativeWholeBloodValues.resize(
                processedPatientValues.size());

            for (size_t i = 0;
                 i < processedPatientValues.size();
                 ++i)
            {
                double pbr = 0.0;
                if (!this->pbrAtTime(
                        processedPatientTimes[i],
                        pbr,
                        errorMessage))
                {
                    return false;
                }

                result.nativeWholeBloodValues[i] =
                    processedPatientValues[i] / pbr;
            }

            result.hasWholeBlood = true;
            result.plasmaIsParent = false;
        }
        else
        {
            // Already metabolite-corrected parent plasma.
            result.nativePlasmaTimesSec =
                processedPatientTimes;
            result.nativePlasmaValues =
                processedPatientValues;
            result.plasmaIsParent = true;

            const bool validWholeBlood =
                this->externalWholeBloodTimesSec.size() >= 2 &&
                this->externalWholeBloodTimesSec.size() ==
                    this->externalWholeBloodConcentrations.size() &&
                this->externalWholeBloodTimesSec.back() + 1e-6 >=
                    petEndTimeSec;

            if (validWholeBlood)
            {
                std::vector<double> companionValues;
                if (!this->convertActivityVector(
                        this->externalWholeBloodConcentrations,
                        this->selectedCompanionWholeBloodActivityUnit(),
                        this->petStoredActivityUnit(),
                        companionValues,
                        errorMessage))
                {
                    return false;
                }

                result.nativeWholeBloodTimesSec =
                    this->externalWholeBloodTimesSec;
                result.nativeWholeBloodValues =
                    std::move(companionValues);
                result.hasWholeBlood = true;
            }
            else if (requireWholeBlood)
            {
                if (errorMessage)
                {
                    *errorMessage =
                        QObject::tr(
                            "The external input already represents "
                            "parent plasma. TCM therefore requires a "
                            "separate total whole-blood CSV covering "
                            "the full PET acquisition.");
                }
                return false;
            }
        }

        if (!fillFrameAverages(
                result.nativePlasmaTimesSec,
                result.nativePlasmaValues,
                inputInterpolation,
                result.framePlasma))
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Could not calculate PET-frame plasma values.");
            }
            return false;
        }

        if (result.hasWholeBlood &&
            !fillFrameAverages(
                result.nativeWholeBloodTimesSec,
                result.nativeWholeBloodValues,
                inputInterpolation,
                result.frameWholeBlood))
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Could not calculate PET-frame whole-blood values.");
            }
            return false;
        }

        result.frameKeep.assign(
            q->durations.size(),
            true);
        result.supportFrameCount = q->durations.size();
    }

    // ------------------------------------------------------------
    // 2. Optional PBIF calibration.
    //    The patient curve calibrates the template; the template itself
    //    becomes the authoritative whole-blood/total-plasma representation.
    // ------------------------------------------------------------

    if (this->PBIFOptionCheckBox->isChecked())
    {
        if (result.plasmaIsParent)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "PBIF calibration is not available when the "
                        "patient input is already parent plasma.");
            }
            return false;
        }

        if (this->pbifTimesSec.size() < 2 ||
            this->pbifTimesSec.size() !=
                this->pbifTemplateValues.size())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Enable PBIF calibration only after loading a "
                        "valid time_s,template_value CSV.");
            }
            return false;
        }

        if (this->pbifTimesSec.back() + 1e-6 <
            petEndTimeSec)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "PBIF coverage is insufficient for the PET "
                        "acquisition. Extrapolation is not performed.");
            }
            return false;
        }

        const double calibrationStart =
            this->PBIFCalibrationStartSpinBox->value();
        const double calibrationEnd =
            this->PBIFCalibrationEndSpinBox->value();

        if (calibrationStart < 0.0 ||
            calibrationEnd <= calibrationStart ||
            calibrationEnd > petEndTimeSec ||
            calibrationStart < this->pbifTimesSec.front() ||
            calibrationEnd > this->pbifTimesSec.back())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Invalid PBIF AUC calibration interval. "
                        "It must lie inside both the patient IF and "
                        "PBIF temporal coverage.");
            }
            return false;
        }

        PBIFTemplateDomain templateChoice =
            this->selectedPBIFTemplateDomain();

        IFCurveDomain calibrationDomain =
            IFCurveDomain::WholeBlood;

        if (templateChoice ==
            PBIFTemplateDomain::WholeBlood)
        {
            calibrationDomain =
                IFCurveDomain::WholeBlood;
        }
        else if (templateChoice ==
                 PBIFTemplateDomain::TotalPlasma)
        {
            calibrationDomain =
                IFCurveDomain::TotalPlasma;
        }
        else
        {
            calibrationDomain =
                sourceDomain == IFCurveDomain::TotalPlasma
                ? IFCurveDomain::TotalPlasma
                : IFCurveDomain::WholeBlood;
        }

        double patientAUC =
            std::numeric_limits<double>::quiet_NaN();

        if (source == 0)
        {
            const std::vector<double>& patientFrames =
                calibrationDomain == IFCurveDomain::WholeBlood
                ? result.frameWholeBlood
                : result.framePlasma;

            patientAUC =
                this->integrateFrameAverageCurveOverInterval(
                    patientFrames,
                    calibrationStart,
                    calibrationEnd);

            result.pbifPatientCalibrationTimesSec.reserve(
                patientFrames.size());
            result.pbifPatientCalibrationValues =
                patientFrames;

            for (size_t i = 0;
                 i < patientFrames.size();
                 ++i)
            {
                result.pbifPatientCalibrationTimesSec.push_back(
                    q->timePoints[i] -
                    0.5 * q->durations[i]);
            }
        }
        else
        {
            const std::vector<double>& patientTimes =
                calibrationDomain == IFCurveDomain::WholeBlood
                ? result.nativeWholeBloodTimesSec
                : result.nativePlasmaTimesSec;

            const std::vector<double>& patientValues =
                calibrationDomain == IFCurveDomain::WholeBlood
                ? result.nativeWholeBloodValues
                : result.nativePlasmaValues;

            patientAUC =
                this->integrateInputFunctionOverInterval(
                    patientTimes,
                    patientValues,
                    calibrationStart,
                    calibrationEnd,
                    inputInterpolation);

            result.pbifPatientCalibrationTimesSec =
                patientTimes;
            result.pbifPatientCalibrationValues =
                patientValues;
        }

        const double pbifAUC =
            this->integrateInputFunctionOverInterval(
                this->pbifTimesSec,
                this->pbifTemplateValues,
                calibrationStart,
                calibrationEnd,
                "linear");

        if (!std::isfinite(patientAUC) ||
            !std::isfinite(pbifAUC) ||
            patientAUC <= 0.0 ||
            pbifAUC <= 0.0)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "PBIF AUC calibration failed because the "
                        "patient or template AUC is not positive and finite.");
            }
            return false;
        }

        result.pbifScale =
            patientAUC / pbifAUC;
        result.pbifApplied = true;
        result.pbifCalibrationDomain =
            calibrationDomain;

        result.pbifScaledValues.resize(
            this->pbifTemplateValues.size());

        for (size_t i = 0;
             i < this->pbifTemplateValues.size();
             ++i)
        {
            result.pbifScaledValues[i] =
                result.pbifScale *
                this->pbifTemplateValues[i];
        }

        result.nativePlasmaTimesSec =
            this->pbifTimesSec;
        result.nativeWholeBloodTimesSec =
            this->pbifTimesSec;
        result.nativePlasmaValues.resize(
            this->pbifTimesSec.size());
        result.nativeWholeBloodValues.resize(
            this->pbifTimesSec.size());

        if (calibrationDomain ==
            IFCurveDomain::WholeBlood)
        {
            result.nativeWholeBloodValues =
                result.pbifScaledValues;

            for (size_t i = 0;
                 i < this->pbifTimesSec.size();
                 ++i)
            {
                double pbr = 0.0;
                if (!this->pbrAtTime(
                        this->pbifTimesSec[i],
                        pbr,
                        errorMessage))
                {
                    return false;
                }

                result.nativePlasmaValues[i] =
                    result.pbifScaledValues[i] * pbr;
            }
        }
        else
        {
            result.nativePlasmaValues =
                result.pbifScaledValues;

            for (size_t i = 0;
                 i < this->pbifTimesSec.size();
                 ++i)
            {
                double pbr = 0.0;
                if (!this->pbrAtTime(
                        this->pbifTimesSec[i],
                        pbr,
                        errorMessage))
                {
                    return false;
                }

                result.nativeWholeBloodValues[i] =
                    result.pbifScaledValues[i] / pbr;
            }
        }

        result.hasWholeBlood = true;
        result.plasmaIsParent = false;
        result.supportFrameCount = q->durations.size();

        if (!fillFrameAverages(
                result.nativePlasmaTimesSec,
                result.nativePlasmaValues,
                "linear",
                result.framePlasma) ||
            !fillFrameAverages(
                result.nativeWholeBloodTimesSec,
                result.nativeWholeBloodValues,
                "linear",
                result.frameWholeBlood))
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Could not calculate PET-frame values from "
                        "the scaled PBIF.");
            }
            return false;
        }
    }

    // ------------------------------------------------------------
    // 3. Optional parent-fraction correction.
    //    Keep the fraction separate for TCM, and integrate Cp,total * fparent
    //    directly for MTGA frame averages.
    // ------------------------------------------------------------

    if (this->MetaboliteCorrectionCheckBox->isChecked())
    {
        if (result.plasmaIsParent)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "The selected external curve already represents "
                        "parent plasma; metabolite correction must be off.");
            }
            return false;
        }

        if (this->parentFractionTimesSec.size() < 2 ||
            this->parentFractionTimesSec.size() !=
                this->parentFractionValues.size())
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Load a valid time_s,parent_fraction CSV "
                        "before enabling parent-fraction correction.");
            }
            return false;
        }

        if (this->parentFractionTimesSec.back() + 1e-6 <
            petEndTimeSec)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Parent-fraction coverage is insufficient for "
                        "the PET acquisition. Extrapolation is not performed.");
            }
            return false;
        }

        result.applyParentFraction = true;
        result.parentFractionTimesSec =
            this->parentFractionTimesSec;
        result.parentFractionValues =
            this->parentFractionValues;

        result.frameModelPlasma.clear();
        result.frameModelPlasma.reserve(
            q->durations.size());

        const bool plasmaIsFrameCurve =
            result.nativePlasmaTimesSec.empty();

        const std::string plasmaInterpolation =
            result.pbifApplied
            ? "linear"
            : inputInterpolation;

        double frameStart = 0.0;

        for (double duration : q->durations)
        {
            const double frameEnd =
                frameStart + duration;

            const std::vector<double>& plasmaValues =
                plasmaIsFrameCurve
                ? result.framePlasma
                : result.nativePlasmaValues;

            const double value =
                this->averagePlasmaTimesParentFractionOverInterval(
                    result.nativePlasmaTimesSec,
                    plasmaValues,
                    plasmaIsFrameCurve,
                    plasmaInterpolation,
                    result.parentFractionTimesSec,
                    result.parentFractionValues,
                    frameStart,
                    frameEnd);

            if (!std::isfinite(value))
            {
                if (errorMessage)
                {
                    *errorMessage =
                        QObject::tr(
                            "Could not calculate parent-plasma "
                            "frame averages.");
                }
                return false;
            }

            result.frameModelPlasma.push_back(value);
            frameStart = frameEnd;
        }
    }
    else
    {
        result.frameModelPlasma =
            result.framePlasma;
    }

    if (result.supportFrameCount == 0)
    {
        result.supportFrameCount = q->durations.size();
    }

    if (result.frameModelPlasma.size() !=
        q->durations.size())
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Final plasma input does not match PET framing.");
        }
        return false;
    }

    // Cache the fully prepared IF. The result itself is independent of whether
    // the immediate caller requires whole blood; that requirement is checked
    // on retrieval as well. This avoids repeating an expensive Feng fit during
    // routine UI/status refreshes.
    this->cachedInputFunction = result;
    this->inputFunctionCacheValid = true;

    if (requireWholeBlood &&
        !result.hasWholeBlood)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "TCM requires total whole blood for its vascular term.");
        }
        return false;
    }

    return true;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
hasValidInputFunction(
    QString* errorMessage,
    bool requireWholeBlood)
{
    InputFunctionResult result;

    return
        this->buildCurrentInputFunction(
            result,
            requireWholeBlood,
            errorMessage);
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
buildCurrentInputFunctionWeights(
    std::vector<double>& weights,
    bool weighted,
    QString* errorMessage,
    bool excludeRemovedInputFrames)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    weights.clear();

    const bool processedInputWithoutUncertainty =
        this->PBIFOptionCheckBox->isChecked() ||
        this->MetaboliteCorrectionCheckBox->isChecked();

    if (processedInputWithoutUncertainty)
    {
        if (weighted)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Weighted voxelwise fitting is not available "
                        "when PBIF calibration or parent-fraction "
                        "correction is active because uncertainty "
                        "propagation is not modeled yet.");
            }
            return false;
        }

        weights.assign(
            q->durations.size(),
            1.0);

        if (excludeRemovedInputFrames &&
            this->IFSourceSelector->currentIndex() == 0 &&
            !q->IFID.empty())
        {
            const auto it = q->segmentTACs.find(q->IFID);
            if (it != q->segmentTACs.end() &&
                it->second.size() == weights.size())
            {
                for (size_t i = 0; i < weights.size(); ++i)
                {
                    if (!it->second[i].keep)
                    {
                        weights[i] = 0.0;
                    }
                }
            }
        }

        return true;
    }

    const int source =
        this->IFSourceSelector->
            currentIndex();

    if (source == 1)
    {
        if (weighted)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Weighted voxelwise fitting is not "
                        "available for the current two-column "
                        "CSV input-function format because the "
                        "file does not provide uncertainty "
                        "information.");
            }

            return false;
        }

        weights.assign(
            q->durations.size(),
            1.0);

        return true;
    }

    if (source != 0 ||
        q->IFID.empty())
    {
        return false;
    }

    const auto it =
        q->segmentTACs.find(
            q->IFID);

    if (it ==
        q->segmentTACs.end())
    {
        return false;
    }

    const std::string statistic =
        this->selectedIFStatistic();

    weights.reserve(
        it->second.size());

    for (const VoxelStatistics& vs :
         it->second)
    {
        double weight = 1.0;

        if (weighted)
        {
            const double sigma =
                statisticDispersionSigma(
                    vs,
                    statistic);
            weight =
                inverseVarianceWeightFromSigma(
                    sigma);
        }

        if (excludeRemovedInputFrames && !vs.keep)
        {
            weight = 0.0;
        }

        weights.push_back(weight);
    }

    if (weighted)
    {
        normalizePositiveWeights(weights);
    }

    return true;
}

void qSlicerDynamicPETModuleWidgetPrivate::populateVOI(std :: string ifID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->VOICheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->VOICheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Clear previous VOI checkboxes
  this->VOICheckContents->blockSignals(true);
  QLayoutItem* child;
  while ((child = this->VOICheckLayout->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }
    delete child;
  }

  const bool externalIFSource =
      this->IFSourceSelector->currentIndex() == 1;

  if (ifID.empty() &&
      !externalIFSource)
  {
      this->VOICheckContents->blockSignals(false);

      q->VOIsegmentIDs.clear();

      q->enableFITbutton();

      this->VOIsegmentSelectAll->
          setEnabled(false);

      return;
  }

  // Get the selected IF segment ID
  q->VOIsegmentIDs.clear();

  // Add checkboxes for all other segments
  for (const std::string& segmentID : this->segmentDisplayOrder)
  {
    if (!ifID.empty() &&
        segmentID == ifID)
    {
        continue;
    }
    auto it = q->segmentTACsnames.find(segmentID);
    if (it == q->segmentTACsnames.end())
      continue;

    const std::string& displayName = it->second;

    QCheckBox* cb = new QCheckBox(QString::fromStdString(displayName));
    cb->setProperty("SegmentID", QString::fromStdString(segmentID));
    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    cb->setChecked(wasSelected);
    this->VOICheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                     q, SLOT(onVOISegmentsChanged()));
    if (wasSelected)
      q->VOIsegmentIDs.push_back(segmentID);
  }
  q->enableFITbutton();

  this->VOICheckLayout->addStretch();
  this->VOIsegmentSelectAll->setEnabled(true);
  this->VOICheckContents->blockSignals(false);

}

void qSlicerDynamicPETModuleWidgetPrivate::populateVOIMTGA(std :: string ifID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->VOIMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->VOIMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->property("SegmentID").toString());
    }
  }

  // Clear previous VOI checkboxes
  this->VOIMTGACheckContents->blockSignals(true);
  QLayoutItem* child;
  while ((child = this->VOIMTGACheckLayout->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }
    delete child;
  }

  const bool externalIFSource =
      this->IFSourceSelector->currentIndex() == 1;

  if (ifID.empty() &&
      !externalIFSource)
  {
      this->VOIMTGACheckContents->
          blockSignals(false);

      q->VOIMTGAsegmentIDs.clear();

      q->enableFITMTGAbutton();

      this->VOIMTGAsegmentSelectAll->
          setEnabled(false);

      return;
  }

  // Get the selected IF segment ID
  q->VOIMTGAsegmentIDs.clear();

  // Add checkboxes for all other segments
  for (const std::string& segmentID : this->segmentDisplayOrder)
  {
    if (!ifID.empty() &&
        segmentID == ifID)
    {
        continue;
    }
    auto it = q->segmentTACsnames.find(segmentID);
    if (it == q->segmentTACsnames.end())
      continue;

    const std::string& displayName = it->second;
    QCheckBox* cb = new QCheckBox(QString::fromStdString(displayName));
    cb->setProperty("SegmentID", QString::fromStdString(segmentID));
    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(segmentID));
    cb->setChecked(wasSelected);
    this->VOIMTGACheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                     q, SLOT(onVOIMTGASegmentsChanged()));
    if (wasSelected)
      q->VOIMTGAsegmentIDs.push_back(segmentID);
  }
  q->enableFITMTGAbutton();

  this->VOIMTGACheckLayout->addStretch();
  this->VOIMTGAsegmentSelectAll->setEnabled(true);
  this->VOIMTGACheckContents->blockSignals(false);

}


void qSlicerDynamicPETModuleWidgetPrivate::populateResultsVOI()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->VOISelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->VOISelector->itemData(currentIndex).toString().toStdString();
  }

  this->VOISelector->blockSignals(true);  // Optional: prevent signal emission
  this->VOISelector->clear();
  this->VOISelector->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentTCM.empty() || q->segmentTACsnames.empty() || q->segmentTACs.empty()) {
    this->TCMResultsButton->setEnabled(false);
    this->VOISelector->blockSignals(false);
    this->populateResultsTable("");
    return;
  }
  this->TCMResultsButton->setEnabled(true);

  int restoredIndex = 0;
  for (const std::string& segmentID : this->segmentDisplayOrder)
  {
    if (q->segmentTCM.find(segmentID) == q->segmentTCM.end())
      continue;
    auto nameIt = q->segmentTACsnames.find(segmentID);
    if (nameIt == q->segmentTACsnames.end())
      continue;
    this->VOISelector->addItem(QString::fromStdString(nameIt->second), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->VOISelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->VOISelector->setCurrentIndex(restoredIndex);
  }

  this->VOISelector->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q->plotTCMVOI = passonID;
  this->populateResultsTable(passonID);
  return;
}

void qSlicerDynamicPETModuleWidgetPrivate::populateResultsVOIMTGA()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->VOISelectorMTGA->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->VOISelectorMTGA->itemData(currentIndex).toString().toStdString();
  }

  this->VOISelectorMTGA->blockSignals(true);  // Optional: prevent signal emission
  this->VOISelectorMTGA->clear();
  this->VOISelectorMTGA->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (q->segmentMTGA.empty() || q->segmentTACsnames.empty() || q->segmentTACs.empty()) {
    this->MTGAResultsButton->setEnabled(false);
    this->VOISelectorMTGA->blockSignals(false);
    this->populateResultsMTGATable("");
    return;
  }
  this->MTGAResultsButton->setEnabled(true);

  int restoredIndex = 0;
  for (const std::string& segmentID : this->segmentDisplayOrder)
  {
    if (q->segmentMTGA.find(segmentID) == q->segmentMTGA.end())
      continue;

    auto nameIt = q->segmentTACsnames.find(segmentID);
    if (nameIt == q->segmentTACsnames.end())
      continue;

    this->VOISelectorMTGA->addItem(QString::fromStdString(nameIt->second), QString::fromStdString(segmentID));
    if (segmentID==currentSelectedID) {
      restoredIndex = this->VOISelectorMTGA->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->VOISelectorMTGA->setCurrentIndex(restoredIndex);
  }

  this->VOISelectorMTGA->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q->plotMTGAVOI = passonID;
  this->populateResultsMTGATable(passonID);
  return;
}

auto makeNumericItem = [](double value, int precision = 6) {
    auto *item = new QTableWidgetItem(QString::number(value));
    item->setData(Qt::EditRole, value);  // ensures numeric sorting
    item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    return item;
};

void qSlicerDynamicPETModuleWidgetPrivate::populateResultsTable(std :: string segmentID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  this->TCMResultsTable->clear();
  this->TCMResultsTable->setRowCount(0);
  this->TCMResultsTable->setColumnCount(0);

  if (segmentID.empty())
  {
    this->populateModelsTCM(segmentID);
    this->populateModelComboTCM(this->TCMModel1, "", "", segmentID);
    this->populateModelComboTCM(this->TCMModel2, "", "", segmentID);
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
    return;
  }

  // Get the parameter map for the selected segment
  const auto& labelMap = q->segmentTCM[segmentID];

  const bool hasLiverDBIF =
      labelMap.find("Liver DBIF") !=
      labelMap.end();

  // Keep ka/fA out of the ordinary TCM table unless the liver model
  // is actually present for this VOI.
  QStringList headers =
      hasLiverDBIF
      ? QStringList{
            "", "K1", "k2", "k3", "k4", "ka", "fA",
            "vb", "td", "Ki", "DV", "AIC", "BIC", "MASE", "chi^2_nu", "Bounds"}
      : QStringList{
            "", "K1", "k2", "k3", "k4",
            "vb", "td", "Ki", "DV", "AIC", "BIC", "MASE", "chi^2_nu", "Bounds"};

  this->TCMResultsTable->setColumnCount(headers.size());
  this->TCMResultsTable->setHorizontalHeaderLabels(headers);
  int row = 0;
  this->TCMResultsTable->setRowCount(labelMap.size());

  int totalRowHeight = 0;
  this->TCMResultsTable->resizeRowsToContents();
  for (const auto& [label, params] : labelMap)
  {
    this->TCMResultsTable->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(label)));
    this->TCMResultsTable->setItem(row, 1, makeNumericItem(params.K1));
    this->TCMResultsTable->setItem(row, 2, makeNumericItem(params.k2));
    this->TCMResultsTable->setItem(row, 3, makeNumericItem(params.k3));
    this->TCMResultsTable->setItem(row, 4, makeNumericItem(params.k4));

    if (hasLiverDBIF)
    {
      this->TCMResultsTable->setItem(row, 5, makeNumericItem(params.ka));
      this->TCMResultsTable->setItem(row, 6, makeNumericItem(params.fa));
      this->TCMResultsTable->setItem(row, 7, makeNumericItem(params.vb));
      this->TCMResultsTable->setItem(row, 8, makeNumericItem(params.td));
      this->TCMResultsTable->setItem(row, 9, makeNumericItem(params.Ki));
      this->TCMResultsTable->setItem(row, 10, makeNumericItem(params.DV));
      this->TCMResultsTable->setItem(row, 11, makeNumericItem(params.AIC));
      this->TCMResultsTable->setItem(row, 12, makeNumericItem(params.BIC));
      this->TCMResultsTable->setItem(row, 13, makeNumericItem(params.MASE));
      this->TCMResultsTable->setItem(row, 14, makeNumericItem(params.chi2));
    }
    else
    {
      this->TCMResultsTable->setItem(row, 5, makeNumericItem(params.vb));
      this->TCMResultsTable->setItem(row, 6, makeNumericItem(params.td));
      this->TCMResultsTable->setItem(row, 7, makeNumericItem(params.Ki));
      this->TCMResultsTable->setItem(row, 8, makeNumericItem(params.DV));
      this->TCMResultsTable->setItem(row, 9, makeNumericItem(params.AIC));
      this->TCMResultsTable->setItem(row, 10, makeNumericItem(params.BIC));
      this->TCMResultsTable->setItem(row, 11, makeNumericItem(params.MASE));
      this->TCMResultsTable->setItem(row, 12, makeNumericItem(params.chi2));
    }

    const QString boundStatus =
        formatTCMBoundStatus(params.boundFlags);
    auto* boundItem =
        new QTableWidgetItem(
            boundStatus.isEmpty()
            ? QString()
            : QString("BOUND: ") + boundStatus);
    boundItem->setToolTip(
        boundStatus.isEmpty()
        ? QObject::tr("No fitted parameter is within 0.1% of an active bound.")
        : QObject::tr("Parameter at/near an active optimization bound. "
             "Consider whether the configured bound is appropriate before interpreting the estimate."));
    this->TCMResultsTable->setItem(
        row,
        headers.size() - 1,
        boundItem);

    totalRowHeight += this->TCMResultsTable->rowHeight(row);
    ++row;
  }
  totalRowHeight += this->TCMResultsTable->horizontalHeader()->height();
  totalRowHeight += 2 * this->TCMResultsTable->frameWidth();
  const int minimumUsableTCMTableHeight =
      this->TCMResultsTable->horizontalHeader()->height() +
      3 * this->TCMResultsTable->verticalHeader()->defaultSectionSize() +
      2 * this->TCMResultsTable->frameWidth();
  this->TCMResultsTable->setMinimumHeight(
      std::max(totalRowHeight, minimumUsableTCMTableHeight));
  // Never collapse a one-model result table to a header-sized strip. Let the
  // surrounding ROI layout allocate additional vertical space normally.
  this->TCMResultsTable->setMaximumHeight(QWIDGETSIZE_MAX);
  // Make table read-only
  this->TCMResultsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  this->TCMResultsTable->resizeColumnsToContents();
  this->populateModelsTCM(segmentID);

  std::string sel1, sel2;
  int idx1 = this->TCMModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = this->TCMModel1->itemData(idx1).toString().toStdString();
  int idx2 = this->TCMModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = this->TCMModel2->itemData(idx2).toString().toStdString();
  this->populateModelComboTCM(this->TCMModel1, sel2, sel1, segmentID);
  this->populateModelComboTCM(this->TCMModel2, sel1, sel2, segmentID);
  if (idx1 > 0 & idx2 >0){
    q->runTCMstat(sel1, sel2, segmentID);
  } else {
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
  }
}

void qSlicerDynamicPETModuleWidgetPrivate::populateResultsMTGATable(std :: string segmentID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  this->MTGAResultsTable->clear();
  this->MTGAResultsTable->setRowCount(0);
  this->MTGAResultsTable->setColumnCount(0);

  if (segmentID.empty())
  {
    this->populateModelsMTGA(segmentID);
    this->populateModelCombo(this->MTGAModel1, "", "", segmentID);
    this->populateModelCombo(this->MTGAModel2, "", "", segmentID);
    this->MTGAVuongP->setText("");
    return;
  }

  // Define column headers
  QStringList headers = { "", "Ki", "DV", "Intercept", "R2", "AIC", "MASE"};
  this->MTGAResultsTable->setColumnCount(headers.size());
  this->MTGAResultsTable->setHorizontalHeaderLabels(headers);

  // Get the parameter map for the selected segment
  const auto& labelMap = q->segmentMTGA[segmentID];
  int row = 0;
  this->MTGAResultsTable->setRowCount(labelMap.size());

  int totalRowHeight = 0;
  this->MTGAResultsTable->resizeRowsToContents();
  for (const auto& [label, params] : labelMap)
  {
    this->MTGAResultsTable->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(label)));
    this->MTGAResultsTable->setItem(row, 1, makeNumericItem(params.Ki));
    this->MTGAResultsTable->setItem(row, 2, makeNumericItem(params.DV));
    this->MTGAResultsTable->setItem(row, 3, makeNumericItem(params.Intercept));
    this->MTGAResultsTable->setItem(row, 4, makeNumericItem(params.R2));
    this->MTGAResultsTable->setItem(row, 5, makeNumericItem(params.AIC));
    this->MTGAResultsTable->setItem(row, 6, makeNumericItem(params.MASE));
    // this->MTGAResultsTable->setItem(row, 7, makeNumericItem(params.chi2));

    totalRowHeight += this->MTGAResultsTable->rowHeight(row);
    ++row;
  }
  totalRowHeight += this->MTGAResultsTable->horizontalHeader()->height();
  totalRowHeight += 4 * this->MTGAResultsTable->frameWidth();
  const int minimumUsableMTGATableHeight =
      this->MTGAResultsTable->horizontalHeader()->height() +
      3 * this->MTGAResultsTable->verticalHeader()->defaultSectionSize() +
      2 * this->MTGAResultsTable->frameWidth();
  this->MTGAResultsTable->setMinimumHeight(
      std::max(totalRowHeight, minimumUsableMTGATableHeight));
  this->MTGAResultsTable->setMaximumHeight(QWIDGETSIZE_MAX);
  // Make table read-only
  this->MTGAResultsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  this->MTGAResultsTable->resizeColumnsToContents();
  this->populateModelsMTGA(segmentID);

  std::string sel1, sel2;
  int idx1 = this->MTGAModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = this->MTGAModel1->itemData(idx1).toString().toStdString();
  int idx2 = this->MTGAModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = this->MTGAModel2->itemData(idx2).toString().toStdString();
  this->populateModelCombo(this->MTGAModel1, sel2, sel1, segmentID);
  this->populateModelCombo(this->MTGAModel2, sel1, sel2, segmentID);
  if (idx1 > 0 & idx2 >0) {
    q->runVuong(sel1, sel2, segmentID);
  } else {
    this->MTGAVuongP->setText("");
  }
}


void qSlicerDynamicPETModuleWidgetPrivate::populateModelsTCM(std :: string segmentID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);
  // Step 1: Save currently selected segment IDs
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < this->ModelsTCMCheckLayout->count(); ++i)
  {
    QLayoutItem* item = this->ModelsTCMCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->text());
    }
  }

  // Clear previous VOI checkboxes
  this->ModelsTCMCheckContents->blockSignals(true);
  QLayoutItem* child;
  while ((child = this->ModelsTCMCheckLayout->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }
    delete child;
  }


  if (segmentID.empty())
  {
    this->ModelsTCMCheckContents->blockSignals(false);
    this->plotTCMButton->setEnabled(false);
    this->ModelsTCMSelectAll->setEnabled(false);
    return;
  }

  // Get the parameter map for the selected segment
  const auto& labelMap = q->segmentTCM[segmentID];
  for (const auto& [label, params] : labelMap) {
    QCheckBox* cb = new QCheckBox(QString::fromStdString(label));
    bool wasSelected = previouslySelectedIDs.contains(QString::fromStdString(label));
    cb->setChecked(wasSelected);
    this->ModelsTCMCheckLayout->addWidget(cb);
  }
  this->plotTCMButton->setEnabled(true);

  this->ModelsTCMCheckLayout->addStretch();
  this->ModelsTCMSelectAll->setEnabled(true);
  this->ModelsTCMCheckContents->blockSignals(false);
  q->onPlotTCMbutton();
}

void qSlicerDynamicPETModuleWidgetPrivate::populateModelsMTGA(std :: string segmentID)
{

  Q_Q(qSlicerDynamicPETModuleWidget);

  std :: string currentSelectedID = "";
  int currentIndex = this->MTGASelector->currentIndex();
  if (currentIndex >= 0)
  {
    currentSelectedID = this->MTGASelector->itemData(currentIndex).toString().toStdString();
  }

  this->MTGASelector->blockSignals(true);  // Optional: prevent signal emission
  this->MTGASelector->clear();
  this->MTGASelector->addItem(QString::fromStdString("None"), QString::fromStdString(""));

  if (segmentID.empty())
  {
    this->MTGASelector->blockSignals(false);
    this->plotMTGAButton->setEnabled(false);
    return;
  }

  const auto& labelMap = q->segmentMTGA[segmentID];
  int restoredIndex = 0;
  for (const auto& [label, params] : labelMap)
  {
    this->MTGASelector->addItem(QString::fromStdString(label), QString::fromStdString(label));
    if (label==currentSelectedID) {
      restoredIndex = this->MTGASelector->count() - 1;
    }
  }

  // Restore previous selection if possible
  if (restoredIndex >= 0)
  {
    this->MTGASelector->setCurrentIndex(restoredIndex);
  }
  this->plotMTGAButton->setEnabled(true);

  this->MTGASelector->blockSignals(false);

  std :: string passonID = restoredIndex>0 ? currentSelectedID : "";
  q->plotMTGAModel = passonID;
  q->onPlotMTGAbutton();
}

void qSlicerDynamicPETModuleWidgetPrivate::populateModelCombo(
    QComboBox* comboToFill,
    const std::string& otherSelectedModel,
    const std::string& currentSelectedModel,
    const std::string& segmentID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (segmentID.empty()) {
    this->MTGAModel1->clear();
    this->MTGAModel2->clear();
    this->MTGAVuongP->setText("");
    this->MTGAStatTestButton->setCollapsed(true);
    this->MTGAStatTestButton->setEnabled(false);
    return;
  }

  auto it = q->segmentMTGA.find(segmentID);
  if (it == q->segmentMTGA.end()) {
    this->MTGAModel1->clear();
    this->MTGAModel2->clear();
    this->MTGAVuongP->setText("");
    this->MTGAStatTestButton->setCollapsed(true);
    this->MTGAStatTestButton->setEnabled(false);
    return;
  }
  this->MTGAStatTestButton->setEnabled(true);
  const auto& modelsForSegment = it->second;

  comboToFill->blockSignals(true);
  comboToFill->clear();
  comboToFill->addItem("", "");  // empty choice

  int restoredIndex = 0;
  for (const auto& [modelName, params] : modelsForSegment)
  {
    if (!otherSelectedModel.empty() && modelName == otherSelectedModel)
      continue;  // skip what’s selected in the other box

    comboToFill->addItem(QString::fromStdString(modelName), QString::fromStdString(modelName));

    if (modelName == currentSelectedModel)
    {
      restoredIndex = comboToFill->count() - 1;
    }
  }

  comboToFill->setCurrentIndex(restoredIndex);
  comboToFill->blockSignals(false);
}

void qSlicerDynamicPETModuleWidgetPrivate::populateModelComboTCM(
    QComboBox* comboToFill,
    const std::string& otherSelectedModel,
    const std::string& currentSelectedModel,
    const std::string& segmentID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (segmentID.empty()) {
    this->TCMModel1->clear();
    this->TCMModel2->clear();
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
    this->TCMStatTestButton->setCollapsed(true);
    this->TCMStatTestButton->setEnabled(false);
    return;
  }

  auto it = q->segmentTCM.find(segmentID);
  if (it == q->segmentTCM.end()) {
    this->TCMModel1->clear();
    this->TCMModel2->clear();
    this->TCMLRTP->setText("");
    this->TCMVuongP->setText("");
    this->TCMStatTestButton->setCollapsed(true);
    this->TCMStatTestButton->setEnabled(false);
    return;
  }
  const auto& modelsForSegment = it->second;

  int comparableModelCount = 0;

  for (const auto& [modelName, params] :
       modelsForSegment)
  {
    if (modelName != "Liver DBIF")
    {
      ++comparableModelCount;
    }
  }

  this->TCMStatTestButton->
      setEnabled(
          comparableModelCount >= 2);

  if (comparableModelCount < 2)
  {
    this->TCMStatTestButton->
        setCollapsed(true);
  }

  comboToFill->blockSignals(true);
  comboToFill->clear();
  comboToFill->addItem("", "");  // empty choice

  int restoredIndex = 0;
  for (const auto& [modelName, params] : modelsForSegment)
  {
    // Liver DBIF has a different parameterization and is intentionally
    // excluded from the generic nested/Vuong TCM comparison machinery.
    if (modelName == "Liver DBIF")
      continue;

    if (!otherSelectedModel.empty() && modelName == otherSelectedModel)
      continue;  // skip what’s selected in the other box

    comboToFill->addItem(QString::fromStdString(modelName), QString::fromStdString(modelName));

    if (modelName == currentSelectedModel)
    {
      restoredIndex = comboToFill->count() - 1;
    }
  }

  comboToFill->setCurrentIndex(restoredIndex);
  comboToFill->blockSignals(false);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
setPostTACEnabled(bool enabled)
{
    const int roiIndex =
        this->PlotsTabWidget->
            indexOf(
                this->ROIModelingWidget);

    if (roiIndex < 0)
    {
        return;
    }

    bool inputReady = false;

    if (enabled)
    {
        QString ignoredError;
        inputReady = this->hasValidInputFunction(
            &ignoredError,
            false);
    }

    this->PlotsTabWidget->setTabEnabled(
        roiIndex,
        enabled && inputReady);

    // Parametric Imaging is intentionally NOT
    // controlled by TAC extraction.
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
ensureParametricPETData(
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (q->petID ==
        vtkMRMLSubjectHierarchyNode::
            INVALID_ITEM_ID)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "No dynamic PET is selected.");
        }

        return false;
    }

    if (q->durations.empty() ||
        q->timePoints.empty() ||
        q->numberOfTimepoints <= 0)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Dynamic PET timing information "
                    "is not available.");
        }

        return false;
    }

    vtkSlicerDynamicPETLogic* logic =
        vtkSlicerDynamicPETLogic::
            SafeDownCast(
                q->logic());

    if (!logic)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "DynamicPET logic is unavailable.");
        }

        return false;
    }

    if (q->PET_flatten_values.empty())
    {
        q->stopRequested = false;

        q->ProgressBar->setVisible(true);
        q->ProgressBar->setValue(0);

        logic->Image2Flatten(
            q->petID,
            q->PET_flatten_values,
            q->PETdims,
            q->numberOfTimepoints,
            q->ProgressBar,
            q->stopButton,
            q->stopRequested);

        q->ProgressBar->setValue(0);
        q->ProgressBar->setVisible(false);

        if (q->stopRequested)
        {
            if (errorMessage)
            {
                *errorMessage =
                    QObject::tr(
                        "Dynamic PET preparation was canceled.");
            }

            return false;
        }
    }

    const vtkIdType expected =
        static_cast<vtkIdType>(
            q->PETdims[0]) *
        static_cast<vtkIdType>(
            q->PETdims[1]) *
        static_cast<vtkIdType>(
            q->PETdims[2]);

    if (expected <= 0 ||
        static_cast<vtkIdType>(
            q->PET_flatten_values.size())
            != expected)
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Dynamic PET voxel dimensions "
                    "are inconsistent.");
        }

        return false;
    }

    if (q->PET_flatten_values.empty() ||
        q->PET_flatten_values.front().size()
            != static_cast<size_t>(
                q->numberOfTimepoints))
    {
        if (errorMessage)
        {
            *errorMessage =
                QObject::tr(
                    "Dynamic PET frame count "
                    "is inconsistent.");
        }

        return false;
    }

    return true;
}

void qSlicerDynamicPETModuleWidgetPrivate::updateMTGAOutputUI()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  const bool saveDICOM =
      this->MTGASaveDICOMCheckBoxImg->isChecked();

  this->MTGADICOMDirectoryLabelImg->setEnabled(saveDICOM);
  this->MTGADICOMDirectoryImg->setEnabled(saveDICOM);

  q->enableFITMTGAImgbutton();

  this->updateMTGAOptimizationUI();
}


void qSlicerDynamicPETModuleWidgetPrivate::updateTCMOutputUI()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  const bool saveDICOM =
      this->TCMSaveDICOMCheckBoxImg->isChecked();

  this->TCMDICOMDirectoryLabelImg->setEnabled(saveDICOM);
  this->TCMDICOMDirectoryImg->setEnabled(saveDICOM);

  q->enableFITTCMImgbutton();

  this->updateTCMOptimizationUI();
}

void qSlicerDynamicPETModuleWidgetPrivate::
updateMTGAOptimizationUI()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  auto hasResult =
      [&](const std::string& modelID)
      {
        auto it =
            q->MTGAImgOutcomes.find(modelID);

        return
            it != q->MTGAImgOutcomes.end() &&
            !it->second.empty();
      };

  const bool hasPatlak =
      hasResult("Patlak");

  const bool hasLogan =
      hasResult("Logan");

  const bool hasRE =
      hasResult("RE");

  const bool available =
      hasPatlak &&
      (hasLogan || hasRE);

  this->MTGAOptimizationCollapsibleButtonImg
      ->setEnabled(available);

  if (!available)
  {
    this->MTGAOptimizationCollapsibleButtonImg
        ->setCollapsed(true);

    this->GenerateMTGAOptimizedImgButton
        ->setEnabled(false);

    this->RefreshMTGARGBButtonImg
        ->setEnabled(false);

    return;
  }

  // Preserve reversible-model selection if possible.
  const QString previousModel =
      this->MTGAReversibleModelComboImg
          ->currentText();

  this->MTGAReversibleModelComboImg
      ->blockSignals(true);

  this->MTGAReversibleModelComboImg
      ->clear();

  if (hasLogan)
  {
    this->MTGAReversibleModelComboImg
        ->addItem(
            "Logan",
            "Logan");
  }

  if (hasRE)
  {
    this->MTGAReversibleModelComboImg
        ->addItem(
            "RE",
            "RE");
  }

  int restoredIndex =
      this->MTGAReversibleModelComboImg
          ->findText(previousModel);

  if (restoredIndex < 0 &&
      this->MTGAReversibleModelComboImg
          ->count() > 0)
  {
    restoredIndex = 0;
  }

  if (restoredIndex >= 0)
  {
    this->MTGAReversibleModelComboImg
        ->setCurrentIndex(restoredIndex);
  }

  this->MTGAReversibleModelComboImg
      ->blockSignals(false);

  const bool useVuong =
      this->MTGAUseVuongCheckBoxImg
          ->isChecked();

  this->MTGAVuongAlphaSpinBoxImg
      ->setEnabled(useVuong);

  const bool show =
      this->MTGAShowInSlicerCheckBoxImg
          ->isChecked();

  const bool save =
      this->MTGASaveDICOMCheckBoxImg
          ->isChecked();

  const bool saveReady =
      !save ||
      !this->MTGADICOMDirectoryImg
           ->currentPath()
           .trimmed()
           .isEmpty();

  this->GenerateMTGAOptimizedImgButton
      ->setEnabled(
          available &&
          (show || save) &&
          saveReady);

  vtkMRMLScene* scene =
      q->mrmlScene();

  const bool rgbReady =
      scene &&
      !this->MTGAOptimizedKiNodeID.empty() &&
      !this->MTGAOptimizedDVNodeID.empty() &&
      scene->GetNodeByID(
          this->MTGAOptimizedKiNodeID.c_str()) &&
      scene->GetNodeByID(
          this->MTGAOptimizedDVNodeID.c_str());

  this->RefreshMTGARGBButtonImg
      ->setEnabled(rgbReady);
}

vtkMRMLScalarVolumeNode*
qSlicerDynamicPETModuleWidgetPrivate::
createMTGAOptimizedScalarVolume(
    const std::vector<double>& values,
    const QString& name,
    vtkMRMLScalarVolumeNode* refPETNode,
    vtkMRMLSubjectHierarchyNode* shNode,
    vtkIdType refPetID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!scene ||
      !refPETNode ||
      !shNode)
  {
    return nullptr;
  }

  const vtkIdType expectedSize =
      static_cast<vtkIdType>(
          q->PETdims[0]) *
      static_cast<vtkIdType>(
          q->PETdims[1]) *
      static_cast<vtkIdType>(
          q->PETdims[2]);

  if (values.size() !=
      static_cast<size_t>(expectedSize))
  {
    return nullptr;
  }

  vtkNew<vtkImageData> image;

  image->SetDimensions(
      q->PETdims[0],
      q->PETdims[1],
      q->PETdims[2]);

  image->AllocateScalars(
      VTK_DOUBLE,
      1);

  double* destination =
      static_cast<double*>(
          image->GetScalarPointer());

  if (!destination)
  {
    return nullptr;
  }

  std::copy(
      values.begin(),
      values.end(),
      destination);

  vtkMRMLScalarVolumeNode* node =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          scene->AddNewNodeByClass(
              "vtkMRMLScalarVolumeNode",
              name.toUtf8().constData()));

  if (!node)
  {
    return nullptr;
  }

  node->SetAndObserveImageData(
      image.GetPointer());

  node->CopyOrientation(
      refPETNode);

  node->SetSpacing(
      refPETNode->GetSpacing());

  node->SetOrigin(
      refPETNode->GetOrigin());

  node->CreateDefaultDisplayNodes();

  vtkMRMLScalarVolumeDisplayNode*
      displayNode =
          vtkMRMLScalarVolumeDisplayNode::
              SafeDownCast(
                  node->GetDisplayNode());

  if (displayNode)
  {
      displayNode->AutoWindowLevelOn();
  }

  const vtkIdType parentItemID =
      shNode->GetItemParent(
          refPetID);

  const vtkIdType newItemID =
      shNode->GetItemByDataNode(
          node);

  shNode->SetItemParent(
      newItemID,
      parentItemID);

  return node;
}

void qSlicerDynamicPETModuleWidgetPrivate::
refreshMTGAOptimizedRGB()
{
  std::cout
      << "[MTGA RGB] refresh START"
      << std::endl;

  Q_Q(qSlicerDynamicPETModuleWidget);

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!scene)
  {
    return;
  }

  vtkMRMLScalarVolumeNode* kiNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          scene->GetNodeByID(
              this->MTGAOptimizedKiNodeID.c_str()));

  vtkMRMLScalarVolumeNode* dvNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          scene->GetNodeByID(
              this->MTGAOptimizedDVNodeID.c_str()));

  if (!kiNode || !dvNode)
  {
    this->RefreshMTGARGBButtonImg
        ->setEnabled(false);

    return;
  }

  std::cout
      << "[MTGA RGB] Ki/DV nodes retrieved."
      << std::endl;

  vtkMRMLScalarVolumeDisplayNode*
      kiDisplay =
          vtkMRMLScalarVolumeDisplayNode::
              SafeDownCast(
                  kiNode->GetDisplayNode());

  vtkMRMLScalarVolumeDisplayNode*
      dvDisplay =
          vtkMRMLScalarVolumeDisplayNode::
              SafeDownCast(
                  dvNode->GetDisplayNode());

  if (!kiDisplay ||
      !dvDisplay)
  {
    return;
  }

  std::cout
      << "[MTGA RGB] Ki/DV display nodes retrieved."
      << std::endl;

  const double kiWindow =
      std::max(
          kiDisplay->GetWindow(),
          1.0e-12);

  const double kiLevel =
      kiDisplay->GetLevel();

  const double kiLow =
      kiLevel -
      0.5 * kiWindow;

  const double kiHigh =
      kiLevel +
      0.5 * kiWindow;


  const double dvWindow =
      std::max(
          dvDisplay->GetWindow(),
          1.0e-12);

  const double dvLevel =
      dvDisplay->GetLevel();

  const double dvLow =
      dvLevel -
      0.5 * dvWindow;

  const double dvHigh =
      dvLevel +
      0.5 * dvWindow;


  const vtkIdType numberOfVoxels =
      static_cast<vtkIdType>(
          q->PETdims[0]) *
      static_cast<vtkIdType>(
          q->PETdims[1]) *
      static_cast<vtkIdType>(
          q->PETdims[2]);

  if (this->MTGAOptimizedSelection.size() !=
          static_cast<size_t>(numberOfVoxels) ||
      this->MTGAOptimizedKiValues.size() !=
          static_cast<size_t>(numberOfVoxels) ||
      this->MTGAOptimizedDVValues.size() !=
          static_cast<size_t>(numberOfVoxels))
  {
    return;
  }

  std::cout
      << "[MTGA RGB] Ki display: "
      << kiLow << " -> " << kiHigh
      << std::endl;

  std::cout
      << "[MTGA RGB] DV display: "
      << dvLow << " -> " << dvHigh
      << std::endl;


  std::cout
      << "[MTGA RGB] Allocating RGB vtkImageData..."
      << std::endl;
  vtkNew<vtkImageData> rgbImage;

  rgbImage->SetDimensions(
      q->PETdims[0],
      q->PETdims[1],
      q->PETdims[2]);

  rgbImage->AllocateScalars(
      VTK_UNSIGNED_CHAR,
      3);

  unsigned char* rgb =
      static_cast<unsigned char*>(
          rgbImage->GetScalarPointer());

  if (!rgb)
  {
    return;
  }

  std::cout
      << "[MTGA RGB] RGB buffer allocated."
      << std::endl;
  auto normalize =
      [](double value,
         double low,
         double high)
      {
        if (!std::isfinite(value))
        {
          return 0.0;
        }

        const double width =
            high - low;

        if (width <= 0.0)
        {
          return 0.0;
        }

        return std::clamp(
            (value - low) / width,
            0.0,
            1.0);
      };


  const unsigned char patlakClass =
      static_cast<unsigned char>(
          MTGAOptimizedClass::Patlak);

  const unsigned char reversibleClass =
      static_cast<unsigned char>(
          MTGAOptimizedClass::Reversible);

  std::cout
      << "[MTGA RGB] Filling RGB buffer, voxels="
      << numberOfVoxels
      << std::endl;

  for (vtkIdType i = 0;
       i < numberOfVoxels;
       ++i)
  {
    unsigned char red = 0;
    unsigned char green = 0;
    unsigned char blue = 0;

    const unsigned char selected =
        this->MTGAOptimizedSelection[
            static_cast<size_t>(i)];

    if (selected == patlakClass)
    {
      const double normalized =
          normalize(
              this->MTGAOptimizedKiValues[
                  static_cast<size_t>(i)],
              kiLow,
              kiHigh);

      red =
          static_cast<unsigned char>(
              std::lround(
                  255.0 *
                  normalized));
    }
    else if (selected ==
             reversibleClass)
    {
      const double normalized =
          normalize(
              this->MTGAOptimizedDVValues[
                  static_cast<size_t>(i)],
              dvLow,
              dvHigh);

      blue =
          static_cast<unsigned char>(
              std::lround(
                  255.0 *
                  normalized));
    }

    rgb[3 * i + 0] = red;
    rgb[3 * i + 1] = green;
    rgb[3 * i + 2] = blue;
  }

  std::cout
      << "[MTGA RGB] RGB buffer filled."
      << std::endl;

  std::cout
      << "[MTGA RGB] Looking for existing RGB node..."
      << std::endl;

  vtkMRMLVolumeNode* rgbNode =
      nullptr;

  if (!this->MTGAOptimizedRGBNodeID.empty())
  {
    rgbNode =
        vtkMRMLVolumeNode::SafeDownCast(
            scene->GetNodeByID(
                this->MTGAOptimizedRGBNodeID.c_str()));
  }

  if (!rgbNode)
  {

    std::cout
        << "[MTGA RGB] Creating vtkMRMLVectorVolumeNode "
           "through MRML factory..."
        << std::endl;
    QString rgbName =
        QString::fromUtf8(
            kiNode->GetName());

    rgbName.replace(
        "MTGA Optimized Ki",
        "MTGA Selection RGB");

    vtkMRMLNode* createdNode =
        scene->AddNewNodeByClass(
            "vtkMRMLVectorVolumeNode",
            rgbName.toUtf8().constData());

    std::cout
        << "[MTGA RGB] AddNewNodeByClass returned: "
        << (createdNode
            ? createdNode->GetClassName()
            : "NULL")
        << std::endl;

    rgbNode =
        vtkMRMLVolumeNode::SafeDownCast(
            createdNode);

    std::cout
        << "[MTGA RGB] vtkMRMLVolumeNode cast: "
        << (rgbNode ? "OK" : "FAILED")
        << std::endl;

    if (!rgbNode)
    {
      if (createdNode)
      {
        scene->RemoveNode(createdNode);
      }

      qWarning()
          << "Could not create "
             "vtkMRMLVectorVolumeNode for "
             "MTGA RGB visualization.";

      return;
    }

    this->MTGAOptimizedRGBNodeID =
        rgbNode->GetID();

    rgbNode->CopyOrientation(
        kiNode);

    rgbNode->SetSpacing(
        kiNode->GetSpacing());

    rgbNode->SetOrigin(
        kiNode->GetOrigin());

    std::cout
        << "[MTGA RGB] Setting VoxelVectorTypeColorRGB..."
        << std::endl;

    rgbNode->SetVoxelVectorType(
        vtkMRMLVolumeNode::
            VoxelVectorTypeColorRGB);

    std::cout
        << "[MTGA RGB] Voxel vector type set."
        << std::endl;

    rgbNode->CreateDefaultDisplayNodes();

    vtkMRMLSubjectHierarchyNode* shNode =
        vtkMRMLSubjectHierarchyNode::
            GetSubjectHierarchyNode(scene);

    if (shNode)
    {
      const vtkIdType kiItemID =
          shNode->GetItemByDataNode(
              kiNode);

      const vtkIdType parentItemID =
          shNode->GetItemParent(
              kiItemID);

      const vtkIdType rgbItemID =
          shNode->GetItemByDataNode(
              rgbNode);

      shNode->SetItemParent(
          rgbItemID,
          parentItemID);
    }
  }

  std::cout
      << "[MTGA RGB] Assigning RGB vtkImageData..."
      << std::endl;

  rgbNode->SetAndObserveImageData(
      rgbImage.GetPointer());

  std::cout
      << "[MTGA RGB] RGB vtkImageData assigned."
      << std::endl;

  rgbNode->Modified();

  this->RefreshMTGARGBButtonImg
      ->setEnabled(true);

  std::cout
      << "[MTGA RGB] refresh END"
      << std::endl;
}

void qSlicerDynamicPETModuleWidgetPrivate::
generateMTGAOptimizedResult()
{
  std::cout
      << "[MTGA OPT] generateMTGAOptimizedResult START"
      << std::endl;
  Q_Q(qSlicerDynamicPETModuleWidget);

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!scene)
  {
    return;
  }

  vtkSlicerDynamicPETLogic* logic =
      vtkSlicerDynamicPETLogic::SafeDownCast(
          q->logic());

  if (!logic)
  {
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode =
      vtkMRMLSubjectHierarchyNode::
          GetSubjectHierarchyNode(scene);

  if (!shNode)
  {
    return;
  }

  vtkMRMLScalarVolumeNode* refPETNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          shNode->GetItemDataNode(
              q->petID));

  if (!refPETNode)
  {
    return;
  }


  // ------------------------------------------------------------------------
  // Required fitted models
  // ------------------------------------------------------------------------

  auto patlakIt =
      q->MTGAImgOutcomes.find(
          "Patlak");

  if (patlakIt ==
          q->MTGAImgOutcomes.end() ||
      patlakIt->second.empty())
  {
    return;
  }


  QString reversibleQString =
      this->MTGAReversibleModelComboImg
          ->currentData()
          .toString();

  if (reversibleQString.isEmpty())
  {
    reversibleQString =
        this->MTGAReversibleModelComboImg
            ->currentText();
  }

  const std::string reversibleModel =
      reversibleQString.toStdString();


  auto reversibleIt =
      q->MTGAImgOutcomes.find(
          reversibleModel);

  if (reversibleIt ==
          q->MTGAImgOutcomes.end() ||
      reversibleIt->second.empty())
  {
    return;
  }


  const auto& patlak =
      patlakIt->second;

  const auto& reversible =
      reversibleIt->second;


  if (patlak.size() !=
          reversible.size() ||
      patlak.size() !=
          q->PET_flatten_values.size())
  {
    QMessageBox::warning(
        q,
        QObject::tr("MTGA model selection"),
        QObject::tr(
            "The fitted MTGA result sizes do not match "
            "the PET voxel count."));

    return;
  }

  std::cout
      << "[MTGA OPT] Cached models retrieved. "
      << "Patlak voxels=" << patlak.size()
      << ", " << reversibleModel
      << " voxels=" << reversible.size()
      << std::endl;

  // ------------------------------------------------------------------------
  // Selection settings
  // ------------------------------------------------------------------------

  const QString criterion =
      this->MTGASelectionCriterionComboImg
          ->currentText();

  const bool useVuong =
      this->MTGAUseVuongCheckBoxImg
          ->isChecked();

  const double alpha =
      this->MTGAVuongAlphaSpinBoxImg
          ->value();

  QString optimizationSuffix =
      reversibleQString +
      " - " +
      criterion;

  if (useVuong)
  {
    optimizationSuffix +=
        QString(
            " - Vuong p<%1")
            .arg(
                alpha,
                0,
                'g',
                3);
  }


  // ------------------------------------------------------------------------
  // Allocate optimized outputs
  // ------------------------------------------------------------------------

  const size_t numberOfVoxels =
      patlak.size();

  this->MTGAOptimizedSelection.assign(
      numberOfVoxels,
      static_cast<unsigned char>(
          MTGAOptimizedClass::Excluded));

  this->MTGAOptimizedKiValues.assign(
      numberOfVoxels,
      0.0);

  this->MTGAOptimizedDVValues.assign(
      numberOfVoxels,
      0.0);


  auto metricValue =
      [&](const MTGAParameters& parameters)
      {
        if (criterion == "R2")
        {
          return parameters.R2;
        }

        if (criterion == "AIC")
        {
          return parameters.AIC;
        }

        return parameters.MASE;
      };


  size_t patlakSelectedCount = 0;
  size_t reversibleSelectedCount = 0;
  size_t nonSignificantCount = 0;
  size_t invalidCount = 0;


  // ------------------------------------------------------------------------
  // ONE voxelwise selection pass.
  //
  // Metric chooses candidate winner.
  // Optional Vuong test only decides whether that choice is accepted.
  // ------------------------------------------------------------------------

  std::cout
      << "[MTGA OPT] Starting voxelwise selection. "
      << "Criterion=" << criterion.toStdString()
      << ", reversible=" << reversibleModel
      << ", Vuong=" << (useVuong ? "ON" : "OFF")
      << ", alpha=" << alpha
      << std::endl;

  for (const int voxelIndex :
       this->parametricFitVoxelIndices)
  {
    if (voxelIndex < 0 ||
        static_cast<size_t>(
            voxelIndex) >= numberOfVoxels)
    {
      continue;
    }

    const size_t v =
        static_cast<size_t>(
            voxelIndex);

    const MTGAParameters& p =
        patlak[v];

    const MTGAParameters& r =
        reversible[v];


    // Empty fitted arrays mean that a valid fit was
    // not produced for this voxel.
    if (p.y.empty() ||
        r.y.empty())
    {
      ++invalidCount;
      continue;
    }


    const double pMetric =
        metricValue(p);

    const double rMetric =
        metricValue(r);

    if (!std::isfinite(pMetric) ||
        !std::isfinite(rMetric))
    {
      ++invalidCount;
      continue;
    }


    bool patlakWins = false;

    if (criterion == "R2")
    {
      if (pMetric == rMetric)
      {
        ++invalidCount;
        continue;
      }

      patlakWins =
          pMetric > rMetric;
    }
    else
    {
      // AIC / MASE: lower is better.
      if (pMetric == rMetric)
      {
        ++invalidCount;
        continue;
      }

      patlakWins =
          pMetric < rMetric;
    }


    // ------------------------------------------------------
    // Optional significance gate.
    //
    // This is the ONLY place where Vuong is calculated.
    // Refreshing RGB later never repeats this.
    // ------------------------------------------------------

    if (useVuong)
    {
      if (p.r.empty() ||
          r.r.empty() ||
          p.r.size() != r.r.size() ||
          p.weights.size() !=
              p.r.size() ||
          r.weights.size() !=
              r.r.size())
      {
        ++invalidCount;
        continue;
      }

      std::vector<double> averageWeights(
          p.weights.size(),
          1.0);

      for (size_t i = 0;
           i < averageWeights.size();
           ++i)
      {
        averageWeights[i] =
            0.5 *
            (
              p.weights[i] +
              r.weights[i]
            );
      }

      const double pValue =
          logic->computeVuongP(
              p.r,
              r.r,
              &averageWeights,
              p.dof,
              r.dof,
              VuongCorrection::BIC,
              Tail::TwoSided);

      if (!std::isfinite(pValue) ||
          pValue >= alpha)
      {
        ++nonSignificantCount;
        continue;
      }
    }


    // ------------------------------------------------------
    // Accept candidate model.
    // ------------------------------------------------------

    if (patlakWins)
    {
      if (!std::isfinite(p.Ki))
      {
        ++invalidCount;
        continue;
      }

      this->MTGAOptimizedSelection[v] =
          static_cast<unsigned char>(
              MTGAOptimizedClass::Patlak);

      this->MTGAOptimizedKiValues[v] =
          p.Ki;

      ++patlakSelectedCount;
    }
    else
    {
      if (!std::isfinite(r.DV))
      {
        ++invalidCount;
        continue;
      }

      this->MTGAOptimizedSelection[v] =
          static_cast<unsigned char>(
              MTGAOptimizedClass::Reversible);

      this->MTGAOptimizedDVValues[v] =
          r.DV;

      ++reversibleSelectedCount;
    }
  }


  const size_t selectedCount =
      patlakSelectedCount +
      reversibleSelectedCount;

  if (selectedCount == 0)
  {
    QMessageBox::warning(
        q,
        QObject::tr("MTGA model selection"),
        QObject::tr(
            "No voxel survived MTGA model selection."));

    return;
  }

  std::cout
      << "[MTGA OPT] Selection pass COMPLETE. "
      << "Patlak=" << patlakSelectedCount
      << ", reversible=" << reversibleSelectedCount
      << ", non-significant=" << nonSignificantCount
      << ", invalid=" << invalidCount
      << std::endl;

  qDebug()
      << "MTGA optimized selection:"
      << "criterion =" << criterion
      << "| reversible ="
      << reversibleQString
      << "| Vuong =" << useVuong
      << "| Patlak =" << patlakSelectedCount
      << "| reversible ="
      << reversibleSelectedCount
      << "| non-significant ="
      << nonSignificantCount
      << "| invalid ="
      << invalidCount;


  // ------------------------------------------------------------------------
  // Replace previous optimized visualization nodes.
  // ------------------------------------------------------------------------

  std::cout
      << "[MTGA OPT] Removing previous optimized scene nodes..."
      << std::endl;

  this->removeMTGAOptimizedSceneNodes();

  std::cout
      << "[MTGA OPT] Previous nodes removed."
      << std::endl;

  // ------------------------------------------------------------------------
  // Create scalar quantitative maps in Slicer.
  // ------------------------------------------------------------------------

  if (this->MTGAShowInSlicerCheckBoxImg
          ->isChecked())
  {
    std::cout
        << "[MTGA OPT] Creating optimized Ki scalar volume..."
        << std::endl;
    vtkMRMLScalarVolumeNode* kiNode =
        this->createMTGAOptimizedScalarVolume(
            this->MTGAOptimizedKiValues,
            "MTGA Optimized Ki - " + optimizationSuffix,
            refPETNode,
            shNode,
            q->petID);

    std::cout
        << "[MTGA OPT] Ki node created: "
        << (kiNode ? kiNode->GetID() : "NULL")
        << std::endl;

    std::cout
        << "[MTGA OPT] Creating optimized DV scalar volume..."
        << std::endl;

    vtkMRMLScalarVolumeNode* dvNode =
        this->createMTGAOptimizedScalarVolume(
            this->MTGAOptimizedDVValues,
            "MTGA Optimized DV - " + optimizationSuffix,
            refPETNode,
            shNode,
            q->petID);

    std::cout
        << "[MTGA OPT] DV node created: "
        << (dvNode ? dvNode->GetID() : "NULL")
        << std::endl;

    if (kiNode)
    {
      this->MTGAOptimizedKiNodeID =
          kiNode->GetID();
    }

    if (dvNode)
    {
      this->MTGAOptimizedDVNodeID =
          dvNode->GetID();
    }


    auto setMetadata =
        [&](vtkMRMLNode* node)
        {
          if (!node)
          {
            return;
          }

          node->SetAttribute(
              "SlicerDynamicPET.MTGA.ReversibleModel",
              reversibleQString
                  .toUtf8()
                  .constData());

          node->SetAttribute(
              "SlicerDynamicPET.MTGA.SelectionCriterion",
              criterion
                  .toUtf8()
                  .constData());

          node->SetAttribute(
              "SlicerDynamicPET.MTGA.UseVuong",
              useVuong
                  ? "1"
                  : "0");

          if (useVuong)
          {
            node->SetAttribute(
                "SlicerDynamicPET.MTGA.VuongAlpha",
                QString::number(
                    alpha,
                    'g',
                    8)
                    .toUtf8()
                    .constData());
          }
        };


    setMetadata(kiNode);
    setMetadata(dvNode);



    if (kiNode && dvNode)
    {
      std::cout
          << "[MTGA OPT] Calling refreshMTGAOptimizedRGB..."
          << std::endl;

      this->refreshMTGAOptimizedRGB();

      std::cout
          << "[MTGA OPT] refreshMTGAOptimizedRGB returned."
          << std::endl;

      vtkMRMLNode* rgbNode =
          scene->GetNodeByID(
              this->MTGAOptimizedRGBNodeID
                  .c_str());

      setMetadata(rgbNode);
    }
  }



  // ------------------------------------------------------------------------
  // DICOM PM export:
  // only the QUANTITATIVE Ki and DV maps.
  // RGB is intentionally not exported.
  // ------------------------------------------------------------------------

  if (this->MTGASaveDICOMCheckBoxImg
          ->isChecked())
  {
    const QString outputDirectory =
        this->MTGADICOMDirectoryImg
            ->currentPath()
            .trimmed();

    const double framingNorm =
        this->framingNormEditImg
            ->text()
            .toDouble();


    QString kiUnitCode;
    QString kiUnitMeaning;

    bool kiUnitValid = true;

    if (std::abs(
            framingNorm - 60.0) <
        1.0e-9)
    {
      kiUnitCode = "/min";
      kiUnitMeaning =
          "per minute";
    }
    else if (std::abs(
                 framingNorm - 1.0) <
             1.0e-9)
    {
      kiUnitCode = "/s";
      kiUnitMeaning =
          "per second";
    }
    else
    {
      kiUnitValid = false;
    }


    if (kiUnitValid)
    {
      q->ProgressBar->setMinimum(0);
      q->ProgressBar->setMaximum(0);

      q->ProgressBar->setFormat(
          "Saving DICOM PMAP: "
          "MTGA Optimized - Ki...");

      q->ProgressBar->setVisible(true);

      QApplication::processEvents();

      this->exportParametricMapDICOM(
          refPETNode,
          this->MTGAOptimizedKiValues,
          "MTGA",
          "MTGAOptimized",
          "Ki",
          outputDirectory,
          7160,
          kiUnitCode,
          kiUnitMeaning);
    }
    else
    {
      QMessageBox::warning(
          q,
          QObject::tr("DICOM PMAP export"),
          QObject::tr(
              "MTGA Optimized - Ki was not exported "
              "because Framing Norm is %1 s.\n\n"
              "Its physical unit cannot be represented "
              "honestly as seconds or minutes without "
              "rescaling the numerical values.")
              .arg(framingNorm));
    }


    q->ProgressBar->setMinimum(0);
    q->ProgressBar->setMaximum(0);

    q->ProgressBar->setFormat(
        "Saving DICOM PMAP: "
        "MTGA Optimized - DV...");

    q->ProgressBar->setVisible(true);

    QApplication::processEvents();

    this->exportParametricMapDICOM(
        refPETNode,
        this->MTGAOptimizedDVValues,
        "MTGA",
        "MTGAOptimized",
        "DV",
        outputDirectory,
        7161,
        "1",
        "1");


    q->ProgressBar->hide();

    q->ProgressBar->setMinimum(0);
    q->ProgressBar->setMaximum(100);
    q->ProgressBar->setValue(0);
    q->ProgressBar->setFormat("%p%");
  }


  this->updateMTGAOptimizationUI();

  std::cout
      << "[MTGA OPT] generateMTGAOptimizedResult END"
      << std::endl;
}

void qSlicerDynamicPETModuleWidgetPrivate::
removeMTGAOptimizedSceneNodes()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  std::cout
      << "[MTGA OPT CLEANUP] START"
      << std::endl;

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!scene)
  {
    std::cout
        << "[MTGA OPT CLEANUP] No MRML scene."
        << std::endl;

    return;
  }

  auto removeNode =
      [&](std::string& nodeID)
      {
        if (nodeID.empty())
        {
          return;
        }

        std::cout
            << "[MTGA OPT CLEANUP] Looking for node: "
            << nodeID
            << std::endl;

        vtkMRMLNode* node =
            scene->GetNodeByID(
                nodeID.c_str());

        if (node)
        {
          std::cout
              << "[MTGA OPT CLEANUP] Removing node: "
              << node->GetName()
              << " ["
              << node->GetID()
              << "]"
              << std::endl;

          scene->RemoveNode(node);
        }
        else
        {
          std::cout
              << "[MTGA OPT CLEANUP] Node no longer exists."
              << std::endl;
        }

        nodeID.clear();
      };

  // Remove RGB first because it depends visually
  // on the Ki/DV optimized volumes.
  removeNode(
      this->MTGAOptimizedRGBNodeID);

  removeNode(
      this->MTGAOptimizedKiNodeID);

  removeNode(
      this->MTGAOptimizedDVNodeID);

  this->RefreshMTGARGBButtonImg
      ->setEnabled(false);

  std::cout
      << "[MTGA OPT CLEANUP] END"
      << std::endl;
}

void qSlicerDynamicPETModuleWidgetPrivate::populateTCMOptimizationModels()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  // Preserve previous states.
  QMap<QString, bool> previousStates;

  for (int i = 0;
       i < this->TCMOptimizationModelsCheckLayoutImg->count();
       ++i)
  {
    QLayoutItem* item =
        this->TCMOptimizationModelsCheckLayoutImg->itemAt(i);

    QCheckBox* cb =
        qobject_cast<QCheckBox*>(item->widget());

    if (cb)
    {
      previousStates[cb->text()] = cb->isChecked();
    }
  }

  // Clear old checkboxes.
  QLayoutItem* child = nullptr;

  while ((child =
      this->TCMOptimizationModelsCheckLayoutImg
          ->takeAt(0)) != nullptr)
  {
    if (child->widget())
    {
      delete child->widget();
    }

    delete child;
  }

  // Keep your established model ordering rather than
  // std::map alphabetical ordering.
  for (const QString& modelName : q->ModelsNamesTCM)
  {
    const std::string modelID =
        modelName.toStdString();

    auto it = q->TCMImgOutcomes.find(modelID);

    if (it == q->TCMImgOutcomes.end() ||
        it->second.empty())
    {
      continue;
    }

    QCheckBox* cb =
        new QCheckBox(
            modelName,
            this->TCMOptimizationModelsCheckContentsImg);

    // Previously existing model -> restore state.
    // Newly fitted model -> selected by default.
    if (previousStates.contains(modelName))
    {
      cb->setChecked(previousStates.value(modelName));
    }
    else
    {
      cb->setChecked(true);
    }

    this->TCMOptimizationModelsCheckLayoutImg
        ->addWidget(cb);

    QObject::connect(
        cb,
        &QCheckBox::toggled,
        q,
        [this](bool)
        {
          this->updateTCMOptimizationUI();
        });
  }

  this->TCMOptimizationModelsCheckLayoutImg
      ->addStretch();

  this->updateTCMOptimizationUI();
}

void qSlicerDynamicPETModuleWidgetPrivate::updateTCMOptimizationUI()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  int fittedModelCount = 0;
  int selectedModelCount = 0;

  for (int i = 0;
       i < this->TCMOptimizationModelsCheckLayoutImg->count();
       ++i)
  {
    QLayoutItem* item =
        this->TCMOptimizationModelsCheckLayoutImg->itemAt(i);

    QCheckBox* cb =
        qobject_cast<QCheckBox*>(item->widget());

    if (!cb)
    {
      continue;
    }

    ++fittedModelCount;

    if (cb->isChecked())
    {
      ++selectedModelCount;
    }
  }

  // Section becomes available when at least one voxelwise
  // TCM result exists.
  this->TCMOptimizationCollapsibleButtonImg
      ->setEnabled(fittedModelCount > 0);

  this->TCMOptimizationModelsSelectAllImg
      ->setEnabled(fittedModelCount > 0);

  if (fittedModelCount == 0)
  {
    this->TCMOptimizationCollapsibleButtonImg
        ->setCollapsed(true);

    this->GenerateTCMOptimizedImgButton
        ->setEnabled(false);

    return;
  }

  const bool useTests =
      this->TCMUseStatTestsCheckBoxImg->isChecked();

  this->TCMStatAlphaSpinBoxImg
      ->setEnabled(useTests);

  const bool anyParameter =
      this->TCMOptK1CheckBoxImg->isChecked() ||
      this->TCMOptk2CheckBoxImg->isChecked() ||
      this->TCMOptk3CheckBoxImg->isChecked() ||
      this->TCMOptk4CheckBoxImg->isChecked() ||
      this->TCMOptvbCheckBoxImg->isChecked() ||
      this->TCMOpttdCheckBoxImg->isChecked() ||
      this->TCMOptKiCheckBoxImg->isChecked() ||
      this->TCMOptDVCheckBoxImg->isChecked();

  const bool show =
      this->TCMShowInSlicerCheckBoxImg->isChecked();

  const bool save =
      this->TCMSaveDICOMCheckBoxImg->isChecked();

  const bool saveReady =
      !save ||
      !this->TCMDICOMDirectoryImg
           ->currentPath()
           .trimmed()
           .isEmpty();

  // Selection itself requires >=2 models.
  this->GenerateTCMOptimizedImgButton
      ->setEnabled(
          selectedModelCount >= 2 &&
          anyParameter &&
          (show || save) &&
          saveReady);
}

void qSlicerDynamicPETModuleWidgetPrivate::
removeTCMOptimizedSceneNodes()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!scene)
  {
    return;
  }

  for (std::string& nodeID :
       this->TCMOptimizedNodeIDs)
  {
    if (nodeID.empty())
    {
      continue;
    }

    vtkMRMLNode* node =
        scene->GetNodeByID(
            nodeID.c_str());

    if (node)
    {
      scene->RemoveNode(node);
    }
  }

  this->TCMOptimizedNodeIDs.clear();

  if (!this->TCMOptimizedModelSelectionNodeID.empty())
  {
    vtkMRMLNode* node =
        scene->GetNodeByID(
            this->
                TCMOptimizedModelSelectionNodeID
                .c_str());

    if (node)
    {
      scene->RemoveNode(node);
    }

    this->TCMOptimizedModelSelectionNodeID.clear();
  }
}

void qSlicerDynamicPETModuleWidgetPrivate::
generateTCMOptimizedResult()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  std::cout
      << "[TCM OPT] START"
      << std::endl;

  vtkSlicerDynamicPETLogic* logic =
      vtkSlicerDynamicPETLogic::SafeDownCast(
          q->logic());

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!logic || !scene)
  {
    return;
  }

  vtkMRMLSubjectHierarchyNode* shNode =
      vtkMRMLSubjectHierarchyNode::
          GetSubjectHierarchyNode(scene);

  if (!shNode)
  {
    return;
  }

  vtkMRMLScalarVolumeNode* refPETNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          shNode->GetItemDataNode(
              q->petID));

  if (!refPETNode)
  {
    return;
  }

  const size_t numberOfVoxels =
      q->PET_flatten_values.size();

  if (numberOfVoxels == 0)
  {
    return;
  }

  // ------------------------------------------------------------------------
  // Selected models
  // ------------------------------------------------------------------------

  std::vector<std::string> selectedModels;

  for (int i = 0;
       i <
       this->TCMOptimizationModelsCheckLayoutImg
           ->count();
       ++i)
  {
    QLayoutItem* item =
        this->TCMOptimizationModelsCheckLayoutImg
            ->itemAt(i);

    QCheckBox* cb =
        qobject_cast<QCheckBox*>(
            item->widget());

    if (cb && cb->isChecked())
    {
      selectedModels.push_back(
          cb->text().toStdString());
    }
  }

  if (selectedModels.size() < 2)
  {
    return;
  }

  // Every selected cached result must correspond to the
  // complete PET voxel grid.
  for (const std::string& modelID :
       selectedModels)
  {
    auto it =
        q->TCMImgOutcomes.find(modelID);

    if (it == q->TCMImgOutcomes.end() ||
        it->second.size() != numberOfVoxels)
    {
      QMessageBox::warning(
          q,
          QObject::tr("TCM model selection"),
          QObject::tr(
              "Cached voxelwise result for %1 "
              "is missing or has an invalid size.")
              .arg(
                  QString::fromStdString(
                      modelID)));

      return;
    }
  }

  // ------------------------------------------------------------------------
  // Requested output parameters
  // ------------------------------------------------------------------------

  std::vector<std::string> outputFields;

  if (this->TCMOptK1CheckBoxImg->isChecked())
    outputFields.push_back("K1");

  if (this->TCMOptk2CheckBoxImg->isChecked())
    outputFields.push_back("k2");

  if (this->TCMOptk3CheckBoxImg->isChecked())
    outputFields.push_back("k3");

  if (this->TCMOptk4CheckBoxImg->isChecked())
    outputFields.push_back("k4");

  if (this->TCMOptvbCheckBoxImg->isChecked())
    outputFields.push_back("vb");

  if (this->TCMOpttdCheckBoxImg->isChecked())
    outputFields.push_back("td");

  if (this->TCMOptKiCheckBoxImg->isChecked())
    outputFields.push_back("Ki");

  if (this->TCMOptDVCheckBoxImg->isChecked())
    outputFields.push_back("DV");

  if (outputFields.empty())
  {
    return;
  }

  const QString criterion =
      this->TCMSelectionCriterionComboImg
          ->currentText();

  const bool useTests =
      this->TCMUseStatTestsCheckBoxImg
          ->isChecked();

  const double alpha =
      this->TCMStatAlphaSpinBoxImg
          ->value();

  // ------------------------------------------------------------------------
  // Helpers
  // ------------------------------------------------------------------------

  auto criterionValue =
      [&](const TCMParameters& p)
      {
        if (criterion == "AIC")
          return p.AIC;

        if (criterion == "BIC")
          return p.BIC;

        if (criterion == "MASE")
          return p.MASE;

        return p.chi2;
      };

  auto modelCode =
      [](const std::string& modelID)
          -> unsigned char
      {
        if (modelID == "1TCM")   return 1;
        if (modelID == "1TdCM")  return 2;
        if (modelID == "1TiCM")  return 3;
        if (modelID == "1TidCM") return 4;
        if (modelID == "2TCM")   return 5;
        if (modelID == "2TdCM")  return 6;
        if (modelID == "2TiCM")  return 7;
        if (modelID == "2TidCM") return 8;

        return 0;
      };

  auto modelHasField =
      [](const std::string& modelID,
         const std::string& field)
      {
        if (field == "K1" ||
            field == "vb" ||
            field == "Ki")
        {
          return true;
        }

        if (field == "k2")
        {
          return
              modelID == "1TCM" ||
              modelID == "1TdCM" ||
              modelID == "2TCM" ||
              modelID == "2TdCM" ||
              modelID == "2TiCM" ||
              modelID == "2TidCM";
        }

        if (field == "k3")
        {
          return
              modelID == "2TCM" ||
              modelID == "2TdCM" ||
              modelID == "2TiCM" ||
              modelID == "2TidCM";
        }

        if (field == "k4")
        {
          return
              modelID == "2TCM" ||
              modelID == "2TdCM";
        }

        if (field == "td")
        {
          return
              modelID == "1TdCM" ||
              modelID == "1TidCM" ||
              modelID == "2TdCM" ||
              modelID == "2TidCM";
        }

        // DV is meaningful only for reversible models.
        if (field == "DV")
        {
          return
              modelID == "1TCM" ||
              modelID == "1TdCM" ||
              modelID == "2TCM" ||
              modelID == "2TdCM";
        }

        return false;
      };

  auto parameterValue =
      [](const TCMParameters& p,
         const std::string& field)
      {
        if (field == "K1") return p.K1;
        if (field == "k2") return p.k2;
        if (field == "k3") return p.k3;
        if (field == "k4") return p.k4;
        if (field == "vb") return p.vb;
        if (field == "td") return p.td;
        if (field == "Ki") return p.Ki;
        if (field == "DV") return p.DV;

        return
            std::numeric_limits<double>::
                quiet_NaN();
      };

  // ------------------------------------------------------------------------
  // Allocate outputs
  // ------------------------------------------------------------------------

  std::vector<unsigned char> selectedModelMap(
      numberOfVoxels,
      static_cast<unsigned char>(0));

  std::map<
      std::string,
      std::vector<double>>
      optimizedValues;

  for (const std::string& field :
       outputFields)
  {
    optimizedValues[field].assign(
        numberOfVoxels,
        0.0);
  }

  std::map<std::string, size_t>
      selectedCounts;

  size_t invalidVoxelCount = 0;
  size_t simplifiedVoxelCount = 0;
  size_t skippedStatComparisonCount = 0;

  size_t lrtComparisonCount = 0;
  size_t vuongComparisonCount = 0;

  // ------------------------------------------------------------------------
  // Progress
  // ------------------------------------------------------------------------

  q->ProgressBar->setMinimum(0);
  q->ProgressBar->setMaximum(100);
  q->ProgressBar->setValue(0);
  q->ProgressBar->setFormat(
      "Selecting TCM model (%p%)");
  q->ProgressBar->setVisible(true);

  const size_t totalEligible =
      this->parametricFitVoxelIndices.size();

  const size_t updateInterval =
      std::max(
          static_cast<size_t>(1),
          totalEligible / 100);

  size_t processed = 0;

  auto updateProgress =
      [&]()
      {
        ++processed;

        if (processed % updateInterval == 0 ||
            processed == totalEligible)
        {
          const int value =
              totalEligible > 0
              ? static_cast<int>(
                    100 * processed /
                    totalEligible)
              : 100;

          q->ProgressBar->setValue(value);

          QApplication::processEvents();
        }
      };

  // ------------------------------------------------------------------------
  // Single voxelwise selection pass
  // ------------------------------------------------------------------------

  for (const int voxelIndex :
       this->parametricFitVoxelIndices)
  {
    if (voxelIndex < 0 ||
        static_cast<size_t>(voxelIndex) >=
            numberOfVoxels)
    {
      updateProgress();
      continue;
    }

    const size_t v =
        static_cast<size_t>(
            voxelIndex);

    // ------------------------------------------------------
    // 1. Criterion-selected model.
    // All four TCM criteria are minimized.
    // ------------------------------------------------------

    const TCMParameters* bestParams =
        nullptr;

    std::string bestModel;

    double bestMetric =
        std::numeric_limits<double>::
            infinity();

    for (const std::string& modelID :
         selectedModels)
    {
      const TCMParameters& p =
          q->TCMImgOutcomes
              .at(modelID)[v];

      const double metric =
          criterionValue(p);

      // This model is invalid for this voxel.
      if (!std::isfinite(metric))
      {
        continue;
      }

      if (!bestParams ||
          metric < bestMetric ||
          (metric == bestMetric &&
           p.dof < bestParams->dof))
      {
        bestParams = &p;
        bestModel = modelID;
        bestMetric = metric;
      }
    }

    // All selected models were invalid.
    if (!bestParams)
    {
      ++invalidVoxelCount;
      updateProgress();
      continue;
    }

    // ------------------------------------------------------
    // 2. Optional statistical simplification.
    //
    // IMPORTANT:
    // all alternatives are compared with the ORIGINAL
    // criterion-selected winner, never with a progressively
    // replaced winner.
    // ------------------------------------------------------

    std::string finalModel =
        bestModel;

    const TCMParameters* finalParams =
        bestParams;

    int finalDof =
        bestParams->dof;

    double finalMetric =
        bestMetric;

    if (useTests)
    {
      for (const std::string& alternativeModel :
           selectedModels)
      {
        if (alternativeModel == bestModel)
        {
          continue;
        }

        const TCMParameters& alternative =
            q->TCMImgOutcomes
                .at(alternativeModel)[v];

        const double alternativeMetric =
            criterionValue(alternative);

        if (!std::isfinite(
                alternativeMetric))
        {
          continue;
        }

        // We only use statistical testing to simplify.
        if (alternative.dof >=
            bestParams->dof)
        {
          continue;
        }

        // compareModels() may need either likelihoods
        // (LRT) or residuals/weights (Vuong).
        // Valid fitted TCM results normally contain all.
        if (!std::isfinite(
                bestParams->loglik) ||
            !std::isfinite(
                alternative.loglik) ||
            bestParams->r.empty() ||
            alternative.r.empty() ||
            bestParams->r.size() !=
                alternative.r.size() ||
            bestParams->weights.size() !=
                alternative.weights.size() ||
            bestParams->weights.size() !=
                bestParams->r.size())
        {
          ++skippedStatComparisonCount;
          continue;
        }

        ModelComparisonResult comparison;

        try
        {
          comparison =
              logic->compareModels(
                  bestModel,
                  alternativeModel,
                  *bestParams,
                  alternative);
          if (comparison.type == "LRT")
          {
            ++lrtComparisonCount;
          }
          else if (comparison.type == "Vuong")
          {
            ++vuongComparisonCount;
          }
        }
        catch (const std::exception& e)
        {
          qWarning()
              << "TCM model comparison failed:"
              << QString::fromStdString(
                     bestModel)
              << "vs"
              << QString::fromStdString(
                     alternativeModel)
              << ":"
              << e.what();

          ++skippedStatComparisonCount;
          continue;
        }

        if (!std::isfinite(
                comparison.p_value))
        {
          ++skippedStatComparisonCount;
          continue;
        }

        // No significant evidence of a difference:
        // the simpler model becomes a candidate.
        if (comparison.p_value >= alpha)
        {
          if (alternative.dof < finalDof ||
              (alternative.dof == finalDof &&
               alternativeMetric <
                   finalMetric))
          {
            finalModel =
                alternativeModel;

            finalParams =
                &alternative;

            finalDof =
                alternative.dof;

            finalMetric =
                alternativeMetric;
          }
        }
      }
    }

    if (finalModel != bestModel)
    {
      ++simplifiedVoxelCount;
    }

    const unsigned char code =
        modelCode(finalModel);

    if (code == 0)
    {
      ++invalidVoxelCount;
      updateProgress();
      continue;
    }

    selectedModelMap[v] =
        code;

    ++selectedCounts[finalModel];

    // ------------------------------------------------------
    // 3. Copy the parameters from the selected model.
    //
    // A parameter absent from that model stays zero.
    // ------------------------------------------------------

    for (const std::string& field :
         outputFields)
    {
      if (!modelHasField(
              finalModel,
              field))
      {
        continue;
      }

      const double value =
          parameterValue(
              *finalParams,
              field);

      if (std::isfinite(value))
      {
        optimizedValues[field][v] =
            value;
      }
    }

    updateProgress();
  }

  q->ProgressBar->setValue(100);

  std::cout
      << "[TCM OPT] Selection COMPLETE"
      << " | criterion="
      << criterion.toStdString()
      << " | tests="
      << (useTests ? "ON" : "OFF")
      << " | invalid="
      << invalidVoxelCount
      << " | simplified="
      << simplifiedVoxelCount
      << " | LRT comparisons="
      << lrtComparisonCount
      << " | Vuong comparisons="
      << vuongComparisonCount
      << " | skipped statistical comparisons="
      << skippedStatComparisonCount
      << std::endl;

  for (const auto& item :
       selectedCounts)
  {
    std::cout
        << "[TCM OPT] "
        << item.first
        << " selected: "
        << item.second
        << std::endl;
  }

  // ------------------------------------------------------------------------
  // Provenance strings
  // ------------------------------------------------------------------------

  QString comparedModels;

  for (size_t i = 0;
       i < selectedModels.size();
       ++i)
  {
    if (i > 0)
    {
      comparedModels += ", ";
    }

    comparedModels +=
        QString::fromStdString(
            selectedModels[i]);
  }

  QString tcmDerivationDetails =
      QString(
          "selectionCriterion=%1; "
          "comparedModels=%2; "
          "statisticalTests=%3")
          .arg(
              criterion,
              comparedModels,
              useTests ? "ON" : "OFF");

  if (useTests)
  {
    tcmDerivationDetails +=
        QString(
            "; alpha=%1; "
            "nestedComparison=LRT; "
            "nonNestedComparison=Vuong; "
            "selectionPolicy=prefer simpler model "
            "when p>=alpha")
            .arg(
                alpha,
                0,
                'g',
                8);
  }

  QString suffix =
      criterion;

  if (useTests)
  {
    suffix +=
        QString(" - Tests alpha=%1")
            .arg(
                alpha,
                0,
                'g',
                3);
  }

  auto setMetadata =
      [&](vtkMRMLNode* node)
      {
        if (!node)
        {
          return;
        }

        node->SetAttribute(
            "SlicerDynamicPET.TCM.SelectionCriterion",
            criterion
                .toUtf8()
                .constData());

        node->SetAttribute(
            "SlicerDynamicPET.TCM.ComparedModels",
            comparedModels
                .toUtf8()
                .constData());

        node->SetAttribute(
            "SlicerDynamicPET.TCM.UseStatTests",
            useTests ? "1" : "0");

        if (useTests)
        {
          node->SetAttribute(
              "SlicerDynamicPET.TCM.StatAlpha",
              QString::number(
                  alpha,
                  'g',
                  8)
                  .toUtf8()
                  .constData());
        }
      };

  // ------------------------------------------------------------------------
  // Show quantitative maps + selection map in Slicer
  // ------------------------------------------------------------------------

  if (this->TCMShowInSlicerCheckBoxImg
          ->isChecked())
  {
    this->removeTCMOptimizedSceneNodes();

    const vtkIdType refItemID =
        shNode->GetItemByDataNode(
            refPETNode);

    const vtkIdType parentItemID =
        shNode->GetItemParent(
            refItemID);

    auto createVolume =
        [&](const std::vector<double>& values,
            const QString& name)
            -> vtkMRMLScalarVolumeNode*
        {
          vtkMRMLScalarVolumeNode* node =
              logic->Flatten2Image(
                  values,
                  q->PETdims,
                  name.toStdString());

          if (!node)
          {
            return nullptr;
          }

          node->CopyOrientation(
              refPETNode);

          node->SetSpacing(
              refPETNode->GetSpacing());

          node->SetOrigin(
              refPETNode->GetOrigin());

          const vtkIdType itemID =
              shNode->GetItemByDataNode(
                  node);

          if (itemID !=
              vtkMRMLSubjectHierarchyNode::
                  INVALID_ITEM_ID)
          {
            shNode->SetItemParent(
                itemID,
                parentItemID);
          }

          return node;
        };

    q->ProgressBar->setMinimum(0);
    q->ProgressBar->setMaximum(0);
    q->ProgressBar->setFormat(
        "Creating optimized TCM maps...");
    QApplication::processEvents();

    for (const std::string& field :
         outputFields)
    {
      const QString name =
          "TCM Optimized - " +
          QString::fromStdString(field) +
          " - " +
          suffix;

      vtkMRMLScalarVolumeNode* node =
          createVolume(
              optimizedValues.at(field),
              name);

      if (!node)
      {
        continue;
      }

      setMetadata(node);

      this->TCMOptimizedNodeIDs
          .push_back(
              node->GetID());
    }

    // ------------------------------------------------------
    // Model-selection visualization.
    // ------------------------------------------------------

    std::vector<double> selectionValues(
        numberOfVoxels,
        0.0);

    for (size_t v = 0;
         v < numberOfVoxels;
         ++v)
    {
      selectionValues[v] =
          static_cast<double>(
              selectedModelMap[v]);
    }

    vtkMRMLScalarVolumeNode*
        selectionNode =
            createVolume(
                selectionValues,
                "TCM Optimized Model Selection - " +
                suffix);

    if (selectionNode)
    {
      setMetadata(selectionNode);

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.0",
          "Excluded");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.1",
          "1TCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.2",
          "1TdCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.3",
          "1TiCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.4",
          "1TidCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.5",
          "2TCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.6",
          "2TdCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.7",
          "2TiCM");

      selectionNode->SetAttribute(
          "SlicerDynamicPET.TCM.ModelCode.8",
          "2TidCM");

      vtkMRMLScalarVolumeDisplayNode*
          displayNode =
              vtkMRMLScalarVolumeDisplayNode::
                  SafeDownCast(
                      selectionNode->
                          GetDisplayNode());

      if (displayNode)
      {
        displayNode->AutoWindowLevelOff();

        displayNode->SetWindowLevelMinMax(
            0.0,
            8.0);

        displayNode->SetAndObserveColorNodeID(
            "vtkMRMLColorTableNodeLabels");
      }

      this->TCMOptimizedModelSelectionNodeID =
          selectionNode->GetID();
    }
  }

  // ------------------------------------------------------------------------
  // DICOM PMAP export:
  // only quantitative parameter maps.
  // Selection map is visualization-only.
  // ------------------------------------------------------------------------

  if (this->TCMSaveDICOMCheckBoxImg
          ->isChecked())
  {
    const QString outputDirectory =
        this->TCMDICOMDirectoryImg
            ->currentPath()
            .trimmed();

    auto fieldIndex =
        [](const std::string& field)
        {
          if (field == "K1") return 0;
          if (field == "k2") return 1;
          if (field == "k3") return 2;
          if (field == "k4") return 3;
          if (field == "vb") return 4;
          if (field == "td") return 5;
          if (field == "Ki") return 6;
          if (field == "DV") return 7;

          return 99;
        };

    for (const std::string& field :
         outputFields)
    {
      QString unitCode = "1";
      QString unitMeaning = "no units";

      if (field == "K1" || field == "Ki")
      {
        unitCode = "mL/(cm3.min)";
        unitMeaning =
            "milliliters per cubic centimeter per minute";
      }
      else if (field == "k2" ||
               field == "k3" ||
               field == "k4")
      {
        unitCode = "/min";
        unitMeaning = "per minute";
      }
      else if (field == "DV")
      {
        unitCode = "mL/cm3";
        unitMeaning = "milliliters per cubic centimeter";
      }
      else if (field == "td")
      {
        unitCode = "s";
        unitMeaning = "seconds";
      }

      q->ProgressBar->setMinimum(0);
      q->ProgressBar->setMaximum(0);

      q->ProgressBar->setFormat(
          "Saving DICOM PMAP: "
          "TCM Optimized - " +
          QString::fromStdString(field) +
          "...");

      QApplication::processEvents();

      const bool ok =
          this->exportParametricMapDICOM(
              refPETNode,
              optimizedValues.at(field),
              "TCM",
              "TCMOptimized",
              field,
              outputDirectory,
              7400 + fieldIndex(field),
              unitCode,
              unitMeaning,
              tcmDerivationDetails);

      if (!ok)
      {
        break;
      }
    }
  }

  q->ProgressBar->hide();
  q->ProgressBar->setMinimum(0);
  q->ProgressBar->setMaximum(100);
  q->ProgressBar->setValue(0);
  q->ProgressBar->setFormat("%p%");

  std::cout
      << "[TCM OPT] END"
      << std::endl;
}

void qSlicerDynamicPETModuleWidgetPrivate::
outputMTGAParametricResult(
    const std::string& modelID,
    vtkSlicerDynamicPETLogic* logic,
    vtkMRMLScalarVolumeNode* refPETNode,
    vtkMRMLSubjectHierarchyNode* shNode,
    vtkIdType refPetID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (!logic || !refPETNode || !shNode)
  {
    return;
  }

  auto resultIt =
      q->MTGAImgOutcomes.find(modelID);

  if (resultIt == q->MTGAImgOutcomes.end() ||
      resultIt->second.empty())
  {
    return;
  }

  // ------------------------------------------------------------------------
  // Show in Slicer
  // ------------------------------------------------------------------------
  std::vector<std::string> fields;

  if (modelID == "Patlak")
  {
    fields =
    {
      "Ki",
      "Intercept",
      "AIC",
      "MASE",
      "R2",
      "chi2"
    };
  }
  else if (modelID == "Logan" ||
           modelID == "RE")
  {
    fields =
    {
      "DV",
      "Intercept",
      "AIC",
      "MASE",
      "R2",
      "chi2"
    };
  }
  else
  {
    qWarning()
        << "Unknown MTGA model:"
        << QString::fromStdString(modelID);
    return;
  }

  // ------------------------------------------------------------------------
  // Show in Slicer
  // ------------------------------------------------------------------------
  if (this->MTGAShowInSlicerCheckBoxImg->isChecked())
  {
    q->ProgressBar->setMinimum(0);
    q->ProgressBar->setMaximum(0);

    q->ProgressBar->setFormat(
        "Creating " +
        QString::fromStdString(modelID) +
        " maps in Slicer...");

    QApplication::processEvents();

    logic->CreateMTGAParametricImages(
        resultIt->second,
        q->PETdims,
        fields,
        modelID,
        refPETNode,
        shNode,
        refPetID);
  }

  if (this->MTGASaveDICOMCheckBoxImg->isChecked())
  {
    const QString outputDirectory =
        this->MTGADICOMDirectoryImg
            ->currentPath()
            .trimmed();

    int modelIndex = 0;

    if (modelID == "Patlak")      modelIndex = 0;
    else if (modelID == "Logan")  modelIndex = 1;
    else if (modelID == "RE")     modelIndex = 2;

    const double framingNorm =
        this->framingNormEditImg
            ->text()
            .toDouble();

    for (int fieldIndex = 0;
         fieldIndex <
             static_cast<int>(fields.size());
         ++fieldIndex)
    {
      const std::string& field =
          fields[fieldIndex];

      QString unitCode = "1";
      QString unitMeaning = "no units";

      bool requiresNormalizedTimeUnit = false;

      // Patlak slope Ki has inverse normalized-time units.
      if (modelID == "Patlak" &&
          field == "Ki")
      {
        requiresNormalizedTimeUnit = true;

        if (std::abs(
                framingNorm - 60.0) < 1e-9)
        {
          unitCode = "/min";
          unitMeaning = "per minute";
        }
        else if (std::abs(
                     framingNorm - 1.0) < 1e-9)
        {
          unitCode = "/s";
          unitMeaning = "per second";
        }
      }

      // Logan/RE regression intercept is a time quantity.
      if ((modelID == "Logan" ||
           modelID == "RE") &&
          field == "Intercept")
      {
        requiresNormalizedTimeUnit = true;

        if (std::abs(
                framingNorm - 60.0) < 1e-9)
        {
          unitCode = "min";
          unitMeaning = "minutes";
        }
        else if (std::abs(
                     framingNorm - 1.0) < 1e-9)
        {
          unitCode = "s";
          unitMeaning = "seconds";
        }
      }

      if (requiresNormalizedTimeUnit &&
          std::abs(framingNorm - 60.0) >= 1e-9 &&
          std::abs(framingNorm - 1.0) >= 1e-9)
      {
        QMessageBox::warning(
            q,
            QObject::tr("DICOM PMAP export"),
            QObject::tr(
                "%1 - %2 was not exported because "
                "Framing Norm is %3 s.\n\n"
                "Its physical unit therefore cannot be "
                "represented honestly as seconds or minutes "
                "without rescaling the numerical values.")
                .arg(
                    QString::fromStdString(modelID),
                    QString::fromStdString(field))
                .arg(framingNorm));

        continue;
      }

      const std::vector<double> values =
          logic->ExtractParameter(
              resultIt->second,
              field);

      const int seriesNumber =
          7100 +
          modelIndex * 20 +
          fieldIndex;

      q->ProgressBar->setMinimum(0);
      q->ProgressBar->setMaximum(0);

      q->ProgressBar->setFormat(
          "Saving DICOM PMAP: " +
          QString::fromStdString(modelID) +
          " - " +
          QString::fromStdString(field) +
          "...");

      QApplication::processEvents();

      q->ProgressBar->setMinimum(0);
      q->ProgressBar->setMaximum(0);

      if (!this->exportParametricMapDICOM(
              refPETNode,
              values,
              "MTGA",
              modelID,
              field,
              outputDirectory,
              seriesNumber,
              unitCode,
              unitMeaning))
      {
        break;
      }
    }
  }
}

bool qSlicerDynamicPETModuleWidgetPrivate::
exportParametricMapDICOM(
    vtkMRMLScalarVolumeNode* refPETNode,
    const std::vector<double>& values,
    const std::string& method,
    const std::string& modelID,
    const std::string& field,
    const QString& outputDirectory,
    int seriesNumber,
    const QString& unitCode,
    const QString& unitMeaning,
    const QString& derivationDetails)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (!refPETNode ||
      values.empty() ||
      outputDirectory.trimmed().isEmpty())
  {
    return false;
  }

  const size_t expectedSize =
      static_cast<size_t>(q->PETdims[0]) *
      static_cast<size_t>(q->PETdims[1]) *
      static_cast<size_t>(q->PETdims[2]);

  if (values.size() != expectedSize)
  {
    QMessageBox::warning(
        q,
        QObject::tr("DICOM PMAP export"),
        QObject::tr(
            "Cannot export %1 - %2:\n"
            "voxel count does not match PET geometry.")
            .arg(
                QString::fromStdString(modelID),
                QString::fromStdString(field)));

    return false;
  }

  // ------------------------------------------------------------------------
  // Find DICOM source instance UIDs.
  //
  // Prefer the first actual PET sequence frame, because the displayed
  // proxy node may or may not carry the original DICOM attributes.
  // ------------------------------------------------------------------------

  QStringList geometryUIDList;
  QStringList allSourceUIDList;

  QSet<QString> geometryUIDSet;
  QSet<QString> allSourceUIDSet;

  // ------------------------------------------------------------------------
  // We need two different source sets:
  //
  // geometryUIDList:
  //   DICOM instances for ONE temporal PET frame only.
  //   Classic PET   -> complete slice stack for one time point.
  //   Enhanced PET  -> one multiframe SOP instance.
  //
  // allSourceUIDList:
  //   every unique SOP Instance UID that contributed to the
  //   dynamic kinetic fit. These are provenance references only.
  // ------------------------------------------------------------------------

  if (q->sequencePETNode)
  {
    const int numberOfFrames =
        q->sequencePETNode->GetNumberOfDataNodes();

    for (int frameIndex = 0;
         frameIndex < numberOfFrames;
         ++frameIndex)
    {
      vtkMRMLNode* frameNode =
          q->sequencePETNode->GetNthDataNode(
              frameIndex);

      if (!frameNode)
      {
        continue;
      }

      const char* attr =
          frameNode->GetAttribute(
              "DICOM.instanceUIDs");

      if (!attr)
      {
        continue;
      }

      const QStringList frameUIDs =
          QString::fromUtf8(attr)
              .split(
                  QRegularExpression("\\s+"),
                  Qt::SkipEmptyParts);

      if (frameUIDs.isEmpty())
      {
        continue;
      }

      // First valid dynamic frame becomes the spatial
      // reference set used by highdicom.
      if (geometryUIDList.isEmpty())
      {
        for (const QString& uid : frameUIDs)
        {
          if (!geometryUIDSet.contains(uid))
          {
            geometryUIDSet.insert(uid);
            geometryUIDList.append(uid);
          }
        }
      }

      // All temporal source SOPs are retained as provenance.
      for (const QString& uid : frameUIDs)
      {
        if (!allSourceUIDSet.contains(uid))
        {
          allSourceUIDSet.insert(uid);
          allSourceUIDList.append(uid);
        }
      }
    }
  }


  // Fallback if sequence frames do not retain source references.
  if (geometryUIDList.isEmpty())
  {
    const char* attr =
        refPETNode->GetAttribute(
            "DICOM.instanceUIDs");

    if (attr)
    {
      const QStringList proxyUIDs =
          QString::fromUtf8(attr)
              .split(
                  QRegularExpression("\\s+"),
                  Qt::SkipEmptyParts);

      for (const QString& uid : proxyUIDs)
      {
        if (!geometryUIDSet.contains(uid))
        {
          geometryUIDSet.insert(uid);
          geometryUIDList.append(uid);
        }

        if (!allSourceUIDSet.contains(uid))
        {
          allSourceUIDSet.insert(uid);
          allSourceUIDList.append(uid);
        }
      }
    }
  }


  const QString geometryInstanceUIDs =
      geometryUIDList.join(" ");

  const QString allInstanceUIDs =
      allSourceUIDList.join(" ");

  if (geometryInstanceUIDs.isEmpty() ||
      allInstanceUIDs.isEmpty())
  {
    QMessageBox::warning(
        q,
        QObject::tr("DICOM PMAP export"),
        QObject::tr(
            "The source PET does not contain "
            "DICOM.instanceUIDs.\n\n"
            "A standards-based DICOM Parametric Map needs "
            "the original DICOM patient/study context, "
            "therefore this map was not exported."));

    return false;
  }

  // ------------------------------------------------------------------------
  // Quantity semantics.
  //
  // These are intentionally LOCAL/private codes under 99SDPET.
  // ------------------------------------------------------------------------

  QString quantityCode;
  QString quantityMeaning;

  if (field == "K1")
  {
    quantityCode = "SDP_K1";
    quantityMeaning = "K1";
  }
  else if (field == "k2")
  {
    quantityCode = "SDP_K2";
    quantityMeaning = "k2";
  }
  else if (field == "k3")
  {
    quantityCode = "SDP_K3";
    quantityMeaning = "k3";
  }
  else if (field == "k4")
  {
    quantityCode = "SDP_K4";
    quantityMeaning = "k4";
  }
  else if (field == "vb")
  {
    quantityCode = "SDP_VB";
    quantityMeaning = "Blood volume fraction";
  }
  else if (field == "td")
  {
    quantityCode = "SDP_TD";
    quantityMeaning = "Time delay";
  }
  else if (field == "Ki")
  {
    quantityCode = "SDP_KI";
    quantityMeaning = "Net influx rate";
  }
  else if (field == "DV")
  {
    quantityCode = "SDP_DV";
    quantityMeaning = "Distribution volume";
  }
  else if (field == "Intercept")
  {
    quantityCode = "SDP_INT";
    quantityMeaning = "Regression intercept";
  }
  else if (field == "AIC")
  {
    quantityCode = "SDP_AIC";
    quantityMeaning = "Akaike information criterion";
  }
  else if (field == "BIC")
  {
    quantityCode = "SDP_BIC";
    quantityMeaning = "Bayesian information criterion";
  }
  else if (field == "MASE")
  {
    quantityCode = "SDP_MASE";
    quantityMeaning = "Mean absolute scaled error";
  }
  else if (field == "R2")
  {
    quantityCode = "SDP_R2";
    quantityMeaning = "Coefficient of determination";
  }
  else if (field == "chi2")
  {
    quantityCode = "SDP_CHI2";
    quantityMeaning = "Chi-square statistic";
  }
  else
  {
    qWarning()
        << "No DICOM PMAP quantity mapping for"
        << QString::fromStdString(field);

    return false;
  }

  // ------------------------------------------------------------------------
  // Measurement method semantics.
  // Also intentionally private/local codes.
  // ------------------------------------------------------------------------

  QString methodCode;

  if (modelID == "Patlak")
    methodCode = "SDP_PATLAK";
  else if (modelID == "Logan")
    methodCode = "SDP_LOGAN";
  else if (modelID == "RE")
    methodCode = "SDP_RE";
  else if (modelID == "MTGAOptimized")
    methodCode = "SDP_MTGAOPT";
  else if (modelID == "1TCM")
    methodCode = "SDP_1TCM";
  else if (modelID == "1TdCM")
    methodCode = "SDP_1TDCM";
  else if (modelID == "1TiCM")
    methodCode = "SDP_1TICM";
  else if (modelID == "1TidCM")
    methodCode = "SDP_1TIDCM";
  else if (modelID == "2TCM")
    methodCode = "SDP_2TCM";
  else if (modelID == "2TdCM")
    methodCode = "SDP_2DTCM";
  else if (modelID == "2TiCM")
    methodCode = "SDP_2TICM";
  else if (modelID == "2TidCM")
    methodCode = "SDP_2TIDCM";
  else if (modelID == "TCMOptimized")
    methodCode = "SDP_TCMOPT";
  else
  {
    qWarning()
        << "No DICOM PMAP method mapping for"
        << QString::fromStdString(modelID);

    return false;
  }

  vtkMRMLScene* scene =
      q->mrmlScene();

  if (!scene)
  {
    return false;
  }

  // ------------------------------------------------------------------------
  // Temporary scalar volume.
  //
  // No display node is created. This exists only long enough to save the
  // NRRD that dcmqi consumes.
  // ------------------------------------------------------------------------

  vtkNew<vtkImageData> image;

  image->SetDimensions(
      q->PETdims[0],
      q->PETdims[1],
      q->PETdims[2]);

  image->AllocateScalars(
      VTK_DOUBLE,
      1);

  double* destination =
      static_cast<double*>(
          image->GetScalarPointer());

  if (!destination)
  {
    return false;
  }

  std::copy(
      values.begin(),
      values.end(),
      destination);

  vtkMRMLScalarVolumeNode* tempNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          scene->AddNewNodeByClass(
              "vtkMRMLScalarVolumeNode"));

  if (!tempNode)
  {
    return false;
  }

  tempNode->SetName(
      "SlicerDynamicPET_PMAP_Temporary");

  tempNode->SetAndObserveImageData(
      image.GetPointer());

  // Exact spatial geometry of the source PET.
  tempNode->CopyOrientation(refPETNode);
  tempNode->SetSpacing(
      refPETNode->GetSpacing());
  tempNode->SetOrigin(
      refPETNode->GetOrigin());

  // ------------------------------------------------------------------------
  // Output directory / filename.
  // ------------------------------------------------------------------------

  if (!QDir(outputDirectory).exists() &&
      !QDir().mkpath(outputDirectory))
  {
    scene->RemoveNode(tempNode);

    QMessageBox::warning(
        q,
        QObject::tr("DICOM PMAP export"),
        QObject::tr(
            "Could not create the DICOM output directory:\n%1")
            .arg(outputDirectory));

    return false;
  }

  QDir outputDir(outputDirectory);

  const QString methodQString =
      QString::fromStdString(method);

  const QString modelQString =
      QString::fromStdString(modelID);

  const QString fieldQString =
      QString::fromStdString(field);

  const QString outputPath =
      outputDir.filePath(
          QString(
              "SlicerDynamicPET_%1_%2_%3.dcm")
              .arg(
                  methodQString,
                  modelQString,
                  fieldQString));

  const QString seriesDescription =
      QString(
          "SlicerDynamicPET %1 %2 %3")
          .arg(
              methodQString,
              modelQString,
              fieldQString);

  const QString methodMeaning =
      QString(
          "SlicerDynamicPET %1")
          .arg(modelQString);

  // ------------------------------------------------------------------------
  // dcmqi conversion.
  // ------------------------------------------------------------------------

  PythonQtObjectPtr mainContext =
      PythonQt::self()->getMainModule();

  QVariant result =
      mainContext.call(
          "DPE_export_parametric_map",
          QVariantList{
              QString::fromUtf8(
                  tempNode->GetID()),
              geometryInstanceUIDs,
              allInstanceUIDs,
              outputPath,
              seriesDescription,
              seriesNumber,
              quantityCode,
              quantityMeaning,
              methodCode,
              methodMeaning,
              unitCode,
              unitMeaning,
              derivationDetails
          });

  // Save-only must not leave anything in the MRML scene.
  scene->RemoveNode(tempNode);

  const QVariantMap resultMap =
      result.toMap();

  const bool ok =
      resultMap.value("ok").toBool();

  if (!ok)
  {
    const QString error =
        resultMap.value("error").toString();

    QMessageBox::warning(
        q,
        QObject::tr("DICOM PMAP export"),
        QObject::tr(
            "Could not export %1 - %2.\n\n%3")
            .arg(
                modelQString,
                fieldQString,
                error));

    return false;
  }

  qDebug()
      << "DICOM PMAP exported:"
      << resultMap.value("path").toString();

  return true;
}


void qSlicerDynamicPETModuleWidgetPrivate::
outputTCMParametricResult(
    const std::string& modelID,
    vtkSlicerDynamicPETLogic* logic,
    vtkMRMLScalarVolumeNode* refPETNode,
    vtkMRMLSubjectHierarchyNode* shNode,
    vtkIdType refPetID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (!logic || !refPETNode || !shNode)
  {
    return;
  }

  auto resultIt =
      q->TCMImgOutcomes.find(modelID);

  if (resultIt == q->TCMImgOutcomes.end() ||
      resultIt->second.empty())
  {
    return;
  }

  std::vector<std::string> fields;

  if (modelID == "1TCM")
  {
    fields =
    {
      "K1", "k2", "vb",
      "Ki", "DV",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "1TdCM")
  {
    fields =
    {
      "K1", "k2", "vb", "td",
      "Ki", "DV",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "1TiCM")
  {
    fields =
    {
      "K1", "vb", "Ki",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "1TidCM")
  {
    fields =
    {
      "K1", "vb", "td", "Ki",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "2TCM")
  {
    fields =
    {
      "K1", "k2", "k3", "k4", "vb",
      "Ki", "DV",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "2TdCM")
  {
    fields =
    {
      "K1", "k2", "k3", "k4", "vb", "td",
      "Ki", "DV",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "2TiCM")
  {
    fields =
    {
      "K1", "k2", "k3", "vb",
      "Ki",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else if (modelID == "2TidCM")
  {
    fields =
    {
      "K1", "k2", "k3", "vb", "td",
      "Ki",
      "AIC", "MASE", "BIC", "chi2"
    };
  }
  else
  {
    qWarning()
        << "Unknown TCM model:"
        << QString::fromStdString(modelID);
    return;
  }

  // ------------------------------------------------------------------------
  // Show in Slicer
  // ------------------------------------------------------------------------
  if (this->TCMShowInSlicerCheckBoxImg->isChecked())
  {
    q->ProgressBar->setMinimum(0);
    q->ProgressBar->setMaximum(0);

    q->ProgressBar->setFormat(
        "Creating " +
        QString::fromStdString(modelID) +
        " maps in Slicer...");

    QApplication::processEvents();

    logic->CreateTCMParametricImages(
        resultIt->second,
        q->PETdims,
        fields,
        modelID,
        refPETNode,
        shNode,
        refPetID);
  }

  if (this->TCMSaveDICOMCheckBoxImg->isChecked())
  {
    const QString outputDirectory =
        this->TCMDICOMDirectoryImg
            ->currentPath()
            .trimmed();

    int modelIndex = 0;

    if (modelID == "1TCM")        modelIndex = 0;
    else if (modelID == "1TdCM")  modelIndex = 1;
    else if (modelID == "1TiCM")  modelIndex = 2;
    else if (modelID == "1TidCM") modelIndex = 3;
    else if (modelID == "2TCM")   modelIndex = 4;
    else if (modelID == "2TdCM")  modelIndex = 5;
    else if (modelID == "2TiCM")  modelIndex = 6;
    else if (modelID == "2TidCM") modelIndex = 7;

    for (int fieldIndex = 0;
         fieldIndex <
             static_cast<int>(fields.size());
         ++fieldIndex)
    {
      const std::string& field =
          fields[fieldIndex];

      QString unitCode = "1";
      QString unitMeaning = "no units";

      // KMAP rate constants are parameterized per minute.  The convolution
      // converts the second-based numerical integration step using td/60.
      if (field == "K1" || field == "Ki")
      {
        unitCode = "mL/(cm3.min)";
        unitMeaning =
            "milliliters per cubic centimeter per minute";
      }
      else if (field == "k2" ||
               field == "k3" ||
               field == "k4")
      {
        unitCode = "/min";
        unitMeaning = "per minute";
      }
      else if (field == "DV")
      {
        unitCode = "mL/cm3";
        unitMeaning = "milliliters per cubic centimeter";
      }
      else if (field == "td")
      {
        unitCode = "s";
        unitMeaning = "seconds";
      }

      const std::vector<double> values =
          logic->ExtractParameter(
              resultIt->second,
              field);

      const int seriesNumber =
          7200 +
          modelIndex * 20 +
          fieldIndex;

      q->ProgressBar->setMinimum(0);
      q->ProgressBar->setMaximum(0);

      q->ProgressBar->setFormat(
          "Saving DICOM PMAP: " +
          QString::fromStdString(modelID) +
          " - " +
          QString::fromStdString(field) +
          "...");

      QApplication::processEvents();

      if (!this->exportParametricMapDICOM(
              refPETNode,
              values,
              "TCM",
              modelID,
              field,
              outputDirectory,
              seriesNumber,
              unitCode,
              unitMeaning))
      {
        // Avoid producing one error dialog for every
        // remaining parameter if dcmqi/source DICOM
        // itself is unavailable.
        break;
      }
    }
  }
}

void qSlicerDynamicPETModuleWidgetPrivate::
invalidateParametricVoxelSelection()
{
    this->parametricVoxelMask.clear();
    this->parametricFitVoxelIndices.clear();

    this->parametricBodySupportImage = nullptr;
    this->parametricFinalFitMaskImage = nullptr;

    this->parametricVoxelSelectionPETID =
        vtkMRMLSubjectHierarchyNode::
            INVALID_ITEM_ID;

    this->parametricVoxelSelectionCTID =
        vtkMRMLSubjectHierarchyNode::
            INVALID_ITEM_ID;

    // A changed support mask means old voxelwise fits no longer represent
    // the current analysis domain.
    this->MTGAImgFitSignatures.clear();
    this->TCMImgFitSignatures.clear();

    Q_Q(qSlicerDynamicPETModuleWidget);

    q->MTGAImgOutcomes.clear();
    q->TCMImgOutcomes.clear();
}

bool qSlicerDynamicPETModuleWidgetPrivate::
ensureParametricVoxelSelection()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (q->petID ==
            vtkMRMLSubjectHierarchyNode::
                INVALID_ITEM_ID)
    {
        return false;
    }

    QString petError;

    if (!this->ensureParametricPETData(
            &petError))
    {
        qWarning()
            << "Parametric PET preparation:"
            << petError;

        return false;
    }

    if (q->PET_flatten_values.empty())
    {
        return false;
    }

    const size_t numberOfVoxels =
        q->PET_flatten_values.size();

    // Cache valid for current PET + CT.
    if (this->parametricVoxelSelectionPETID ==
            q->petID &&
        this->parametricVoxelSelectionCTID ==
            q->ctID &&
        this->parametricVoxelMask.size() ==
            numberOfVoxels &&
        this->parametricFinalFitMaskImage)
    {
        return true;
    }

    vtkMRMLScene* scene =
        q->mrmlScene();

    if (!scene)
    {
        return false;
    }

    vtkMRMLSubjectHierarchyNode* shNode =
        vtkMRMLSubjectHierarchyNode::
            GetSubjectHierarchyNode(scene);

    if (!shNode)
    {
        return false;
    }

    vtkMRMLScalarVolumeNode* petNode =
        vtkMRMLScalarVolumeNode::SafeDownCast(
            shNode->GetItemDataNode(
                q->petID));

    vtkMRMLScalarVolumeNode* ctNode =
        nullptr;

    if (q->ctID !=
        vtkMRMLSubjectHierarchyNode::
            INVALID_ITEM_ID)
    {
        ctNode =
            vtkMRMLScalarVolumeNode::SafeDownCast(
                shNode->GetItemDataNode(
                    q->ctID));
    }

    if (!petNode)
    {
        return false;
    }

    vtkSlicerDynamicPETLogic* logic =
        vtkSlicerDynamicPETLogic::SafeDownCast(
            q->logic());

    if (!logic)
    {
        return false;
    }

    // For now these are fixed defaults.
    const double ctThresholdHU =
        this->BodySupportCTThresholdImg->
            value();

    const double ctBodyMarginMm =
        this->BodySupportMarginImg->
            value();

    const bool fillHoles =
        this->BodySupportFillHolesCheckBoxImg->
            isChecked();

    const bool supportEnabled =
        this->BodySupportEnabledCheckBoxImg->
            isChecked();

    const BodySupportSource requestedSource =
        static_cast<BodySupportSource>(
            this->BodySupportSourceImg->
                currentIndex());

    BodySupportSource effectiveSource =
        requestedSource;

    if (requestedSource ==
        BodySupportSource::Auto)
    {
        effectiveSource =
            ctNode
            ? BodySupportSource::CT
            : BodySupportSource::PET;
    }

    // PET composite choice.
    const PETCompositeMode compositeMode =
        this->BodySupportPETCompositeImg->
            currentIndex() == 1
            ? PETCompositeMode::DurationWeightedSum
            : PETCompositeMode::UnweightedSum;

    vtkSmartPointer<vtkOrientedImageData> ctMask;
    vtkSmartPointer<vtkOrientedImageData> petMask;


    // ------------------------------------------------------------------
    // Body-support enabled
    // ------------------------------------------------------------------

    if (supportEnabled)
    {
        // --------------------------------------------------------------
        // CT mask if required.
        // --------------------------------------------------------------

        const bool needCT =
            effectiveSource ==
                BodySupportSource::CT ||
            effectiveSource ==
                BodySupportSource::Union ||
            effectiveSource ==
                BodySupportSource::Intersection;

        if (needCT)
        {
            if (!ctNode)
            {
                std::cerr
                    << "Patient support requires CT, "
                       "but no CT is selected."
                    << std::endl;

                return false;
            }

            ctMask =
                logic->CreateCTBodySupportMask(
                    ctNode,
                    petNode,
                    ctThresholdHU,
                    ctBodyMarginMm,
                    fillHoles);

            if (!ctMask)
            {
                std::cerr
                    << "Could not create CT body-support mask."
                    << std::endl;

                return false;
            }
        }


        // --------------------------------------------------------------
        // PET mask if required.
        // --------------------------------------------------------------

        const bool needPET =
            effectiveSource ==
                BodySupportSource::PET ||
            effectiveSource ==
                BodySupportSource::Union ||
            effectiveSource ==
                BodySupportSource::Intersection;

        if (needPET)
        {
            double petThreshold =
                std::numeric_limits<double>::
                    quiet_NaN();

            bool usedOtsuFallback =
                false;

            petMask =
                logic->CreatePETBodySupportMask(
                    q->PET_flatten_values,
                    q->PETdims,
                    petNode,
                    q->durations,
                    compositeMode,

                    // PET is already in PET geometry.
                    // Keep no extra PET dilation for now.
                    0.0,

                    fillHoles,

                    // Keep components >= 1% of largest.
                    // Make this Advanced later.
                    0.01,

                    &petThreshold,
                    &usedOtsuFallback);

            if (!petMask)
            {
                std::cerr
                    << "Could not create PET body-support mask."
                    << std::endl;

                return false;
            }

            qDebug()
                << "PET patient support:"
                << "threshold ="
                << petThreshold
                << ", threshold method ="
                << (usedOtsuFallback
                        ? "log-Otsu fallback"
                        : "multiscale log-histogram")
                << ", composite ="
                << (compositeMode ==
                        PETCompositeMode::
                            DurationWeightedSum
                        ? "duration-weighted"
                        : "unweighted");
        }


        // --------------------------------------------------------------
        // Select / combine final support.
        // --------------------------------------------------------------

        switch (effectiveSource)
        {
            case BodySupportSource::CT:
            {
                this->parametricBodySupportImage =
                    ctMask;

                break;
            }

            case BodySupportSource::PET:
            {
                this->parametricBodySupportImage =
                    petMask;

                break;
            }

            case BodySupportSource::Union:
            {
                this->parametricBodySupportImage =
                    logic->CombineBodySupportMasks(
                        ctMask,
                        petMask,
                        BodySupportCombination::Union);

                break;
            }

            case BodySupportSource::Intersection:
            {
                this->parametricBodySupportImage =
                    logic->CombineBodySupportMasks(
                        ctMask,
                        petMask,
                        BodySupportCombination::Intersection);

                break;
            }

            default:
            {
                return false;
            }
        }
    }
    else
    {
        // No anatomical patient-support restriction.
        // Numerical finite/non-zero TAC filtering below still applies.
        this->parametricBodySupportImage =
            logic->CreateFullPETSupportMask(
                petNode);
    }

    if (!this->parametricBodySupportImage)
    {
        return false;
    }

    vtkDataArray* supportArray =
        this->parametricBodySupportImage
            ->GetPointData()
            ->GetScalars();

    if (!supportArray)
    {
        return false;
    }

    if (static_cast<size_t>(
            supportArray->GetNumberOfTuples()) !=
        numberOfVoxels)
    {
        std::cerr
            << "Body support mask size mismatch: "
            << supportArray->GetNumberOfTuples()
            << " vs "
            << numberOfVoxels
            << std::endl;

        return false;
    }

    this->parametricVoxelMask.assign(
        numberOfVoxels,
        static_cast<unsigned char>(0));

    this->parametricFitVoxelIndices.clear();

    this->parametricFitVoxelIndices.reserve(
        numberOfVoxels);

    // Start from the anatomical support image and turn off voxels that
    // fail the numerical TAC sanity checks.
    this->parametricFinalFitMaskImage =
        vtkSmartPointer<vtkOrientedImageData>::New();

    this->parametricFinalFitMaskImage->DeepCopy(
        this->parametricBodySupportImage);

    vtkDataArray* finalMaskArray =
        this->parametricFinalFitMaskImage
            ->GetPointData()
            ->GetScalars();

    if (!finalMaskArray)
    {
        return false;
    }

    size_t outsideBodyCount = 0;
    size_t nonFiniteCount = 0;
    size_t zeroTacCount = 0;

    for (size_t v = 0;
         v < numberOfVoxels;
         ++v)
    {
        const bool insideBody =
            supportArray->GetComponent(
                static_cast<vtkIdType>(v),
                0) != 0.0;

        if (!insideBody)
        {
            ++outsideBodyCount;

            finalMaskArray->SetComponent(
                static_cast<vtkIdType>(v),
                0,
                0.0);

            continue;
        }

        const auto& tac =
            q->PET_flatten_values[v];

        if (tac.empty())
        {
            ++zeroTacCount;

            finalMaskArray->SetComponent(
                static_cast<vtkIdType>(v),
                0,
                0.0);

            continue;
        }

        const bool allFinite =
            std::all_of(
                tac.begin(),
                tac.end(),
                [](double value)
                {
                    return std::isfinite(value);
                });

        if (!allFinite)
        {
            ++nonFiniteCount;

            finalMaskArray->SetComponent(
                static_cast<vtkIdType>(v),
                0,
                0.0);

            continue;
        }

        const bool allZero =
            std::all_of(
                tac.begin(),
                tac.end(),
                [](double value)
                {
                    return value == 0.0;
                });

        if (allZero)
        {
            ++zeroTacCount;

            finalMaskArray->SetComponent(
                static_cast<vtkIdType>(v),
                0,
                0.0);

            continue;
        }

        this->parametricVoxelMask[v] =
            static_cast<unsigned char>(1);

        this->parametricFitVoxelIndices.push_back(
            static_cast<int>(v));

        finalMaskArray->SetComponent(
            static_cast<vtkIdType>(v),
            0,
            1.0);
    }

    this->parametricFinalFitMaskImage->Modified();

    this->parametricVoxelSelectionPETID =
        q->petID;

    this->parametricVoxelSelectionCTID =
        q->ctID;

    qDebug()
        << "Parametric fitting support:"
        << this->parametricFitVoxelIndices.size()
        << "/"
        << numberOfVoxels
        << "voxels eligible;"
        << outsideBodyCount
        << "outside body,"
        << nonFiniteCount
        << "non-finite,"
        << zeroTacCount
        << "zero TAC.";

    return true;
}

void qSlicerDynamicPETModuleWidgetPrivate::
resetParametricImagingSelections()
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  // ------------------------------------------------------------------------
  // MTGA model selection
  // ------------------------------------------------------------------------
  q->modelsMTGAImgID.clear();

  for (int i = 0;
       i < this->ModelsCheckLayoutMTGAImg->count();
       ++i)
  {
    QLayoutItem* item =
        this->ModelsCheckLayoutMTGAImg->itemAt(i);

    QCheckBox* cb =
        qobject_cast<QCheckBox*>(item->widget());

    if (!cb)
    {
      continue;
    }

    cb->blockSignals(true);
    cb->setChecked(false);
    cb->blockSignals(false);
  }

  // ------------------------------------------------------------------------
  // TCM model selection
  // ------------------------------------------------------------------------
  q->modelsTCMImgID.clear();

  for (int i = 0;
       i < this->ModelsCheckLayoutTCMImg->count();
       ++i)
  {
    QLayoutItem* item =
        this->ModelsCheckLayoutTCMImg->itemAt(i);

    QCheckBox* cb =
        qobject_cast<QCheckBox*>(item->widget());

    if (!cb)
    {
      continue;
    }

    cb->blockSignals(true);
    cb->setChecked(false);
    cb->blockSignals(false);
  }

  // Model-selection results.
  q->TCMImgOutcomes.clear();
  q->MTGAImgOutcomes.clear();

  this->TCMImgFitSignatures.clear();
  this->MTGAImgFitSignatures.clear();

  this->TCMOptimizedNodeIDs.clear();
  this->TCMOptimizedModelSelectionNodeID.clear();

  this->MTGAOptimizedSelection.clear();
  this->MTGAOptimizedKiValues.clear();
  this->MTGAOptimizedDVValues.clear();

  this->MTGAOptimizedKiNodeID.clear();
  this->MTGAOptimizedDVNodeID.clear();
  this->MTGAOptimizedRGBNodeID.clear();

  this->RefreshMTGARGBButtonImg
      ->setEnabled(false);

  this->populateTCMOptimizationModels();
  this->updateMTGAOptimizationUI();

  // Explicitly disable them now.
  this->FITbuttonTCMImg->setEnabled(false);
  this->FITbuttonMTGAImg->setEnabled(false);
}

void qSlicerDynamicPETModuleWidgetPrivate::
setPETItemID(vtkIdType newPetID)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (newPetID == q->petID)
  {
    return;
  }

  // Everything below belongs to the previous PET.
  this->invalidateParametricVoxelSelection();

  q->PET_flatten_values.clear();

  q->sequencePETNode = nullptr;
  q->sequenceBrowserPETNode = nullptr;
  q->segSequenceNode = nullptr;

  if (q->SegWatcher)
  {
    q->SegWatcher->browser = nullptr;
  }

  q->durations.clear();
  q->timePoints.clear();
  q->suvFactors.clear();

  q->numberOfTimepoints = 0;

  q->petID = newPetID;

  this->updateBodySupportUI();
}


//-----------------------------------------------------------------------------
// qSlicerDynamicPETModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerDynamicPETModuleWidget::qSlicerDynamicPETModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerDynamicPETModuleWidgetPrivate(*this) )
{
  Q_D(qSlicerDynamicPETModuleWidget);
  this->SubjectHierarchyNode = nullptr;
  this->patID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  this->sequencePETNode = nullptr;
  this->segSequenceNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  this->numberOfTimepoints = 0;
  QMainWindow* mainWindow = qobject_cast<QMainWindow*>(qSlicerApplication::application()->mainWindow());
  this->ProgressBar = new QProgressBar(mainWindow);
  this->ProgressBar->setObjectName(QString::fromUtf8("ProgressBar"));
  this->ProgressBar->setMaximum(100);
  this->ProgressBar->setValue(0);
  this->ProgressBar->resize(300, 30);
  QRect parentGeometry = this->ProgressBar->parentWidget()->geometry();
  QRect barGeometry = this->ProgressBar->frameGeometry();
  this->ProgressBar->move(
      (parentGeometry.width() - barGeometry.width()) / 2,
      (parentGeometry.height() - barGeometry.height()) / 2
  );
  this->stopButton = new QPushButton("Stop", mainWindow);
  stopButton->move(this->ProgressBar->x(), this->ProgressBar->y() + 40);

  this->stopButton->setStyleSheet(
      "QPushButton {"
      "  background-color: red;"
      "  color: white;"        // Text color
      "  font-weight: bold;"
      "}"
      "QPushButton:pressed {"
      "  background-color: darkred;"
      "}"
  );
  this->stopRequested = false;
  QObject::connect(
      this->stopButton,
      &QPushButton::clicked,
      this,
      [this]()
      {
        this->stopRequested.store(true);
        this->stopButton->setEnabled(false);
        this->stopButton->setText("Stopping...");
      });

  this->checkboxNames = QStringList{
    "Mean", "Median", "Peak", "Min", "Max"//, "VoxelCount", "Volume(cc)"
  };
  this->ModelsNamesMTGA = QStringList{
    "Patlak", "Logan", "RE"
  };
  this->ModelsNamesTCM = QStringList{
    "1TCM", "1TdCM", "1TiCM", "1TidCM", "2TCM", "2TdCM", "2TiCM", "2TidCM"
  };
  this->StatsNames = QStringList{
    "Mean", "Median", "Peak"
  };
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI = "";
  // Install global key watcher
  if (mainWindow)
  {
    this->keyWatcher = new KeyPressWatcher(mainWindow);
    mainWindow->installEventFilter(this->keyWatcher);
    QObject::connect(this->keyWatcher, SIGNAL(deletePressed()), this, SLOT(onDeleteKeyPressed()));
  }

  // Optional: start with watcher disabled
  this->keyWatcher->setActive(true);
  this->PET_flatten_values.clear();
  this->PETdims[0] = 0;
  this->PETdims[1] = 0;
  this->PETdims[2] = 0;
  // qRegisterMetaType<std::vector<TCMParameters>>("std::vector<TCMParameters>");
  d->init();
}

//-----------------------------------------------------------------------------
qSlicerDynamicPETModuleWidget::~qSlicerDynamicPETModuleWidget()
{
}

// void qSlicerDynamicPETModuleWidget::setNodeSelectorEnabled(qMRMLNodeComboBox* selector, bool enabled)
// {
//   selector->setEnabled(enabled);
//   const auto children = selector->findChildren<QWidget*>();
//   for (QWidget* child : children)
//   {
//     child->setEnabled(enabled);
//   }
// }

// std::map<std::string, vtkIdType> qSlicerDynamicPETModuleWidget::GetStudyAndPatientAncestors(
//   vtkMRMLSubjectHierarchyNode* shNode,
//   vtkIdType itemID)
// {
//   std::map<std::string, vtkIdType> result;
//
//   if (shNode==nullptr)
//     return result;
//
//   if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
//     return result;
//
//   vtkIdType currentID = shNode->GetItemParent(itemID);  // skip self
//   while (currentID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
//   {
//     // Check and retrieve "Level" attribute
//     if (shNode->HasItemAttribute(currentID, "Level"))
//     {
//       std :: string levelStr = shNode->GetItemAttribute(currentID, "Level");
//       if (levelStr == "Study" || levelStr == "Patient")
//       {
//         // Insert if not already in map (closer ancestor takes precedence)
//         if (result.find(levelStr) == result.end())
//         {
//           result[levelStr] = currentID;
//         }
//
//         // Optional: stop once both found
//         if (result.size() == 2)
//           break;
//       }
//     }
//
//     currentID = shNode->GetItemParent(currentID);
//   }
//
//   return result;
// }


// void qSlicerDynamicPETModuleWidget::showProgressBar()
// {
//   this->ui->progressBar->setVisible(true);
//   this->ui->progressBar->setValue(0);
// }
//
// void qSlicerDynamicPETModuleWidget::updateProgress(double progress)
// {
//   this->ui->progressBar->setValue(static_cast<int>(progress * 100.0));
//   qApp->processEvents(); // To force GUI update
// }
//
// void qSlicerDynamicPETModuleWidget::hideProgressBar()
// {
//   this->ui->progressBar->setVisible(false);
// }

void qSlicerDynamicPETModuleWidget::enter()
{
  this->IsActive = true;
  this->Superclass::enter();

  // Optional: force refresh when the user enters the module
  this->onSubjectHierarchyChanged();
}

void qSlicerDynamicPETModuleWidget::exit()
{
  this->IsActive = false;
  this->Superclass::exit();
}


void qSlicerDynamicPETModuleWidget::onSubjectHierarchyChanged() {
  if (!this->IsActive)
  {
    return;  // Don't do anything if the module is not active
  }
  Q_D(qSlicerDynamicPETModuleWidget);
  d->populatePatientComboBox();
}

void qSlicerDynamicPETModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  this->Superclass::setMRMLScene(scene);
  Q_D(qSlicerDynamicPETModuleWidget);

  this->qvtkDisconnectAll();

  this->SubjectHierarchyNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (this->SubjectHierarchyNode)
  {
    this->qvtkConnect(this->SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemAddedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
    this->qvtkConnect(this->SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemRemovedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
    this->qvtkConnect(this->SubjectHierarchyNode, vtkMRMLSubjectHierarchyNode::SubjectHierarchyItemModifiedEvent,
                      this, SLOT(onSubjectHierarchyChanged()));
  }
  this->SegWatcher = vtkSmartPointer<SegmentationChangeWatcher>::New();
  this->SegWatcher->GetSequencePET = [this]() { return this->sequencePETNode; };
  this->SegWatcher->GetSequenceSeg = [this]() { return this->segSequenceNode; };
  this->SegWatcher->GetLogic = [this]() { return vtkSlicerDynamicPETLogic::SafeDownCast(this->logic()); };
  this->SegWatcher->GetsegmentTACs = [this]() { return &this->segmentTACs; };
  this->SegWatcher->GetSegEditCorr = [d]() { return d->PlotLiveSegEdit->isChecked(); };
  this->SegWatcher->RunPlot = [this]() { this->onPlotbutton(); };
  this->SegWatcher->GetCurrentSegID = [this]() { return this->SubjectHierarchyNode->GetItemName(this->segID); };
  this->SegWatcher->OnSegmentStructureChanged = [this, d]()
  {
    // Defer UI/cache changes until Slicer finishes its own segment reparent
    // callback. This avoids modifying subject-hierarchy-dependent state from
    // inside vtkSegmentation::SegmentAdded/SegmentRemoved processing.
    QTimer::singleShot(0, this, [this, d]()
    {
      // Segment add/remove changes the structure of every dynamic segmentation
      // frame. Cached TACs and any segment-derived IF are no longer trustworthy.
      this->clearTACdata();
      d->populateSegmentCheckboxes(this->segID);
      d->updateInputFunctionStatus();
      this->enableTACbutton();
    });
  };
  this->SegWatcher->OnSegmentTACChanged = [this, d](const std::string& segmentID)
  {
    // Any edited tissue TAC invalidates ROI fits. If the edited segment is
    // the active IDIF, voxelwise fits and IF previews are stale as well.
    this->clearFITdata();
    this->clearFITMTGAdata();

    if (segmentID == this->IFID)
    {
      d->invalidateInputFunctionResults();
      d->updateInputFunctionStatus();
    }
    else
    {
      this->enableFITbutton();
      this->enableFITMTGAbutton();
    }
  };

}




void qSlicerDynamicPETModuleWidget::onPatChanged (int index) {
  Q_D(qSlicerDynamicPETModuleWidget);
  this->patID = d->PatSelector->itemData(index).value<vtkIdType>();
  d->resetInputFunctionData(true);
  d->populateStudyComboBox(this->patID);
  if (this->patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // "None" selected — ignore or reset state
    d->fileexcel->setText(QString::fromStdString(".xlsx"));
    return;
  }
  std :: string name = this->SubjectHierarchyNode->GetItemName(this->patID);
  std :: string excelfile = name + "_TAC.xlsx";
  d->fileexcel->setText(QString::fromStdString(excelfile));
  std :: string excelfiletcm = name + "_TCMparameters.xlsx";
  d->fileexceltcm->setText(QString::fromStdString(excelfiletcm));
  std :: string excelfiletcmfitted = name + "_TCMfitted.xlsx";
  d->fileexceltcmfitted->setText(QString::fromStdString(excelfiletcmfitted));
  std :: string excelfilemtga = name + "_MTGAparameters.xlsx";
  d->fileexcelmtga->setText(QString::fromStdString(excelfilemtga));
  std :: string excelfilemtgafitted = name + "_MTGAfitted.xlsx";
  d->fileexcelmtgafitted->setText(QString::fromStdString(excelfilemtgafitted));
}


void qSlicerDynamicPETModuleWidget::onStuChanged (int index) {
  Q_D(qSlicerDynamicPETModuleWidget);
  d->resetInputFunctionData(true);
  this->stuID = d->StuSelector->itemData(index).value<vtkIdType>();
  d->populateNodeComboBox(d->CTSelector,
                          this->stuID,
                          "vtkMRMLScalarVolumeNode",
                          "CT"
                          );
  d->populateNodeComboBox(d->PETSelector,
                          this->stuID,
                          "vtkMRMLScalarVolumeNode",
                          "PT"
                          );
  d->populateNodeComboBox(d->SegSelector,
                          this->stuID,
                          "vtkMRMLSegmentationNode",
                          ""
                          );
}


void qSlicerDynamicPETModuleWidget::onCTChanged(int index)
{
    Q_D(qSlicerDynamicPETModuleWidget);

    this->ctID =
        d->CTSelector
            ->itemData(index)
            .value<vtkIdType>();

    d->updateBodySupportUI();
    d->invalidateParametricVoxelSelection();

    this->enableTACbutton();
}

void qSlicerDynamicPETModuleWidget::resetPETSelection()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->resetInputFunctionData(true);
  this->sequencePETNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  this->SegWatcher->browser = nullptr;

  this->durations.clear();
  this->timePoints.clear();
  this->suvFactors.clear();
  d->suvbwFactorValidated = false;

  int noneIndex = d->PETSelector->findData(
    QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));

  if (noneIndex >= 0)
  {
    d->PETSelector->blockSignals(true);
    d->PETSelector->setCurrentIndex(noneIndex);
    d->PETSelector->blockSignals(false);
  }

  this->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
  d->updateQuantitativeUnitUI();
  d->populateNodeComboBox(d->SegSelector,
                          this->stuID,
                          "vtkMRMLSegmentationNode",
                          ""
                          );
  d->invalidateParametricVoxelSelection();
}

void qSlicerDynamicPETModuleWidget::onPETChanged (int index) {
  Q_D(qSlicerDynamicPETModuleWidget);
  d->resetInputFunctionData(true);
  const vtkIdType newPetID =
      d->PETSelector
          ->itemData(index)
          .value<vtkIdType>();

  d->setPETItemID(newPetID);

  this->sequencePETNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  this->SegWatcher->browser = nullptr;

  vtkMRMLScene* scene = this->mrmlScene();
  if (scene==nullptr) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode) {
    return;
  }
  // Fetch PET
  if (this->petID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    this->resetPETSelection();
    return;
  }
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(this->petID));
  if (!petNode) {
    this->resetPETSelection();
    return;
  }

  const char* loadedBy = petNode->GetAttribute("dPETImporter.LoadedBy");
  bool proxyIsValid = (loadedBy && std::string(loadedBy) == "dPETImporterPlugin");

  // Collect the sequence for the dynamic PET
  vtkMRMLSequenceNode* foundSeqNode = nullptr;
  vtkMRMLSequenceBrowserNode* foundBrowser = nullptr;

  for (int i = 0; i < scene->GetNumberOfNodesByClass("vtkMRMLSequenceBrowserNode"); ++i)
  {
    vtkMRMLSequenceBrowserNode* browser =
      vtkMRMLSequenceBrowserNode::SafeDownCast(
        scene->GetNthNodeByClass(i, "vtkMRMLSequenceBrowserNode"));

    if (!browser)
      continue;

    vtkMRMLSequenceNode* seqNode = browser->GetSequenceNode(petNode);
    if (seqNode)
    {
      foundSeqNode = seqNode;
      foundBrowser = browser;

      // If proxy was not valid, check sequence
      if (!proxyIsValid)
      {
        const char* seqLoadedBy =
          seqNode->GetAttribute("dPETImporter.LoadedBy");

        if (!(seqLoadedBy && std::string(seqLoadedBy) == "dPETImporterPlugin"))
        {
          QMessageBox::warning(nullptr,
                               tr("Invalid PET"),
                               tr("Selected PET was not loaded using dPETImporter."));
          this->resetPETSelection();
          return;
        }
      }

      break;
    }
  }

  // If no sequence found → invalid
  if (!foundSeqNode || !foundBrowser)
  {
    QMessageBox::warning(nullptr,
                         tr("Missing node"),
                         tr("Could not find sequence or browser node for PET."));
    this->resetPETSelection();
    return;
  }

  // Assign AFTER validation
  this->sequencePETNode = foundSeqNode;
  this->sequenceBrowserPETNode = foundBrowser;
  this->SegWatcher->browser = foundBrowser;
  this->numberOfTimepoints = this->sequencePETNode->GetNumberOfDataNodes();


  this->durations.clear();
  this->timePoints.clear();
  this->suvFactors.clear();
  d->suvbwFactorValidated = true;

  double cumulativeTime = 0.0;
  std::vector<std::string> valueTypes;

  this->durations.reserve(this->numberOfTimepoints);
  this->timePoints.reserve(this->numberOfTimepoints);
  this->suvFactors.reserve(this->numberOfTimepoints);
  valueTypes.reserve(this->numberOfTimepoints);

  for (int i = 0; i < this->numberOfTimepoints; ++i)
  {
    vtkMRMLNode* frameNode = this->sequencePETNode->GetNthDataNode(i);
    vtkMRMLScalarVolumeNode* volNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(frameNode);

    if (!volNode)
    {
      QMessageBox::critical(nullptr,
                            tr("Invalid PET"),
                            tr("Invalid frame node at index %1").arg(i));
      this->resetPETSelection();
      return;
    }

    const char* durStr = volNode->GetAttribute("dPET.Duration");
    const char* suvStr = volNode->GetAttribute("dPET.SUVbwFactor");
    const char* suvValidStr = volNode->GetAttribute("dPET.SUVbwFactorValid");
    const char* typeStr = volNode->GetAttribute("dPET.ValueType");

    if (!durStr)
    {
      QMessageBox::critical(nullptr,
                            tr("Invalid PET"),
                            tr("Missing dPET.Duration on frame %1").arg(i));
      this->resetPETSelection();
      return;
    }

    double duration = QString(durStr).toDouble();

    cumulativeTime += duration;

    this->durations.push_back(duration);
    this->timePoints.push_back(cumulativeTime);

    // SUVbw conversion factor. New dPETImporter versions explicitly mark
    // whether the factor is physiologically validated. For legacy sequences,
    // a positive non-unity factor is accepted as the old importer only stored
    // a real factor when it had converted BQML to SUVbw.
    const double suvFactor =
        suvStr ? QString(suvStr).toDouble() : 0.0;
    this->suvFactors.push_back(suvFactor);

    bool frameFactorValid =
        suvValidStr && QString(suvValidStr) == "1";

    if (!frameFactorValid &&
        std::isfinite(suvFactor) &&
        suvFactor > 0.0 &&
        std::abs(suvFactor - 1.0) > 1e-12)
    {
        frameFactorValid = true;
    }

    if (!frameFactorValid)
    {
        d->suvbwFactorValidated = false;
    }

    // ValueType (safe)
    if (typeStr)
      valueTypes.emplace_back(typeStr);
    else
    {
      QMessageBox::critical(nullptr,
                            tr("Invalid PET"),
                            tr("Missing dPET.ValueType on frame %1").arg(i));
      this->resetPETSelection();
      return;
    }
  }
  std::set<std::string> uniqueTypes(valueTypes.begin(), valueTypes.end());
  if (uniqueTypes.size() != 1)
  {
    QMessageBox::critical(nullptr,
                          tr("Invalid PET"),
                          tr("Inconsistent dPET.ValueType across frames."));
    this->resetPETSelection();
    return;
  }
  this->dPETvalueType = *uniqueTypes.begin();
  d->updateQuantitativeUnitUI();
  d->populateTimeBarMTGA();
  d->populateTimeBarMTGAImg();

  d->populateNodeComboBox(d->SegSelector,
                          this->stuID,
                          "vtkMRMLSegmentationNode",
                          ""
                          );

  this->enableTACbutton();

  // Get proxy node for current time/frame
  this->sequenceBrowserPETNode->SetSelectedItemNumber(this->numberOfTimepoints - 1);
  vtkMRMLScalarVolumeNode* proxyVolume =
      vtkMRMLScalarVolumeNode::SafeDownCast(this->sequenceBrowserPETNode->GetProxyNode(this->sequencePETNode));
  if (!proxyVolume)
  {
      qCritical() << "Cannot get proxy volume for PET frame";
      return;
  }

  vtkMRMLProceduralColorNode* petLUTNode =
      vtkMRMLProceduralColorNode::SafeDownCast(
          scene->GetFirstNodeByName("PET-DICOM"));
  if (!petLUTNode)
  {
      qCritical() << "Could not find PET-DICOM procedural color node in the scene";
      return;
  }

  // Apply LUT and auto window/level
  vtkMRMLScalarVolumeDisplayNode* displayNode =
      vtkMRMLScalarVolumeDisplayNode::SafeDownCast(proxyVolume->GetDisplayNode());

  if (displayNode)
  {
      displayNode->SetAndObserveColorNodeID(petLUTNode->GetID());
      displayNode->AutoWindowLevelOn();
      displayNode->SetAutoWindowLevel(1);
      displayNode->SetVisibility(true);
  }

  vtkSlicerApplicationLogic* appLogic = qSlicerApplication::application()->applicationLogic();
  if (appLogic)
  {
      vtkMRMLSliceCompositeNode* compositeNode;

      // Iterate over all slice views
      for (int i = 0; i < appLogic->GetMRMLScene()->GetNumberOfNodesByClass("vtkMRMLSliceCompositeNode"); ++i)
      {
          compositeNode = vtkMRMLSliceCompositeNode::SafeDownCast(
              appLogic->GetMRMLScene()->GetNthNodeByClass(i, "vtkMRMLSliceCompositeNode"));
          if (compositeNode)
          {
              compositeNode->SetBackgroundVolumeID(proxyVolume->GetID());
          }
      }
  }
}

void qSlicerDynamicPETModuleWidget::onSegChanged (int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);

  // A segment-derived IDIF belongs to the current
  // segmentation/TAC set. An external CSV does not.
  if (d->IFSourceSelector->currentIndex() == 0)
  {
    d->resetInputFunctionData(false);
  }

  this->segID =
      d->SegSelector
          ->itemData(index)
          .value<vtkIdType>();

  this->segSequenceNode = nullptr;

  vtkMRMLScene* scene = this->mrmlScene();
  if (scene==nullptr) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode) {
    return;
  }
  // Fetch PET
  if (this->petID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->populateNodeComboBox(d->SegSelector,
                            this->stuID,
                            "vtkMRMLSegmentationNode",
                            ""
                            );
    return;
  }
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(this->petID));
  if (!petNode || this->sequencePETNode == nullptr || this->sequenceBrowserPETNode == nullptr) {
    d->populateNodeComboBox(d->SegSelector,
                            this->stuID,
                            "vtkMRMLSegmentationNode",
                            ""
                            );
    return;
  }
  // Fetch Segmentation
  if (this->segID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->populateSegmentCheckboxes(this->segID);
    return;
  }
  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(this->segID));
  if (!segNode) {
    d->populateSegmentCheckboxes(this->segID);
    return;
  }
  // Make sure the source of the segmentation is a binary label map, alongside a created closed surface
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    return;
  }

  this->ProgressBar->setMinimum(0);
  this->ProgressBar->setMaximum(100);
  this->ProgressBar->setValue(0);
  this->ProgressBar->setFormat("Preparing segmentation...");
  this->ProgressBar->setVisible(true);
  this->ProgressBar->show();

  this->stopButton->setVisible(false);
  qApp->processEvents();

  this->ProgressBar->setValue(5);
  this->ProgressBar->setFormat(
    "Preparing segmentation representations...");

  qApp->processEvents();

  logic->setupSeg(segNode);

  this->ProgressBar->setValue(80);
  this->ProgressBar->setFormat(
    "Preparing frame-by-frame segmentation...");

  qApp->processEvents();

  vtkMRMLSequenceNode* seqNode =
  this->sequenceBrowserPETNode->GetSequenceNode(segNode);
  if (!seqNode)
  {
    vtkSmartPointer<vtkMRMLSequenceNode> newSeqNode =
      vtkSmartPointer<vtkMRMLSequenceNode>::New();

    newSeqNode->SetName(
      shNode->GetItemName(segID).c_str());

    scene->AddNode(newSeqNode);

    std::cout
        << "[DynamicPET] About to synchronize segmentation sequence\n"
        << "  PET index name = '" << this->sequencePETNode->GetIndexName() << "'\n"
        << "  PET index unit = '" << this->sequencePETNode->GetIndexUnit() << "'\n"
        << "  PET index type = " << this->sequencePETNode->GetIndexType() << "\n"
        << "  SEG index name = '" << newSeqNode->GetIndexName() << "'\n"
        << "  SEG index unit = '" << newSeqNode->GetIndexUnit() << "'\n"
        << "  SEG index type = " << newSeqNode->GetIndexType()
        << std::endl;

    newSeqNode->CopySequenceIndex(
        this->sequencePETNode);

    std::string indexValue;

    for (int i = 0; i < this->numberOfTimepoints; ++i)
    {
      indexValue =
        this->sequencePETNode->GetNthIndexValue(i);

      if (!newSeqNode->GetDataNodeAtValue(indexValue))
      {
        newSeqNode->SetDataNodeAtValue(
          segNode,
          indexValue);
      }

      // Remaining 20% belongs to sequence creation.
      const int progress =
        80 +
        static_cast<int>(
          20.0 *
          static_cast<double>(i + 1) /
          static_cast<double>(this->numberOfTimepoints));

      this->ProgressBar->setValue(progress);

      this->ProgressBar->setFormat(
        QString("Preparing segmentation frames %1/%2 (%p%)")
          .arg(i + 1)
          .arg(this->numberOfTimepoints));

      // Sequence creation consists of many individual operations,
      // so here the progress bar can genuinely update.
      qApp->processEvents();
    }

    // Add the proxy only after the sequence has a data-node class.
    this->sequenceBrowserPETNode->AddProxyNode(
      segNode,
      newSeqNode,
      false);

    this->sequenceBrowserPETNode->SetSaveChanges(
      newSeqNode,
      true);

    this->SegWatcher->ObserveSegmentationNode(segNode);

    seqNode = newSeqNode;
  }

  this->segSequenceNode = seqNode;

  d->populateSegmentCheckboxes(this->segID);
  this->enableTACbutton();

  // Finished
  this->ProgressBar->setValue(100);
  this->ProgressBar->setFormat("Segmentation ready");
  qApp->processEvents(QEventLoop::ExcludeUserInputEvents);

  // Hide and reset
  this->ProgressBar->hide();
  this->ProgressBar->setValue(0);
  this->ProgressBar->setFormat("%p%");

  qApp->processEvents(QEventLoop::ExcludeUserInputEvents);

}

void qSlicerDynamicPETModuleWidget::onSegmentsChanged()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  // std::vector<QString> selectedSegmentIDs;
  //
  // for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
  // {
  //   QLayoutItem* item = d->segmentCheckLayout->itemAt(i);
  //   QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
  //   if (checkbox && checkbox->isChecked())
  //   {
  //     QString segmentID = checkbox->property("SegmentID").toString();
  //     selectedSegmentIDs.push_back(segmentID);
  //   }
  // }
  // this->segmentIDs = selectedSegmentIDs;
  d->populateSegmentCheckboxes(this->segID);
}

void qSlicerDynamicPETModuleWidget::clearTACdata()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  const bool segmentInputSource =
      d->IFSourceSelector->currentIndex() == 0;

  // Disable all workflows that depend on TAC computation.
  d->setPostTACEnabled(false);

  // Segment-derived TAC data.
  //
  // PET_flatten_values is PET data, not TAC data.
  // It is cleared when the selected PET changes.
  this->segmentTACs.clear();
  this->segmentTACsnames.clear();

  // Current TAC-dependent selections
  this->IFID.clear();
  this->VOIsegmentIDs.clear();
  this->VOIMTGAsegmentIDs.clear();
  this->plotTCMVOI.clear();
  this->plotMTGAVOI.clear();

  // ROI modeling results
  this->segmentTCM.clear();
  this->segmentMTGA.clear();

  // TCM plotting/fitted data
  this->segmentTAC4TCMfits.clear();
  this->segmentkeep4TCMfits.clear();
  this->segmentTCMfits.clear();

  // Parametric imaging depends on extracted segment TACs
  // only when the input function is segment-derived.
  //
  // With an external CSV IF, changing/clearing the
  // segmentation must not invalidate voxelwise PET fits.
  if (segmentInputSource)
  {
    d->resetParametricImagingSelections();
  }

  // Plot nodes
  this->RemoveExistingPlotChartAndTable();

  // Refresh TAC-dependent UI
  d->populatePlotSegmentCheckboxes();
  d->populateIF();

  d->TACCollapsibleButton->setCollapsed(true);

  for (int i = 0; i < d->PlotStatsCheckLayout->count(); ++i)
  {
    QWidget* widget =
        d->PlotStatsCheckLayout->itemAt(i)->widget();

    QCheckBox* cb =
        qobject_cast<QCheckBox*>(widget);

    if (cb)
    {
      cb->setChecked(false);
    }
  }

  d->direxcel->setCurrentPath(QString());
  d->saveExcelButton->setEnabled(false);
}

void qSlicerDynamicPETModuleWidget::enableTACbutton() {
  Q_D(qSlicerDynamicPETModuleWidget);
  if (this->petID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->segID==vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->segmentIDs.empty()) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->sequencePETNode==nullptr || this->sequenceBrowserPETNode==nullptr) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  if (this->segSequenceNode==nullptr) {
    d->TACbutton->setEnabled(false);
    this->clearTACdata();
    return;
  }
  d->TACbutton->setEnabled(true);
}


vtkMRMLTableNode* qSlicerDynamicPETModuleWidget::GetOrCreatePlotTable()
{
  vtkMRMLTableNode* tableNode = vtkMRMLTableNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("DynamicPET.PlotTable"));
  if (!tableNode)
  {
    tableNode = vtkMRMLTableNode::New();
    tableNode->SetName("DynamicPET.PlotTable");
    this->mrmlScene()->AddNode(tableNode);
    tableNode->Delete();
  }
  else
  {
    tableNode->RemoveAllColumns();
  }
  return tableNode;
}

vtkMRMLPlotChartNode* qSlicerDynamicPETModuleWidget::GetOrCreatePlotChart()
{
  vtkMRMLPlotChartNode* chartNode = vtkMRMLPlotChartNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("DynamicPET.PlotChart"));
  if (!chartNode)
  {
    chartNode = vtkMRMLPlotChartNode::New();
    chartNode->SetName("DynamicPET.PlotChart");
    this->mrmlScene()->AddNode(chartNode);
    chartNode->Delete();
  }
  else
  {
    chartNode->RemoveAllPlotSeriesNodeIDs();
  }
  return chartNode;
}

void qSlicerDynamicPETModuleWidget::RemoveExistingPlotChartAndTable()
{
  this->MapPlotSeriesNodeIDToPlot.clear();
  this->ColNameToSegmentID.clear();
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI.clear();
  this->lastSelection.clear();

  vtkMRMLPlotChartNode* chartNode = vtkMRMLPlotChartNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("DynamicPET.PlotChart"));
  if (chartNode)
  {
    vtkCollection* viewNodes = this->mrmlScene()->GetNodesByClass("vtkMRMLPlotViewNode");
    if (viewNodes)
    {
      for (int i = 0; i < viewNodes->GetNumberOfItems(); ++i)
      {
        vtkMRMLPlotViewNode* viewNode = vtkMRMLPlotViewNode::SafeDownCast(
          viewNodes->GetItemAsObject(i));
        if (viewNode && viewNode->GetPlotChartNodeID() &&
            std::string(viewNode->GetPlotChartNodeID()) == chartNode->GetID())
        {
          viewNode->SetPlotChartNodeID(nullptr);
        }
      }
      viewNodes->Delete();
    }

    std::vector<std::string> seriesIDs;
    chartNode->GetPlotSeriesNodeIDs(seriesIDs);
    this->mrmlScene()->RemoveNode(chartNode);
    for (const std::string& id : seriesIDs)
    {
      vtkMRMLNode* node = this->mrmlScene()->GetNodeByID(id);
      if (node)
        this->mrmlScene()->RemoveNode(node);
    }
  }

  vtkMRMLTableNode* tableNode = vtkMRMLTableNode::SafeDownCast(
    this->mrmlScene()->GetFirstNodeByName("DynamicPET.PlotTable"));
  if (tableNode)
    this->mrmlScene()->RemoveNode(tableNode);
}

void qSlicerDynamicPETModuleWidget::onTACbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  vtkMRMLScene* scene = this->mrmlScene();
  if (!scene) {
    return;
  }
  // Run TAC computation
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    return;
  }

  if (this->durations.empty() || this->timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  std::vector<QString> segmentsToCompute;

  for (const QString& segmentID_qt : this->segmentIDs)
  {
    std::string segmentID = segmentID_qt.toStdString();
    if (this->segmentTACsnames.find(segmentID) == this->segmentTACsnames.end())
    {
      segmentsToCompute.push_back(segmentID_qt);
    }
    // else if (this->segmentHasChanged(segmentID))  // Implement this!
    // {
    //   // Segment changed — needs recomputing
    //   segmentsToCompute.push_back(segmentID_qt);
    // }
  }

  this->ProgressBar->setVisible(true);
  this->ProgressBar->show();
  // logic->computeTAC(this->ctID, this->petID, this->segID, segmentsToCompute, this->segmentTACs, this->segmentTACsnames, this->ProgressBar);
  this->stopRequested = false;
  if (this->PET_flatten_values.empty()) {
    logic->Image2Flatten(this->petID, this->PET_flatten_values, this->PETdims, this->numberOfTimepoints, this->ProgressBar, this->stopButton, this->stopRequested);
  }
  if (this->stopRequested) {
    return;
  }
  logic->TAC(this->sequencePETNode, this->segSequenceNode, segmentsToCompute, this->segmentTACs, this->segmentTACsnames, this->ProgressBar, this->stopButton, this->stopRequested);
  this->ProgressBar->setValue(0);
  this->ProgressBar->setVisible(false);
  if (this->stopRequested) {
    return;
  }
  // qApp->processEvents();
  d->populatePlotSegmentCheckboxes();
  d->populateIF();
  d->populateTimeBarMTGA();
  d->populateTimeBarMTGAImg();

  const bool tacReady =
      !this->segmentTACs.empty() &&
      !this->segmentTACsnames.empty();

  d->setPostTACEnabled(tacReady);
  d->updateParametricImagingAvailability();
  return;
}

void qSlicerDynamicPETModuleWidget::onSelectAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->SegmentCheckContents->blockSignals(true);
  if (this->segmentIDs.size()==(d->segmentCheckLayout->count()-1)) {
    for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
    {
      QWidget* widget = d->segmentCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->segmentCheckLayout->count(); ++i)
    {
      QWidget* widget = d->segmentCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->SegmentCheckContents->blockSignals(false);
  d->populateSegmentCheckboxes(this->segID);
}

void qSlicerDynamicPETModuleWidget::onExcelPathChanged(const QString& path)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->saveExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerDynamicPETModuleWidget::onExcelTCMPathChanged(const QString& path)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->saveTCMExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerDynamicPETModuleWidget::onExcelMTGAPathChanged(const QString& path)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->saveMTGAExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerDynamicPETModuleWidget::onExcelTCMfittedPathChanged(const QString& path)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->saveTCMfittedExcelButton->setEnabled(!path.trimmed().isEmpty());
}

void qSlicerDynamicPETModuleWidget::onExcelMTGAfittedPathChanged(const QString& path)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->saveMTGAfittedExcelButton->setEnabled(!path.trimmed().isEmpty());
}

QVariantMap qSlicerDynamicPETModuleWidget::TACtoPythonDict()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QVariantMap out;

  for (const auto& [segmentName, statsVec] : this->segmentTACs)
  {
    QVariantList voxelStatsList;
    const size_t N = statsVec.size();
    if (N != this->durations.size() || N != this->timePoints.size())
    {
      std::cerr << "Mismatch in vector sizes: statsVec (" << N
                << "), durations (" << durations.size()
                << "), timePoints (" << timePoints.size() << ")" << std::endl;
    }
    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        const auto& vs = statsVec[i];
        QVariantMap vsMap;
        // PET measurements are frame averages.  Keep the familiar frame-end time
        // in the primary Time(s) column while preserving exact start/mid/end bounds.
        const double frameEndSec = timePoints[i];
        const double frameStartSec = frameEndSec - durations[i];
        const double frameMidSec =
            frameStartSec + 0.5 * durations[i];

        vsMap["Time(s)"] = frameEndSec;
        vsMap["FrameStart_s"] = frameStartSec;
        vsMap["FrameMid_s"] = frameMidSec;
        vsMap["FrameEnd_s"] = frameEndSec;
        vsMap["Duration_s"] = durations[i];
        // Add stats
        vsMap["VoxelCount"] = vs.count;
        vsMap["Mean"] = vs.mean;
        vsMap["Median"] = vs.median;
        vsMap["Min"] = vs.min;
        vsMap["Max"] = vs.max;
        vsMap["StDev"] = vs.stddev;
        vsMap["Q1"] = vs.q1;
        vsMap["Q3"] = vs.q3;
        vsMap["IQR"] = vs.iqr;
        vsMap["Peak"] = vs.peak;
        vsMap["PeakStDev"] = vs.peakStddev;
        vsMap["PeakVoxelCount"] = vs.peakCount;
        vsMap["Volume(mm3)"] = vs.volume_mm3;
        vsMap["Volume(cm3)"] = vs.volume_cm3;

        voxelStatsList.append(vsMap);
      }
    }

    std :: string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }
    out[QString::fromStdString(sheetName)] = voxelStatsList;
  }

  return out;
}

QVariantMap qSlicerDynamicPETModuleWidget::TCMParamsToPythonDict()
{
  QVariantMap out;

  for (const auto& [segmentName, modelParamsMap] : this->segmentTCM)
  {
    QVariantList rows;

    for (const auto& [modelName, params] : modelParamsMap)
    {
      QVariantMap row;
      row["Model"] = QString::fromStdString(modelName);
      row["K1"]    = params.K1;
      row["k2"]    = params.k2;
      row["k3"]    = params.k3;
      row["k4"]    = params.k4;
      row["ka"]    = params.ka;
      row["fA"]    = params.fa;
      row["vb"]    = params.vb;
      row["td"]    = params.td;
      row["Ki"]    = params.Ki;
      row["DV"]    = params.DV;
      row["AIC"]   = params.AIC;
      row["BIC"]   = params.BIC;
      row["MASE"]  = params.MASE;
      row["chi^2_nu"]  = params.chi2;
      row["BoundHits"] = formatTCMBoundStatus(params.boundFlags);

      rows.append(row);
    }

    // Ensure Excel-friendly sheet name
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rows;
  }

  return out;
}

QVariantMap qSlicerDynamicPETModuleWidget::MTGAParamsToPythonDict()
{
  QVariantMap out;

  for (const auto& [segmentName, modelParamsMap] : this->segmentMTGA)
  {
    QVariantList rows;

    for (const auto& [modelName, params] : modelParamsMap)
    {
      QVariantMap row;
      row["Model"] = QString::fromStdString(modelName);
      row["Ki"]    = params.Ki;
      row["DV"]    = params.DV;
      row["Intercept"] = params.Intercept;
      row["R2"]   = params.R2;
      row["AIC"]   = params.AIC;
      row["MASE"]  = params.MASE;
      // row["chi^2_nu"]  = params.chi2;

      rows.append(row);
    }

    // Ensure Excel-friendly sheet name
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rows;
  }

  return out;
}

QVariantMap qSlicerDynamicPETModuleWidget::fittedTCMtoPythonDict()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QVariantMap out;

  int statIDQString = d->StatSelector->currentIndex();
  std::string currentSelectedStatID = d->StatSelector->itemData(statIDQString).toString().toStdString();

  for (const auto& [segmentName, tacvoi] : this->segmentTAC4TCMfits)
  {
    QVariantList rowList; // Will hold rows for this VOI
    const size_t N = this->timePoints.size();
    auto fitMapIt = this->segmentTCMfits.find(segmentName);
    if (fitMapIt == this->segmentTCMfits.end())
      continue;
    const auto& fits = fitMapIt->second;
    for (size_t i = 0; i < N; ++i)
    {
      QVariantMap row;
      const double frameEndSec = this->timePoints[i];
      const double frameStartSec =
          frameEndSec - this->durations[i];
      const double frameMidSec =
          frameStartSec + 0.5 * this->durations[i];
      row["Time(s)"] = frameEndSec;
      row["FrameStart_s"] = frameStartSec;
      row["FrameMid_s"] = frameMidSec;
      row["FrameEnd_s"] = frameEndSec;
      row["Duration_s"] = this->durations[i];

      // Add TAC VOI
      if (!tacvoi.empty() && tacvoi.size() == N)
      {
        row[QString::fromStdString("TAC (" + currentSelectedStatID + ")")] = tacvoi[i][0];
      }

      // Add each TCM fit value for this time point
      for (const auto& [modelName, fitPtr] : fits)
      {
        if (fitPtr != nullptr)
        {
          row[QString::fromStdString(modelName)] = fitPtr[i];
        }
        else
        {
          row[QString::fromStdString(modelName)] = QVariant(); // blank cell
        }
      }
      rowList.append(row);
    }

    // Shorten sheet name for Excel if needed
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rowList;
  }

  return out;
}

QVariantMap qSlicerDynamicPETModuleWidget::fittedMTGAtoPythonDict()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QVariantMap out;

  int statIDQString = d->StatSelectorMTGA->currentIndex();
  std::string currentSelectedStatID = d->StatSelectorMTGA->itemData(statIDQString).toString().toStdString();

  for (const auto& [segmentName, modelParamsMap] : this->segmentMTGA)
  {
    QVariantList rowList;

    if (!modelParamsMap.empty())
    {
        const auto& firstParams = modelParamsMap.begin()->second;
        const size_t N = firstParams.x.size();

        for (size_t i = 0; i < N; ++i)
        {
          QVariantMap row;

          for (const auto& [modelName, params] : modelParamsMap)
          {

            row[QString::fromStdString(modelName + "_x")] = params.x[i];
            row[QString::fromStdString(modelName + "_y")] = params.y[i];
            row[QString::fromStdString(modelName + "_fitted")] = params.fitted[i];

          }
          rowList.append(row);
        }
    }

    // Shorten sheet name for Excel if needed
    std::string sheetName = this->segmentTACsnames[segmentName];
    if (sheetName.length() > 30)
    {
      sheetName = sheetName.substr(0, 30);
    }

    out[QString::fromStdString(sheetName)] = rowList;

  }

  return out;
}

void qSlicerDynamicPETModuleWidget::onSaveExcelbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QString path = d->direxcel->currentPath();
  QString filename = d->fileexcel->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->TACtoPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_save_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerDynamicPETModuleWidget::onSaveTCMExcelbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QString path = d->direxceltcm->currentPath();
  QString filename = d->fileexceltcm->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->TCMParamsToPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_saveTCM_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerDynamicPETModuleWidget::onSaveMTGAExcelbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QString path = d->direxcelmtga->currentPath();
  QString filename = d->fileexcelmtga->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->MTGAParamsToPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_saveMTGA_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerDynamicPETModuleWidget::onSaveTCMfittedExcelbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QString path = d->direxceltcmfitted->currentPath();
  QString filename = d->fileexceltcmfitted->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->fittedTCMtoPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_generic_save_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerDynamicPETModuleWidget::onSaveMTGAfittedExcelbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  QString path = d->direxcelmtgafitted->currentPath();
  QString filename = d->fileexcelmtgafitted->text();
  QString fullPath = QDir(path).filePath(filename);

  QVariantMap segmentDict = this->fittedMTGAtoPythonDict();
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  PythonQtObjectPtr result = mainContext.call("DPE_genericMTGA_save_multisheet_excel", QVariantList{ fullPath, segmentDict });
}

void qSlicerDynamicPETModuleWidget::onSelectedPoint(vtkStringArray* mrmlPlotSeriesIDs, vtkCollection* selectionCol)
{
  if (!this->checkdisplayedDynamicPET())
    return;

  vtkMRMLScene* scene = this->mrmlScene();
  if (scene==nullptr)
    return;

  if (!mrmlPlotSeriesIDs || !selectionCol || this->MapPlotSeriesNodeIDToPlot.empty())
    return;

  QSet<QPair<QString, vtkIdType>> newSelection;

  // Loop over each series that has selected points
  int psf_value = -1;
  std :: string psv_value = "";
  std :: string lastseriesID = "";
  for (vtkIdType i = 0; i < mrmlPlotSeriesIDs->GetNumberOfValues(); ++i)
  {
    QString seriesID = QString::fromStdString(mrmlPlotSeriesIDs->GetValue(i));
    vtkIdTypeArray* selectedPoints = vtkIdTypeArray::SafeDownCast(
        selectionCol->GetItemAsObject(i));

    if (!selectedPoints)
      continue;

    vtkIdType pointIndex = selectedPoints->GetValue(selectedPoints->GetNumberOfTuples()-1);
    QPair<QString, vtkIdType> candidate(seriesID, pointIndex);
    vtkMRMLPlotSeriesNode* seriesNode = vtkMRMLPlotSeriesNode::SafeDownCast(
        scene->GetNodeByID(seriesID.toStdString()));
    if (!seriesNode)
        return;
    vtkMRMLTableNode* tableNode = seriesNode->GetTableNode();
    if (!tableNode)
        return;
    vtkTable* table = tableNode->GetTable();
    if (!table)
        return;
    std::string labelColName = seriesNode->GetLabelColumnName();
    vtkAbstractArray* labelArray = table->GetColumnByName(labelColName.c_str());
    if (!labelArray)
        return;

    vtkStringArray* strArray = vtkStringArray::SafeDownCast(labelArray);
    if (!strArray)
        return;

    QString labelValue = QString::fromStdString(strArray->GetValue(pointIndex));
    QStringList parts = labelValue.split(',');
    if (!parts.isEmpty())
    {
      QString framePart = parts[0].trimmed();
      framePart.remove("Frame:");
      psf_value = framePart.trimmed().toInt();
      psv_value = this->ColNameToSegmentID[seriesNode->GetName()];
    }
    if (!newSelection.contains(candidate))
      newSelection.insert(candidate);
    if (!this->lastSelection.contains(candidate))
    {
      this->PlotSelectedFrame = psf_value;
      this->PlotSelectedVOI = psv_value;
    } else {
      if (this->PlotSelectedFrame == psf_value && this->PlotSelectedVOI == psv_value) {
        lastseriesID = seriesID.toStdString();
        continue;
      }
      if (!this->MapPlotSeriesNodeIDToPlot.contains(seriesID)) {
        vtkGenericWarningMacro("MapPlotSeriesNodeIDToPlot does not contain seriesID=" << seriesID.toStdString());
      }
      vtkPlot* vtkplot = this->MapPlotSeriesNodeIDToPlot.value(seriesID);
      vtkSmartPointer<vtkIdTypeArray> emptySelection = vtkSmartPointer<vtkIdTypeArray>::New();
      vtkplot->SetSelection(emptySelection);
    }
  }
  if (newSelection.size()==0) {
    this->PlotSelectedFrame = -1;
    this->PlotSelectedVOI = "";
  } else if (newSelection.size()==1 && lastseriesID!="") {
  } else {
    if (!this->MapPlotSeriesNodeIDToPlot.contains(QString::fromStdString(lastseriesID))) {
        vtkGenericWarningMacro("MapPlotSeriesNodeIDToPlot does not contain seriesID=" << lastseriesID);
    }
    vtkPlot* vtkplot = this->MapPlotSeriesNodeIDToPlot.value(QString::fromStdString(lastseriesID));
    vtkSmartPointer<vtkIdTypeArray> emptySelection = vtkSmartPointer<vtkIdTypeArray>::New();
    vtkplot->SetSelection(emptySelection);
  }
  this->lastSelection = newSelection;
}

void qSlicerDynamicPETModuleWidget::onPlotbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  this->ColNameToSegmentID.clear();
  this->MapPlotSeriesNodeIDToPlot.clear();
  vtkMRMLScene* scene = this->mrmlScene();
  // Get selected segments
  std::vector<std::string> PlotSelectedIDs;
  for (int i = 0; i < d->PlotsegmentCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->PlotsegmentCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      std::string segmentid = checkbox->property("SegmentID").toString().toStdString();
      PlotSelectedIDs.push_back(segmentid);
    }
  }

  // Get selected stats
  std::vector<std::string> PlotSelectedStats;
  for (int i = 0; i < d->PlotStatsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->PlotStatsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      PlotSelectedStats.push_back(checkbox->text().toStdString());
    }
  }

  if (PlotSelectedIDs.empty() || PlotSelectedStats.empty())
    return;

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  // Add time column
  vtkNew<vtkDoubleArray> timeArray;
  timeArray->SetName("Time (min)");
  vtkNew<vtkStringArray> labelArray;
  labelArray->SetName("ToolTipLabelTAC");
  for (int i = 0; i < this->timePoints.size(); ++i)
  {
    const double frameEndSec = this->timePoints[i];
    const double frameStartSec =
        frameEndSec - this->durations[i];
    const double frameMidSec =
        frameStartSec + 0.5 * this->durations[i];

    timeArray->InsertNextValue(frameEndSec / 60.0);

    std::ostringstream oss;
    oss << "Frame: " << i
        << ", start(s): " << frameStartSec
        << ", midpoint(s): " << frameMidSec
        << ", end(s): " << frameEndSec
        << ", duration(s): " << this->durations[i];
    labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(timeArray);
  tableNode->AddColumn(labelArray);

  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();
  chartNode->SetTitle("Time Activity Curve");
  chartNode->SetXAxisTitle("Time (min)");
  const ActivityUnit displayUnit = d->selectedDisplayActivityUnit();
  chartNode->SetYAxisTitle(
      d->activityUnitLabel(displayUnit).toStdString().c_str());

  auto isActivityStatistic = [](const std::string& name)
  {
    return name == "Mean" || name == "Median" || name == "Peak" ||
           name == "Min" || name == "Max";
  };

  auto convertPlotActivity = [&](double nativeValue, const std::string& statName)
  {
    if (!isActivityStatistic(statName) || !std::isfinite(nativeValue))
    {
      return nativeValue;
    }

    double converted = nativeValue;
    if (!d->convertActivityValue(
            nativeValue,
            d->petStoredActivityUnit(),
            displayUnit,
            converted,
            nullptr))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return converted;
  };

  std::unordered_map<std::string, std::string> LabelToSeriesID;
  for (const std::string& segmentID : PlotSelectedIDs)
  {
    std :: string segmentName = this->segmentTACsnames[segmentID];
    for (const std::string& statName : PlotSelectedStats)
    {
      std::string colName = segmentName + " - " + statName;
      this->ColNameToSegmentID[colName] = segmentID;
      vtkNew<vtkDoubleArray> statArray;
      statArray->SetName(colName.c_str());

      vtkNew<vtkDoubleArray> statErrArray;
      std::string colErrName = colName + " Error";
      statErrArray->SetName(colErrName.c_str());

      vtkNew<vtkDoubleArray> statArrayLine;
      std::string statArrayLineName = colName + " Line";
      statArrayLine->SetName(statArrayLineName.c_str());

      for (int ivs=0; ivs<this->segmentTACs[segmentID].size(); ++ivs)
      {
        const VoxelStatistics& vs = this->segmentTACs[segmentID][ivs];
        double value = std::numeric_limits<double>::quiet_NaN();
        if (vs.keep) {
          if (statName == "Mean")
          {
            value = vs.mean;
            // if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            //   statErrArray->InsertNextValue(vs.stddev);
          }
          else if (statName == "Median")
          {
            value = vs.median;
            // if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            //   statErrArray->InsertNextValue(vs.iqr);
          }
          else if (statName == "Peak")
          {
            value = vs.peak;
            // if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
            //   statErrArray->InsertNextValue(vs.iqr);
          }
          else if (statName == "VoxelCount") value = vs.count;
          else if (statName == "Min")        value = vs.min;
          else if (statName == "Max")        value = vs.max;
          else if (statName == "Volume(cc)") value = vs.volume_cm3;
          else vtkGenericWarningMacro("Unknown stat name: " << statName);

          value = convertPlotActivity(value, statName);
        }

        statArray->InsertNextValue(value);
        // Line points
        if (!std::isnan(value)) {
          statArrayLine->InsertNextValue(value);
        } else {
          double nextValue = std::numeric_limits<double>::quiet_NaN();
          double x1 = std::numeric_limits<double>::quiet_NaN();
          for (int next_ivs = ivs+1; next_ivs<this->timePoints.size(); ++next_ivs) {
            const VoxelStatistics& vsNext = this->segmentTACs[segmentID][next_ivs];
            if (vsNext.keep)
            {
                x1 = this->timePoints[next_ivs];
                if (statName == "Mean") nextValue = vsNext.mean;
                else if (statName == "Median") nextValue = vsNext.median;
                else if (statName == "Peak") nextValue = vsNext.peak;
                else if (statName == "VoxelCount") nextValue = vsNext.count;
                else if (statName == "Min") nextValue = vsNext.min;
                else if (statName == "Max") nextValue = vsNext.max;
                else if (statName == "Volume(cc)") nextValue = vsNext.volume_cm3;
                nextValue = convertPlotActivity(nextValue, statName);
                break;
            }
          }
          if (std::isnan(nextValue)) {
            statArrayLine->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
            continue;
          }
          double prevValue = std::numeric_limits<double>::quiet_NaN();
          double x0 = std::numeric_limits<double>::quiet_NaN();
          for (int prev_ivs = ivs-1; prev_ivs>=0; --prev_ivs) {
            const VoxelStatistics& vsPrev = this->segmentTACs[segmentID][prev_ivs];
            if (vsPrev.keep)
            {
                x0 = this->timePoints[prev_ivs];
                if (statName == "Mean") prevValue = vsPrev.mean;
                else if (statName == "Median") prevValue = vsPrev.median;
                else if (statName == "Peak") prevValue = vsPrev.peak;
                else if (statName == "VoxelCount") prevValue = vsPrev.count;
                else if (statName == "Min") prevValue = vsPrev.min;
                else if (statName == "Max") prevValue = vsPrev.max;
                else if (statName == "Volume(cc)") prevValue = vsPrev.volume_cm3;
                prevValue = convertPlotActivity(prevValue, statName);
                break;
            }
          }
          if (std::isnan(prevValue)) {
            statArrayLine->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
            continue;
          }

          double x  = this->timePoints[ivs];
          // Proper linear interpolation
          value = prevValue + ((x - x0) / (x1 - x0)) * (nextValue - prevValue);
          statArrayLine->InsertNextValue(value);
        }

      }

      tableNode->AddColumn(statArray);
      tableNode->AddColumn(statArrayLine);
      // if (statErrArray->GetNumberOfTuples() > 0)
      //   tableNode->AddColumn(statErrArray);

      vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(lineSeries);
      lineSeries->SetName("");
      lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
      lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
      lineSeries->SetXColumnName("Time (min)");
      lineSeries->SetYColumnName(statArrayLineName.c_str());
      lineSeries->SetLabelColumnName("ToolTipLabelTAC");
      lineSeries->SetUniqueColor();
      lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
      chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());

      vtkSmartPointer<vtkMRMLPlotSeriesNode> series = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(series);
      series->SetName(colName.c_str());
      series->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
      series->SetAndObserveTableNodeID(tableNode->GetID());
      series->SetXColumnName("Time (min)");
      series->SetYColumnName(colName.c_str());
      series->SetLabelColumnName("ToolTipLabelTAC");
      series->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
      series->SetColor(lineSeries->GetColor());
      chartNode->AddAndObservePlotSeriesNodeID(series->GetID());
      LabelToSeriesID[colName] = series->GetID();
    }
  }

  // Show plot view
  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
    scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode)
  {
    plotViewNode->SetPlotChartNodeID(chartNode->GetID());
    qMRMLPlotWidget* plotWidget = nullptr;
    if (qSlicerApplication::application())
    {
      qSlicerLayoutManager* layoutManager =
          qSlicerApplication::application()->layoutManager();
      qMRMLPlotWidget* plotWidget = nullptr;
      plotWidget = layoutManager->plotWidget(0);
      qMRMLPlotView* plotView = plotWidget->plotView();
      if (plotView)
      {
        QObject::connect(plotView, SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
                         this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)));

        vtkSmartPointer<vtkChartXY> chart = plotView->chart();
        for (int i = 0; i < chart->GetNumberOfPlots(); ++i)
        {
           vtkPlot* plot = chart->GetPlot(i);
           std :: string PlotLabel = plot->GetLabel();
           QString seriesNodeID = QString::fromStdString(LabelToSeriesID[PlotLabel]);
           this->MapPlotSeriesNodeIDToPlot[seriesNodeID] = plot;
        }
      }
    }
  }
}

void
qSlicerDynamicPETModuleWidget::
onIFSelectionChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);

  if (index < 0)
  {
    this->IFID.clear();
  }
  else
  {
    this->IFID =
        d->IFSelector
            ->itemData(index)
            .toString()
            .toStdString();
  }

  // Selecting None must immediately invalidate all IF-dependent QC/results.
  // The source selector itself remains available so a new IDIF can be chosen.
  d->invalidateInputFunctionResults();

  d->populateVOI(this->IFID);
  d->populateVOIMTGA(this->IFID);

  d->updateInputFunctionStatus();

  // Any result using the old IF is now stale.
  this->clearFITdata();
  this->clearFITMTGAdata();

  d->MTGAImgFitSignatures.clear();
  d->TCMImgFitSignatures.clear();

  this->MTGAImgOutcomes.clear();
  this->TCMImgOutcomes.clear();

  this->enableFITbutton();
  this->enableFITMTGAbutton();
  this->enableFITMTGAImgbutton();
  this->enableFITTCMImgbutton();
}

void qSlicerDynamicPETModuleWidget::onVOISelectionChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std :: string segmentID = d->VOISelector->itemData(index).toString().toStdString();
  this->plotTCMVOI=segmentID;
  d->populateResultsTable(segmentID);
  if (segmentID == "")
  {
    return;
  }
}

void qSlicerDynamicPETModuleWidget::onVOIMTGASelectionChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std :: string segmentID = d->VOISelectorMTGA->itemData(index).toString().toStdString();
  this->plotMTGAVOI=segmentID;
  d->populateResultsMTGATable(segmentID);
  if (segmentID == "")
  {
    return;
  }
}

void qSlicerDynamicPETModuleWidget::onVOISelectAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->VOICheckContents->blockSignals(true);
  if (this->VOIsegmentIDs.size()==(d->VOICheckLayout->count()-1)) {
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOICheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOICheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->VOICheckContents->blockSignals(false);
  d->populateVOI(this->IFID);
}

void qSlicerDynamicPETModuleWidget::onVOIMTGASelectAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->VOIMTGACheckContents->blockSignals(true);
  if (this->VOIMTGAsegmentIDs.size()==(d->VOIMTGACheckLayout->count()-1)) {
    for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOIMTGACheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
    {
      QWidget* widget = d->VOIMTGACheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->VOIMTGACheckContents->blockSignals(false);
  d->populateVOIMTGA(this->IFID);
}

void qSlicerDynamicPETModuleWidget::onOLSclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->weightedFitCheckBox->blockSignals(true);
  d->robustFitCheckBox->blockSignals(true);
  d->weightedFitCheckBox->setChecked(false);
  d->robustFitCheckBox->setChecked(false);
  d->weightedFitCheckBox->blockSignals(false);
  d->robustFitCheckBox->blockSignals(false);
  d->robustParamsWidget->setVisible(false);
  if (!d->olsFitCheckBox->isChecked())
    d->olsFitCheckBox->setChecked(true);
}

void qSlicerDynamicPETModuleWidget::onOLSImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->weightedFitCheckBoxImg->blockSignals(true);
  d->robustFitCheckBoxImg->blockSignals(true);
  d->weightedFitCheckBoxImg->setChecked(false);
  d->robustFitCheckBoxImg->setChecked(false);
  d->weightedFitCheckBoxImg->blockSignals(false);
  d->robustFitCheckBoxImg->blockSignals(false);
  d->robustParamsWidgetImg->setVisible(false);
  if (!d->olsFitCheckBoxImg->isChecked())
    d->olsFitCheckBoxImg->setChecked(true);
}

void qSlicerDynamicPETModuleWidget::onWLSclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->olsFitCheckBox->blockSignals(true);
  d->robustFitCheckBox->blockSignals(true);
  d->olsFitCheckBox->setChecked(false);
  d->robustFitCheckBox->setChecked(false);
  d->olsFitCheckBox->blockSignals(false);
  d->robustFitCheckBox->blockSignals(false);
  d->robustParamsWidget->setVisible(false);
  if (!d->weightedFitCheckBox->isChecked())
    d->olsFitCheckBox->setChecked(true);
}

void qSlicerDynamicPETModuleWidget::onWLSImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->olsFitCheckBoxImg->blockSignals(true);
  d->robustFitCheckBoxImg->blockSignals(true);
  d->olsFitCheckBoxImg->setChecked(false);
  d->robustFitCheckBoxImg->setChecked(false);
  d->olsFitCheckBoxImg->blockSignals(false);
  d->robustFitCheckBoxImg->blockSignals(false);
  d->robustParamsWidgetImg->setVisible(false);
  if (!d->weightedFitCheckBoxImg->isChecked())
    d->olsFitCheckBoxImg->setChecked(true);
}

void qSlicerDynamicPETModuleWidget::onRLSclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->olsFitCheckBox->blockSignals(true);
  d->weightedFitCheckBox->blockSignals(true);
  d->olsFitCheckBox->setChecked(false);
  d->weightedFitCheckBox->setChecked(false);
  d->olsFitCheckBox->blockSignals(false);
  d->weightedFitCheckBox->blockSignals(false);
  if (!d->robustFitCheckBox->isChecked()) {
    d->olsFitCheckBox->setChecked(true);
    d->robustParamsWidget->setVisible(false);
  } else {
    d->robustParamsWidget->setVisible(true);
  }
}

void qSlicerDynamicPETModuleWidget::onRLSImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->olsFitCheckBoxImg->blockSignals(true);
  d->weightedFitCheckBoxImg->blockSignals(true);
  d->olsFitCheckBoxImg->setChecked(false);
  d->weightedFitCheckBoxImg->setChecked(false);
  d->olsFitCheckBoxImg->blockSignals(false);
  d->weightedFitCheckBoxImg->blockSignals(false);
  if (!d->robustFitCheckBoxImg->isChecked()) {
    d->olsFitCheckBoxImg->setChecked(true);
    d->robustParamsWidgetImg->setVisible(false);
  } else {
    d->robustParamsWidgetImg->setVisible(true);
  }
}

void qSlicerDynamicPETModuleWidget::onStdFitclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->weightFitCheckBox->blockSignals(true);
  d->weightFitCheckBox->setChecked(false);
  d->weightFitCheckBox->blockSignals(false);
  if (!d->standardFitCheckBox->isChecked()) {
    d->standardFitCheckBox->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onStdFitImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->weightFitCheckBoxImg->blockSignals(true);
  d->weightFitCheckBoxImg->setChecked(false);
  d->weightFitCheckBoxImg->blockSignals(false);
  if (!d->standardFitCheckBoxImg->isChecked()) {
    d->standardFitCheckBoxImg->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onWFitclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->standardFitCheckBox->blockSignals(true);
  d->standardFitCheckBox->setChecked(false);
  d->standardFitCheckBox->blockSignals(false);
  if (!d->weightFitCheckBox->isChecked()) {
    d->standardFitCheckBox->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onWFitImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->standardFitCheckBoxImg->blockSignals(true);
  d->standardFitCheckBoxImg->setChecked(false);
  d->standardFitCheckBoxImg->blockSignals(false);
  if (!d->weightFitCheckBoxImg->isChecked()) {
    d->standardFitCheckBoxImg->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onSliderChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  const double timeSec =
      this->timePoints[index - 1];
  const double timeMin = timeSec / 60.0;

  if (d->timeEndSlider->value() < index)
  {
    d->timeEndSlider->setValue(index);
  }

  d->frameEdit->setText(QString::number(index));
  d->timeSecEdit->setText(QString::number(timeSec, 'f', 2));
  d->timeMinEdit->setText(QString::number(timeMin, 'f', 2));
}

void qSlicerDynamicPETModuleWidget::onSliderImgChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  const double timeSec =
      this->timePoints[index - 1];
  const double timeMin = timeSec / 60.0;

  if (d->timeEndSliderImg->value() < index)
  {
    d->timeEndSliderImg->setValue(index);
  }

  d->frameEditImg->setText(QString::number(index));
  d->timeSecEditImg->setText(QString::number(timeSec, 'f', 2));
  d->timeMinEditImg->setText(QString::number(timeMin, 'f', 2));
}

void qSlicerDynamicPETModuleWidget::runVuong(std::string sel1,
                                       std::string sel2,
                                       std::string segmentID
                                     )
{
  Q_D(qSlicerDynamicPETModuleWidget);
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  if (segmentID.empty()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }

  auto it = this->segmentMTGA.find(segmentID);
  if (it == this->segmentMTGA.end()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }
  auto& modelsForSegment = it->second;

  const int N1 = modelsForSegment[sel1].y.size();
  const int N2 = modelsForSegment[sel2].y.size();
  if (N1 != N2) {
    throw std::runtime_error(
        sel1 + " has not been fitted with the same number of datapoints (" + std::to_string(N1) +
        ") of " + sel2 + " (" + std::to_string(N2) + ")."
    );
  }
  std::vector<double> w1 = modelsForSegment[sel1].weights;
  std::vector<double> w2 = modelsForSegment[sel2].weights;
  // Check they are the same length
  if (w1.size() != w2.size()) {
      throw std::invalid_argument("Weight vectors must have the same length");
  }
  // Compute average weights
  std::vector<double> wgt_avg(w1.size());
  for (size_t i = 0; i < w1.size(); ++i) {
      wgt_avg[i] = 0.5 * (w1[i] + w2[i]);
  }
  const std::vector<double>* wgt = &wgt_avg;
  double p = logic->computeVuongP(modelsForSegment[sel1].r,
                                  modelsForSegment[sel2].r,
                                  wgt,
                                  modelsForSegment[sel1].dof,
                                  modelsForSegment[sel2].dof,
                                  VuongCorrection::BIC,
                                  Tail::TwoSided
                                );
  d->MTGAVuongP->setText(QString::number(p, 'g', 4));
  // d->MTGAVuongP->adjustSize();
  return;
}

void qSlicerDynamicPETModuleWidget::runTCMstat(std::string sel1,
                                         std::string sel2,
                                         std::string segmentID
                                        )
{
  Q_D(qSlicerDynamicPETModuleWidget);
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  auto clearStats = [&]()
  {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
  };

  if (segmentID.empty()) {
    clearStats();
    return;
  }

  auto it = this->segmentTCM.find(segmentID);
  if (it == this->segmentTCM.end()) {
    clearStats();
    return;
  }
  auto& modelsForSegment = it->second;

  if (!modelsForSegment.count(sel1) || !modelsForSegment.count(sel2))
  {
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }

  const TCMParameters& m1 = modelsForSegment.at(sel1);
  const TCMParameters& m2 = modelsForSegment.at(sel2);

  if (m1.weights.size() != m2.weights.size())
  {
    throw std::runtime_error(
        sel1 + " has not been fitted with the same number of datapoints (" +
        std::to_string(m1.weights.size()) + ") of " + sel2 + " (" +
        std::to_string(m2.weights.size()) + ")."
    );
  }

  if (m1.r.size() != m2.r.size())
  {
    throw std::runtime_error(
        sel1 + " and " + sel2 + " have different residual vector sizes."
    );
  }

  ModelComparisonResult res =
      logic->compareModels(sel1, sel2, m1, m2);

  d->TCMLRTP->setText("");
  d->TCMVuongP->setText("");

  if (res.type == "LRT")
  {
    d->TCMLRTP->setText(QString::number(res.p_value, 'g', 4));
  }
  else if (res.type == "Vuong")
  {
    d->TCMVuongP->setText(QString::number(res.p_value, 'g', 4));
  }

  return;
}

void qSlicerDynamicPETModuleWidget::onMTGAModelBox(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);

  std::string selectedVOI = this->plotMTGAVOI;
  if (selectedVOI.empty()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }

  auto it = this->segmentMTGA.find(selectedVOI);
  if (it == this->segmentMTGA.end()) {
    d->MTGAModel1->clear();
    d->MTGAModel2->clear();
    d->MTGAVuongP->setText("");
    return;
  }
  const auto& modelsForSegment = it->second;
  std::string sel1, sel2;
  int idx1 = d->MTGAModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = d->MTGAModel1->itemData(idx1).toString().toStdString();
  int idx2 = d->MTGAModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = d->MTGAModel2->itemData(idx2).toString().toStdString();
  d->populateModelCombo(d->MTGAModel1, sel2, sel1, selectedVOI);
  d->populateModelCombo(d->MTGAModel2, sel1, sel2, selectedVOI);
  if (idx1>0 & idx2>0) {
    this->runVuong(sel1, sel2, selectedVOI);
  } else {
    d->MTGAVuongP->setText("");
  }
  return;
}

void qSlicerDynamicPETModuleWidget::onTCMModelBox(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);

  std::string selectedVOI = this->plotTCMVOI;
  if (selectedVOI.empty()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }

  auto it = this->segmentTCM.find(selectedVOI);
  if (it == this->segmentTCM.end()) {
    d->TCMModel1->clear();
    d->TCMModel2->clear();
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
    return;
  }
  const auto& modelsForSegment = it->second;
  std::string sel1, sel2;
  int idx1 = d->TCMModel1->currentIndex();
  if (idx1 >= 0)
    sel1 = d->TCMModel1->itemData(idx1).toString().toStdString();
  int idx2 = d->TCMModel2->currentIndex();
  if (idx2 >= 0)
    sel2 = d->TCMModel2->itemData(idx2).toString().toStdString();
  d->populateModelComboTCM(d->TCMModel1, sel2, sel1, selectedVOI);
  d->populateModelComboTCM(d->TCMModel2, sel1, sel2, selectedVOI);
  if (idx1>0 & idx2>0) {
    this->runTCMstat(sel1, sel2, selectedVOI);
  } else {
    d->TCMLRTP->setText("");
    d->TCMVuongP->setText("");
  }
  return;
}

void qSlicerDynamicPETModuleWidget::clearFITdata() {
  Q_D(qSlicerDynamicPETModuleWidget);
  this->segmentTCM.clear();
  d->populateResultsVOI();
  d->TCMResultsButton->setCollapsed(true);
  return;

}

void qSlicerDynamicPETModuleWidget::clearFITMTGAdata() {
  Q_D(qSlicerDynamicPETModuleWidget);
  this->segmentMTGA.clear();
  d->populateResultsVOIMTGA();
  d->MTGAResultsButton->setCollapsed(true);
  return;
}



void qSlicerDynamicPETModuleWidget::enableFITbutton()
{
    Q_D(qSlicerDynamicPETModuleWidget);

    QString ifError;

    if (!d->hasValidInputFunction(&ifError, true) ||
        this->VOIsegmentIDs.empty() ||
        this->modelsID.empty())
    {
        d->FITbutton->setEnabled(false);
        this->clearFITdata();
        return;
    }

    d->FITbutton->setEnabled(true);
}

void qSlicerDynamicPETModuleWidget::enableFITTCMImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  if (d->parametricFitRunning)
  {
    d->FITbuttonTCMImg->setEnabled(false);
    return;
  }

  if (this->petID ==
          vtkMRMLSubjectHierarchyNode::
              INVALID_ITEM_ID)
  {
      d->FITbuttonTCMImg->
          setEnabled(false);
      return;
  }

  QString ifError;

  if (!d->hasValidInputFunction(
          &ifError,
          true))
  {
      d->FITbuttonTCMImg->
          setEnabled(false);
      return;
  }

  if (this->modelsTCMImgID.empty())
  {
    d->FITbuttonTCMImg->setEnabled(false);
    return;
  }

  const bool show =
      d->TCMShowInSlicerCheckBoxImg->isChecked();

  const bool save =
      d->TCMSaveDICOMCheckBoxImg->isChecked();

  if (!show && !save)
  {
    d->FITbuttonTCMImg->setEnabled(false);
    return;
  }

  if (save &&
      d->TCMDICOMDirectoryImg
          ->currentPath()
          .trimmed()
          .isEmpty())
  {
    d->FITbuttonTCMImg->setEnabled(false);
    return;
  }

  d->FITbuttonTCMImg->setEnabled(true);
}

void qSlicerDynamicPETModuleWidget::enableFITMTGAbutton()
{
    Q_D(qSlicerDynamicPETModuleWidget);

    QString ifError;

    if (!d->hasValidInputFunction(&ifError) ||
        this->VOIMTGAsegmentIDs.empty() ||
        this->modelsMTGAID.empty())
    {
        d->FITMTGAbutton->setEnabled(false);
        this->clearFITMTGAdata();
        return;
    }

    d->FITMTGAbutton->setEnabled(true);
}

void qSlicerDynamicPETModuleWidget::enableFITMTGAImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  if (d->parametricFitRunning)
  {
    d->FITbuttonMTGAImg->setEnabled(false);
    return;
  }

  if (this->petID ==
          vtkMRMLSubjectHierarchyNode::
              INVALID_ITEM_ID)
  {
      d->FITbuttonMTGAImg->
          setEnabled(false);
      return;
  }

  QString ifError;

  if (!d->hasValidInputFunction(
          &ifError))
  {
      d->FITbuttonMTGAImg->
          setEnabled(false);
      return;
  }

  if (this->modelsMTGAImgID.empty())
  {
    d->FITbuttonMTGAImg->setEnabled(false);
    return;
  }

  const bool show =
      d->MTGAShowInSlicerCheckBoxImg->isChecked();

  const bool save =
      d->MTGASaveDICOMCheckBoxImg->isChecked();

  if (!show && !save)
  {
    d->FITbuttonMTGAImg->setEnabled(false);
    return;
  }

  if (save &&
      d->MTGADICOMDirectoryImg
          ->currentPath()
          .trimmed()
          .isEmpty())
  {
    d->FITbuttonMTGAImg->setEnabled(false);
    return;
  }

  d->FITbuttonMTGAImg->setEnabled(true);
}



void qSlicerDynamicPETModuleWidget::onVOISegmentsChanged()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std::vector<std::string> VOIselectedSegmentIDs;

  for (int i = 0; i < d->VOICheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->VOICheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      std :: string segmentID = checkbox->property("SegmentID").toString().toStdString();
      VOIselectedSegmentIDs.push_back(segmentID);
    }
  }
  this->VOIsegmentIDs = VOIselectedSegmentIDs;

  this->enableFITbutton();
}

void qSlicerDynamicPETModuleWidget::onVOIMTGASegmentsChanged()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std::vector<std::string> VOIselectedSegmentIDs;

  for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->VOIMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      std :: string segmentID = checkbox->property("SegmentID").toString().toStdString();
      VOIselectedSegmentIDs.push_back(segmentID);
    }
  }
  this->VOIMTGAsegmentIDs = VOIselectedSegmentIDs;

  this->enableFITMTGAbutton();
}

std::vector<double> qSlicerDynamicPETModuleWidget::extractColumn(const std::vector<std::vector<double>>& mat, const int index)
{
    std::vector<double> col;
    col.reserve(mat.size());
    for (const auto& row : mat)
    {
        col.push_back(row[index]); // assumes at least one column
    }
    return col;
}

void qSlicerDynamicPETModuleWidget::onFITbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  if (segmentTACsnames.empty() || segmentTACs.empty()) {
    std::cerr << "Missing TACs!" << std::endl;
    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  int statIDQString = d->StatSelector->currentIndex();
  if (statIDQString<0) {
    std::cerr << "Missing stat choice!" << std::endl;
    return;
  }
  std::string currentSelectedStatID = d->StatSelector->itemData(statIDQString).toString().toStdString();


  QString ifError;

  if (!d->hasValidInputFunction(
          &ifError))
  {
      QMessageBox::warning(
          this,
          tr("Input Function"),
          ifError);

      return;
  }

  if (VOIsegmentIDs.empty()) {
    std::cerr << "Missing VOIs to fit!" << std::endl;
    return;
  }

  if (modelsID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  // Collect parameters using a lambda for brevity
  auto getParamTriplet = [&](QLineEdit* init, QLineEdit* lb, QLineEdit* ub) {
    return std::tuple<double, double, double>{
      init->text().toDouble(), lb->text().toDouble(), ub->text().toDouble()
    };
  };

  double k1Init, k1Lower, k1Upper;
  std::tie(k1Init, k1Lower, k1Upper) =
      getParamTriplet(
          d->k1Initial,
          d->k1Lower,
          d->k1Upper);

  double k2Init, k2Lower, k2Upper;
  std::tie(k2Init, k2Lower, k2Upper) =
      getParamTriplet(
          d->k2Initial,
          d->k2Lower,
          d->k2Upper);

  double k3Init, k3Lower, k3Upper;
  std::tie(k3Init, k3Lower, k3Upper) =
      getParamTriplet(
          d->k3Initial,
          d->k3Lower,
          d->k3Upper);

  double k4Init, k4Lower, k4Upper;
  std::tie(k4Init, k4Lower, k4Upper) =
      getParamTriplet(
          d->k4Initial,
          d->k4Lower,
          d->k4Upper);

  double vbInit, vbLower, vbUpper;
  std::tie(vbInit, vbLower, vbUpper) =
      getParamTriplet(
          d->vbInitial,
          d->vbLower,
          d->vbUpper);

  double tdInit, tdLower, tdUpper;
  std::tie(tdInit, tdLower, tdUpper) =
      getParamTriplet(
          d->tdInitial,
          d->tdLower,
          d->tdUpper);

  double liverKaInit, liverKaLower, liverKaUpper;
  std::tie(
      liverKaInit,
      liverKaLower,
      liverKaUpper) =
      getParamTriplet(
          d->liverKaInitial,
          d->liverKaLower,
          d->liverKaUpper);

  double liverFaInit, liverFaLower, liverFaUpper;
  std::tie(
      liverFaInit,
      liverFaLower,
      liverFaUpper) =
      getParamTriplet(
          d->liverFaInitial,
          d->liverFaLower,
          d->liverFaUpper);

  const long Nframe = timePoints.size();
  const long Nvox = 1;

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }

  const std::vector<double>* wgt = nullptr;
  std::map< std::string, std::vector<double>> wgtVec;

  std::map<std::string, std::vector<std::vector<double>>> tac;
  std::map<std::string, std::vector<bool>> keeptacvec;
  for (const auto& [segmentName, statsVec] : segmentTACs)
  {
    if (statsVec.size() != static_cast<size_t>(Nframe))
    {
      std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
      return;
    }

    tac[segmentName].reserve(Nframe);
    keeptacvec[segmentName].reserve(Nframe);
    for (int ivs=0; ivs<statsVec.size(); ++ivs)
    {
      const auto& vs = statsVec[ivs];
      double value;
      if (currentSelectedStatID == "Mean") {
        value = vs.mean;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(vs.stddev));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Median") {
        value = vs.median;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(vs.iqr / 1.3489795003921634));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Peak") {
        value = vs.peak;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(vs.peakStddev));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else
      {
        vtkGenericWarningMacro("Unknown stat: " << currentSelectedStatID);
        return;
      }
      if (!vs.keep) {
        wgtVec[segmentName][ivs] = 0.;
      }
      tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
      keeptacvec[segmentName].push_back(vs.keep);  // Adds one-element row (column vector)
    }
  }
  if (d->weightFitCheckBox->isChecked())
  {
    for (auto& [segmentName, weights] : wgtVec)
    {
      normalizePositiveWeights(weights);
    }
  }

  this->segmentTAC4TCMfits = tac;
  this->segmentkeep4TCMfits = keeptacvec;

  const double dk = d->decayConstEdit->text().toDouble();
  const double timestep = d->timeStepEdit->text().toDouble();

  const int maxiter = d->maxIterTCMEdit->text().toInt();

  InputFunctionResult ifResult;
  QString inputFunctionError;

  if (!d->buildCurrentInputFunction(
          ifResult,
          true,
          &inputFunctionError))
  {
    QMessageBox::warning(
        this,
        tr("Input Function"),
        inputFunctionError);
    return;
  }

  const size_t requestedTCMEndCount =
      std::clamp<size_t>(
          static_cast<size_t>(d->TCMEndSlider->value()),
          1,
          static_cast<size_t>(Nframe));

  const size_t ifSupportCount =
      std::min(
          ifResult.supportFrameCount > 0
              ? ifResult.supportFrameCount
              : static_cast<size_t>(Nframe),
          static_cast<size_t>(Nframe));

  std::vector<std::vector<double>> Cp;
  std::vector<std::vector<double>> Cwb;

  Cp.reserve(ifResult.framePlasma.size());
  Cwb.reserve(ifResult.frameWholeBlood.size());

  for (double value : ifResult.framePlasma)
  {
    Cp.emplace_back(1, value);
  }

  for (double value : ifResult.frameWholeBlood)
  {
    Cwb.emplace_back(1, value);
  }

  const std::string interpolationType =
      d->selectedIFInterpolation();

  for (const std::string& segmentID : VOIsegmentIDs)
  {
    const auto& tacVOI = tac[segmentID];

    const auto statsIt = segmentTACs.find(segmentID);
    if (statsIt == segmentTACs.end())
    {
      continue;
    }

    const size_t maximumFitCount =
        std::min(requestedTCMEndCount, ifSupportCount);

    size_t tissueSupportCount = 0;
    for (size_t i = maximumFitCount; i > 0; --i)
    {
      if (statsIt->second[i - 1].keep)
      {
        tissueSupportCount = i;
        break;
      }
    }

    const size_t fitFrameCount =
        std::min(maximumFitCount, tissueSupportCount);

    if (fitFrameCount < 2)
    {
      QMessageBox::warning(
          this,
          tr("TCM fit range"),
          tr("At least two retained tissue/IF frames are required for TCM fitting."));
      continue;
    }

    std::vector<std::vector<double>> tacFit(
        tacVOI.begin(),
        tacVOI.begin() + fitFrameCount);
    std::vector<std::vector<double>> cpFit(
        Cp.begin(),
        Cp.begin() + fitFrameCount);
    std::vector<std::vector<double>> cwbFit(
        Cwb.begin(),
        Cwb.begin() + fitFrameCount);
    std::vector<std::vector<double>> framingFit(
        framing.begin(),
        framing.begin() + fitFrameCount);
    std::vector<double> weightFit(
        wgtVec[segmentID].begin(),
        wgtVec[segmentID].begin() + fitFrameCount);

    wgt = &weightFit;
    const long NframeFit = static_cast<long>(fitFrameCount);

    for (const std::string& modelID : modelsID)
    {
      if (modelID == "1TCM") {
        bool sens[] = {true, true, true, false};
        double lb_1tcm[]   = {vbLower, k1Lower, k2Lower, 0.};
        double ub_1tcm[]   = {vbUpper, k1Upper, k2Upper, 0.};
        double init_1tcm[] = {vbInit,  k1Init,  k2Init,  0.};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, maxiter, 1,
                       this->segmentTCM[segmentID]["1TCM"],
                       this->segmentTCMfits[segmentID]["1TCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TCM"]);
      }
      else if (modelID == "1TdCM") {
        bool sens[] = {true, true, true, true};
        double lb_1tcm[]   = {vbLower, k1Lower, k2Lower, tdLower};
        double ub_1tcm[]   = {vbUpper, k1Upper, k2Upper, tdUpper};
        double init_1tcm[] = {vbInit,  k1Init,  k2Init,  tdInit};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, maxiter, 1,
                       this->segmentTCM[segmentID]["1TdCM"],
                       this->segmentTCMfits[segmentID]["1TdCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TdCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TdCM"]);
      }
      else if (modelID == "1TiCM") {
        bool sens[] = {true, true, false, false};
        double lb_1tcm[]   = {vbLower, k1Lower, 0., 0.};
        double ub_1tcm[]   = {vbUpper, k1Upper, 0., 0.};
        double init_1tcm[] = {vbInit,  k1Init,  0.,  0.};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, maxiter, 1,
                       this->segmentTCM[segmentID]["1TiCM"],
                       this->segmentTCMfits[segmentID]["1TiCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TiCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TiCM"]);
      }
      else if (modelID == "1TidCM") {
        bool sens[] = {true, true, false, true};
        double lb_1tcm[]   = {vbLower, k1Lower, 0., tdLower};
        double ub_1tcm[]   = {vbUpper, k1Upper, 0., tdUpper};
        double init_1tcm[] = {vbInit,  k1Init,  0.,  tdInit};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_1tcm, lb_1tcm, ub_1tcm, sens,
                       dk, timestep, maxiter, 1,
                       this->segmentTCM[segmentID]["1TidCM"],
                       this->segmentTCMfits[segmentID]["1TidCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["1TidCM"],
        //                     Cp, framing, Nframe, Nvox, init_1tcm, lb_1tcm,
        //                     ub_1tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["1TidCM"]);
      }
      else if (modelID == "2TCM") {
        bool sens[] = {true, true, true, true, true, false};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, k4Lower, 0.};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, k4Upper, 0.};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  k4Init,  0.};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, maxiter, 2,
                       this->segmentTCM[segmentID]["2TCM"],
                       this->segmentTCMfits[segmentID]["2TCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TCM"]);
      }
      else if (modelID == "2TdCM") {
        bool sens[] = {true, true, true, true, true, true};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, k4Lower, tdLower};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, k4Upper, tdUpper};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  k4Init,  tdInit};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, maxiter, 2,
                       this->segmentTCM[segmentID]["2TdCM"],
                       this->segmentTCMfits[segmentID]["2TdCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TdCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TdCM"]);
      }
      else if (modelID == "2TiCM") {
        bool sens[] = {true, true, true, true, false, false};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, 0., 0.};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, 0., 0.};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  0., 0.};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, maxiter, 2,
                       this->segmentTCM[segmentID]["2TiCM"],
                       this->segmentTCMfits[segmentID]["2TiCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TiCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TiCM"]);
      }
      else if (modelID == "2TidCM") {
        bool sens[] = {true, true, true, true, false, true};
        double lb_2tcm[]   = {vbLower, k1Lower, k2Lower, k3Lower, 0., tdLower};
        double ub_2tcm[]   = {vbUpper, k1Upper, k2Upper, k3Upper, 0., tdUpper};
        double init_2tcm[] = {vbInit,  k1Init,  k2Init,  k3Init,  0.,  tdInit};
        logic->callTCM(tacFit, cpFit, cwbFit, framingFit, NframeFit, Nvox,
                       init_2tcm, lb_2tcm, ub_2tcm, sens,
                       dk, timestep, maxiter, 2,
                       this->segmentTCM[segmentID]["2TidCM"],
                       this->segmentTCMfits[segmentID]["2TidCM"],
                       wgt,
                       interpolationType,
                       ifResult.nativePlasmaTimesSec.empty()
                           ? nullptr : &ifResult.nativePlasmaTimesSec,
                       ifResult.nativePlasmaValues.empty()
                           ? nullptr : &ifResult.nativePlasmaValues,
                       ifResult.nativeWholeBloodTimesSec.empty()
                           ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                       ifResult.nativeWholeBloodValues.empty()
                           ? nullptr : &ifResult.nativeWholeBloodValues,
                       ifResult.parentFractionTimesSec.empty()
                           ? nullptr : &ifResult.parentFractionTimesSec,
                       ifResult.parentFractionValues.empty()
                           ? nullptr : &ifResult.parentFractionValues,
                       ifResult.plasmaIsParent
                       );
        // logic->getFittedTCM(this->segmentTCMfits[segmentID]["2TidCM"],
        //                     Cp, framing, Nframe, Nvox, init_2tcm, lb_2tcm,
        //                     ub_2tcm, sens, dk, timestep, maxiter,
        //                     1, this->segmentTCM[segmentID]["2TidCM"]);
      }
      else if (modelID == "Liver DBIF")
      {
        // Published optimization-derived liver DBIF parameterization:
        // [vb, K1, k2, k3, k4, ka, fA].
        //
        // The KMAP kernel contains td as an eighth slot. Keep td fixed
        // to zero so this first Slicer implementation remains the
        // established seven-parameter model.
        bool sens_liver[] =
        {
          true,  // vb
          true,  // K1
          true,  // k2
          true,  // k3
          true,  // k4
          true,  // ka
          true,  // fA
          false  // td fixed to 0
        };

        double lb_liver[] =
        {
          vbLower,
          k1Lower,
          k2Lower,
          k3Lower,
          k4Lower,
          liverKaLower,
          liverFaLower,
          0.0
        };

        double ub_liver[] =
        {
          vbUpper,
          k1Upper,
          k2Upper,
          k3Upper,
          k4Upper,
          liverKaUpper,
          liverFaUpper,
          0.0
        };

        double init_liver[] =
        {
          vbInit,
          k1Init,
          k2Init,
          k3Init,
          k4Init,
          liverKaInit,
          liverFaInit,
          0.0
        };

        logic->callLiverTCM(
            tacFit,
            cwbFit,
            framingFit,
            NframeFit,
            init_liver,
            lb_liver,
            ub_liver,
            sens_liver,
            dk,
            timestep,
            maxiter,
            this->segmentTCM[
                segmentID]["Liver DBIF"],
            this->segmentTCMfits[
                segmentID]["Liver DBIF"],
            wgt,
            interpolationType,
            ifResult.nativeWholeBloodTimesSec.empty()
                ? nullptr
                : &ifResult.nativeWholeBloodTimesSec,
            ifResult.nativeWholeBloodValues.empty()
                ? nullptr
                : &ifResult.nativeWholeBloodValues);
      }
      else
      {
        std::cerr << "Unknown model ID: " << modelID << std::endl;
        return;
      }

      padFittedCurveToFullFrames(
          this->segmentTCMfits[segmentID][modelID],
          fitFrameCount,
          static_cast<size_t>(Nframe));
    }
  }
  d->populateResultsVOI();
}

void qSlicerDynamicPETModuleWidget::onFITTCMImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  QString petError;

  if (!d->ensureParametricPETData(
          &petError))
  {
      QMessageBox::warning(
          this,
          tr("Parametric imaging"),
          petError);

      return;
  }

  int numVoxels = this->PETdims[0]*this->PETdims[1]*this->PETdims[2];
  if (this->PET_flatten_values.size()!=numVoxels) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of dynamic PET is not as expected."));
    return;
  }

  if (this->PET_flatten_values[0].size()!=this->numberOfTimepoints) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of timepoints is not as expected."));
    return;
  }

  if (!d->ensureParametricVoxelSelection())
  {
    QMessageBox::warning(
        nullptr,
        tr("Parametric imaging"),
        tr("Could not determine eligible PET voxels."));

    return;
  }

  if (d->parametricFitVoxelIndices.empty())
  {
    QMessageBox::warning(
        nullptr,
        tr("Parametric imaging"),
        tr("No PET voxels are eligible for fitting."));

    return;
  }


  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }


  if (modelsTCMImgID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  // Get PET reference node from the subject hierarchy item.
  // This is the same mechanism already used by Image2Flatten() and computeTAC().
  vtkMRMLScene* scene = logic->GetMRMLScene();

  if (!scene)
  {
      std::cerr << "Missing MRML scene!" << std::endl;
      return;
  }

  vtkMRMLSubjectHierarchyNode* shNode =
      vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);

  if (!shNode)
  {
      std::cerr << "Missing subject hierarchy!" << std::endl;
      return;
  }

  if (this->petID ==
      vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
      std::cerr << "Invalid PET subject hierarchy item!" << std::endl;
      return;
  }

  vtkMRMLScalarVolumeNode* refPETNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          shNode->GetItemDataNode(this->petID));

  if (!refPETNode)
  {
      std::cerr
          << "Could not retrieve PET scalar volume from petID = "
          << this->petID
          << std::endl;
      return;
  }

  // Collect parameters using a lambda for brevity
  auto getParamTriplet = [&](QLineEdit* init, QLineEdit* lb, QLineEdit* ub) {
    return std::tuple<double, double, double>{
      init->text().toDouble(), lb->text().toDouble(), ub->text().toDouble()
    };
  };

  double k1Init, k1Lower, k1Upper;
  std::tie(k1Init, k1Lower, k1Upper) =
      getParamTriplet(
          d->k1InitialImg,
          d->k1LowerImg,
          d->k1UpperImg);

  double k2Init, k2Lower, k2Upper;
  std::tie(k2Init, k2Lower, k2Upper) =
      getParamTriplet(
          d->k2InitialImg,
          d->k2LowerImg,
          d->k2UpperImg);

  double k3Init, k3Lower, k3Upper;
  std::tie(k3Init, k3Lower, k3Upper) =
      getParamTriplet(
          d->k3InitialImg,
          d->k3LowerImg,
          d->k3UpperImg);

  double k4Init, k4Lower, k4Upper;
  std::tie(k4Init, k4Lower, k4Upper) =
      getParamTriplet(
          d->k4InitialImg,
          d->k4LowerImg,
          d->k4UpperImg);

  double vbInit, vbLower, vbUpper;
  std::tie(vbInit, vbLower, vbUpper) =
      getParamTriplet(
          d->vbInitialImg,
          d->vbLowerImg,
          d->vbUpperImg);

  double tdInit, tdLower, tdUpper;
  std::tie(tdInit, tdLower, tdUpper) =
      getParamTriplet(
          d->tdInitialImg,
          d->tdLowerImg,
          d->tdUpperImg);

  const long Nframe = timePoints.size();
  const long Nvox = 1;

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }

  std::vector<double> framing_flatten =
      extractColumn(framing);

  InputFunctionResult ifResult;
  QString ifError;

  if (!d->buildCurrentInputFunction(
          ifResult,
          true,
          &ifError))
  {
      QMessageBox::warning(
          this,
          tr("Input Function"),
          ifError);
      return;
  }

  const size_t requestedTCMEndCount =
      std::clamp<size_t>(
          static_cast<size_t>(d->TCMEndSliderImg->value()),
          1,
          static_cast<size_t>(Nframe));

  const size_t ifSupportCount =
      std::min(
          ifResult.supportFrameCount > 0
              ? ifResult.supportFrameCount
              : static_cast<size_t>(Nframe),
          static_cast<size_t>(Nframe));

  const size_t fitFrameCount =
      std::min(requestedTCMEndCount, ifSupportCount);

  if (fitFrameCount < 2)
  {
    QMessageBox::warning(
        this,
        tr("TCM fit range"),
        tr("At least two IF/PET frames are required for voxelwise TCM fitting."));
    return;
  }

  std::vector<double> Cp_flatten(
      ifResult.framePlasma.begin(),
      ifResult.framePlasma.begin() + fitFrameCount);
  std::vector<double> Cwb_flatten(
      ifResult.frameWholeBlood.begin(),
      ifResult.frameWholeBlood.begin() + fitFrameCount);
  framing_flatten.resize(fitFrameCount);

  std::vector<double> inputWeights;

  if (!d->buildCurrentInputFunctionWeights(
          inputWeights,
          d->weightFitCheckBoxImg->
              isChecked(),
          &ifError,
          false))
  {
      QMessageBox::warning(
          this,
          tr("Input Function"),
          ifError);

      return;
  }

  inputWeights.resize(fitFrameCount);

  const std::vector<double>* wgt =
      &inputWeights;

  const std::string interpolationType =
      d->selectedIFInterpolation();

  const double dk =
      d->decayConstEditImg->
          text().toDouble();

  const double timestep =
      d->timeStepEdit->
          text().toDouble();

  const int maxiter =
      d->maxIterTCMEditImg->
          text().toInt();

  const int numThreads =
      d->numThreadsTCM->
          text().toInt();

  auto appendDouble =
      [](QString& key, double value)
      {
        key += "|" + QString::number(value, 'g', 17);
      };

  auto appendTriplet =
      [&](QString& key,
          double init,
          double lower,
          double upper)
      {
        appendDouble(key, init);
        appendDouble(key, lower);
        appendDouble(key, upper);
      };

  auto appendVector =
      [&](QString& key,
          const std::vector<double>& values)
      {
        for (double value : values)
        {
          appendDouble(key, value);
        }
      };

  auto makeTCMFitSignature =
      [&](const std::string& modelID)
      {
        QString key =
            QString::fromStdString(modelID)
            + "|IFsource="
            + QString::number(
                d->IFSourceSelector->currentIndex())
            + "|IFID="
            + QString::fromStdString(this->IFID)
            + "|IFfile="
            + d->externalIFPath
            + "|WBfile="
            + d->externalWholeBloodPath
            + "|PBIFfile="
            + d->pbifPath
            + "|ParentFractionFile="
            + d->parentFractionPath
            + "|curveType="
            + QString::number(
                d->IFCurveTypeSelector->currentIndex())
            + "|interp="
            + QString::fromStdString(
                interpolationType)
            + "|weighted="
            + QString::number(
                d->weightFitCheckBoxImg->isChecked());

        // Common TCM settings
        appendDouble(key, dk);
        appendDouble(key, timestep);
        key += "|fitFrameCount=" + QString::number(fitFrameCount);

        key += "|" + QString::number(maxiter);

        // All models use vb and K1.
        appendTriplet(
            key, vbInit, vbLower, vbUpper);

        appendTriplet(
            key, k1Init, k1Lower, k1Upper);

        const bool usesK2 =
            modelID == "1TCM"  ||
            modelID == "1TdCM" ||
            modelID == "2TCM"  ||
            modelID == "2TdCM" ||
            modelID == "2TiCM" ||
            modelID == "2TidCM";

        const bool usesK3 =
            modelID == "2TCM"  ||
            modelID == "2TdCM" ||
            modelID == "2TiCM" ||
            modelID == "2TidCM";

        const bool usesK4 =
            modelID == "2TCM" ||
            modelID == "2TdCM";

        const bool usesTd =
            modelID == "1TdCM"  ||
            modelID == "1TidCM" ||
            modelID == "2TdCM"  ||
            modelID == "2TidCM";

        if (usesK2)
        {
          appendTriplet(
              key, k2Init, k2Lower, k2Upper);
        }

        if (usesK3)
        {
          appendTriplet(
              key, k3Init, k3Lower, k3Upper);
        }

        if (usesK4)
        {
          appendTriplet(
              key, k4Init, k4Lower, k4Upper);
        }

        if (usesTd)
        {
          appendTriplet(
              key, tdInit, tdLower, tdUpper);
        }

        // Effective biological input data + weights.
        appendVector(key, Cp_flatten);
        appendVector(key, Cwb_flatten);
        appendVector(key, *wgt);
        appendVector(key, ifResult.nativePlasmaTimesSec);
        appendVector(key, ifResult.nativePlasmaValues);
        appendVector(key, ifResult.nativeWholeBloodTimesSec);
        appendVector(key, ifResult.nativeWholeBloodValues);
        appendVector(key, ifResult.parentFractionTimesSec);
        appendVector(key, ifResult.parentFractionValues);
        key += "|plasmaIsParent=" +
            QString::number(ifResult.plasmaIsParent);

        return key;
      };

  std::vector<std::string> modelsToFit;
  std::map<std::string, QString> pendingFitSignatures;

  for (const std::string& modelID :
       this->modelsTCMImgID)
  {
    const QString signature =
        makeTCMFitSignature(modelID);

    const auto resultIt =
        this->TCMImgOutcomes.find(modelID);

    const auto signatureIt =
        d->TCMImgFitSignatures.find(modelID);

    const bool alreadyValid =
        resultIt != this->TCMImgOutcomes.end() &&
        !resultIt->second.empty() &&
        signatureIt !=
            d->TCMImgFitSignatures.end() &&
        signatureIt->second == signature;

    if (alreadyValid)
    {
      qDebug()
          << "Reusing existing TCM voxelwise fit:"
          << QString::fromStdString(modelID);

      this->ProgressBar->setMinimum(0);
      this->ProgressBar->setMaximum(0);

      this->ProgressBar->setFormat(
          "Creating " +
          QString::fromStdString(modelID) +
          " maps in Slicer...");

      QApplication::processEvents();

      d->outputTCMParametricResult(
          modelID,
          logic,
          refPETNode,
          shNode,
          this->petID);

      continue;
    }

    // The old result no longer represents current settings.
    this->TCMImgOutcomes.erase(modelID);
    d->TCMImgFitSignatures.erase(modelID);

    modelsToFit.push_back(modelID);
    pendingFitSignatures[modelID] =
        signature;
  }

  d->populateTCMOptimizationModels();

  if (modelsToFit.empty())
  {
    qDebug()
        << "All requested TCM models already have valid fits.";
    return;
  }

  this->stopRequested = false;

  d->parametricFitRunning = true;
  d->FITbuttonMTGAImg->setEnabled(false);
  d->FITbuttonTCMImg->setEnabled(false);

  TCMWorker* worker = new TCMWorker(
    logic,
    this->PET_flatten_values,
    Cp_flatten,
    Cwb_flatten,
    framing_flatten,
    modelsToFit,
    d->parametricFitVoxelIndices,
    vbInit, vbLower, vbUpper,
    k1Init, k1Lower, k1Upper,
    k2Init, k2Lower, k2Upper,
    k3Init, k3Lower, k3Upper,
    k4Init, k4Lower, k4Upper,
    tdInit, tdLower, tdUpper,
    dk,
    timestep,
    maxiter,
    this->stopRequested,
    wgt,
    numThreads,
    interpolationType,
    ifResult.nativePlasmaTimesSec.empty()
        ? nullptr : &ifResult.nativePlasmaTimesSec,
    ifResult.nativePlasmaValues.empty()
        ? nullptr : &ifResult.nativePlasmaValues,
    ifResult.nativeWholeBloodTimesSec.empty()
        ? nullptr : &ifResult.nativeWholeBloodTimesSec,
    ifResult.nativeWholeBloodValues.empty()
        ? nullptr : &ifResult.nativeWholeBloodValues,
    ifResult.parentFractionTimesSec.empty()
        ? nullptr : &ifResult.parentFractionTimesSec,
    ifResult.parentFractionValues.empty()
        ? nullptr : &ifResult.parentFractionValues,
    ifResult.plasmaIsParent
  );
  this->ProgressBar->setMinimum(0);
  this->ProgressBar->setMaximum(100);
  this->ProgressBar->setValue(0);

  QObject::connect(
      worker,
      &TCMWorker::modelStarted,
      this,
      [this](const QString& modelID)
      {
        this->ProgressBar->setMinimum(0);
        this->ProgressBar->setMaximum(100);
        this->ProgressBar->setValue(0);

        this->ProgressBar->setFormat(
            "Fitting " +
            modelID +
            " (%p%)");

        this->ProgressBar->setVisible(true);

        this->stopButton->setEnabled(true);
        this->stopButton->setText("Stop");
        this->stopButton->setVisible(true);
        this->stopButton->show();
      });

  QObject::connect(worker, &TCMWorker::progressChanged, this, [this](int value){
      this->ProgressBar->setValue(value);
  });

  QObject::connect(worker, &TCMWorker::canceled, this, [this, worker](const QString& modelID){
      this->ProgressBar->setVisible(false);
      this->stopButton->setVisible(false);
  });

  vtkWeakPointer<vtkMRMLScalarVolumeNode> refPETNodeWeak =
      refPETNode;

  vtkWeakPointer<vtkMRMLSubjectHierarchyNode> shNodeWeak =
      shNode;

  const vtkIdType refPetID = this->petID;

  QObject::connect(worker, &TCMWorker::finishedProcessing,
    this, [this, logic, worker, refPETNodeWeak, shNodeWeak, refPetID, pendingFitSignatures](const QString& modelID){//, const std::vector<TCMParameters>& results) {

      if (!refPETNodeWeak)
      {
          std::cerr
              << "Reference PET node no longer exists."
              << std::endl;
          return;
      }

      if (!shNodeWeak)
      {
          std::cerr
              << "Subject hierarchy node no longer exists."
              << std::endl;
          return;
      }

      Q_D(qSlicerDynamicPETModuleWidget);

      const std::string id =
          modelID.toStdString();

      this->TCMImgOutcomes[id] =
          std::move(worker->results);

      auto signatureIt =
          pendingFitSignatures.find(id);

      if (signatureIt !=
          pendingFitSignatures.end())
      {
        d->TCMImgFitSignatures[id] =
            signatureIt->second;
      }

      d->populateTCMOptimizationModels();

      this->ProgressBar->setMinimum(0);
      this->ProgressBar->setMaximum(0);

      this->ProgressBar->setFormat(
          "Preparing " +
          modelID +
          " parametric outputs...");

      this->ProgressBar->setVisible(true);

      this->stopButton->setEnabled(false);
      this->stopButton->setText("Finalizing");

      QApplication::processEvents();

      d->outputTCMParametricResult(
          id,
          logic,
          refPETNodeWeak.GetPointer(),
          shNodeWeak.GetPointer(),
          refPetID);
  }, Qt::BlockingQueuedConnection);

  // connect(this->stopButton, &QPushButton::clicked, this, [this](){
  //     this->stopRequested = true;
  // });

  QObject::connect(
      worker,
      &TCMWorker::finishedAll,
      this,
      [this]()
      {
        Q_D(qSlicerDynamicPETModuleWidget);

        this->ProgressBar->setVisible(false);

        this->stopButton->setVisible(false);
        this->stopButton->setEnabled(true);
        this->stopButton->setText("Stop");

        d->parametricFitRunning = false;

        this->enableFITMTGAImgbutton();
        this->enableFITTCMImgbutton();

      });

  QObject::connect(
      worker,
      &QThread::finished,
      worker,
      &QObject::deleteLater);

  worker->start();

}

void qSlicerDynamicPETModuleWidget::onFITMTGAbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  if (segmentTACsnames.empty() || segmentTACs.empty()) {
    std::cerr << "Missing TACs!" << std::endl;
    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  int statIDQString = d->StatSelectorMTGA->currentIndex();
  if (statIDQString<0) {
    std::cerr << "Missing stat choice!" << std::endl;
    return;
  }
  std::string currentSelectedStatID = d->StatSelectorMTGA->itemData(statIDQString).toString().toStdString();


  QString ifError;

  if (!d->hasValidInputFunction(
          &ifError))
  {
      QMessageBox::warning(
          this,
          tr("Input Function"),
          ifError);

      return;
  }

  if (this->VOIMTGAsegmentIDs.empty()) {
    std::cerr << "Missing VOIs to fit!" << std::endl;
    return;
  }

  if (this->modelsMTGAID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  const long Nframe = timePoints.size();

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }
  std::vector<double> framing_flatten = extractColumn(framing);

  std::map<std::string, std::vector<std::vector<double>>> tac;
  const std::vector<double>* wgt = nullptr;
  std::map< std::string, std::vector<double>> wgtVec;

  for (const auto& [segmentName, statsVec] : segmentTACs)
  {
    if (statsVec.size() != static_cast<size_t>(Nframe))
    {
      std::cerr << "Mismatch in TAC frame size for segment " << segmentName << std::endl;
      return;
    }

    tac[segmentName].reserve(Nframe);
    for (int ivs=0; ivs<statsVec.size(); ++ivs)
    {
      const auto& vs = statsVec[ivs];
      double value;
      if (currentSelectedStatID == "Mean") {
        value = vs.mean;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(vs.stddev));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Median") {
        value = vs.median;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(vs.iqr / 1.3489795003921634));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Peak") {
        value = vs.peak;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(vs.peakStddev));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      } else {
        std::cerr << "Unknown stat: " << currentSelectedStatID << std::endl;
        return;
      }
      if (!vs.keep) {
        wgtVec[segmentName][ivs] = 0.;
      }
      tac[segmentName].emplace_back(1, value);  // Adds one-element row (column vector)
    }
  }

  if (d->weightedFitCheckBox->isChecked())
  {
    for (auto& [segmentName, weights] : wgtVec)
    {
      normalizePositiveWeights(weights);
    }
  }

  // const double timeOffset =  d->timeOffsetEdit->text().toDouble();
  const double framingNorm = d->framingNormEdit->text().toDouble();
  const size_t startFrameIndex =
      static_cast<size_t>(d->timeOffsetSlider->value() - 1);
  const double timeOffset =
      (this->timePoints[startFrameIndex]
       - 0.5 * this->durations[startFrameIndex])
      / framingNorm;
  const bool robust = d->robustFitCheckBox->isChecked();
  const bool std = d->standardizationCheckBox->isChecked();
  const double huber_tune = d->huberTuneEdit->text().toDouble();
  const double tol = d->tolEdit->text().toDouble();
  const int max_iter = d->maxIterEdit->text().toInt();


  InputFunctionResult ifResult;
  QString inputFunctionError;

  if (!d->buildCurrentInputFunction(
          ifResult,
          false,
          &inputFunctionError))
  {
    QMessageBox::warning(
        this,
        tr("Input Function"),
        inputFunctionError);
    return;
  }

  const size_t requestedMTGAEndCount =
      std::clamp<size_t>(
          static_cast<size_t>(d->timeEndSlider->value()),
          1,
          static_cast<size_t>(Nframe));
  const size_t ifSupportCount =
      std::min(
          ifResult.supportFrameCount > 0
              ? ifResult.supportFrameCount
              : static_cast<size_t>(Nframe),
          static_cast<size_t>(Nframe));

  if (ifResult.frameKeep.size() == static_cast<size_t>(Nframe))
  {
    for (auto& [segmentID, weights] : wgtVec)
    {
      for (size_t i = 0; i < weights.size(); ++i)
      {
        if (!ifResult.frameKeep[i])
        {
          weights[i] = 0.0;
        }
      }
    }
  }

  auto reconstructRemovedTissueFrames =
      [&](const std::string& segmentID,
          std::vector<double>& values) -> bool
      {
        const auto statsIt = segmentTACs.find(segmentID);
        if (statsIt == segmentTACs.end() ||
            statsIt->second.size() < values.size())
        {
          return false;
        }

        const auto& stats = statsIt->second;
        const bool hasRemoved =
            std::any_of(
                stats.begin(),
                stats.begin() + values.size(),
                [](const VoxelStatistics& vs)
                {
                  return !vs.keep;
                });

        if (!hasRemoved)
        {
          return true;
        }

        if (values.size() < 2 || !stats.front().keep)
        {
          QMessageBox::warning(
              this,
              tr("MTGA point exclusion"),
              tr("The first tissue TAC observation must remain available. "
                 "A removed terminal observation is handled as a shorter acquisition."));
          return false;
        }

        std::vector<double> retainedTimes;
        std::vector<double> retainedValues;

        for (size_t i = 0; i < values.size(); ++i)
        {
          if (!stats[i].keep)
          {
            continue;
          }

          retainedTimes.push_back(
              timePoints[i] - 0.5 * durations[i]);
          retainedValues.push_back(values[i]);
        }

        if (retainedTimes.size() < 2)
        {
          QMessageBox::warning(
              this,
              tr("MTGA point exclusion"),
              tr("At least two tissue TAC observations must remain."));
          return false;
        }

        for (size_t i = 0; i < values.size(); ++i)
        {
          if (stats[i].keep)
          {
            continue;
          }

          const double t =
              timePoints[i] - 0.5 * durations[i];

          values[i] = d->interpolateInputFunction(
              retainedTimes,
              retainedValues,
              t,
              "linear");
        }

        return true;
      };

  for (const std::string& segmentID : this->VOIMTGAsegmentIDs)
  {
    const auto& tacVOI = tac[segmentID];
    const auto statsIt = segmentTACs.find(segmentID);
    if (statsIt == segmentTACs.end())
    {
      continue;
    }

    const size_t maximumFitCount =
        std::min(requestedMTGAEndCount, ifSupportCount);

    size_t tissueSupportCount = 0;
    for (size_t i = maximumFitCount; i > 0; --i)
    {
      if (statsIt->second[i - 1].keep)
      {
        tissueSupportCount = i;
        break;
      }
    }

    const size_t fitFrameCount =
        std::min(maximumFitCount, tissueSupportCount);

    if (fitFrameCount < 2 ||
        startFrameIndex >= fitFrameCount ||
        fitFrameCount - startFrameIndex < 2)
    {
      QMessageBox::warning(
          this,
          tr("MTGA fit range"),
          tr("The selected MTGA start is outside the usable IF/tissue acquisition range."));
      continue;
    }

    auto tac_flatten = extractColumn(tacVOI);
    tac_flatten.resize(fitFrameCount);

    std::vector<double> cpFit(
        ifResult.frameModelPlasma.begin(),
        ifResult.frameModelPlasma.begin() + fitFrameCount);
    std::vector<double> framingFit(
        framing_flatten.begin(),
        framing_flatten.begin() + fitFrameCount);
    std::vector<double> weightFit(
        wgtVec[segmentID].begin(),
        wgtVec[segmentID].begin() + fitFrameCount);

    if (!reconstructRemovedTissueFrames(
            segmentID,
            tac_flatten))
    {
      continue;
    }

    wgt = &weightFit;

    for (const std::string& modelID : this->modelsMTGAID)
    {
      if (modelID == "Patlak") {
        logic->Patlak(tac_flatten,
                      cpFit,
                      framingFit,
                      this->segmentMTGA[segmentID]["Patlak"],
                      wgt,
                      timeOffset,
                      framingNorm,
                      robust,
                      std,
                      huber_tune,
                      tol,
                      max_iter
                      );
      }
      else if (modelID == "Logan") {
        logic->Logan(tac_flatten,
                     cpFit,
                     framingFit,
                     this->segmentMTGA[segmentID]["Logan"],
                     wgt,
                     timeOffset,
                     framingNorm,
                     robust,
                     std,
                     huber_tune,
                     tol,
                     max_iter
                     );
      }
      else if (modelID == "RE") {
        logic->RE(tac_flatten,
                  cpFit,
                  framingFit,
                  this->segmentMTGA[segmentID]["RE"],
                  wgt,
                  timeOffset,
                  framingNorm,
                  robust,
                  std,
                  huber_tune,
                  tol,
                  max_iter
                 );
      } else {
        std::cerr << "Unknown model ID: " << modelID << std::endl;
        return;
      }
    }
  }
  d->populateResultsVOIMTGA();
}

void qSlicerDynamicPETModuleWidget::onFITMTGAImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  QString petError;

  if (!d->ensureParametricPETData(
          &petError))
  {
      QMessageBox::warning(
          this,
          tr("Parametric imaging"),
          petError);

      return;
  }

  int numVoxels = this->PETdims[0]*this->PETdims[1]*this->PETdims[2];
  if (this->PET_flatten_values.size()!=numVoxels) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of voxels is not as expected."));
    return;
  }

  if (this->PET_flatten_values[0].size()!=this->numberOfTimepoints) {
    QMessageBox::warning(nullptr,
                         tr("Dynamic PET: mismatch values"),
                         tr("Number of timepoints is not as expected."));
    return;
  }

  if (!d->ensureParametricVoxelSelection())
  {
    QMessageBox::warning(
        nullptr,
        tr("Parametric imaging"),
        tr("Could not determine eligible PET voxels."));

    return;
  }

  if (d->parametricFitVoxelIndices.empty())
  {
    QMessageBox::warning(
        nullptr,
        tr("Parametric imaging"),
        tr("No PET voxels are eligible for fitting."));

    return;
  }

  if (durations.empty() || timePoints.empty()) {
    std::cerr << "Missing frame time information!" << std::endl;
    return;
  }

  if (this->modelsMTGAImgID.empty()) {
    std::cerr << "Missing Models to fit!" << std::endl;
    return;
  }

  // Run TAC computation
  vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(this->logic());
  if (!logic) {
    std::cerr << "Missing Logic!" << std::endl;
    return;
  }

  // Get PET reference node from the subject hierarchy item.
  // This is the same mechanism already used by Image2Flatten() and computeTAC().
  vtkMRMLScene* scene = logic->GetMRMLScene();

  if (!scene)
  {
      std::cerr << "Missing MRML scene!" << std::endl;
      return;
  }

  vtkMRMLSubjectHierarchyNode* shNode =
      vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);

  if (!shNode)
  {
      std::cerr << "Missing subject hierarchy!" << std::endl;
      return;
  }

  if (this->petID ==
      vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
      std::cerr << "Invalid PET subject hierarchy item!" << std::endl;
      return;
  }

  vtkMRMLScalarVolumeNode* refPETNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          shNode->GetItemDataNode(this->petID));

  if (!refPETNode)
  {
      std::cerr
          << "Could not retrieve PET scalar volume from petID = "
          << this->petID
          << std::endl;
      return;
  }

  const bool showInSlicer =
      d->MTGAShowInSlicerCheckBoxImg
          ->isChecked();

  const bool saveDICOM =
      d->MTGASaveDICOMCheckBoxImg
          ->isChecked();

  const QString dicomOutputDirectory =
      d->MTGADICOMDirectoryImg
          ->currentPath();

  const long Nframe = timePoints.size();

  std::vector<std::vector<double>> framing;
  framing.reserve(Nframe);
  for (double d : durations)
  {
    framing.emplace_back(1, d);  // Adds a vector with 1 element (column vector)
  }
  std::vector<double> framing_flatten = extractColumn(framing);
  InputFunctionResult ifResult;
  QString ifError;

  if (!d->buildCurrentInputFunction(
          ifResult,
          false,
          &ifError))
  {
      QMessageBox::warning(
          this,
          tr("Input Function"),
          ifError);
      return;
  }

  const size_t requestedMTGAEndCount =
      std::clamp<size_t>(
          static_cast<size_t>(d->timeEndSliderImg->value()),
          1,
          static_cast<size_t>(Nframe));
  const size_t ifSupportCount =
      std::min(
          ifResult.supportFrameCount > 0
              ? ifResult.supportFrameCount
              : static_cast<size_t>(Nframe),
          static_cast<size_t>(Nframe));
  const size_t fitFrameCount =
      std::min(requestedMTGAEndCount, ifSupportCount);

  const size_t startFrameIndex =
      static_cast<size_t>(d->timeOffsetSliderImg->value() - 1);

  if (fitFrameCount < 2 ||
      startFrameIndex >= fitFrameCount ||
      fitFrameCount - startFrameIndex < 2)
  {
    QMessageBox::warning(
        this,
        tr("MTGA fit range"),
        tr("The selected MTGA start is outside the usable IF/PET acquisition range."));
    return;
  }

  std::vector<double> Cp_flatten(
      ifResult.frameModelPlasma.begin(),
      ifResult.frameModelPlasma.begin() + fitFrameCount);
  framing_flatten.resize(fitFrameCount);

  std::vector<double> inputWeights;

  if (!d->buildCurrentInputFunctionWeights(
          inputWeights,
          d->weightedFitCheckBoxImg->
              isChecked(),
          &ifError))
  {
      QMessageBox::warning(
          this,
          tr("Input Function"),
          ifError);

      return;
  }

  inputWeights.resize(fitFrameCount);

  const std::vector<double>* wgt =
      &inputWeights;
  const double framingNorm = d->framingNormEditImg->text().toDouble();
  const double timeOffset =
      (this->timePoints[startFrameIndex]
       - 0.5 * this->durations[startFrameIndex])
      / framingNorm;
  const bool robust = d->robustFitCheckBoxImg->isChecked();
  const bool std = d->standardizationCheckBoxImg->isChecked();
  const double huber_tune = d->huberTuneEditImg->text().toDouble();
  const double tol = d->tolEditImg->text().toDouble();
  const int max_iter = d->maxIterEditImg->text().toInt();
  const int numThreads = d->numThreadsMTGA->text().toInt();


  // --------------------------------------------------------------------------
  // Fit signature used to determine whether this MTGA model must be refitted.
  // --------------------------------------------------------------------------
  auto appendDouble =
      [](QString& key, double value)
      {
        key += "|" + QString::number(value, 'g', 17);
      };

  auto appendVector =
      [&](QString& key,
          const std::vector<double>& values)
      {
        for (double value : values)
        {
          appendDouble(key, value);
        }
      };

  auto makeMTGAFitSignature =
      [&](const std::string& modelID)
      {
        QString key =
            QString::fromStdString(modelID)
            + "|IFsource="
            + QString::number(
                d->IFSourceSelector->currentIndex())
            + "|IFID="
            + QString::fromStdString(this->IFID)
            + "|IFfile="
            + d->externalIFPath
            + "|PBIFfile="
            + d->pbifPath
            + "|ParentFractionFile="
            + d->parentFractionPath
            + "|curveType="
            + QString::number(
                d->IFCurveTypeSelector->currentIndex())
            + "|interp="
            + QString::fromStdString(
                d->selectedIFInterpolation())
            + "|weighted="
            + QString::number(
                d->weightedFitCheckBoxImg->isChecked())
            + "|robust="
            + QString::number(robust)
            + "|standardize="
            + QString::number(std);

        appendDouble(key, framingNorm);
        appendDouble(key, timeOffset);
        key += "|fitFrameCount=" + QString::number(fitFrameCount);

        if (robust)
        {
          appendDouble(key, huber_tune);
          appendDouble(key, tol);

          key += "|maxiter="
              + QString::number(max_iter);
        }

        // Include actual effective input data.
        appendVector(key, Cp_flatten);
        appendVector(key, *wgt);

        return key;
      };

  std::vector<std::string> modelsToFit;

  std::map<std::string, QString>
      pendingFitSignatures;

  for (const std::string& modelID :
       this->modelsMTGAImgID)
  {
    const QString fitSignature =
        makeMTGAFitSignature(modelID);

    const auto resultIt =
        this->MTGAImgOutcomes.find(
            modelID);

    const auto signatureIt =
        d->MTGAImgFitSignatures.find(
            modelID);

    const bool alreadyValid =
        resultIt !=
            this->MTGAImgOutcomes.end() &&
        !resultIt->second.empty() &&
        signatureIt !=
            d->MTGAImgFitSignatures.end() &&
        signatureIt->second ==
            fitSignature;

    if (alreadyValid)
    {
      qDebug()
          << "Reusing existing MTGA voxelwise fit:"
          << QString::fromStdString(modelID);

      d->outputMTGAParametricResult(
          modelID,
          logic,
          refPETNode,
          shNode,
          this->petID);

      continue;
    }

    this->MTGAImgOutcomes.erase(
        modelID);

    d->MTGAImgFitSignatures.erase(
        modelID);

    modelsToFit.push_back(
        modelID);

    pendingFitSignatures[modelID] =
        fitSignature;
  }

  d->updateMTGAOptimizationUI();

  if (modelsToFit.empty())
  {
    qDebug()
        << "All requested MTGA models "
           "already have valid fits.";

    return;
  }

  this->stopRequested.store(false);

  this->stopButton->setEnabled(true);
  this->stopButton->setText("Stop");

  d->parametricFitRunning = true;
  d->FITbuttonMTGAImg->setEnabled(false);
  d->FITbuttonTCMImg->setEnabled(false);

  MTGAWorker* worker =
      new MTGAWorker(
          logic,
          this->PET_flatten_values,
          Cp_flatten,
          framing_flatten,
          modelsToFit,
          d->parametricFitVoxelIndices,
          wgt,
          timeOffset,
          framingNorm,
          robust,
          std,
          huber_tune,
          tol,
          max_iter,
          this->stopRequested,
          numThreads);

  QObject::connect(
      worker,
      &MTGAWorker::modelStarted,
      this,
      [this](const QString& modelID)
      {
        this->ProgressBar->setFormat(
            "Fitting " +
            modelID +
            " (%p%)");

        this->ProgressBar->setMinimum(0);
        this->ProgressBar->setMaximum(100);
        this->ProgressBar->setValue(0);
        this->ProgressBar->setVisible(true);

        this->stopButton->setEnabled(true);
        this->stopButton->setText("Stop");
        this->stopButton->setVisible(true);
      });

  QObject::connect(
      worker,
      &MTGAWorker::progressChanged,
      this,
      [this](int value)
      {
        this->ProgressBar->setValue(
            value);
      });

  QObject::connect(
      worker,
      &MTGAWorker::canceled,
      this,
      [this](const QString&)
      {
        this->ProgressBar->setVisible(
            false);

        this->stopButton->setVisible(
            false);
      });

  vtkWeakPointer<
      vtkMRMLScalarVolumeNode>
      refPETNodeWeak =
          refPETNode;

  vtkWeakPointer<
      vtkMRMLSubjectHierarchyNode>
      shNodeWeak =
          shNode;

  const vtkIdType refPetID =
      this->petID;

  QObject::connect(
      worker,
      &MTGAWorker::finishedProcessing,
      this,
      [this,
       logic,
       worker,
       refPETNodeWeak,
       shNodeWeak,
       refPetID,
       pendingFitSignatures]
      (const QString& modelID)
      {
        if (!refPETNodeWeak ||
            !shNodeWeak)
        {
          return;
        }

        Q_D(qSlicerDynamicPETModuleWidget);

        const std::string id =
            modelID.toStdString();

        this->MTGAImgOutcomes[id] =
            std::move(worker->results);

        auto signatureIt =
            pendingFitSignatures.find(id);

        if (signatureIt !=
            pendingFitSignatures.end())
        {
          d->MTGAImgFitSignatures[id] =
              signatureIt->second;
        }

        // Fitting is complete.  Output creation/export is a
        // separate synchronous phase and may take noticeable time.
        this->ProgressBar->setMinimum(0);
        this->ProgressBar->setMaximum(0);

        this->ProgressBar->setFormat(
            "Preparing " +
            modelID +
            " parametric outputs...");

        this->ProgressBar->setVisible(true);

        this->stopButton->setEnabled(false);
        this->stopButton->setText("Finalizing");

        QApplication::processEvents();

        d->outputMTGAParametricResult(
            id,
            logic,
            refPETNodeWeak.GetPointer(),
            shNodeWeak.GetPointer(),
            refPetID);

        d->updateMTGAOptimizationUI();
      },
      Qt::BlockingQueuedConnection);

  QObject::connect(
      worker,
      &MTGAWorker::finishedAll,
      this,
      [this]()
      {
        Q_D(qSlicerDynamicPETModuleWidget);

        this->ProgressBar->setVisible(
            false);

        this->stopButton->setVisible(
            false);

        this->stopButton->setEnabled(
            true);

        this->stopButton->setText(
            "Stop");

        d->parametricFitRunning = false;

        d->updateMTGAOptimizationUI();

        this->enableFITMTGAImgbutton();
        this->enableFITTCMImgbutton();

      });

  QObject::connect(
      worker,
      &QThread::finished,
      worker,
      &QObject::deleteLater);

  worker->start();

}


void qSlicerDynamicPETModuleWidget::onModelsAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  d->ModelsCheckContents->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsCheckLayout->count())) {
    this->modelsID.clear();
    for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsID.clear();
    for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsCheckContents->blockSignals(false);
  d->updateLiverParameterUI();
  this->enableFITbutton();

}

void qSlicerDynamicPETModuleWidget::onModelsMTGAAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  d->ModelsMTGACheckContents->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsMTGACheckLayout->count())) {
    this->modelsMTGAID.clear();
    for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsMTGAID.clear();
    for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsMTGAID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsMTGACheckContents->blockSignals(false);
  this->enableFITMTGAbutton();

}

void qSlicerDynamicPETModuleWidget::onModelsSelectAllMTGAImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  d->ModelsCheckContentsMTGAImg->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsCheckLayoutMTGAImg->count())) {
    this->modelsMTGAImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsMTGAImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsMTGAImgID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsCheckContentsMTGAImg->blockSignals(false);
  this->enableFITMTGAImgbutton();

}

void qSlicerDynamicPETModuleWidget::onModelsChanged()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  this->modelsID.clear();
  for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsID.push_back(checkbox->text().toStdString());
    }
  }
  d->updateLiverParameterUI();
  this->enableFITbutton();
  return ;
}

void qSlicerDynamicPETModuleWidget::onModelsTCMImgChanged()
{
  Q_D(const qSlicerDynamicPETModuleWidget);

  this->modelsTCMImgID.clear();
  for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsTCMImgID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITTCMImgbutton();
  return ;
}

void qSlicerDynamicPETModuleWidget::onModelsMTGAChanged()
{
  Q_D(const qSlicerDynamicPETModuleWidget);

  this->modelsMTGAID.clear();
  for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsMTGAID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITMTGAbutton();
  return ;
}

void qSlicerDynamicPETModuleWidget::onModelsMTGAImgChanged()
{
  Q_D(const qSlicerDynamicPETModuleWidget);

  this->modelsMTGAImgID.clear();
  for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutMTGAImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      this->modelsMTGAImgID.push_back(checkbox->text().toStdString());
    }
  }
  this->enableFITMTGAImgbutton();
  return ;
}

void qSlicerDynamicPETModuleWidget::onModelsTCMSelectAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->ModelsTCMCheckContents->blockSignals(true);
  QSet<QString> previouslySelectedIDs;
  for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsTCMCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedIDs.insert(checkbox->text());
    }
  }

  if (previouslySelectedIDs.size()==(d->ModelsTCMCheckLayout->count()-1)) {
    for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
    {
      QWidget* widget = d->ModelsTCMCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(false);
        cb->blockSignals(false);
      }
    }
  } else {
    for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
    {
      QWidget* widget = d->ModelsTCMCheckLayout->itemAt(i)->widget();
      QCheckBox* cb = qobject_cast<QCheckBox*>(widget);
      if (cb)
      {
        cb->blockSignals(true);
        cb->setChecked(true);
        cb->blockSignals(false);
      }
    }
  }
  d->ModelsTCMCheckContents->blockSignals(false);
  d->populateModelsTCM(this->plotTCMVOI);
}

void qSlicerDynamicPETModuleWidget::onModelsSelectAllTCMImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  d->ModelsCheckContentsTCMImg->blockSignals(true);
  std :: vector < std :: string > previouslySelectedModels;
  for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
  {
    QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      previouslySelectedModels.push_back(checkbox->text().toStdString());
    }
  }

  if (previouslySelectedModels.size()==(d->ModelsCheckLayoutTCMImg->count())) {
    this->modelsTCMImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(false);
        checkbox->blockSignals(false);
      }
    }
  } else {
    this->modelsTCMImgID.clear();
    for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
    {
      QLayoutItem* item = d->ModelsCheckLayoutTCMImg->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (checkbox)
      {
        checkbox->blockSignals(true);
        checkbox->setChecked(true);
        checkbox->blockSignals(false);
        this->modelsTCMImgID.push_back(checkbox->text().toStdString());
      }
    }
  }
  d->ModelsCheckContentsTCMImg->blockSignals(false);
  this->enableFITTCMImgbutton();

}


void qSlicerDynamicPETModuleWidget::onPlotTCMbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  this->ColNameToSegmentID.clear();
  this->MapPlotSeriesNodeIDToPlot.clear();
  if (this->plotTCMVOI.empty()) {
    return;
  }

  vtkMRMLScene* scene = this->mrmlScene();

  // Get selected VOI to plot TCM fit
  std :: string selectedVOI = this->plotTCMVOI;

  // Get TAC to plot
  std::vector<std::vector<double>> tacvoi = this->segmentTAC4TCMfits[selectedVOI];
  std::vector<bool> keepvoi = this->segmentkeep4TCMfits[selectedVOI];

  // Get TCM models to plot
  std::vector<std::string> PlotSelectedTCMs;
  for (int i = 0; i < d->ModelsTCMCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->ModelsTCMCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      PlotSelectedTCMs.push_back(checkbox->text().toStdString());
    }
  }

  std::map<std::string, double*> TCMfits = this->segmentTCMfits[selectedVOI];

  if (selectedVOI.empty() || segmentTAC4TCMfits.empty() || PlotSelectedTCMs.empty() || TCMfits.empty())
    return;

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  // Add time column (convert to minutes)
  vtkNew<vtkDoubleArray> timeArray;
  timeArray->SetName("Time (min)");
  for (size_t i = 0; i < this->timePoints.size(); ++i)
  {
    timeArray->InsertNextValue(this->timePoints[i] / 60.0);
  }
  tableNode->AddColumn(timeArray);

  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();
  const ActivityUnit displayUnit = d->selectedDisplayActivityUnit();

  auto toDisplayValue = [&](double nativeValue)
  {
    if (!std::isfinite(nativeValue))
    {
      return nativeValue;
    }

    double converted = nativeValue;
    if (!d->convertActivityValue(
            nativeValue,
            d->petStoredActivityUnit(),
            displayUnit,
            converted,
            nullptr))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return converted;
  };

  // TAC as scatterpoints
  vtkNew<vtkDoubleArray> tacArray;
  tacArray->SetName("TAC");
  vtkNew<vtkStringArray> labelArray;
  labelArray->SetName("ToolTipLabelTAC");
  for (size_t i = 0; i < tacvoi.size(); ++i)
  {
      // Assuming tacvoi[0] holds measured values for VOI (adapt if structured differently)
      tacArray->InsertNextValue(
          keepvoi[i]
          ? toDisplayValue(tacvoi[i][0])
          : std::numeric_limits<double>::quiet_NaN());
      std::ostringstream oss;
      oss << "Frame: " << i
          << ", Time(s): " << this->timePoints[i]
          << ", Time(min): " << this->timePoints[i]/60.0;
      labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(tacArray);
  tableNode->AddColumn(labelArray);

  // TCM fits as line plots
  for (const std::string& modelName : PlotSelectedTCMs)
  {
      auto it = TCMfits.find(modelName);
      if (it == TCMfits.end() || !it->second)
        continue;
      TCMParameters params = this->segmentTCM[selectedVOI][modelName];

      double* fitArrayPtr = it->second;

      vtkNew<vtkDoubleArray> fitArray;
      fitArray->SetName(modelName.c_str());
      for (size_t ivs = 0; ivs < this->timePoints.size(); ++ivs)
      {
        double fv = fitArrayPtr[ivs];
        fitArray->InsertNextValue(toDisplayValue(fv));
        // if (params.keep[ivs]) {
        //   double fv = fitArrayPtr[ivs];
        //   fitArray->InsertNextValue(fv);
        // } else {
        //   double nextValue = std::numeric_limits<double>::quiet_NaN();
        //   double x1 = std::numeric_limits<double>::quiet_NaN();
        //   for (int next_ivs = ivs+1; next_ivs<this->timePoints.size(); ++next_ivs) {
        //     if (params.keep[next_ivs])
        //     {
        //         x1 = this->timePoints[next_ivs];
        //         nextValue = fitArrayPtr[next_ivs];
        //     }
        //   }
        //   if (std::isnan(nextValue)) {
        //     fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        //     continue;
        //   }
        //   double prevValue = std::numeric_limits<double>::quiet_NaN();
        //   double x0 = std::numeric_limits<double>::quiet_NaN();
        //   for (int prev_ivs = ivs-1; prev_ivs>=0; --prev_ivs) {
        //     if (params.keep[prev_ivs])
        //     {
        //         x0 = this->timePoints[prev_ivs];
        //         prevValue = fitArrayPtr[prev_ivs];
        //     }
        //   }
        //   if (std::isnan(prevValue)) {
        //     fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        //     continue;
        //   }
        //
        //   double x  = this->timePoints[ivs];
        //   // Proper linear interpolation
        //   double value = prevValue + ((x - x0) / (x1 - x0)) * (nextValue - prevValue);
        //   fitArray->InsertNextValue(value);
        // }
      }

      tableNode->AddColumn(fitArray);

      vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
      scene->AddNode(lineSeries);
      lineSeries->SetName(modelName.c_str());
      lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);  // line plot
      lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
      lineSeries->SetXColumnName("Time (min)");
      lineSeries->SetYColumnName(modelName.c_str());
      lineSeries->SetLabelColumnName("ToolTipLabelTAC");
      lineSeries->SetUniqueColor();
      lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
      chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());
  }

  this->ColNameToSegmentID["TAC"]=selectedVOI;
  vtkSmartPointer<vtkMRMLPlotSeriesNode> scatterSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(scatterSeries);
  scatterSeries->SetName("TAC");
  scatterSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);  // scatter points
  scatterSeries->SetAndObserveTableNodeID(tableNode->GetID());
  scatterSeries->SetXColumnName("Time (min)");
  scatterSeries->SetYColumnName("TAC");
  scatterSeries->SetLabelColumnName("ToolTipLabelTAC");
  scatterSeries->SetUniqueColor();
  scatterSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(scatterSeries->GetID());
  chartNode->SetTitle(this->segmentTACsnames[selectedVOI].c_str());
  chartNode->SetXAxisTitle("Time (min)");
  chartNode->SetYAxisTitle(
      d->activityUnitLabel(displayUnit).toStdString().c_str());

  // Show plot view
  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
      scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode)
  {
      plotViewNode->SetPlotChartNodeID(chartNode->GetID());
      // qMRMLPlotWidget* plotWidget = nullptr;
      // if (qSlicerApplication::application())
      // {
      //     qSlicerLayoutManager* layoutManager =
      //         qSlicerApplication::application()->layoutManager();
      //     qMRMLPlotWidget* plotWidget = nullptr;
      //     plotWidget = layoutManager->plotWidget(0);
      //     qMRMLPlotView* plotView = plotWidget->plotView();
      //     if (plotView)
      //     {
      //       QObject::connect(plotView, SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
      //                        this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)));
      //     }
      // }
  }

}

void qSlicerDynamicPETModuleWidget::onPlotMTGAbutton() {
  Q_D(qSlicerDynamicPETModuleWidget);
  this->ColNameToSegmentID.clear();
  this->MapPlotSeriesNodeIDToPlot.clear();

  if (this->plotMTGAVOI.empty()) {
    return;
  }

  vtkMRMLScene* scene = this->mrmlScene();

  // Get selected VOI for MTGA plot
  std::string selectedVOI = this->plotMTGAVOI;

  // Get selected MTGA model from combo box
  QString modelName_qstr = d->MTGASelector->currentText();
  if (modelName_qstr.toStdString().empty())
    return;
  std::string modelName = modelName_qstr.toStdString();

  // Check data availability
  auto voiIt = this->segmentMTGA.find(selectedVOI);
  if (voiIt == this->segmentMTGA.end())
    return;

  auto modelIt = voiIt->second.find(modelName);
  if (modelIt == voiIt->second.end())
    return;

  const MTGAParameters& params = modelIt->second;
  if (params.x.empty() || params.y.empty() || params.fitted.empty())
    return;

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  vtkNew<vtkDoubleArray> xArray;
  vtkNew<vtkDoubleArray> yArray;
  vtkNew<vtkStringArray> labelArray;
  xArray->SetName("X");
  yArray->SetName("Data");
  labelArray->SetName("ToolTipData");
  for (int i=0;  i < params.x.size(); ++i) {
    double xv = params.x[i];
    double yv = params.y[i];
    xArray->InsertNextValue(xv);
    yArray->InsertNextValue(params.keep[i] ? yv : std::numeric_limits<double>::quiet_NaN());
    const int frameIndex = params.frame[i] - 1;
    std::ostringstream oss;
    oss << "Frame: " << frameIndex
        << ", Time(s): " << this->timePoints[frameIndex]
        << ", Time(min): " << this->timePoints[frameIndex]/60.0;
    labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(xArray);
  tableNode->AddColumn(yArray);
  tableNode->AddColumn(labelArray);


  // Create plot chart
  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();

  // Fitted values as line plot
  vtkNew<vtkDoubleArray> fitArray;
  fitArray->SetName(modelName.c_str());
  for (int ivs=0; ivs<params.fitted.size(); ++ivs) {
    if (params.keep[ivs]) {
      double fv = params.fitted[ivs];
      fitArray->InsertNextValue(fv);
    } else {
      double nextValue = std::numeric_limits<double>::quiet_NaN();
      double x1 = std::numeric_limits<double>::quiet_NaN();
      for (int next_ivs = ivs+1; next_ivs<params.fitted.size(); ++next_ivs) {
        if (params.keep[next_ivs])
        {
            x1 = params.x[next_ivs];
            nextValue = params.fitted[next_ivs];
        }
      }
      if (std::isnan(nextValue)) {
        fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        continue;
      }
      double prevValue = std::numeric_limits<double>::quiet_NaN();
      double x0 = std::numeric_limits<double>::quiet_NaN();
      for (int prev_ivs = ivs-1; prev_ivs>=0; --prev_ivs) {
        if (params.keep[prev_ivs])
        {
            x0 = params.x[prev_ivs];
            prevValue = params.fitted[prev_ivs];
        }
      }
      if (std::isnan(prevValue)) {
        fitArray->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        continue;
      }

      double x  = params.x[ivs];
      // Proper linear interpolation
      double value = prevValue + ((x - x0) / (x1 - x0)) * (nextValue - prevValue);
      fitArray->InsertNextValue(value);
    }
  }
  tableNode->AddColumn(fitArray);

  vtkSmartPointer<vtkMRMLPlotSeriesNode> lineSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(lineSeries);
  lineSeries->SetName(modelName.c_str());
  lineSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
  lineSeries->SetAndObserveTableNodeID(tableNode->GetID());
  lineSeries->SetXColumnName("X");
  lineSeries->SetYColumnName(modelName.c_str());
  lineSeries->SetLabelColumnName("ToolTipData");
  lineSeries->SetUniqueColor();
  lineSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(lineSeries->GetID());

  // Scatter series for measured values
  this->ColNameToSegmentID["Data"] = selectedVOI;
  vtkSmartPointer<vtkMRMLPlotSeriesNode> scatterSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(scatterSeries);
  scatterSeries->SetName("Data");
  scatterSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
  scatterSeries->SetAndObserveTableNodeID(tableNode->GetID());
  scatterSeries->SetXColumnName("X");
  scatterSeries->SetYColumnName("Data");
  scatterSeries->SetLabelColumnName("ToolTipData");
  scatterSeries->SetUniqueColor();
  scatterSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(scatterSeries->GetID());
  chartNode->SetTitle((modelName + " - " + this->segmentTACsnames[selectedVOI]).c_str());
  if (modelName=="Patlak") {
    chartNode->SetXAxisTitle("intCp/Cp");
    chartNode->SetYAxisTitle("Ct/Cp");
  } else if (modelName=="Logan") {
    chartNode->SetXAxisTitle("intCp/Ct");
    chartNode->SetYAxisTitle("intCt/Ct");
  } else if (modelName=="RE") {
    chartNode->SetXAxisTitle("intCp/Cp");
    chartNode->SetYAxisTitle("intCt/Cp");
  } else {
    std::cerr << "Unknown model: " << modelName << std::endl;
  }

  // Show plot view
  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
      scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode) {
    plotViewNode->SetPlotChartNodeID(chartNode->GetID());
    // qMRMLPlotWidget* plotWidget = nullptr;
    // if (qSlicerApplication::application())
    // {
    //     qSlicerLayoutManager* layoutManager =
    //         qSlicerApplication::application()->layoutManager();
    //     qMRMLPlotWidget* plotWidget = nullptr;
    //     plotWidget = layoutManager->plotWidget(0);
    //     qMRMLPlotView* plotView = plotWidget->plotView();
    //     if (plotView)
    //     {
    //       QObject::connect(plotView, SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
    //                        this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)));
    //     }
    // }
  }
}

bool qSlicerDynamicPETModuleWidget::checkdisplayedDynamicPET() {
  vtkMRMLPlotViewNode* plotViewNode = nullptr;
  vtkCollection* viewNodes = this->mrmlScene()->GetNodesByClass("vtkMRMLPlotViewNode");
  if (!viewNodes || viewNodes->GetNumberOfItems() == 0)
  {
      return false; // no plot view node
  }

  plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(viewNodes->GetItemAsObject(0));
  if (!plotViewNode)
  {
      return false;
  }

  // Get the currently displayed plot chart
  vtkMRMLPlotChartNode* currentPlot = vtkMRMLPlotChartNode::SafeDownCast(
      this->mrmlScene()->GetNodeByID(plotViewNode->GetPlotChartNodeID())
  );
  if (!currentPlot)
  {
      return false;
  }

  // Check if the currently displayed plot is your "DynamicPET.PlotChart"
  if (currentPlot->GetName() == nullptr || std::string(currentPlot->GetName()) != "DynamicPET.PlotChart")
  {
      return false; // not the DynamicPET plot, do nothing
  }
  return true;
}

void qSlicerDynamicPETModuleWidget::onDeleteKeyPressed() {
  Q_D(qSlicerDynamicPETModuleWidget);
  if (!this->checkdisplayedDynamicPET())
    return;

  if (this->PlotSelectedFrame == -1)
    return;

  if (this->PlotSelectedVOI == "")
    return;

  if (this->segmentTACs.empty())
    return;

  vtkGenericWarningMacro("Segment " << this->segmentTACsnames[this->PlotSelectedVOI] << " removed at frame " << this->PlotSelectedFrame);
  this->segmentTACs[this->PlotSelectedVOI][this->PlotSelectedFrame].keep = false;
  this->segmentTACs[this->PlotSelectedVOI][this->PlotSelectedFrame].empty = false;

  if (this->PlotSelectedVOI == this->IFID)
  {
    d->invalidateInputFunctionResults();
    d->updateInputFunctionStatus();
    d->updateROIModelingAvailability();
    d->updateParametricImagingAvailability();
  }

  this->onPlotbutton();
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI = "";
  return;
}

void qSlicerDynamicPETModuleWidget::onResetbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  vtkGenericWarningMacro("Restoring removed points ");

  if (this->segmentTACs.empty())
    return;

  bool currentInputWasModified = false;
  const auto ifIt = this->segmentTACs.find(this->IFID);
  if (ifIt != this->segmentTACs.end())
  {
    currentInputWasModified =
        std::any_of(
            ifIt->second.begin(),
            ifIt->second.end(),
            [](const VoxelStatistics& vs)
            {
              return !vs.keep && !vs.empty;
            });
  }

  for (auto& segmentPair : this->segmentTACs)
  {
    const std::string& segmentID = segmentPair.first;
    std::vector<VoxelStatistics>& tacVector = segmentPair.second;

    for (size_t frameID = 0; frameID < tacVector.size(); ++frameID)
    {
      VoxelStatistics& vs = tacVector[frameID];
      if (!vs.empty)
      {
          vs.keep = true;
      }
    }
  }
  if (currentInputWasModified)
  {
    d->invalidateInputFunctionResults();
    d->updateInputFunctionStatus();
    d->updateROIModelingAvailability();
    d->updateParametricImagingAvailability();
  }

  if (this->checkdisplayedDynamicPET()) {
    this->onPlotbutton();
  }
  return;
}
