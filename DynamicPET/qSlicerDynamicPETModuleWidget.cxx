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

#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif
#include <PythonQt.h>


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
#include <QGroupBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QBoxLayout>
#include <QLineEdit>
#include <QLabel>
#include <QSlider>
#include <QFileInfo>
#include <QSignalBlocker>
#include <QSizePolicy>
#include <QDateTime>
#include <QDir>
#include <QTableWidget>
#include <QHeaderView>
#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>

#include <vtkWeakPointer.h>
#include <vtkDataArray.h>
#include <vtkSlicerSegmentationsModuleLogic.h>
#include <vtkSegmentationConverter.h>

#include <cmath>
#include <limits>
#include <algorithm>
#include <array>
#include <chrono>
#include <map>
#include <set>
#include <tuple>
#include <exception>

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
  WholeBlood = 0,
  TotalPlasma = 1
};

enum class ActivityUnit
{
  BqPerMl = 0,
  KBqPerMl = 1,
  MBqPerMl = 2,
  SUVbw = 3
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
  size_t supportFrameStartIndex{0};
  size_t supportFrameCount{0};

  bool plasmaIsParent{false};
  bool applyParentFraction{false};
  bool hasWholeBlood{false};

  // Internal support classification for the forthcoming relative-model gate.
  // This is deliberately not exported as workbook metadata.
  bool inputCoversFromInjection{false};
  bool inputCoverageReconstructedByPBIF{false};
  double earliestAvailableInputTimeSec{
      std::numeric_limits<double>::quiet_NaN()};

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

  // Explicit model-support provenance. These flags distinguish measured input
  // support from intervals supplied only by a fitted analytic model.
  bool fengExtrapolationApplied{false};
  double sourceMeasuredEndTimeSec{std::numeric_limits<double>::quiet_NaN()};
  double sourceModeledEndTimeSec{std::numeric_limits<double>::quiet_NaN()};
  bool parentFractionExtrapolationApplied{false};
  double parentFractionMeasuredEndTimeSec{std::numeric_limits<double>::quiet_NaN()};
  double parentFractionModeledEndTimeSec{std::numeric_limits<double>::quiet_NaN()};
};

struct PreviewCurve
{
  QString name;
  std::vector<double> times;
  std::vector<double> values;
  bool pointsOnly{false};
  bool sourceObservations{false};
  std::string observationRole;
};

struct TACModeState
{
  bool valid{false};
  std::map<std::string, std::vector<VoxelStatistics>> segmentTACs;
  std::map<std::string, std::string> segmentTACsnames;
  std::vector<double> timePoints;
  std::vector<double> durations;
  int numberOfTimepoints{0};
  std::vector<std::string> segmentDisplayOrder;
};

struct AcquisitionTimingContext
{
  bool timingAvailable{false};
  bool delayedAcquisition{false};
  bool tableTimesAlreadyPostInjection{false};
  double rawInjectionToAcquisitionOffsetSec{
      std::numeric_limits<double>::quiet_NaN()};
  double acquisitionStartPostInjectionSec{0.0};
  QString injectionDateTime;
  QString firstFrameDateTime;
  QString source;
};

struct MultiTimepointCandidate
{
  vtkIdType studyItemID{vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID};
  vtkIdType petItemID{vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID};
  QString studyName;
  QString petName;
  QString acquisitionType;
  QString timingText;
  QDateTime acquisitionStart;
  QDateTime acquisitionEnd;
  double durationSec{std::numeric_limits<double>::quiet_NaN()};
  bool dynamic{false};
  bool kineticMetadataReady{false};
  QString metadataNodeID;
};

struct MultiTimepointFrameInterval
{
  QDateTime start;
  QDateTime end;
  double durationSec{std::numeric_limits<double>::quiet_NaN()};
};

struct PreparedMultiTimepointAcquisition
{
  int sourceRow{-1};
  bool dynamic{false};
  bool segmentationTemporal{false};
  QString petName;
  QString petNodeID;
  QString segmentationNodeID;
  QString metadataNodeID;
  QString petSequenceNodeID;
  QString petBrowserNodeID;
  QString segmentationSequenceNodeID;
  QDateTime start;
  QDateTime end;
  double durationSec{std::numeric_limits<double>::quiet_NaN()};
  QString valueType;
  QString decayCorrection;
  double sourceSUVbwFactor{std::numeric_limits<double>::quiet_NaN()};
};

struct PreparedMultiTimepointObservation
{
  int acquisitionIndex{-1};
  int frameIndex{-1};
  bool dynamic{false};
  QString acquisitionName;
  QString sequenceIndex;
  QString petNodeID;
  QString segmentationNodeID;
  QString petSequenceNodeID;
  QString segmentationSequenceNodeID;
  QString metadataNodeID;
  QDateTime start;
  QDateTime end;
  double durationSec{std::numeric_limits<double>::quiet_NaN()};
  std::map<std::string, std::string> segmentIDsByName;
};

struct PendingMultiSegStructureChange
{
  int acquisitionIndex{-1};
  int frameIndex{-1};
  vtkWeakPointer<vtkMRMLSegmentationNode> sourceNode;
  std::string segmentID;
  SegmentationChangeWatcher::StructureChangeType changeType{
      SegmentationChangeWatcher::StructureChangeType::Added};
};

struct TACStatisticOption
{
  QString label;
  std::string id;
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

static std::vector<double> gaussianKernelPredict(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& evalX,
    double sigmaSec)
{
    if (x.empty() || x.size() != y.size() || !(sigmaSec > 0.0))
    {
        return {};
    }

    std::vector<double> out;
    out.reserve(evalX.size());
    const double invTwoSigma2 = 1.0 / (2.0 * sigmaSec * sigmaSec);
    const double cutoff = 4.0 * sigmaSec;

    for (double t : evalX)
    {
        double sw = 0.0;
        double swy = 0.0;
        size_t nearest = 0;
        double nearestDistance = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < x.size(); ++i)
        {
            const double distance = std::abs(x[i] - t);
            if (distance < nearestDistance)
            {
                nearestDistance = distance;
                nearest = i;
            }
            if (distance > cutoff)
            {
                continue;
            }
            const double weight = std::exp(-distance * distance * invTwoSigma2);
            sw += weight;
            swy += weight * y[i];
        }
        out.push_back(sw > 1e-14 ? swy / sw : y[nearest]);
    }
    return out;
}

static bool pchipSlopes(
    const std::vector<double>& x,
    const std::vector<double>& y,
    std::vector<double>& slopes)
{
    const size_t n = x.size();
    if (n < 2 || n != y.size())
    {
        return false;
    }
    for (size_t i = 1; i < n; ++i)
    {
        if (!(x[i] > x[i - 1]))
        {
            return false;
        }
    }

    slopes.assign(n, 0.0);
    if (n == 2)
    {
        const double d = (y[1] - y[0]) / (x[1] - x[0]);
        slopes[0] = slopes[1] = d;
        return true;
    }

    std::vector<double> h(n - 1);
    std::vector<double> delta(n - 1);
    for (size_t i = 0; i + 1 < n; ++i)
    {
        h[i] = x[i + 1] - x[i];
        delta[i] = (y[i + 1] - y[i]) / h[i];
    }

    for (size_t i = 1; i + 1 < n; ++i)
    {
        if (delta[i - 1] == 0.0 || delta[i] == 0.0 ||
            delta[i - 1] * delta[i] <= 0.0)
        {
            slopes[i] = 0.0;
        }
        else
        {
            const double w1 = 2.0 * h[i] + h[i - 1];
            const double w2 = h[i] + 2.0 * h[i - 1];
            slopes[i] = (w1 + w2) /
                (w1 / delta[i - 1] + w2 / delta[i]);
        }
    }

    auto endpointSlope = [](double h0, double h1, double d0, double d1)
    {
        double m = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (m * d0 <= 0.0)
        {
            return 0.0;
        }
        if (d0 * d1 < 0.0 && std::abs(m) > 3.0 * std::abs(d0))
        {
            return 3.0 * d0;
        }
        return m;
    };

    slopes.front() = endpointSlope(h[0], h[1], delta[0], delta[1]);
    slopes.back() = endpointSlope(
        h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
    return true;
}

static double pchipInterpolate(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double target)
{
    if (x.size() < 2 || x.size() != y.size() ||
        target < x.front() || target > x.back())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (target == x.back())
    {
        return y.back();
    }

    std::vector<double> m;
    if (!pchipSlopes(x, y, m))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const auto upper = std::upper_bound(x.begin(), x.end(), target);
    const size_t right = static_cast<size_t>(std::distance(x.begin(), upper));
    const size_t left = right - 1;
    const double h = x[right] - x[left];
    const double u = (target - x[left]) / h;
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    const double h10 = u3 - 2.0 * u2 + u;
    const double h01 = -2.0 * u3 + 3.0 * u2;
    const double h11 = u3 - u2;
    return std::max(0.0,
        h00 * y[left] + h10 * h * m[left] +
        h01 * y[right] + h11 * h * m[right]);
}

QDateTime parseDICOMDateTimeText(const QString& value)
{
  QString text = value.trimmed();
  if (text.size() < 14)
  {
    return QDateTime();
  }

  QDateTime dt = QDateTime::fromString(
      text.left(14),
      QStringLiteral("yyyyMMddHHmmss"));
  if (!dt.isValid())
  {
    return QDateTime();
  }

  // Preserve the fractional part when present. DICOM DT permits arbitrary
  // fractional-second precision; QDateTime stores milliseconds, which is more
  // than sufficient for PET frame/acquisition chronology here.
  if (text.size() > 15 && text.at(14) == QLatin1Char('.'))
  {
    QString fraction;
    for (int i = 15; i < text.size() && text.at(i).isDigit(); ++i)
    {
      fraction.append(text.at(i));
    }
    if (!fraction.isEmpty())
    {
      const QString millisecondsText = fraction.left(3).leftJustified(3, QLatin1Char('0'));
      bool ok = false;
      const int milliseconds = millisecondsText.toInt(&ok);
      if (ok)
      {
        dt = dt.addMSecs(milliseconds);
      }
    }
  }
  return dt;
}

QString nodeAttributeText(vtkMRMLNode* node, const char* name)
{
  if (!node || !name)
  {
    return QString();
  }
  const char* value = node->GetAttribute(name);
  return value ? QString::fromUtf8(value).trimmed() : QString();
}

QString kineticMetadataCommonText(vtkMRMLNode* node, const QString& key)
{
  if (!node || key.isEmpty())
  {
    return QString();
  }

  const QString metadataText = nodeAttributeText(node, "dPET.KineticMetadata");
  if (metadataText.isEmpty())
  {
    return QString();
  }

  QJsonParseError parseError;
  const QJsonDocument document = QJsonDocument::fromJson(metadataText.toUtf8(), &parseError);
  if (parseError.error != QJsonParseError::NoError || !document.isObject())
  {
    return QString();
  }

  const QJsonObject root = document.object();
  const QJsonValue commonValue = root.value(QStringLiteral("common"));
  if (!commonValue.isObject())
  {
    return QString();
  }

  const QJsonValue value = commonValue.toObject().value(key);
  if (value.isString())
  {
    return value.toString().trimmed();
  }
  if (value.isDouble())
  {
    return QString::number(value.toDouble(), 'g', 16);
  }
  if (value.isBool())
  {
    return value.toBool() ? QStringLiteral("1") : QStringLiteral("0");
  }
  return QString();
}

QString nodeOrKineticMetadataText(
    vtkMRMLNode* node,
    const char* attributeName,
    const QString& commonKey)
{
  const QString direct = nodeAttributeText(node, attributeName);
  if (!direct.isEmpty())
  {
    return direct;
  }
  return kineticMetadataCommonText(node, commonKey);
}

bool readPersistedKineticTiming(
    vtkMRMLNode* node,
    QDateTime& startDT,
    QDateTime& endDT,
    double& durationSec)
{
  startDT = parseDICOMDateTimeText(nodeAttributeText(node, "dPET.AcquisitionStartDateTime"));
  endDT = parseDICOMDateTimeText(nodeAttributeText(node, "dPET.AcquisitionEndDateTime"));
  durationSec = std::numeric_limits<double>::quiet_NaN();

  const QString metadataText = nodeAttributeText(node, "dPET.KineticMetadata");
  if (!metadataText.isEmpty())
  {
    QJsonParseError parseError;
    const QJsonDocument document = QJsonDocument::fromJson(metadataText.toUtf8(), &parseError);
    if (parseError.error == QJsonParseError::NoError && document.isObject())
    {
      const QJsonObject root = document.object();
      const QDateTime jsonStart = parseDICOMDateTimeText(
          root.value(QStringLiteral("acquisitionStartDateTime")).toString());
      const QDateTime jsonEnd = parseDICOMDateTimeText(
          root.value(QStringLiteral("acquisitionEndDateTime")).toString());

      if (!startDT.isValid() && jsonStart.isValid())
      {
        startDT = jsonStart;
      }
      if ((!endDT.isValid() || (startDT.isValid() && endDT <= startDT)) && jsonEnd.isValid())
      {
        endDT = jsonEnd;
      }

      // The top-level end time should already be authoritative.  For robustness,
      // recover the full range from per-frame/per-slice records if a stale or
      // incomplete scalar attribute was encountered.
      QDateTime recordStart;
      QDateTime recordEnd;
      const auto accumulateRecords = [&](const QJsonArray& records)
      {
        for (const QJsonValue& value : records)
        {
          if (!value.isObject())
          {
            continue;
          }
          const QJsonObject record = value.toObject();
          const QDateTime rs = parseDICOMDateTimeText(
              record.value(QStringLiteral("acquisitionStartDateTime")).toString());
          QDateTime re = parseDICOMDateTimeText(
              record.value(QStringLiteral("acquisitionEndDateTime")).toString());
          const double d = record.value(QStringLiteral("durationSec")).toDouble(
              std::numeric_limits<double>::quiet_NaN());
          if (!re.isValid() && rs.isValid() && std::isfinite(d) && d > 0.0)
          {
            re = rs.addMSecs(static_cast<qint64>(std::llround(d * 1000.0)));
          }
          if (rs.isValid() && (!recordStart.isValid() || rs < recordStart))
          {
            recordStart = rs;
          }
          if (re.isValid() && (!recordEnd.isValid() || re > recordEnd))
          {
            recordEnd = re;
          }
        }
      };

      if (root.value(QStringLiteral("frames")).isArray())
      {
        accumulateRecords(root.value(QStringLiteral("frames")).toArray());
      }
      if (root.value(QStringLiteral("spatialTiming")).isArray())
      {
        accumulateRecords(root.value(QStringLiteral("spatialTiming")).toArray());
      }

      if (!startDT.isValid() && recordStart.isValid())
      {
        startDT = recordStart;
      }
      if ((!endDT.isValid() || (startDT.isValid() && endDT <= startDT)) && recordEnd.isValid())
      {
        endDT = recordEnd;
      }
    }
  }

  if (!startDT.isValid())
  {
    startDT = parseDICOMDateTimeText(
        nodeAttributeText(node, "dPET.FirstFrameAcquisitionDateTime"));
  }

  if (startDT.isValid() && endDT.isValid() && endDT > startDT)
  {
    durationSec = static_cast<double>(startDT.msecsTo(endDT)) / 1000.0;
  }

  return startDT.isValid();
}


bool readPersistedDynamicFrameIntervals(
    vtkMRMLNode* metadataNode,
    std::vector<MultiTimepointFrameInterval>& intervals)
{
  intervals.clear();
  if (!metadataNode)
  {
    return false;
  }

  const QString metadataText = nodeAttributeText(metadataNode, "dPET.KineticMetadata");
  if (metadataText.isEmpty())
  {
    return false;
  }

  QJsonParseError parseError;
  const QJsonDocument document = QJsonDocument::fromJson(metadataText.toUtf8(), &parseError);
  if (parseError.error != QJsonParseError::NoError || !document.isObject())
  {
    return false;
  }

  const QJsonValue framesValue = document.object().value(QStringLiteral("frames"));
  if (!framesValue.isArray())
  {
    return false;
  }

  const QJsonArray frames = framesValue.toArray();
  intervals.reserve(static_cast<size_t>(frames.size()));
  for (const QJsonValue& value : frames)
  {
    if (!value.isObject())
    {
      intervals.clear();
      return false;
    }

    const QJsonObject frame = value.toObject();
    MultiTimepointFrameInterval interval;
    interval.start = parseDICOMDateTimeText(
        frame.value(QStringLiteral("acquisitionStartDateTime")).toString());
    interval.end = parseDICOMDateTimeText(
        frame.value(QStringLiteral("acquisitionEndDateTime")).toString());
    interval.durationSec = frame.value(QStringLiteral("durationSec")).toDouble(
        std::numeric_limits<double>::quiet_NaN());

    if (!interval.start.isValid() || !std::isfinite(interval.durationSec) ||
        !(interval.durationSec > 0.0))
    {
      intervals.clear();
      return false;
    }
    if (!interval.end.isValid() || interval.end <= interval.start)
    {
      interval.end = interval.start.addMSecs(
          static_cast<qint64>(std::llround(interval.durationSec * 1000.0)));
    }
    intervals.push_back(interval);
  }

  return !intervals.empty();
}

vtkMRMLSequenceNode* findSequenceForProxy(
    vtkMRMLScene* scene,
    vtkMRMLScalarVolumeNode* proxyVolume,
    vtkMRMLSequenceBrowserNode** browserOut = nullptr)
{
  if (browserOut)
  {
    *browserOut = nullptr;
  }
  if (!scene || !proxyVolume)
  {
    return nullptr;
  }

  for (int i = 0; i < scene->GetNumberOfNodesByClass("vtkMRMLSequenceBrowserNode"); ++i)
  {
    vtkMRMLSequenceBrowserNode* browser = vtkMRMLSequenceBrowserNode::SafeDownCast(
        scene->GetNthNodeByClass(i, "vtkMRMLSequenceBrowserNode"));
    if (!browser)
    {
      continue;
    }
    vtkMRMLSequenceNode* sequence = browser->GetMasterSequenceNode();
    if (!sequence || browser->GetProxyNode(sequence) != proxyVolume)
    {
      continue;
    }
    if (browserOut)
    {
      *browserOut = browser;
    }
    return sequence;
  }
  return nullptr;
}

std::string exactSegmentIDForName(
    vtkMRMLSegmentationNode* segmentationNode,
    const QString& exactName)
{
  if (!segmentationNode || !segmentationNode->GetSegmentation())
  {
    return std::string();
  }

  vtkSegmentation* segmentation = segmentationNode->GetSegmentation();
  const std::vector<std::string> ids = segmentation->GetSegmentIDs();
  for (const std::string& id : ids)
  {
    vtkSegment* segment = segmentation->GetSegment(id);
    if (!segment || !segment->GetName())
    {
      continue;
    }
    if (QString::fromStdString(segment->GetName()).trimmed() == exactName)
    {
      return id;
    }
  }
  return std::string();
}

double nodeSUVbwFactor(vtkMRMLNode* node, bool& valid)
{
  valid = false;
  if (!node)
  {
    return 0.0;
  }
  bool ok = false;
  const double factor = nodeAttributeText(node, "dPET.SUVbwFactor").toDouble(&ok);
  const QString validText = nodeAttributeText(node, "dPET.SUVbwFactorValid");
  valid = ok && std::isfinite(factor) && factor > 0.0 &&
      (validText == QStringLiteral("1") || std::abs(factor - 1.0) > 1e-12);
  return ok && std::isfinite(factor) ? factor : 0.0;
}

QString normalizedDecayCorrection(vtkMRMLNode* node)
{
  return nodeOrKineticMetadataText(
      node, "DecayCorrection", QStringLiteral("DecayCorrection")).trimmed().toUpper();
}

bool multiAcquisitionSUVbwFactor(
    vtkMRMLNode* petNode,
    vtkMRMLNode* metadataNode,
    const QString& decayCorrection,
    double& factor,
    QString* errorMessage = nullptr)
{
  factor = std::numeric_limits<double>::quiet_NaN();

  bool valid = false;
  factor = nodeSUVbwFactor(petNode, valid);
  if (!valid && metadataNode && metadataNode != petNode)
  {
    factor = nodeSUVbwFactor(metadataNode, valid);
  }
  if (valid && std::isfinite(factor) && factor > 0.0)
  {
    return true;
  }

  // ADMIN does not require the administration datetime for quantitative
  // normalization: an ADMIN-referenced activity concentration is paired with
  // the administered (non-decayed) dose.  This intentionally avoids any
  // back-extrapolation through a potentially unreliable administration time.
  if (decayCorrection == QStringLiteral("ADMIN"))
  {
    bool weightOK = false;
    bool doseOK = false;
    const double weightKg = nodeOrKineticMetadataText(
        metadataNode, "PatientWeight", QStringLiteral("PatientWeight")).toDouble(&weightOK);
    const double doseBq = nodeOrKineticMetadataText(
        metadataNode, "RadionuclideTotalDose", QStringLiteral("RadionuclideTotalDose")).toDouble(&doseOK);
    if (weightOK && doseOK && weightKg > 0.0 && doseBq > 0.0)
    {
      factor = weightKg / (doseBq * 0.001); // kg/kBq == g/Bq numerically
      return std::isfinite(factor) && factor > 0.0;
    }
  }

  if (errorMessage)
  {
    *errorMessage = decayCorrection == QStringLiteral("START")
        ? QObject::tr("A START-corrected PET lacks the validated dPETImporter SUVbw factor required for Multi-timepoint normalization. Re-import it with the current dPETImporter.")
        : QObject::tr("Could not establish a valid SUVbw factor for this PET acquisition.");
  }
  return false;
}

void scaleVoxelStatistics(VoxelStatistics& stats, double scale)
{
  if (stats.empty || !std::isfinite(scale))
  {
    return;
  }
  stats.mean *= scale;
  stats.median *= scale;
  stats.min *= scale;
  stats.max *= scale;
  stats.stddev *= std::abs(scale);
  stats.q1 *= scale;
  stats.q3 *= scale;
  stats.iqr *= std::abs(scale);
  stats.peak *= scale;
  stats.peakStddev *= std::abs(scale);
}

QString formatAcquisitionTiming(
    const QDateTime& startDT,
    const QDateTime& endDT,
    double durationSec)
{
  if (!startDT.isValid())
  {
    return QString();
  }

  QString result = startDT.toString(QStringLiteral("yyyy-MM-dd HH:mm:ss"));
  if (endDT.isValid() && endDT > startDT)
  {
    const QString endFormat = startDT.date() == endDT.date()
        ? QStringLiteral("HH:mm:ss")
        : QStringLiteral("yyyy-MM-dd HH:mm:ss");
    result += QStringLiteral(" -> ") + endDT.toString(endFormat);
  }

  if (std::isfinite(durationSec) && durationSec > 0.0)
  {
    if (durationSec >= 3600.0)
    {
      result += QObject::tr(" (%1 h)").arg(durationSec / 3600.0, 0, 'f', 2);
    }
    else if (durationSec >= 60.0)
    {
      result += QObject::tr(" (%1 min)").arg(durationSec / 60.0, 0, 'f', 1);
    }
    else
    {
      result += QObject::tr(" (%1 s)").arg(durationSec, 0, 'f', 1);
    }
  }
  return result;
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

double medianValidSigma(const std::vector<double>& sigmas)
{
  std::vector<double> valid;
  valid.reserve(sigmas.size());
  for (double sigma : sigmas)
  {
    if (std::isfinite(sigma) && sigma > 1e-12)
    {
      valid.push_back(sigma);
    }
  }

  if (valid.empty())
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::sort(valid.begin(), valid.end());
  const size_t n = valid.size();
  return (n % 2 == 1)
      ? valid[n / 2]
      : 0.5 * (valid[n / 2 - 1] + valid[n / 2]);
}


double inverseVarianceWeightFromSigma(
    double sigma,
    double fallbackSigma = std::numeric_limits<double>::quiet_NaN())
{
  if (!std::isfinite(sigma) || sigma <= 1e-12)
  {
    sigma = fallbackSigma;
  }

  if (!std::isfinite(sigma) || sigma <= 1e-12)
  {
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

void replaceFittedCurveForPlot(
    double*& fittedCurve,
    size_t fittedFrameCount,
    size_t fullFrameCount,
    const std::vector<double>* extendedPrediction = nullptr)
{
  if (!fittedCurve || fullFrameCount == 0)
  {
    return;
  }

  double* fullCurve = new double[fullFrameCount];
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::fill(fullCurve, fullCurve + fullFrameCount, nan);

  if (extendedPrediction && !extendedPrediction->empty())
  {
    const size_t n = std::min(fullFrameCount, extendedPrediction->size());
    for (size_t i = 0; i < n; ++i)
    {
      fullCurve[i] = (*extendedPrediction)[i];
    }
  }
  else
  {
    const size_t n = std::min(fittedFrameCount, fullFrameCount);
    for (size_t i = 0; i < n; ++i)
    {
      fullCurve[i] = fittedCurve[i];
    }
  }

  delete[] fittedCurve;
  fittedCurve = fullCurve;
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
  bool previewGroupExists(const std::string& groupName);
  void refreshInputFunctionPreviewIfVisible();
  void showCurvePreview(
      const std::string& groupName,
      const QString& title,
      const QString& yAxisTitle,
      const std::vector<PreviewCurve>& curves);
  std::vector<unsigned char>
      MTGAOptimizedSelection;

  // Optional Designer widget. Keep an explicit pointer so this feature also
  // compiles against an older generated ui_qSlicerDynamicPETModuleWidget.h.
  QCheckBox* fengExtrapolationCheckBox{nullptr};

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
  void updateSegmentationAdvancedUI();
  void updateROIModelingAvailability();
  void updateParametricImagingAvailability();
  void updateKineticModelAvailability();
  std::string selectedIFStatistic();
  bool buildCurrentSegmentInputFunction(
      std::vector<double>& values,
      QString* errorMessage = nullptr,
      std::vector<bool>* keepMask = nullptr);

  void updateInputFunctionUI();

  // Image-based / table-based analysis mode.  Table mode intentionally
  // reuses the same ROI fitting, plotting and external-IF pipeline after
  // normalizing workbook data into the active TAC representation.
  void initializeTableBasedUI();
  void setTableBasedMode(bool enabled);
  bool isTableBasedMode() const { return this->tableBasedMode; }
  void initializeMultiTimepointUI();
  void setMultiTimepointMode(bool enabled);
  void scheduleSingleModeAcquisitionRefresh();
  bool isMultiTimepointMode() const { return this->multiTimepointMode; }
  void populateMultiTimepointAcquisitionTable();
  void updateMultiTimepointSelectionStatus();
  void populateMultiTimepointCommonSegmentCheckboxes();
  void syncMultiTimepointSelectedSegments();
  bool prepareMultiTimepointAcquisitions(QString* errorMessage = nullptr);
  bool prepareDynamicMultiTimepointAcquisition(
      PreparedMultiTimepointAcquisition& acquisition,
      std::vector<PreparedMultiTimepointObservation>& observations,
      QString* errorMessage = nullptr);
  bool prepareStaticMultiTimepointAcquisition(
      PreparedMultiTimepointAcquisition& acquisition,
      std::vector<PreparedMultiTimepointObservation>& observations,
      QString* errorMessage = nullptr);
  bool ensureMultiTimepointBinaryRepresentation(
      vtkMRMLSegmentationNode* segmentationNode,
      vtkMRMLScalarVolumeNode* referencePET,
      QString* errorMessage = nullptr,
      double* elapsedMs = nullptr,
      bool* converted = nullptr);
  bool computeMultiTimepointTAC(QString* errorMessage = nullptr);
  void showMultiTimepointSelectionDialog();
  void invalidateMultiTimepointDerivedState();
  void setImageSetupVisible(bool visible);
  void captureActiveTACState(TACModeState& state);
  void restoreActiveTACState(const TACModeState& state);
  void clearActiveTACState();
  bool loadTableWorkbook(const QString& filePath, QString* errorMessage = nullptr);
  void clearTableWorkbook();
  void rebuildTACStatisticUI();
  void updateTableWeightingAvailability();
  void updateTableUnitUI();
  void updateAcquisitionTimingContext(bool logMessage = true);
  double frameTimeShiftForInputSec() const;
  double frameEndForInputSec(size_t frameIndex) const;
  double frameStartForInputSec(size_t frameIndex) const;
  double frameMidForInputSec(size_t frameIndex) const;
  double currentObservedInputStartSec() const;
  void propagateOutputDirectory(const QString& path);
  ActivityUnit selectedTableActivityUnit() const;
  QString selectedTableTimeMode() const;
  double tissueSigmaForWeighting(
      const std::string& segmentID,
      size_t frameIndex,
      const std::string& statistic,
      const VoxelStatistics& stats) const;
  bool plotDistributionSelected() const;
  bool plotStatisticSelected(const QString& statisticID) const;
  void updatePlotDispersionAvailability();
  void enforceDistributionSelection();
  void updateDistributionFrameUI(bool resetRange = false);
  void updateDistributionFrameInfo();
  void refreshDistributionPlotIfActive();
  bool plotROIDistribution(
      const std::string& segmentID,
      QString* errorMessage = nullptr);
  void updateSegmentationFrameUI(bool resetRange = false);
  void updateSegmentationFrameInfo();
  void displaySelectedSegmentationFrame();

  void clearMultiTimepointSegmentationWatchers();
  void setupMultiTimepointSegmentationWatchers();
  void queueMultiTimepointSegmentEdit(
      int acquisitionIndex,
      int frameIndex,
      const std::string& segmentID);
  void processQueuedMultiTimepointSegmentEdits();
  void queueMultiTimepointStructureChange(
      int acquisitionIndex,
      int frameIndex,
      vtkMRMLSegmentationNode* sourceNode,
      const std::string& segmentID,
      SegmentationChangeWatcher::StructureChangeType changeType);
  void processQueuedMultiTimepointStructureChanges();
  int preparedObservationIndexForAcquisitionFrame(
      int acquisitionIndex,
      int frameIndex) const;
  bool recomputePreparedMultiTimepointSegmentObservation(
      int observationIndex,
      const std::string& commonSegmentName,
      QString* errorMessage = nullptr);

  void invalidateInputFunctionResults();

  std::vector<bool>& activeExternalIFKeep();
  const std::vector<bool>& activeExternalIFKeep() const;
  ParentFractionModel selectedParentFractionModel() const;
  void logToPythonConsole(const QString& message) const;
  bool exportFinalInputFunctionCSV(
      const QString& filePath,
      QString* errorMessage = nullptr);
  bool buildProcessedParentFraction(
      double requiredEndTimeSec,
      std::vector<double>& processedTimesSec,
      std::vector<double>& processedValues,
      ParentFractionFitParameters* fitParameters = nullptr,
      QString* processingLabel = nullptr,
      QString* fitSummary = nullptr,
      QString* errorMessage = nullptr);

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
  double initialModelPlasmaIntegralSec(
      const InputFunctionResult& result,
      double endTimeSec);
  double integrateModelPlasmaOverInterval(
      const InputFunctionResult& result,
      double startTimeSec,
      double endTimeSec);
  bool pbrAtTime(
      double timeSec,
      double& pbr,
      QString* errorMessage = nullptr) const;

  bool tableBasedMode{false};
  bool multiTimepointMode{false};
  bool updatingMultiTimepointTable{false};
  bool multiTimepointSelectionValidated{false};
  bool multiTimepointPreparationValid{false};
  bool multiTimepointPreparationRunning{false};
  bool multiTimepointExtractionRunning{false};
  // Prevent synchronous Subject Hierarchy callbacks from rebuilding Single
  // selectors while Multi mode is being entered/exited. Multi preparation and
  // plotting create/remove MRML nodes, so mode teardown must be atomic.
  bool multiTimepointModeTransitionRunning{false};
  // Subject Hierarchy item creation/reparenting can emit several synchronous
  // events before a newly added node has reached its final study parent.
  // Coalesce those events and rebuild the Single-mode selectors once on the
  // next event-loop turn.
  bool subjectHierarchyRefreshQueued{false};
  QSet<QString> multiTimepointCommonSegmentNames;
  QString multiTimepointReferenceMetadataNodeID;
  QString lastMultiTimepointValidationLog;
  std::vector<PreparedMultiTimepointAcquisition> preparedMultiTimepointAcquisitions;
  std::vector<PreparedMultiTimepointObservation> preparedMultiTimepointObservations;
  QPointer<QDialog> multiTimepointSelectionDialog;
  bool tableDataLoaded{false};
  TACModeState imageTACState;
  TACModeState tableTACState;
  int imageIFSourceIndex{0};
  std::string imageIFID;
  int tableIFSourceIndex{1};
  std::string tableIFID;

  QString tableWorkbookPath;
  bool tableFramingExact{false};
  QString tableTimingSummary;
  QVariantMap tableWorkbookMetadata;
  AcquisitionTimingContext acquisitionTiming;
  QString sharedOutputDirectory;
  bool propagatingOutputDirectory{false};
  std::vector<double> tablePlotTimesSec;
  std::vector<TACStatisticOption> tableFitStatistics;
  std::vector<TACStatisticOption> tablePlotStatistics;
  std::map<std::string, std::map<std::string, std::vector<double>>> tableSigma;
  std::string lastPlotSegmentID;
  QLabel* distributionFrameLabel{nullptr};
  QWidget* distributionFrameWidget{nullptr};
  QSlider* distributionFrameSlider{nullptr};
  QLineEdit* distributionFrameInfoEdit{nullptr};
  QLabel* segmentationFrameLabel{nullptr};
  QWidget* segmentationFrameWidget{nullptr};
  QSlider* segmentationFrameSlider{nullptr};
  QLineEdit* segmentationFrameInfoEdit{nullptr};
  bool updatingSegmentationFrameSlider{false};
  // A Single <-> Multi transition represents a new displayed temporal context.
  // Keep the Plot/Distribution frame selector pinned to the first observation
  // the next time TAC-backed distribution data become available.
  bool resetDistributionFrameToFirstPending{false};
  bool syncingVOICheckSelection{false};

  std::vector<vtkSmartPointer<SegmentationChangeWatcher>> multiSegWatchers;
  QTimer* multiSegEditTimer{nullptr};
  std::set<std::tuple<int, int, std::string>> multiDirtySegEdits;
  std::vector<PendingMultiSegStructureChange> pendingMultiSegStructureChanges;
  bool multiSegStructureUpdateQueued{false};
  bool processingMultiSegmentationChanges{false};

  QString externalIFPath;
  std::vector<double> externalIFTimesSec;
  std::vector<double> externalIFConcentrations;
  std::vector<bool> imageExternalIFKeep;
  std::vector<bool> tableExternalIFKeep;
  std::vector<size_t> externalIFPreviewIndexMap;
  std::vector<double> externalIFPreviewTimesSec;
  std::vector<double> externalIFPreviewDisplayValues;
  int externalIFPreviewSelectedIndex{-1};
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
  double multiTimepointReferenceSUVbwFactor{std::numeric_limits<double>::quiet_NaN()};
  QString multiTimepointReferenceDecayCorrection;

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
  void populateTimeBarMTGA(bool resetRange = false);
  void populateTimeBarMTGAImg(bool resetRange = false);
  void updateFitRangeSliders(const InputFunctionResult& result);
  void refreshFitRangeSliderLabels();
  void resetAcquisitionTimingDisplay();
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
      PET = 0,
      CT = 1,
      Union = 2,
      Intersection = 3
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
  bool syncingResultVOISelection{false};

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
    // safely return to the PET default.
    if (!hasCT && sourceMode != 0)
    {
        this->BodySupportSourceImg->blockSignals(true);
        this->BodySupportSourceImg->setCurrentIndex(0);
        this->BodySupportSourceImg->blockSignals(false);
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
            model->item(0)->setEnabled(hasPET); // PET
        if (model->item(1))
            model->item(1)->setEnabled(hasPET && hasCT); // CT
        if (model->item(2))
            model->item(2)->setEnabled(hasPET && hasCT); // Union
        if (model->item(3))
            model->item(3)->setEnabled(hasPET && hasCT); // Intersection
    }

    const bool sourceUsesCT =
        supportEnabled &&
        hasCT &&
        (
            sourceMode == 1 || // CT
            sourceMode == 2 || // Union
            sourceMode == 3    // Intersection
        );

    const bool sourceUsesPET =
        supportEnabled &&
        hasPET &&
        (
            sourceMode == 0 || // PET
            sourceMode == 2 || // Union
            sourceMode == 3    // Intersection
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

bool
qSlicerDynamicPETModuleWidgetPrivate::
previewGroupExists(const std::string& groupName)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    vtkMRMLScene* scene = q->mrmlScene();
    if (!scene)
    {
        return false;
    }
    for (int i = 0; i < scene->GetNumberOfNodes(); ++i)
    {
        vtkMRMLNode* node = scene->GetNthNode(i);
        if (!node)
        {
            continue;
        }
        const char* group = node->GetAttribute("SlicerDynamicPET.IFPreviewGroup");
        if (group && groupName == group &&
            vtkMRMLPlotChartNode::SafeDownCast(node))
        {
            return true;
        }
    }
    return false;
}

void
qSlicerDynamicPETModuleWidgetPrivate::
refreshInputFunctionPreviewIfVisible()
{
    if (!this->previewGroupExists("InputFunctionPreview"))
    {
        return;
    }

    // Replace the already-open preview in place. This is intentionally called
    // only after the user has opened Preview once, so routine parameter edits
    // do not force a plot layout change.
    this->previewInputFunction();
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

    // Draw derived curves first and measured/selectable observations last.
    // VTK therefore renders the source markers on top instead of hiding them
    // underneath continuous plasma/IF curves.
    std::vector<PreviewCurve> orderedCurves = curves;
    std::stable_sort(
        orderedCurves.begin(),
        orderedCurves.end(),
        [](const PreviewCurve& a, const PreviewCurve& b)
        {
            return static_cast<int>(a.sourceObservations) <
                   static_cast<int>(b.sourceObservations);
        });

    for (const PreviewCurve& curve : orderedCurves)
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
        if (curve.sourceObservations)
        {
            seriesNode->SetAttribute(
                "SlicerDynamicPET.SourceObservations",
                "1");
            if (!curve.observationRole.empty())
            {
                seriesNode->SetAttribute(
                    "SlicerDynamicPET.SourceObservationRole",
                    curve.observationRole.c_str());
            }
        }
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
            seriesNode->SetMarkerSize(14.0f);
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

        if (qSlicerApplication::application())
        {
            qSlicerLayoutManager* layoutManager =
                qSlicerApplication::application()->layoutManager();
            qMRMLPlotWidget* plotWidget =
                layoutManager ? layoutManager->plotWidget(0) : nullptr;
            qMRMLPlotView* plotView =
                plotWidget ? plotWidget->plotView() : nullptr;
            if (plotView)
            {
                QObject::connect(
                    plotView,
                    SIGNAL(dataSelected(vtkStringArray*, vtkCollection*)),
                    q,
                    SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)),
                    Qt::UniqueConnection);
            }
        }
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

    this->externalIFPreviewIndexMap.clear();
    this->externalIFPreviewTimesSec.clear();
    this->externalIFPreviewDisplayValues.clear();
    this->externalIFPreviewSelectedIndex = -1;

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
                this->frameEndForInputSec(i));
            sourceValues.push_back(allSourceValues[i]);
        }
    }
    else
    {
        const std::vector<bool>& keep = this->activeExternalIFKeep();
        for (size_t i = 0; i < this->externalIFTimesSec.size(); ++i)
        {
            const bool retained =
                keep.size() == this->externalIFTimesSec.size()
                ? keep[i]
                : true;
            if (!retained)
            {
                continue;
            }

            // An automatically inserted (0 s, 0) anchor belongs to the
            // interpolation definition, not to the measured observations.
            // Keep it in the fitting source, but omit it from the selectable
            // raw-point preview entirely so point indices remain unambiguous.
            if (this->externalIFZeroAnchorAdded && i == 0)
            {
                continue;
            }

            sourceTimes.push_back(this->externalIFTimesSec[i]);
            sourceValues.push_back(this->externalIFConcentrations[i]);
            this->externalIFPreviewIndexMap.push_back(i);
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

    if (!segmentSource)
    {
        this->externalIFPreviewTimesSec = sourceTimes;
        this->externalIFPreviewDisplayValues = sourceDisplayValues;
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
            this->frameEndForInputSec(i));
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

    const std::string previewInterpolation = this->selectedIFInterpolation();
    const double previewEndSec =
        previewFrameCount > 0 ? this->frameEndForInputSec(previewFrameCount - 1) : 0.0;
    const double requestedStep = this->timeStepEdit->text().toDouble();
    const double densePreviewStepSec =
        std::max(0.25, std::min(1.0, requestedStep > 0.0 ? requestedStep : 1.0));

    auto densifyNativeForPreview =
        [&](const std::vector<double>& times,
            const std::vector<double>& values,
            std::vector<double>& denseTimes,
            std::vector<double>& denseValues)
        {
            denseTimes.clear();
            denseValues.clear();
            if (previewInterpolation != "pchip" ||
                times.size() < 2 || times.size() != values.size())
            {
                return false;
            }
            const double endSec = std::min(previewEndSec, times.back());
            if (endSec < times.front())
            {
                return false;
            }
            for (double t = times.front(); t < endSec; t += densePreviewStepSec)
            {
                const double value = this->interpolateInputFunction(
                    times, values, t, previewInterpolation);
                if (std::isfinite(value))
                {
                    denseTimes.push_back(t);
                    denseValues.push_back(value);
                }
            }
            const double lastValue = this->interpolateInputFunction(
                times, values, endSec, previewInterpolation);
            if (std::isfinite(lastValue))
            {
                denseTimes.push_back(endSec);
                denseValues.push_back(lastValue);
            }
            return denseTimes.size() >= 2;
        };

    auto densifyFramesForPreview =
        [&](const std::vector<double>& frameValues,
            std::vector<double>& denseTimes,
            std::vector<double>& denseValues)
        {
            denseTimes.clear();
            denseValues.clear();
            if (previewInterpolation != "pchip" ||
                frameValues.size() < previewFrameCount || previewFrameCount < 2)
            {
                return false;
            }
            const double previewStartSec = this->frameStartForInputSec(0);
            for (double t = previewStartSec; t < previewEndSec; t += densePreviewStepSec)
            {
                const double value = this->evaluateFrameCurve(
                    frameValues, t, previewInterpolation);
                if (std::isfinite(value))
                {
                    denseTimes.push_back(t);
                    denseValues.push_back(value);
                }
            }
            const double lastValue = this->evaluateFrameCurve(
                frameValues, previewEndSec, previewInterpolation);
            if (std::isfinite(lastValue))
            {
                denseTimes.push_back(previewEndSec);
                denseValues.push_back(lastValue);
            }
            return denseTimes.size() >= 2;
        };

    std::vector<PreviewCurve> curves;
    curves.push_back(
        {QObject::tr("Original source"),
         sourceTimes,
         sourceDisplayValues,
         true,
         !segmentSource,
         !segmentSource ? "ExternalIF" : ""});

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
        std::vector<double> processedPreviewTimes =
            result.processedSourcePreviewTimesSec;
        if (segmentSource &&
            processedDisplayValues.size() == previewFrameCount &&
            result.sourceProcessingLabel != QObject::tr("Feng model"))
        {
            processedPreviewTimes = frameTimes;
        }
        QString processedCurveLabel = result.sourceProcessingLabel;
        if (result.fengExtrapolationApplied)
        {
            processedCurveLabel += QObject::tr(" (extrapolated)");
        }
        curves.push_back(
            {processedCurveLabel,
             processedPreviewTimes,
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

        std::vector<double> curveTimes =
            useNative ? result.nativeWholeBloodTimesSec : frameTimes;
        std::vector<double> curveValues = displayValues;
        std::vector<double> denseTimes;
        std::vector<double> denseNativeValues;
        if (useNative
                ? densifyNativeForPreview(result.nativeWholeBloodTimesSec, nativeValues, denseTimes, denseNativeValues)
                : densifyFramesForPreview(frameValues, denseTimes, denseNativeValues))
        {
            std::vector<double> denseDisplayValues;
            if (!toDisplay(denseNativeValues, denseDisplayValues))
            {
                QMessageBox::warning(q, QObject::tr("Input Function"), error);
                return;
            }
            curveTimes = std::move(denseTimes);
            curveValues = std::move(denseDisplayValues);
        }

        curves.push_back(
            {QObject::tr("Total whole blood"),
             curveTimes,
             curveValues,
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

        std::vector<double> curveTimes =
            useNative ? result.nativePlasmaTimesSec : frameTimes;
        std::vector<double> curveValues = displayValues;
        std::vector<double> denseTimes;
        std::vector<double> denseNativeValues;
        if (useNative
                ? densifyNativeForPreview(result.nativePlasmaTimesSec, nativeValues, denseTimes, denseNativeValues)
                : densifyFramesForPreview(frameValues, denseTimes, denseNativeValues))
        {
            std::vector<double> denseDisplayValues;
            if (!toDisplay(denseNativeValues, denseDisplayValues))
            {
                QMessageBox::warning(q, QObject::tr("Input Function"), error);
                return;
            }
            curveTimes = std::move(denseTimes);
            curveValues = std::move(denseDisplayValues);
        }

        curves.push_back(
            {QObject::tr("Total plasma"),
             curveTimes,
             curveValues,
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
            {QObject::tr("Parent plasma"),
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

    std::vector<double> measuredTimes;
    std::vector<double> measuredValues;

    for (size_t i = 0; i < this->parentFractionTimesSec.size(); ++i)
    {
        if (this->parentFractionZeroAnchorAdded && i == 0)
            continue;

        measuredTimes.push_back(this->parentFractionTimesSec[i]);
        measuredValues.push_back(this->parentFractionValues[i]);
    }

    const double requiredEndTimeSec =
        !q->timePoints.empty()
        ? this->frameEndForInputSec(q->timePoints.size() - 1)
        : this->parentFractionTimesSec.back();

    std::vector<double> processedTimes;
    std::vector<double> processedValues;
    ParentFractionFitParameters fitParameters;
    QString processingLabel;
    QString fitSummary;
    QString error;

    if (!this->buildProcessedParentFraction(
            requiredEndTimeSec,
            processedTimes,
            processedValues,
            &fitParameters,
            &processingLabel,
            &fitSummary,
            &error))
    {
        QMessageBox::warning(q, QObject::tr("Parent fraction"), error);
        return;
    }

    QString title = QObject::tr("Parent Fraction");
    if (!processingLabel.isEmpty())
    {
        title += QObject::tr(" - %1").arg(processingLabel);
    }
    if (this->selectedParentFractionModel() != ParentFractionModel::Linear &&
        !measuredTimes.empty() &&
        !processedTimes.empty() &&
        processedTimes.back() > measuredTimes.back() + 1e-6)
    {
        title += QObject::tr(" (extrapolated %1 -> %2 s)")
            .arg(measuredTimes.back(), 0, 'g', 7)
            .arg(processedTimes.back(), 0, 'g', 7);
    }

    std::vector<PreviewCurve> curves;
    curves.push_back(
        {processingLabel,
         processedTimes,
         processedValues,
         false,
         false,
         ""});
    curves.push_back(
        {QObject::tr("Measured parent fraction"),
         measuredTimes,
         measuredValues,
         true,
         false,
         ""});

    this->showCurvePreview(
        "ParentFractionPreview",
        title,
        QObject::tr("Parent fraction"),
        curves);
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
        this->imageExternalIFKeep.clear();
        this->tableExternalIFKeep.clear();
        this->externalIFPreviewIndexMap.clear();
        this->externalIFPreviewTimesSec.clear();
        this->externalIFPreviewDisplayValues.clear();
        this->externalIFPreviewSelectedIndex = -1;
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
    const int index = this->IFInterpolationSelector->currentIndex();
    if (index == 1)
    {
        return "const";
    }
    if (index == 2)
    {
        return "pchip";
    }
    return "linear";
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
    return this->PBIFDomainSelector->currentIndex() == 1
        ? PBIFTemplateDomain::TotalPlasma
        : PBIFTemplateDomain::WholeBlood;
}

ActivityUnit
qSlicerDynamicPETModuleWidgetPrivate::
selectedDisplayActivityUnit() const
{
    switch (this->QuantitativeDisplayUnitSelector->currentIndex())
    {
      case 1: return ActivityUnit::BqPerMl;
      case 2: return ActivityUnit::KBqPerMl;
      case 3: return ActivityUnit::MBqPerMl;
      default: return ActivityUnit::SUVbw;
    }
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

    if (this->tableBasedMode)
    {
        return this->selectedTableActivityUnit();
    }

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
    if (unit == ActivityUnit::MBqPerMl)
    {
        return QObject::tr("MBq/mL");
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

    if (this->multiTimepointMode && this->multiTimepointPreparationValid &&
        std::isfinite(this->multiTimepointReferenceSUVbwFactor) &&
        this->multiTimepointReferenceSUVbwFactor > 0.0)
    {
        factor = this->multiTimepointReferenceSUVbwFactor;
        return true;
    }

    if (this->tableBasedMode)
    {
        if (!this->TableSUVbwFactorEdit)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "A SUVbw conversion factor is required for this table-based unit conversion.");
            }
            return false;
        }

        bool ok = false;
        const double tableFactor =
            this->TableSUVbwFactorEdit->text().trimmed().toDouble(&ok);

        if (!ok || !std::isfinite(tableFactor) || tableFactor <= 0.0)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "Enter a positive SUVbw conversion factor in Table Setup. "
                    "It is interpreted using the same convention as the image-based PET SUVbw factor.");
            }
            return false;
        }

        factor = tableFactor;
        return true;
    }

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
                *errorMessage = this->multiTimepointMode
                    ? QObject::tr(
                        "The selected acquisitions use different validated SUVbw factors. "
                        "For Multi-timepoint analysis, supply an external input function already in SUVbw for now; "
                        "time-dependent Bq/mL <-> SUVbw conversion across separated acquisitions is not implemented yet.")
                    : QObject::tr(
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
    else if (from == ActivityUnit::MBqPerMl)
    {
        bqPerMl = value * 1000000.0;
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
    if (to == ActivityUnit::MBqPerMl)
    {
        converted = bqPerMl / 1000000.0;
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

    const bool hasQuantitativeSource =
        this->tableBasedMode
        ? this->tableDataLoaded
        : (this->multiTimepointMode
            ? (this->multiTimepointPreparationValid &&
               std::isfinite(this->multiTimepointReferenceSUVbwFactor) &&
               this->multiTimepointReferenceSUVbwFactor > 0.0)
            : (q->petID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID &&
               q->sequencePETNode != nullptr));

    this->QuantitativeDisplayUnitLabel->setEnabled(hasQuantitativeSource);
    this->QuantitativeDisplayUnitSelector->setEnabled(hasQuantitativeSource);

    if (!hasQuantitativeSource)
    {
        return;
    }

    const ActivityUnit nativeUnit = this->petStoredActivityUnit();

    // For a newly prepared Multi problem SUVbw is the canonical initial
    // display, but do not overwrite a user's later Bq/kBq/MBq choice every
    // time the UI refreshes. Single/Table retain their established behavior.
    if (!this->multiTimepointMode)
    {
        this->QuantitativeDisplayUnitSelector->blockSignals(true);
        this->QuantitativeDisplayUnitSelector->setCurrentIndex(
            nativeUnit == ActivityUnit::SUVbw ? 0 : 1);
        this->QuantitativeDisplayUnitSelector->blockSignals(false);
    }

    QStandardItemModel* model = qobject_cast<QStandardItemModel*>(
        this->QuantitativeDisplayUnitSelector->model());

    if (model)
    {
        QStandardItem* suvItem = model->item(0);
        QStandardItem* bqItem = model->item(1);
        QStandardItem* kbqItem = model->item(2);
        QStandardItem* mbqItem = model->item(3);

        if (suvItem)
        {
            double factor = 0.0;
            const bool hasSUVFactor = this->getCommonSUVbwFactor(factor, nullptr);
            suvItem->setEnabled(
                nativeUnit == ActivityUnit::SUVbw || hasSUVFactor);
        }
        double factor = 0.0;
        const bool hasSUVFactor = this->getCommonSUVbwFactor(factor, nullptr);
        const bool activityEnabled =
            nativeUnit == ActivityUnit::BqPerMl ||
            nativeUnit == ActivityUnit::KBqPerMl ||
            nativeUnit == ActivityUnit::MBqPerMl ||
            hasSUVFactor;
        if (bqItem) bqItem->setEnabled(activityEnabled);
        if (kbqItem) kbqItem->setEnabled(activityEnabled);
        if (mbqItem) mbqItem->setEnabled(activityEnabled);
    }
}

//-----------------------------------------------------------------------------
ActivityUnit
qSlicerDynamicPETModuleWidgetPrivate::
selectedTableActivityUnit() const
{
    if (!this->TableActivityUnitSelector)
    {
        return ActivityUnit::BqPerMl;
    }

    switch (this->TableActivityUnitSelector->currentIndex())
    {
      case 1:
        return ActivityUnit::KBqPerMl;
      case 2:
        return ActivityUnit::SUVbw;
      default:
        return ActivityUnit::BqPerMl;
    }
}

//-----------------------------------------------------------------------------
QString
qSlicerDynamicPETModuleWidgetPrivate::
selectedTableTimeMode() const
{
    if (!this->TableTimeModeSelector)
    {
        return QStringLiteral("end");
    }

    switch (this->TableTimeModeSelector->currentIndex())
    {
      case 1:
        return QStringLiteral("midpoint");
      case 0:
        return QStringLiteral("auto");
      case 2:
      default:
        return QStringLiteral("end");
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
captureActiveTACState(TACModeState& state)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    state.segmentTACs = q->segmentTACs;
    state.segmentTACsnames = q->segmentTACsnames;
    state.timePoints = q->timePoints;
    state.durations = q->durations;
    state.numberOfTimepoints = q->numberOfTimepoints;
    state.segmentDisplayOrder = this->segmentDisplayOrder;
    state.valid = true;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
restoreActiveTACState(const TACModeState& state)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!state.valid)
    {
        this->clearActiveTACState();
        return;
    }

    q->segmentTACs = state.segmentTACs;
    q->segmentTACsnames = state.segmentTACsnames;
    q->timePoints = state.timePoints;
    q->durations = state.durations;
    q->numberOfTimepoints = state.numberOfTimepoints;
    this->segmentDisplayOrder = state.segmentDisplayOrder;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
clearActiveTACState()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    q->segmentTACs.clear();
    q->segmentTACsnames.clear();
    q->timePoints.clear();
    q->durations.clear();
    q->numberOfTimepoints = 0;
    this->segmentDisplayOrder.clear();
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
initializeMultiTimepointUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    this->multiTimepointMode = false;
    this->MultiTimepointAnalysisCheckBox->setChecked(false);
    this->MultiTimepointSelectionLabel->setVisible(false);
    this->MultiTimepointSelectionRow->setVisible(false);

    // The dialog contents are designed in qSlicerDynamicPETModuleWidget.ui so
    // they remain fully editable in Qt Designer. Only the QDialog shell is
    // created here, then the Designer-defined group is reparented into it.
    this->multiTimepointSelectionDialog = new QDialog(q);
    this->multiTimepointSelectionDialog->setWindowTitle(
        QObject::tr("Select Multi-timepoint Acquisitions"));
    this->multiTimepointSelectionDialog->setModal(true);
    this->multiTimepointSelectionDialog->resize(980, 520);
    QVBoxLayout* dialogLayout = new QVBoxLayout(this->multiTimepointSelectionDialog);
    this->MultiTimepointGroupBox->setParent(this->multiTimepointSelectionDialog);
    this->MultiTimepointGroupBox->setTitle(QObject::tr("Acquisitions"));
    dialogLayout->addWidget(this->MultiTimepointGroupBox);

    QHeaderView* header = this->MultiTimepointAcquisitionTable->horizontalHeader();
    if (header)
    {
        header->setSectionResizeMode(0, QHeaderView::ResizeToContents);
        header->setSectionResizeMode(1, QHeaderView::ResizeToContents);
        header->setSectionResizeMode(2, QHeaderView::Stretch);
        header->setSectionResizeMode(3, QHeaderView::Stretch);
        header->setSectionResizeMode(4, QHeaderView::ResizeToContents);
        header->setSectionResizeMode(5, QHeaderView::ResizeToContents);
    }
    this->MultiTimepointAcquisitionTable->verticalHeader()->setVisible(false);

    QObject::connect(
        this->MultiTimepointAnalysisCheckBox,
        &QCheckBox::toggled,
        q,
        [this](bool checked)
        {
            // The three top-level states are intentionally exclusive:
            // neither checked = Single Image, left = Multi-timepoint Image,
            // right = Table mode.
            if (checked && this->tableBasedMode)
            {
                {
                    QSignalBlocker blocker(this->TableBasedAnalysisCheckBox);
                    this->TableBasedAnalysisCheckBox->setChecked(false);
                }
                this->setTableBasedMode(false);
            }
            this->setMultiTimepointMode(checked);
        });

    QObject::connect(
        this->MultiTimepointSelectionButton,
        &QPushButton::clicked,
        q,
        [this]()
        {
            this->showMultiTimepointSelectionDialog();
        });

    QObject::connect(
        this->MultiTimepointDialogButtonBox,
        &QDialogButtonBox::rejected,
        this->multiTimepointSelectionDialog,
        &QDialog::accept);

    QObject::connect(
        this->MultiTimepointAcquisitionTable,
        &QTableWidget::itemChanged,
        q,
        [this](QTableWidgetItem* item)
        {
            if (this->updatingMultiTimepointTable || !item || item->column() != 0)
            {
                return;
            }
            this->invalidateMultiTimepointDerivedState();
            this->updateMultiTimepointSelectionStatus();
        });
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
showMultiTimepointSelectionDialog()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->multiTimepointMode || !this->multiTimepointSelectionDialog)
    {
        return;
    }

    if (q->patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    {
        QMessageBox::information(
            q,
            QObject::tr("Multi-timepoint acquisitions"),
            QObject::tr("Select a patient first."));
        return;
    }

    this->populateMultiTimepointAcquisitionTable();
    this->multiTimepointSelectionDialog->exec();
    this->updateMultiTimepointSelectionStatus();

    if (this->multiTimepointSelectionValidated)
    {
        QString preparationError;
        if (!this->prepareMultiTimepointAcquisitions(&preparationError))
        {
            this->multiTimepointSelectionValidated = false;
            this->multiTimepointPreparationValid = false;
            this->MultiTimepointStatusLabel->setText(
                QObject::tr("Acquisition preparation failed: %1").arg(preparationError));
            this->logToPythonConsole(
                QObject::tr("[SlicerDynamicPET multi-timepoint PREP] FAILED: %1")
                .arg(preparationError));
            q->enableTACbutton();
            QMessageBox::warning(
                q,
                QObject::tr("Multi-timepoint preparation"),
                preparationError);
        }
        else
        {
            // Preparation creates sequence/representation MRML nodes.  Common
            // ROI widgets are intentionally rebuilt only once preparation has
            // finished, instead of reacting to every Subject Hierarchy event.
            this->populateMultiTimepointCommonSegmentCheckboxes();
            q->enableTACbutton();
        }
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
invalidateMultiTimepointDerivedState()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    // Changing the acquisition set changes the ROI temporal problem, but an
    // external CSV input function is independent patient data and is therefore
    // preserved. Segment-derived IF/TAC/model results cannot be preserved.
    q->clearTACdata();
    q->clearFITdata();
    q->clearFITMTGAdata();

    q->timePoints.clear();
    q->durations.clear();
    q->suvFactors.clear();
    q->numberOfTimepoints = 0;
    q->PET_flatten_values.clear();
    q->MTGAImgOutcomes.clear();
    q->TCMImgOutcomes.clear();
    this->clearMultiTimepointSegmentationWatchers();
    this->multiTimepointSelectionValidated = false;
    this->multiTimepointPreparationValid = false;
    this->preparedMultiTimepointAcquisitions.clear();
    this->preparedMultiTimepointObservations.clear();
    this->multiTimepointReferenceMetadataNodeID.clear();

    if (this->IFSourceSelector->currentIndex() == 0)
    {
        this->resetInputFunctionData(false);
    }
    else
    {
        this->invalidateInputFunctionResults();
        this->updateInputFunctionStatus();
    }

    this->setPostTACEnabled(false);
    this->updateParametricImagingAvailability();
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
scheduleSingleModeAcquisitionRefresh()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    // Selector repopulation restores values with signals blocked. Re-run the
    // Single acquisition callbacks on the next event-loop turn so PET display,
    // segmentation display, and the legacy watcher all match those values.
    QTimer::singleShot(0, q, [this, q]()
    {
        if (this->multiTimepointMode || this->tableBasedMode)
        {
            return;
        }

        const int ctIndex = this->CTSelector ? this->CTSelector->currentIndex() : -1;
        if (ctIndex >= 0 &&
            this->CTSelector->itemData(ctIndex).value<vtkIdType>() !=
                vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        {
            q->onCTChanged(ctIndex);
        }

        const int petIndex = this->PETSelector ? this->PETSelector->currentIndex() : -1;
        if (petIndex >= 0 &&
            this->PETSelector->itemData(petIndex).value<vtkIdType>() !=
                vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        {
            q->onPETChanged(petIndex);

            // The Segmentation frame slider is a display/navigation control.
            // Keep its restored frame 1 synchronized with the PET browser.
            if (q->sequenceBrowserPETNode && q->numberOfTimepoints > 0)
            {
                q->sequenceBrowserPETNode->SetSelectedItemNumber(0);
            }
        }

        const int segIndex = this->SegSelector ? this->SegSelector->currentIndex() : -1;
        if (segIndex >= 0 &&
            this->SegSelector->itemData(segIndex).value<vtkIdType>() !=
                vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        {
            q->onSegChanged(segIndex);
        }
        else
        {
            this->updateSegmentationAdvancedUI();
        }
    });
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
setMultiTimepointMode(bool enabled)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (enabled == this->multiTimepointMode)
    {
        return;
    }

    // clearTACdata()/plot cleanup and selector rebuilding modify MRML/Subject
    // Hierarchy nodes.  Keep the whole mode transition atomic so synchronous
    // hierarchy callbacks cannot execute the opposite mode's refresh logic on
    // half-reset state.
    this->multiTimepointModeTransitionRunning = true;
    this->multiTimepointMode = enabled;
    this->multiTimepointCommonSegmentNames.clear();
    this->invalidateMultiTimepointDerivedState();

    // The single-acquisition PET/segmentation state must never leak into a
    // multi-timepoint analysis. Keep the selected patient, but clear all
    // acquisition-dependent handles and timing state.
    q->sequencePETNode = nullptr;
    q->sequenceBrowserPETNode = nullptr;
    q->segSequenceNode = nullptr;
    if (q->SegWatcher)
    {
        // Do not leave the legacy Single watcher attached while Multi owns
        // segmentation editing. Stale observers were able to receive edits
        // from the old Single context and later survive the return to Single.
        q->SegWatcher->Clear();
        q->SegWatcher->browser = nullptr;
    }
    q->petID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    q->segID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    q->ctID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    q->segmentIDs.clear();
    if (enabled)
    {
        // Do not inherit checked Single-mode segment boxes into the new
        // patient-wide common-ROI selection. The validated common list will
        // be rebuilt after the acquisition dialog is configured.
        this->populateMultiTimepointCommonSegmentCheckboxes();
    }
    q->dPETvalueType.clear();
    this->suvbwFactorValidated = false;
    this->resetAcquisitionTimingDisplay();

    // Single and Multi are distinct temporal display contexts. Do not carry
    // the Plot/Distribution or Advanced Segmentation observation selection
    // across the mode boundary. Distribution may not have TAC data yet, so
    // remember that its first usable refresh must also start at observation 1.
    q->PlotSelectedFrame = -1;
    q->PlotSelectedVOI.clear();
    this->resetDistributionFrameToFirstPending = true;

    if (this->distributionFrameSlider)
    {
        QSignalBlocker blocker(this->distributionFrameSlider);
        this->distributionFrameSlider->setRange(1, 1);
        this->distributionFrameSlider->setValue(1);
        this->distributionFrameSlider->setEnabled(false);
    }
    if (this->distributionFrameInfoEdit)
    {
        this->distributionFrameInfoEdit->clear();
    }

    if (this->segmentationFrameSlider)
    {
        this->updatingSegmentationFrameSlider = true;
        {
            QSignalBlocker blocker(this->segmentationFrameSlider);
            this->segmentationFrameSlider->setRange(1, 1);
            this->segmentationFrameSlider->setValue(1);
            this->segmentationFrameSlider->setEnabled(false);
        }
        this->updatingSegmentationFrameSlider = false;
    }
    if (this->segmentationFrameInfoEdit)
    {
        this->segmentationFrameInfoEdit->clear();
    }

    for (QSlider* slider : {this->timeOffsetSlider, this->timeEndSlider, this->TCMEndSlider,
                            this->timeOffsetSliderImg, this->timeEndSliderImg, this->TCMEndSliderImg})
    {
        if (!slider) continue;
        QSignalBlocker blocker(slider);
        slider->setRange(1, 1);
        slider->setValue(1);
    }
    for (QLineEdit* edit : {this->frameEdit, this->timeSecEdit, this->timeMinEdit,
                            this->timeEndInfoEdit, this->TCMEndInfoEdit,
                            this->frameEditImg, this->timeSecEditImg, this->timeMinEditImg,
                            this->timeEndInfoEditImg, this->TCMEndInfoEditImg})
    {
        if (edit) edit->clear();
    }

    if (enabled)
    {
        // Multi-timepoint selection is patient-wide. There is deliberately no
        // active Study selection in the main workflow.
        q->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
        this->populateMultiTimepointAcquisitionTable();
        if (q->patID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        {
            QTimer::singleShot(0, q, [this]()
            {
                this->showMultiTimepointSelectionDialog();
            });
        }
    }
    else
    {
        if (this->multiTimepointSelectionDialog)
        {
            this->multiTimepointSelectionDialog->close();
        }
        this->updatingMultiTimepointTable = true;
        this->MultiTimepointAcquisitionTable->setRowCount(0);
        this->updatingMultiTimepointTable = false;
        this->MultiTimepointStatusLabel->setText(QObject::tr("No acquisitions selected."));
        this->MultiTimepointSelectionSummaryLabel->setText(QObject::tr("None selected"));
        this->MultiTimepointSelectionButton->setText(QObject::tr("Select acquisitions..."));

        // Rebuild the standard hidden selectors from the retained patient so
        // returning to Single always starts from a clean, predictable state.
        this->populateStudyComboBox(q->patID);
    }

    this->setImageSetupVisible(!this->tableBasedMode);

    // Segmentation editing and voxelwise imaging are intentionally unavailable
    // for multi-timepoint ROI analysis.
    this->SegmentationAdvancedCollapsibleButton->setVisible(
        !this->tableBasedMode && !this->multiTimepointMode);

    const int imagingIndex = this->PlotsTabWidget->indexOf(this->ImagingWidget);
    if (imagingIndex >= 0)
    {
#if QT_VERSION >= QT_VERSION_CHECK(5, 15, 0)
        this->PlotsTabWidget->setTabVisible(
            imagingIndex,
            !this->tableBasedMode && !this->multiTimepointMode);
#else
        this->PlotsTabWidget->setTabEnabled(
            imagingIndex,
            !this->tableBasedMode && !this->multiTimepointMode);
#endif
    }

    this->updateQuantitativeUnitUI();
    this->updateInputFunctionStatus();
    this->updateSegmentationAdvancedUI();
    this->updateParametricImagingAvailability();
    q->enableTACbutton();
    this->multiTimepointModeTransitionRunning = false;

    if (!enabled)
    {
        this->scheduleSingleModeAcquisitionRefresh();
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
populateMultiTimepointAcquisitionTable()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->multiTimepointMode || this->tableBasedMode ||
        this->multiTimepointPreparationRunning || this->multiTimepointExtractionRunning)
    {
        return;
    }

    // Preserve the user's current include/segmentation choices across subject
    // hierarchy refreshes whenever the corresponding PET node still exists.
    struct PreviousChoice
    {
        bool checked{false};
        vtkIdType segmentationItemID{vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID};
    };
    std::map<vtkIdType, PreviousChoice> previousChoices;
    for (int row = 0; row < this->MultiTimepointAcquisitionTable->rowCount(); ++row)
    {
        QTableWidgetItem* petItem = this->MultiTimepointAcquisitionTable->item(row, 2);
        QTableWidgetItem* useItem = this->MultiTimepointAcquisitionTable->item(row, 0);
        if (!petItem || !useItem)
            continue;
        const vtkIdType petItemID = petItem->data(Qt::UserRole).value<vtkIdType>();
        PreviousChoice choice;
        choice.checked = useItem->checkState() == Qt::Checked;
        if (QComboBox* segCombo = qobject_cast<QComboBox*>(
                this->MultiTimepointAcquisitionTable->cellWidget(row, 3)))
        {
            choice.segmentationItemID = segCombo->currentData().value<vtkIdType>();
        }
        previousChoices[petItemID] = choice;
    }

    this->updatingMultiTimepointTable = true;
    this->MultiTimepointAcquisitionTable->setRowCount(0);

    vtkMRMLScene* scene = q->mrmlScene();
    vtkMRMLSubjectHierarchyNode* shNode = scene
        ? vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene)
        : nullptr;

    if (!scene || !shNode ||
        q->patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    {
        this->updatingMultiTimepointTable = false;
        this->MultiTimepointStatusLabel->setText(
            QObject::tr("Choose a patient to list PET acquisitions."));
        return;
    }

    std::vector<MultiTimepointCandidate> candidates;
    std::vector<vtkIdType> studies;
    shNode->GetItemChildren(q->patID, studies);

    for (vtkIdType studyID : studies)
    {
        if (!shNode->HasItemAttribute(studyID, "Level") ||
            shNode->GetItemAttribute(studyID, "Level") != "Study")
        {
            continue;
        }

        const QString studyName = QString::fromStdString(shNode->GetItemName(studyID));
        std::function<void(vtkIdType)> collectPET;
        collectPET = [&](vtkIdType itemID)
        {
            vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(
                shNode->GetItemDataNode(itemID));
            if (petNode)
            {
                const std::string modality = shNode->HasItemAttribute(itemID, "DICOM.Modality")
                    ? shNode->GetItemAttribute(itemID, "DICOM.Modality")
                    : std::string();
                const char* internalAttr = petNode->GetAttribute("SlicerDynamicPET.InternalNode");
                const bool internalNode = internalAttr && std::string(internalAttr) == "1";

                if (modality == "PT" && !internalNode)
                {
                    MultiTimepointCandidate candidate;
                    candidate.studyItemID = studyID;
                    candidate.petItemID = itemID;
                    candidate.studyName = studyName;
                    candidate.petName = QString::fromStdString(shNode->GetItemName(itemID));

                    vtkMRMLSequenceNode* dynamicSequence = nullptr;
                    for (int browserIndex = 0;
                         browserIndex < scene->GetNumberOfNodesByClass("vtkMRMLSequenceBrowserNode");
                         ++browserIndex)
                    {
                        vtkMRMLSequenceBrowserNode* browser = vtkMRMLSequenceBrowserNode::SafeDownCast(
                            scene->GetNthNodeByClass(browserIndex, "vtkMRMLSequenceBrowserNode"));
                        if (!browser) continue;
                        vtkMRMLSequenceNode* master = browser->GetMasterSequenceNode();
                        if (!master || browser->GetProxyNode(master) != petNode) continue;
                        const char* proxyLoadedBy = petNode->GetAttribute("dPETImporter.LoadedBy");
                        const char* seqLoadedBy = master->GetAttribute("dPETImporter.LoadedBy");
                        if ((proxyLoadedBy && std::string(proxyLoadedBy) == "dPETImporterPlugin") ||
                            (seqLoadedBy && std::string(seqLoadedBy) == "dPETImporterPlugin"))
                        {
                            dynamicSequence = master;
                            candidate.dynamic = true;
                            if (master->GetName() && std::string(master->GetName()).size() > 0)
                            {
                                candidate.petName = QString::fromUtf8(master->GetName());
                            }
                            break;
                        }
                    }

                    vtkMRMLNode* metadataNode = dynamicSequence
                        ? static_cast<vtkMRMLNode*>(dynamicSequence)
                        : static_cast<vtkMRMLNode*>(petNode);
                    if (metadataNode && metadataNode->GetID())
                    {
                        candidate.metadataNodeID = QString::fromUtf8(metadataNode->GetID());
                    }

                    const char* kindAttr = metadataNode->GetAttribute("dPET.AcquisitionKind");
                    const QString kind = kindAttr ? QString::fromUtf8(kindAttr).trimmed() : QString();
                    if (kind.compare(QStringLiteral("DYNAMIC"), Qt::CaseInsensitive) == 0)
                    {
                        candidate.dynamic = true;
                    }

                    const char* metadataJson = metadataNode->GetAttribute("dPET.KineticMetadata");
                    const char* metadataSchema = metadataNode->GetAttribute("dPET.KineticMetadataSchemaVersion");
                    candidate.kineticMetadataReady =
                        (metadataJson && std::string(metadataJson).size() > 0) ||
                        (metadataSchema && std::string(metadataSchema).size() > 0);

                    const char* wholeBodyAttr = metadataNode->GetAttribute("dPET.WholeBody");
                    const bool wholeBody = wholeBodyAttr && std::string(wholeBodyAttr) == "1";
                    if (candidate.dynamic)
                        candidate.acquisitionType = QObject::tr("Dynamic");
                    else if (wholeBody)
                        candidate.acquisitionType = QObject::tr("Static / whole-body");
                    else
                        candidate.acquisitionType = QObject::tr("Static");
                    if (!candidate.kineticMetadataReady)
                        candidate.acquisitionType += QObject::tr(" (fallback)");

                    readPersistedKineticTiming(
                        metadataNode,
                        candidate.acquisitionStart,
                        candidate.acquisitionEnd,
                        candidate.durationSec);

                    if (candidate.acquisitionStart.isValid())
                    {
                        candidate.timingText = formatAcquisitionTiming(
                            candidate.acquisitionStart,
                            candidate.acquisitionEnd,
                            candidate.durationSec);
                    }
                    else
                    {
                        candidate.timingText = candidate.kineticMetadataReady
                            ? QObject::tr("Timing incomplete")
                            : QObject::tr("DICOM fallback required");
                    }

                    const char* petLoadedBy = petNode->GetAttribute("dPETImporter.LoadedBy");
                    const char* metadataLoadedBy = metadataNode
                        ? metadataNode->GetAttribute("dPETImporter.LoadedBy") : nullptr;
                    const bool loadedByDPET =
                        (petLoadedBy && std::string(petLoadedBy) == "dPETImporterPlugin") ||
                        (metadataLoadedBy && std::string(metadataLoadedBy) == "dPETImporterPlugin");
                    if (loadedByDPET)
                    {
                        candidates.push_back(candidate);
                    }
                }
            }

            std::vector<vtkIdType> children;
            shNode->GetItemChildren(itemID, children);
            for (vtkIdType childID : children)
            {
                collectPET(childID);
            }
        };
        collectPET(studyID);
    }

    std::stable_sort(
        candidates.begin(), candidates.end(),
        [](const MultiTimepointCandidate& a, const MultiTimepointCandidate& b)
        {
            if (a.acquisitionStart.isValid() != b.acquisitionStart.isValid())
                return a.acquisitionStart.isValid();
            if (a.acquisitionStart.isValid() && b.acquisitionStart.isValid() &&
                a.acquisitionStart != b.acquisitionStart)
                return a.acquisitionStart < b.acquisitionStart;
            const int studyCompare = QString::compare(a.studyName, b.studyName, Qt::CaseInsensitive);
            if (studyCompare != 0) return studyCompare < 0;
            return QString::compare(a.petName, b.petName, Qt::CaseInsensitive) < 0;
        });

    for (const MultiTimepointCandidate& candidate : candidates)
    {
        const int row = this->MultiTimepointAcquisitionTable->rowCount();
        this->MultiTimepointAcquisitionTable->insertRow(row);

        QTableWidgetItem* useItem = new QTableWidgetItem();
        useItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
        useItem->setCheckState(Qt::Unchecked);
        useItem->setData(Qt::UserRole, QVariant::fromValue(candidate.petItemID));
        useItem->setData(Qt::UserRole + 1, QVariant::fromValue(candidate.studyItemID));
        useItem->setData(Qt::UserRole + 2, candidate.dynamic);
        useItem->setData(Qt::UserRole + 3, candidate.kineticMetadataReady);
        useItem->setData(Qt::UserRole + 4, candidate.metadataNodeID);
        this->MultiTimepointAcquisitionTable->setItem(row, 0, useItem);

        QTableWidgetItem* studyItem = new QTableWidgetItem(candidate.studyName);
        this->MultiTimepointAcquisitionTable->setItem(row, 1, studyItem);

        QTableWidgetItem* petItem = new QTableWidgetItem(candidate.petName);
        petItem->setData(Qt::UserRole, QVariant::fromValue(candidate.petItemID));
        this->MultiTimepointAcquisitionTable->setItem(row, 2, petItem);

        QComboBox* segCombo = new QComboBox(this->MultiTimepointAcquisitionTable);
        segCombo->addItem(
            QObject::tr("None"),
            QVariant::fromValue(vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID));

        std::function<void(vtkIdType)> collectSeg;
        collectSeg = [&](vtkIdType itemID)
        {
            vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(
                shNode->GetItemDataNode(itemID));
            if (segNode)
            {
                segCombo->addItem(
                    QString::fromStdString(shNode->GetItemName(itemID)),
                    QVariant::fromValue(itemID));
            }
            std::vector<vtkIdType> children;
            shNode->GetItemChildren(itemID, children);
            for (vtkIdType childID : children)
                collectSeg(childID);
        };
        collectSeg(candidate.studyItemID);

        const auto previousIt = previousChoices.find(candidate.petItemID);
        if (previousIt != previousChoices.end())
        {
            useItem->setCheckState(previousIt->second.checked ? Qt::Checked : Qt::Unchecked);
            const int previousSegIndex = segCombo->findData(
                QVariant::fromValue(previousIt->second.segmentationItemID));
            if (previousSegIndex >= 0)
                segCombo->setCurrentIndex(previousSegIndex);
        }

        QObject::connect(
            segCombo,
            QOverload<int>::of(&QComboBox::currentIndexChanged),
            q,
            [this](int)
            {
                if (this->updatingMultiTimepointTable)
                    return;
                this->invalidateMultiTimepointDerivedState();
                this->updateMultiTimepointSelectionStatus();
            });
        this->MultiTimepointAcquisitionTable->setCellWidget(row, 3, segCombo);

        this->MultiTimepointAcquisitionTable->setItem(
            row, 4, new QTableWidgetItem(candidate.acquisitionType));
        this->MultiTimepointAcquisitionTable->setItem(
            row, 5, new QTableWidgetItem(candidate.timingText));
    }

    this->updatingMultiTimepointTable = false;
    this->updateMultiTimepointSelectionStatus();
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateMultiTimepointSelectionStatus()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->multiTimepointMode)
    {
        return;
    }

    vtkMRMLScene* scene = q->mrmlScene();
    vtkMRMLSubjectHierarchyNode* shNode = scene
        ? vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene)
        : nullptr;

    int selectedCount = 0;
    int missingSegmentationCount = 0;
    int fallbackCount = 0;
    int missingTimingCount = 0;
    int earliestSelectedRow = -1;

    bool duplicateSegmentNames = false;
    QString duplicateSegmentDescription;
    bool haveCommonSegments = false;
    QSet<QString> commonSegmentNames;

    bool missingAdministrationMetadata = false;
    bool injectionMismatch = false;
    bool tracerMismatch = false;
    bool unsupportedQuantitativeType = false;
    bool unsupportedDecayCorrection = false;
    bool missingDecayCorrectionMetadata = false;
    bool invalidSUVFactor = false;
    bool spatialStaticTimingUnsupported = false;

    QDateTime referenceInjection;
    QString referenceRadionuclideCode;
    QString referenceRadiopharmaceuticalCode;
    QString referenceRadiopharmaceuticalName;
    QString referenceDecayCorrection;
    double referenceSUVbwFactor = std::numeric_limits<double>::quiet_NaN();
    double referenceHalfLife = std::numeric_limits<double>::quiet_NaN();
    QStringList validationDiagnostics;

    auto firstNonEmptyReference = [](QString& reference, const QString& value)
    {
        if (reference.isEmpty() && !value.isEmpty())
        {
            reference = value;
        }
    };

    for (int row = 0; row < this->MultiTimepointAcquisitionTable->rowCount(); ++row)
    {
        QTableWidgetItem* useItem = this->MultiTimepointAcquisitionTable->item(row, 0);
        if (!useItem || useItem->checkState() != Qt::Checked)
            continue;

        ++selectedCount;
        if (earliestSelectedRow < 0)
            earliestSelectedRow = row; // table is chronological when timing is known

        const bool metadataReady = useItem->data(Qt::UserRole + 3).toBool();
        if (!metadataReady)
            ++fallbackCount;

        QTableWidgetItem* timingItem = this->MultiTimepointAcquisitionTable->item(row, 5);
        if (!timingItem || timingItem->text().contains("required", Qt::CaseInsensitive) ||
            timingItem->text().contains("incomplete", Qt::CaseInsensitive))
        {
            ++missingTimingCount;
        }

        QComboBox* segCombo = qobject_cast<QComboBox*>(
            this->MultiTimepointAcquisitionTable->cellWidget(row, 3));
        const vtkIdType segItemID = segCombo
            ? segCombo->currentData().value<vtkIdType>()
            : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
        if (segItemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        {
            ++missingSegmentationCount;
        }
        else if (shNode)
        {
            vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(
                shNode->GetItemDataNode(segItemID));
            vtkSegmentation* segmentation = segNode ? segNode->GetSegmentation() : nullptr;
            if (segmentation)
            {
                QSet<QString> namesInThisSegmentation;
                const std::vector<std::string> ids = segmentation->GetSegmentIDs();
                for (const std::string& id : ids)
                {
                    vtkSegment* segment = segmentation->GetSegment(id);
                    if (!segment)
                    {
                        continue;
                    }
                    const QString name = QString::fromStdString(segment->GetName()).trimmed();
                    if (name.isEmpty())
                    {
                        continue;
                    }
                    if (namesInThisSegmentation.contains(name))
                    {
                        duplicateSegmentNames = true;
                        if (duplicateSegmentDescription.isEmpty())
                        {
                            duplicateSegmentDescription = QObject::tr(
                                "Segmentation '%1' contains the duplicate segment name '%2'.")
                                .arg(QString::fromStdString(shNode->GetItemName(segItemID)), name);
                        }
                    }
                    namesInThisSegmentation.insert(name);
                }

                if (!haveCommonSegments)
                {
                    commonSegmentNames = namesInThisSegmentation;
                    haveCommonSegments = true;
                }
                else
                {
                    commonSegmentNames.intersect(namesInThisSegmentation);
                }
            }
        }

        if (!shNode)
        {
            continue;
        }
        const vtkIdType petItemID = useItem->data(Qt::UserRole).value<vtkIdType>();
        vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(
            shNode->GetItemDataNode(petItemID));
        if (!petNode || !metadataReady)
        {
            continue;
        }

        vtkMRMLNode* metadataNode = petNode;
        const QString metadataNodeID = useItem->data(Qt::UserRole + 4).toString();
        if (!metadataNodeID.isEmpty() && scene)
        {
            if (vtkMRMLNode* candidateMetadataNode = scene->GetNodeByID(metadataNodeID.toUtf8().constData()))
            {
                metadataNode = candidateMetadataNode;
            }
        }

        const bool selectedDynamic = useItem->data(Qt::UserRole + 2).toBool();
        if (!selectedDynamic &&
            nodeAttributeText(metadataNode, "dPET.SpatiallyVaryingTiming") == QStringLiteral("1"))
        {
            spatialStaticTimingUnsupported = true;
        }

        // Same-administration validation. Prefer the authoritative metadata node
        // selected while the acquisition table was populated (the sequence node
        // for dPET dynamic PET, the scalar volume for static PET). Fall back to
        // the compact persisted JSON if a convenience MRML attribute is absent.
        const QString injectionText = nodeOrKineticMetadataText(
            metadataNode, "RadionuclideStartDateTime", QStringLiteral("RadionuclideStartDateTime"));
        const QString injectionSource = nodeOrKineticMetadataText(
            metadataNode, "dPET.InjectionDateTimeSource", QStringLiteral("InjectionDateTimeSource"));
        const QString rawStartDT = nodeOrKineticMetadataText(
            metadataNode, "RadiopharmaceuticalStartDateTime", QStringLiteral("RadiopharmaceuticalStartDateTime"));
        const QString rawStartTime = nodeOrKineticMetadataText(
            metadataNode, "RadiopharmaceuticalStartTime", QStringLiteral("RadiopharmaceuticalStartTime"));
        const QString injectionOffset = nodeOrKineticMetadataText(
            metadataNode, "dPET.InjectionToAcquisitionOffsetSec", QStringLiteral("InjectionToAcquisitionOffsetSec"));
        const QString acquisitionStartText = nodeAttributeText(metadataNode, "dPET.AcquisitionStartDateTime");
        const QDateTime injectionDT = parseDICOMDateTimeText(injectionText);
        if (!injectionDT.isValid())
        {
            missingAdministrationMetadata = true;
        }
        else if (!referenceInjection.isValid())
        {
            referenceInjection = injectionDT;
        }
        else if (std::abs(static_cast<double>(referenceInjection.secsTo(injectionDT))) > 120.0)
        {
            injectionMismatch = true;
        }

        const QString petName = this->MultiTimepointAcquisitionTable->item(row, 2)
            ? this->MultiTimepointAcquisitionTable->item(row, 2)->text()
            : QStringLiteral("<unnamed PET>");
        const QString studyName = this->MultiTimepointAcquisitionTable->item(row, 1)
            ? this->MultiTimepointAcquisitionTable->item(row, 1)->text()
            : QStringLiteral("<unnamed study>");
        validationDiagnostics << QObject::tr(
            "  [%1] Study='%2' PET='%3' metadataNode='%4' kind='%5' acquisitionStart='%6' "
            "RadionuclideStartDateTime='%7' source='%8' rawStartDT='%9' rawStartTime='%10' offsetSec='%11'")
            .arg(row + 1)
            .arg(studyName)
            .arg(petName)
            .arg(metadataNode && metadataNode->GetName() ? QString::fromUtf8(metadataNode->GetName()) : QStringLiteral("<none>"))
            .arg(nodeAttributeText(metadataNode, "dPET.AcquisitionKind"))
            .arg(acquisitionStartText)
            .arg(injectionText)
            .arg(injectionSource)
            .arg(rawStartDT)
            .arg(rawStartTime)
            .arg(injectionOffset);

        const QString radionuclideCode = nodeOrKineticMetadataText(
            metadataNode, "dPET.RadionuclideCode", QStringLiteral("RadionuclideCode"));
        const QString radiopharmaceuticalCode = nodeOrKineticMetadataText(
            metadataNode, "dPET.RadiopharmaceuticalCode", QStringLiteral("RadiopharmaceuticalCode"));
        const QString radiopharmaceuticalName = nodeOrKineticMetadataText(
            metadataNode, "RadiopharmaceuticalName", QStringLiteral("RadiopharmaceuticalName"));

        if (!referenceRadionuclideCode.isEmpty() && !radionuclideCode.isEmpty() &&
            referenceRadionuclideCode != radionuclideCode)
        {
            tracerMismatch = true;
        }
        if (!referenceRadiopharmaceuticalCode.isEmpty() && !radiopharmaceuticalCode.isEmpty() &&
            referenceRadiopharmaceuticalCode != radiopharmaceuticalCode)
        {
            tracerMismatch = true;
        }
        if (referenceRadiopharmaceuticalCode.isEmpty() && radiopharmaceuticalCode.isEmpty() &&
            !referenceRadiopharmaceuticalName.isEmpty() && !radiopharmaceuticalName.isEmpty() &&
            QString::compare(referenceRadiopharmaceuticalName, radiopharmaceuticalName,
                             Qt::CaseInsensitive) != 0)
        {
            tracerMismatch = true;
        }
        firstNonEmptyReference(referenceRadionuclideCode, radionuclideCode);
        firstNonEmptyReference(referenceRadiopharmaceuticalCode, radiopharmaceuticalCode);
        firstNonEmptyReference(referenceRadiopharmaceuticalName, radiopharmaceuticalName);

        bool halfLifeOK = false;
        const double halfLife = nodeOrKineticMetadataText(
            metadataNode, "RadionuclideHalfLife", QStringLiteral("RadionuclideHalfLife")).toDouble(&halfLifeOK);
        if (halfLifeOK && halfLife > 0.0)
        {
            if (std::isfinite(referenceHalfLife) && referenceHalfLife > 0.0 &&
                std::abs(halfLife - referenceHalfLife) /
                    std::max(referenceHalfLife, halfLife) > 0.01)
            {
                tracerMismatch = true;
            }
            else if (!std::isfinite(referenceHalfLife))
            {
                referenceHalfLife = halfLife;
            }
        }

        QString valueType = nodeAttributeText(petNode, "dPET.ValueType").toUpper();
        if (valueType.isEmpty())
        {
            valueType = nodeAttributeText(metadataNode, "dPET.ValueType").toUpper();
        }
        if (valueType.isEmpty())
        {
            const QString units = nodeOrKineticMetadataText(
                metadataNode, "Units", QStringLiteral("Units")).toUpper();
            if (units == QStringLiteral("GML"))
                valueType = QStringLiteral("SUVBW");
            else if (units == QStringLiteral("BQML"))
                valueType = QStringLiteral("BQML");
        }
        if (valueType != QStringLiteral("BQML") &&
            valueType != QStringLiteral("SUVBW") &&
            valueType != QStringLiteral("SUV"))
        {
            unsupportedQuantitativeType = true;
        }
        const QString decayCorrection = normalizedDecayCorrection(metadataNode);
        if (decayCorrection.isEmpty())
        {
            missingDecayCorrectionMetadata = true;
        }
        else if (decayCorrection != QStringLiteral("START") &&
                 decayCorrection != QStringLiteral("ADMIN"))
        {
            unsupportedDecayCorrection = true;
        }
        double sourceFactor = std::numeric_limits<double>::quiet_NaN();
        QString factorError;
        if ((decayCorrection == QStringLiteral("START") ||
             decayCorrection == QStringLiteral("ADMIN")) &&
            !multiAcquisitionSUVbwFactor(
                petNode, metadataNode, decayCorrection, sourceFactor, &factorError))
        {
            invalidSUVFactor = true;
        }
        if (!std::isfinite(referenceSUVbwFactor) &&
            std::isfinite(sourceFactor) && sourceFactor > 0.0)
        {
            referenceSUVbwFactor = sourceFactor;
            referenceDecayCorrection = decayCorrection;
        }

        validationDiagnostics.last() += QObject::tr(
            " DecayCorrection='%1' sourceSUVfactor='%2'")
            .arg(decayCorrection.isEmpty() ? QStringLiteral("<missing>") : decayCorrection)
            .arg(std::isfinite(sourceFactor) ? QString::number(sourceFactor, 'g', 12) : QStringLiteral("<invalid>"));
    }

    QString status;
    if (this->MultiTimepointAcquisitionTable->rowCount() == 0)
    {
        status = QObject::tr("No PET acquisitions were found for the selected patient.");
    }
    else if (selectedCount == 0)
    {
        status = QObject::tr("No acquisitions selected.");
    }
    else if (selectedCount < 2)
    {
        status = QObject::tr("Select at least two PET acquisitions for Multi-timepoint analysis.");
    }
    else if (missingTimingCount > 0)
    {
        status = QObject::tr(
            "%1 selected acquisition(s) need timing metadata recovery before chronological validation.")
            .arg(missingTimingCount);
    }
    else if (earliestSelectedRow >= 0)
    {
        QTableWidgetItem* earliestItem = this->MultiTimepointAcquisitionTable->item(
            earliestSelectedRow, 0);
        if (earliestItem && !earliestItem->data(Qt::UserRole + 2).toBool())
        {
            status = QObject::tr("The earliest selected acquisition must be dynamic.");
        }
    }

    if (status.isEmpty() && missingSegmentationCount > 0)
    {
        status = QObject::tr("Assign one segmentation to every selected PET acquisition.");
    }
    if (status.isEmpty() && fallbackCount > 0)
    {
        status = QObject::tr(
            "%1 selected acquisition(s) require the one-time source-DICOM metadata fallback before validation.")
            .arg(fallbackCount);
    }
    if (status.isEmpty() && injectionMismatch)
    {
        status = QObject::tr(
            "The selected PET acquisitions do not appear to belong to the same radiotracer administration (administration times differ by more than 2 minutes).");
    }
    if (status.isEmpty() && tracerMismatch)
    {
        status = QObject::tr(
            "Radiopharmaceutical/radionuclide metadata are inconsistent across the selected acquisitions.");
    }
    if (status.isEmpty() && missingAdministrationMetadata)
    {
        status = QObject::tr(
            "Administration timing is incomplete for at least one selected acquisition; same-injection validation is required before TAC merging.");
    }
    if (status.isEmpty() && unsupportedQuantitativeType)
    {
        status = QObject::tr(
            "Multi-timepoint PET requires dPETImporter quantitative values in BQML or SUVbw.");
    }
    if (status.isEmpty() && missingDecayCorrectionMetadata)
    {
        status = QObject::tr(
            "Decay Correction metadata is missing for at least one selected PET. Multi-timepoint mode requires START or ADMIN.");
    }
    if (status.isEmpty() && unsupportedDecayCorrection)
    {
        status = QObject::tr(
            "Multi-timepoint mode supports only PET with Decay Correction START or ADMIN. NONE, missing, and unknown values are rejected.");
    }
    if (status.isEmpty() && invalidSUVFactor)
    {
        status = QObject::tr(
            "A selected PET lacks a valid SUVbw normalization factor for its DICOM decay convention. Re-import START data with the current dPETImporter; ADMIN requires valid weight and administered dose.");
    }
    if (status.isEmpty() && spatialStaticTimingUnsupported)
    {
        status = QObject::tr(
            "A selected static/whole-body PET has spatially varying acquisition timing. ROI-specific timing must be derived from the stored slice/bed timing map before TAC merging; a single global time will not be invented.");
    }
    if (status.isEmpty() && duplicateSegmentNames)
    {
        status = duplicateSegmentDescription + QObject::tr(
            " Multi-timepoint matching requires unique exact segment names within every segmentation.");
    }
    if (status.isEmpty() && (!haveCommonSegments || commonSegmentNames.isEmpty()))
    {
        status = QObject::tr(
            "No exact segment name is present in every selected segmentation.");
    }
    if (status.isEmpty())
    {
        status = QObject::tr(
            "%1 acquisitions validated: dPETImporter provenance, same administration/tracer, START/ADMIN decay normalization available, and %2 exact common segment(s).")
            .arg(selectedCount)
            .arg(commonSegmentNames.size());
    }

    this->MultiTimepointStatusLabel->setText(status);

    const bool validationPassed =
        !status.isEmpty() && status.contains(QStringLiteral("acquisitions validated"), Qt::CaseInsensitive);
    this->multiTimepointSelectionValidated = validationPassed;
    if (!this->multiTimepointPreparationRunning && !this->multiTimepointExtractionRunning)
    {
        this->multiTimepointPreparationValid = false;
        this->preparedMultiTimepointAcquisitions.clear();
        this->preparedMultiTimepointObservations.clear();
    }
    if (validationPassed)
    {
        this->multiTimepointCommonSegmentNames = commonSegmentNames;
        // Multi uses SUVbw as its canonical scalar bridge regardless of the
        // native per-acquisition storage units.
        q->dPETvalueType = "SUVbw";
        this->multiTimepointReferenceSUVbwFactor = referenceSUVbwFactor;
        this->multiTimepointReferenceDecayCorrection = referenceDecayCorrection;
        this->suvbwFactorValidated =
            std::isfinite(referenceSUVbwFactor) && referenceSUVbwFactor > 0.0;
        QTableWidgetItem* referenceItem = earliestSelectedRow >= 0
            ? this->MultiTimepointAcquisitionTable->item(earliestSelectedRow, 0)
            : nullptr;
        this->multiTimepointReferenceMetadataNodeID = referenceItem
            ? referenceItem->data(Qt::UserRole + 4).toString()
            : QString();
    }
    else
    {
        this->multiTimepointCommonSegmentNames.clear();
        this->multiTimepointReferenceMetadataNodeID.clear();
        q->dPETvalueType.clear();
        this->multiTimepointReferenceSUVbwFactor = std::numeric_limits<double>::quiet_NaN();
        this->multiTimepointReferenceDecayCorrection.clear();
    }

    this->updateQuantitativeUnitUI();
    this->populateMultiTimepointCommonSegmentCheckboxes();
    this->setImageSetupVisible(!this->tableBasedMode);
    q->enableTACbutton();

    if (selectedCount >= 2)
    {
        QString diagnosticMessage;
        if (!validationPassed)
        {
            diagnosticMessage = QObject::tr(
                "[SlicerDynamicPET multi-timepoint] Validation failed: %1\n%2")
                .arg(status, validationDiagnostics.join(QStringLiteral("\n")));
        }
        else
        {
            diagnosticMessage = QObject::tr(
                "[SlicerDynamicPET multi-timepoint] Validation passed. Reference administration='%1'; "
                "reference dynamic decay='%2'; reference SUV factor=%3.")
                .arg(referenceInjection.isValid()
                    ? referenceInjection.toString(QStringLiteral("yyyy-MM-dd HH:mm:ss"))
                    : QStringLiteral("<unavailable>"))
                .arg(referenceDecayCorrection.isEmpty() ? QStringLiteral("<unavailable>") : referenceDecayCorrection)
                .arg(std::isfinite(referenceSUVbwFactor)
                    ? QString::number(referenceSUVbwFactor, 'g', 12)
                    : QStringLiteral("<invalid>"));
        }

        if (!diagnosticMessage.isEmpty() && diagnosticMessage != this->lastMultiTimepointValidationLog)
        {
            this->logToPythonConsole(diagnosticMessage);
            this->lastMultiTimepointValidationLog = diagnosticMessage;
        }
    }
    else
    {
        this->lastMultiTimepointValidationLog.clear();
    }

    if (selectedCount == 0)
    {
        this->MultiTimepointSelectionSummaryLabel->setText(QObject::tr("None selected"));
        this->MultiTimepointSelectionButton->setText(QObject::tr("Select acquisitions..."));
    }
    else
    {
        QString summary = QObject::tr("%1 selected").arg(selectedCount);
        if (missingSegmentationCount > 0)
        {
            summary += QObject::tr(" · %1 segmentation(s) missing")
                .arg(missingSegmentationCount);
        }
        else if (haveCommonSegments)
        {
            summary += QObject::tr(" · %1 common ROI(s)").arg(commonSegmentNames.size());
        }
        this->MultiTimepointSelectionSummaryLabel->setText(summary);
        this->MultiTimepointSelectionButton->setText(QObject::tr("Edit acquisitions..."));
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
populateMultiTimepointCommonSegmentCheckboxes()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->multiTimepointMode)
    {
        return;
    }

    QSet<QString> previouslySelectedNames;
    for (int i = 0; i < this->segmentCheckLayout->count(); ++i)
    {
        QCheckBox* checkbox = qobject_cast<QCheckBox*>(
            this->segmentCheckLayout->itemAt(i)->widget());
        if (checkbox && checkbox->isChecked())
        {
            previouslySelectedNames.insert(
                checkbox->property("SegmentID").toString());
        }
    }

    this->SegmentCheckContents->blockSignals(true);
    QLayoutItem* item = nullptr;
    while ((item = this->segmentCheckLayout->takeAt(0)) != nullptr)
    {
        if (QWidget* widget = item->widget())
        {
            widget->deleteLater();
        }
        delete item;
    }

    q->segmentIDs.clear();
    this->segmentDisplayOrder.clear();

    if (!this->multiTimepointSelectionValidated ||
        this->multiTimepointCommonSegmentNames.isEmpty())
    {
        this->segmentSelectAll->setEnabled(false);
        this->segmentCheckLayout->addStretch();
        this->SegmentCheckContents->blockSignals(false);
        return;
    }

    QStringList names = this->multiTimepointCommonSegmentNames.values();
    std::sort(names.begin(), names.end(), [](const QString& a, const QString& b)
    {
        return QString::compare(a, b, Qt::CaseInsensitive) < 0;
    });

    for (const QString& name : names)
    {
        QCheckBox* checkbox = new QCheckBox(name, this->SegmentCheckContents);
        // In Multi-timepoint mode the stable cross-acquisition key is the exact
        // segment name. Each acquisition resolves this name back to its own
        // local segment ID during extraction.
        checkbox->setProperty("SegmentID", name);
        const bool wasSelected = previouslySelectedNames.contains(name);
        checkbox->setChecked(wasSelected);
        this->segmentCheckLayout->addWidget(checkbox);
        QObject::connect(
            checkbox, SIGNAL(stateChanged(int)),
            q, SLOT(onSegmentsChanged()));
        if (wasSelected)
        {
            q->segmentIDs.push_back(name);
        }
        this->segmentDisplayOrder.push_back(name.toStdString());
    }

    this->segmentCheckLayout->addStretch();
    this->segmentSelectAll->setEnabled(!names.isEmpty());
    this->SegmentCheckContents->blockSignals(false);
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
syncMultiTimepointSelectedSegments()
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!this->multiTimepointMode)
    {
        return;
    }

    q->segmentIDs.clear();
    for (int i = 0; i < this->segmentCheckLayout->count(); ++i)
    {
        QCheckBox* checkbox = qobject_cast<QCheckBox*>(
            this->segmentCheckLayout->itemAt(i)->widget());
        if (checkbox && checkbox->isChecked())
        {
            const QString key = checkbox->property("SegmentID").toString();
            if (!key.isEmpty())
            {
                q->segmentIDs.push_back(key);
            }
        }
    }
    q->enableTACbutton();
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
ensureMultiTimepointBinaryRepresentation(
    vtkMRMLSegmentationNode* segmentationNode,
    vtkMRMLScalarVolumeNode* referencePET,
    QString* errorMessage,
    double* elapsedMs,
    bool* converted)
{
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();

    if (elapsedMs)
    {
        *elapsedMs = 0.0;
    }
    if (converted)
    {
        *converted = false;
    }

    if (!segmentationNode || !segmentationNode->GetSegmentation() ||
        !referencePET || !referencePET->GetImageData())
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Cannot prepare Binary labelmap: segmentation or PET reference is unavailable.");
        }
        return false;
    }

    vtkSegmentation* segmentation = segmentationNode->GetSegmentation();
    const std::string binaryRep =
        vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();

    // Binary is a derived computational representation.  If it is already
    // present then preparation is complete: never ask the converter for a
    // Binary->Binary path.  This also safely accepts legacy scenes in which
    // older Single-mode code made Binary the source representation.
    if (segmentation->ContainsRepresentation(binaryRep))
    {
        if (elapsedMs)
        {
            *elapsedMs = std::chrono::duration<double, std::milli>(
                Clock::now() - start).count();
        }
        return true;
    }

    const std::string sourceRep = segmentation->GetSourceRepresentationName();
    if (sourceRep == binaryRep)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Segmentation '%1' reports Binary labelmap as its source representation, "
                "but no Binary representation is stored. Refusing an invalid Binary-to-Binary conversion.")
                .arg(segmentationNode->GetName()
                    ? QString::fromUtf8(segmentationNode->GetName())
                    : QObject::tr("<unnamed segmentation>"));
        }
        return false;
    }

    segmentationNode->SetReferenceImageGeometryParameterFromVolumeNode(referencePET);
    if (!segmentation->CreateRepresentation(binaryRep) ||
        !segmentation->ContainsRepresentation(binaryRep))
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Could not create the derived Binary labelmap representation for segmentation '%1' "
                "using PET '%2' as reference geometry. Source representation remains '%3'.")
                .arg(segmentationNode->GetName()
                    ? QString::fromUtf8(segmentationNode->GetName())
                    : QObject::tr("<unnamed segmentation>"))
                .arg(referencePET->GetName()
                    ? QString::fromUtf8(referencePET->GetName())
                    : QObject::tr("<unnamed PET>"))
                .arg(QString::fromStdString(sourceRep));
        }
        return false;
    }

    if (converted)
    {
        *converted = true;
    }
    if (elapsedMs)
    {
        *elapsedMs = std::chrono::duration<double, std::milli>(
            Clock::now() - start).count();
    }
    return true;
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
prepareDynamicMultiTimepointAcquisition(
    PreparedMultiTimepointAcquisition& acquisition,
    std::vector<PreparedMultiTimepointObservation>& observations,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    vtkMRMLScene* scene = q->mrmlScene();
    if (!scene)
    {
        if (errorMessage) *errorMessage = QObject::tr("MRML scene is unavailable.");
        return false;
    }

    vtkMRMLScalarVolumeNode* petProxy = vtkMRMLScalarVolumeNode::SafeDownCast(
        scene->GetNodeByID(acquisition.petNodeID.toUtf8().constData()));
    vtkMRMLSegmentationNode* segmentationProxy = vtkMRMLSegmentationNode::SafeDownCast(
        scene->GetNodeByID(acquisition.segmentationNodeID.toUtf8().constData()));
    vtkMRMLNode* metadataNode = scene->GetNodeByID(
        acquisition.metadataNodeID.toUtf8().constData());
    if (!petProxy || !segmentationProxy || !metadataNode)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Dynamic preparation could not recover PET, segmentation, or metadata for '%1'.")
                .arg(acquisition.petName);
        }
        return false;
    }

    vtkMRMLSequenceBrowserNode* browser = nullptr;
    vtkMRMLSequenceNode* petSequence = findSequenceForProxy(scene, petProxy, &browser);
    if (!petSequence || !browser)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Could not resolve the PET sequence/browser for dynamic acquisition '%1'.")
                .arg(acquisition.petName);
        }
        return false;
    }

    acquisition.petSequenceNodeID = QString::fromUtf8(petSequence->GetID());
    acquisition.petBrowserNodeID = QString::fromUtf8(browser->GetID());

    const int frameCount = petSequence->GetNumberOfDataNodes();
    if (frameCount <= 0)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr("Dynamic acquisition '%1' contains no PET frames.")
                .arg(acquisition.petName);
        }
        return false;
    }

    std::vector<MultiTimepointFrameInterval> intervals;
    bool exactTiming = readPersistedDynamicFrameIntervals(metadataNode, intervals) &&
        static_cast<int>(intervals.size()) == frameCount;
    if (!exactTiming)
    {
        intervals.clear();
        QDateTime cursor = acquisition.start;
        for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex)
        {
            vtkMRMLNode* frameNode = petSequence->GetNthDataNode(frameIndex);
            bool durationOK = false;
            const double durationSec = nodeAttributeText(
                frameNode, "dPET.Duration").toDouble(&durationOK);
            if (!durationOK || !std::isfinite(durationSec) || !(durationSec > 0.0))
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "Dynamic acquisition '%1' has incomplete per-frame timing metadata.")
                        .arg(acquisition.petName);
                }
                return false;
            }
            MultiTimepointFrameInterval interval;
            interval.start = cursor;
            interval.durationSec = durationSec;
            interval.end = cursor.addMSecs(
                static_cast<qint64>(std::llround(durationSec * 1000.0)));
            intervals.push_back(interval);
            cursor = interval.end;
        }
    }

    vtkMRMLSequenceNode* segmentationSequence = browser->GetSequenceNode(segmentationProxy);

    // First try the same dRTImporter adoption used by Single mode when the
    // selected segmentation explicitly identifies itself as a dynamic RTSTRUCT.
    const char* importedDynamicAttribute =
        segmentationProxy->GetAttribute("dRTImporter.DynamicRTStruct");
    const bool importedDynamicSegmentation =
        importedDynamicAttribute &&
        QString::fromUtf8(importedDynamicAttribute) == QStringLiteral("1");

    if (!segmentationSequence && importedDynamicSegmentation)
    {
        PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
        const QVariant resultVariant = mainContext.call(
            "DPE_adopt_dynamic_rtstruct",
            QVariantList{
                QString::fromUtf8(segmentationProxy->GetID()),
                QString::fromUtf8(petSequence->GetID()),
                QString::fromUtf8(browser->GetID())});
        const QVariantMap result = resultVariant.toMap();
        if (!result.value("ok").toBool())
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The dynamic segmentation for '%1' could not be synchronized to its PET sequence: %2")
                    .arg(acquisition.petName, result.value("error").toString());
            }
            return false;
        }

        const QByteArray sequenceNodeID =
            result.value("sequence_node_id").toString().toUtf8();
        segmentationSequence = vtkMRMLSequenceNode::SafeDownCast(
            scene->GetNodeByID(sequenceNodeID.constData()));
        if (!segmentationSequence)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The dynamic segmentation for '%1' was synchronized, but its sequence could not be recovered.")
                    .arg(acquisition.petName);
            }
            return false;
        }
    }

    // Match Single-mode preparation for an ordinary/static segmentation
    // assigned to a dynamic PET: prepare Binary once on the proxy, attach a
    // segmentation sequence to the PET browser, then clone that prepared
    // segmentation at every PET index.  Crucially, Binary remains derived;
    // the source representation is not changed.
    if (!segmentationSequence)
    {
        vtkMRMLScalarVolumeNode* firstPET = vtkMRMLScalarVolumeNode::SafeDownCast(
            petSequence->GetNthDataNode(0));
        double preparationMs = 0.0;
        bool converted = false;
        if (!this->ensureMultiTimepointBinaryRepresentation(
                segmentationProxy, firstPET, errorMessage,
                &preparationMs, &converted))
        {
            return false;
        }

        const std::string stableName =
            segmentationProxy->GetName() ? segmentationProxy->GetName() : "Segmentation";

        vtkSmartPointer<vtkMRMLSequenceNode> newSequence =
            vtkSmartPointer<vtkMRMLSequenceNode>::New();
        newSequence->SetName(stableName.c_str());
        newSequence->SetAttribute("SlicerDynamicPET.MultiPreparedStaticSegmentation", "1");
        scene->AddNode(newSequence);

        browser->AddProxyNode(segmentationProxy, newSequence, false);
        browser->SetSaveChanges(newSequence, true);

        for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex)
        {
            const std::string indexValue = petSequence->GetNthIndexValue(frameIndex);
            if (!newSequence->GetDataNodeAtValue(indexValue))
            {
                newSequence->SetDataNodeAtValue(segmentationProxy, indexValue);
            }
        }
        browser->SetOverwriteProxyName(newSequence, false);
        segmentationSequence = newSequence;
        acquisition.segmentationTemporal = false;

        this->logToPythonConsole(QObject::tr(
            "[SlicerDynamicPET multi-timepoint PREP][DYNAMIC] '%1': no PET-linked segmentation sequence existed; "
            "created a %2-frame synchronized sequence from the assigned segmentation. Binary preparation=%3 ms (%4).")
            .arg(acquisition.petName)
            .arg(frameCount)
            .arg(preparationMs, 0, 'f', 3)
            .arg(converted ? QObject::tr("created") : QObject::tr("reused")));
    }
    else
    {
        const char* staticPrepared = segmentationSequence->GetAttribute(
            "SlicerDynamicPET.MultiPreparedStaticSegmentation");
        acquisition.segmentationTemporal = !(staticPrepared && std::string(staticPrepared) == "1");
    }

    acquisition.segmentationSequenceNodeID =
        QString::fromUtf8(segmentationSequence->GetID());

    double binaryPreparationMs = 0.0;
    int binaryConversions = 0;
    int binaryReuses = 0;

    for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex)
    {
        if (q->stopRequested)
        {
            if (errorMessage) *errorMessage = QObject::tr("Multi-timepoint preparation was cancelled.");
            return false;
        }

        const std::string indexValue = petSequence->GetNthIndexValue(frameIndex);
        vtkMRMLScalarVolumeNode* framePET = vtkMRMLScalarVolumeNode::SafeDownCast(
            petSequence->GetDataNodeAtValue(indexValue));
        vtkMRMLSegmentationNode* frameSegmentation = vtkMRMLSegmentationNode::SafeDownCast(
            segmentationSequence->GetDataNodeAtValue(indexValue));
        if (!framePET || !frameSegmentation || !frameSegmentation->GetSegmentation())
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "Dynamic preparation found no matching PET/segmentation item at frame %1 of '%2'.")
                    .arg(frameIndex + 1)
                    .arg(acquisition.petName);
            }
            return false;
        }

        if (q->ProgressBar)
        {
            q->ProgressBar->setFormat(QObject::tr(
                "Preparing frame %1/%2 (%p%)")
                .arg(frameIndex + 1)
                .arg(frameCount));
            q->ProgressBar->setMaximum(frameCount);
            q->ProgressBar->setValue(frameIndex);
        }

        double framePreparationMs = 0.0;
        bool converted = false;
        if (!this->ensureMultiTimepointBinaryRepresentation(
                frameSegmentation, framePET, errorMessage,
                &framePreparationMs, &converted))
        {
            return false;
        }
        binaryPreparationMs += framePreparationMs;
        converted ? ++binaryConversions : ++binaryReuses;

        PreparedMultiTimepointObservation observation;
        observation.frameIndex = frameIndex;
        observation.dynamic = true;
        observation.acquisitionName = acquisition.petName;
        observation.sequenceIndex = QString::fromStdString(indexValue);
        observation.petSequenceNodeID = acquisition.petSequenceNodeID;
        observation.segmentationSequenceNodeID = acquisition.segmentationSequenceNodeID;
        observation.segmentationNodeID = acquisition.segmentationNodeID;
        observation.metadataNodeID = acquisition.metadataNodeID;
        observation.start = intervals[static_cast<size_t>(frameIndex)].start;
        observation.end = intervals[static_cast<size_t>(frameIndex)].end;
        observation.durationSec = intervals[static_cast<size_t>(frameIndex)].durationSec;

        for (const QString& commonName : this->multiTimepointCommonSegmentNames)
        {
            const std::string localID = exactSegmentIDForName(frameSegmentation, commonName);
            if (!localID.empty())
            {
                observation.segmentIDsByName[commonName.toStdString()] = localID;
            }
        }
        observations.push_back(std::move(observation));

        if (frameIndex == 0 || frameIndex == frameCount / 2 || frameIndex + 1 == frameCount)
        {
            vtkSegmentation* segmentation = frameSegmentation->GetSegmentation();
            const std::string binaryRep =
                vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();
            this->logToPythonConsole(QObject::tr(
                "[SlicerDynamicPET DIAG][PREP] dynamic='%1' frame=%2/%3 index='%4' "
                "SEG='%5' source='%6' BinaryAvailable=%7 cachedCommonROIs=%8.")
                .arg(acquisition.petName)
                .arg(frameIndex + 1)
                .arg(frameCount)
                .arg(QString::fromStdString(indexValue))
                .arg(frameSegmentation->GetID()
                    ? QString::fromUtf8(frameSegmentation->GetID()) : QStringLiteral("<none>"))
                .arg(QString::fromStdString(segmentation->GetSourceRepresentationName()))
                .arg(segmentation->ContainsRepresentation(binaryRep) ? 1 : 0)
                .arg(static_cast<int>(observations.back().segmentIDsByName.size())));
        }

        qApp->processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    if (q->ProgressBar)
    {
        q->ProgressBar->setValue(frameCount);
    }

    this->logToPythonConsole(QObject::tr(
        "[SlicerDynamicPET PERF][MULTI PREP] Dynamic '%1': frames=%2; Binary preparation=%3 ms; "
        "conversions=%4; reuses=%5; segmentationSequence='%6'; mode=%7.")
        .arg(acquisition.petName)
        .arg(frameCount)
        .arg(binaryPreparationMs, 0, 'f', 3)
        .arg(binaryConversions)
        .arg(binaryReuses)
        .arg(acquisition.segmentationSequenceNodeID)
        .arg(acquisition.segmentationTemporal
            ? QObject::tr("temporal") : QObject::tr("static-reused")));

    return true;
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
prepareStaticMultiTimepointAcquisition(
    PreparedMultiTimepointAcquisition& acquisition,
    std::vector<PreparedMultiTimepointObservation>& observations,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    vtkMRMLScene* scene = q->mrmlScene();
    vtkMRMLScalarVolumeNode* petNode = scene
        ? vtkMRMLScalarVolumeNode::SafeDownCast(
            scene->GetNodeByID(acquisition.petNodeID.toUtf8().constData()))
        : nullptr;
    vtkMRMLSegmentationNode* segmentationNode = scene
        ? vtkMRMLSegmentationNode::SafeDownCast(
            scene->GetNodeByID(acquisition.segmentationNodeID.toUtf8().constData()))
        : nullptr;
    if (!petNode || !segmentationNode || !segmentationNode->GetSegmentation())
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Static preparation could not recover PET or segmentation for '%1'.")
                .arg(acquisition.petName);
        }
        return false;
    }

    double preparationMs = 0.0;
    bool converted = false;
    if (!this->ensureMultiTimepointBinaryRepresentation(
            segmentationNode, petNode, errorMessage,
            &preparationMs, &converted))
    {
        return false;
    }

    PreparedMultiTimepointObservation observation;
    observation.frameIndex = -1;
    observation.dynamic = false;
    observation.acquisitionName = acquisition.petName;
    observation.petNodeID = acquisition.petNodeID;
    observation.segmentationNodeID = acquisition.segmentationNodeID;
    observation.metadataNodeID = acquisition.metadataNodeID;
    observation.start = acquisition.start;
    observation.end = acquisition.end;
    observation.durationSec = acquisition.durationSec;

    for (const QString& commonName : this->multiTimepointCommonSegmentNames)
    {
        const std::string localID = exactSegmentIDForName(segmentationNode, commonName);
        if (!localID.empty())
        {
            observation.segmentIDsByName[commonName.toStdString()] = localID;
        }
    }
    observations.push_back(std::move(observation));

    vtkSegmentation* segmentation = segmentationNode->GetSegmentation();
    const std::string binaryRep =
        vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();
    this->logToPythonConsole(QObject::tr(
        "[SlicerDynamicPET PERF][MULTI PREP] Static '%1': Binary preparation=%2 ms (%3); "
        "source='%4'; BinaryAvailable=%5; cachedCommonROIs=%6.")
        .arg(acquisition.petName)
        .arg(preparationMs, 0, 'f', 3)
        .arg(converted ? QObject::tr("created") : QObject::tr("reused"))
        .arg(QString::fromStdString(segmentation->GetSourceRepresentationName()))
        .arg(segmentation->ContainsRepresentation(binaryRep) ? 1 : 0)
        .arg(static_cast<int>(observations.back().segmentIDsByName.size())));

    return true;
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
prepareMultiTimepointAcquisitions(QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    using Clock = std::chrono::steady_clock;

    if (!this->multiTimepointMode || !this->multiTimepointSelectionValidated)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Multi-timepoint acquisition validation must pass before preparation.");
        }
        return false;
    }

    vtkMRMLScene* scene = q->mrmlScene();
    vtkMRMLSubjectHierarchyNode* shNode = scene
        ? vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene)
        : nullptr;
    if (!scene || !shNode)
    {
        if (errorMessage) *errorMessage = QObject::tr("MRML scene is unavailable.");
        return false;
    }

    this->clearMultiTimepointSegmentationWatchers();
    this->multiTimepointPreparationValid = false;
    this->preparedMultiTimepointAcquisitions.clear();
    this->preparedMultiTimepointObservations.clear();
    this->multiTimepointPreparationRunning = true;
    q->stopRequested = false;

    auto finishUI = [&]()
    {
        this->multiTimepointPreparationRunning = false;
        if (q->ProgressBar)
        {
            q->ProgressBar->setVisible(false);
            q->ProgressBar->setValue(0);
            q->ProgressBar->setMinimum(0);
            q->ProgressBar->setMaximum(100);
            q->ProgressBar->setFormat("%p%");
        }
        if (q->stopButton)
        {
            q->stopButton->setVisible(false);
        }
    };

    if (q->ProgressBar)
    {
        q->ProgressBar->setMinimum(0);
        q->ProgressBar->setMaximum(100);
        q->ProgressBar->setValue(0);
        q->ProgressBar->setFormat(QObject::tr("Preparing Multi-timepoint acquisitions (%p%)"));
        q->ProgressBar->setVisible(true);
        q->ProgressBar->show();
    }
    if (q->stopButton)
    {
        q->stopButton->setVisible(true);
        q->stopButton->show();
    }

    std::vector<PreparedMultiTimepointAcquisition> acquisitions;
    for (int row = 0; row < this->MultiTimepointAcquisitionTable->rowCount(); ++row)
    {
        QTableWidgetItem* useItem = this->MultiTimepointAcquisitionTable->item(row, 0);
        if (!useItem || useItem->checkState() != Qt::Checked)
        {
            continue;
        }

        const vtkIdType petItemID = useItem->data(Qt::UserRole).value<vtkIdType>();
        QComboBox* segCombo = qobject_cast<QComboBox*>(
            this->MultiTimepointAcquisitionTable->cellWidget(row, 3));
        const vtkIdType segmentationItemID = segCombo
            ? segCombo->currentData().value<vtkIdType>()
            : vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;

        vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(
            shNode->GetItemDataNode(petItemID));
        vtkMRMLSegmentationNode* segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
            shNode->GetItemDataNode(segmentationItemID));
        if (!petNode || !segmentationNode)
        {
            finishUI();
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "A selected acquisition lost its PET or segmentation node during preparation.");
            }
            return false;
        }

        PreparedMultiTimepointAcquisition acquisition;
        acquisition.sourceRow = row;
        acquisition.dynamic = useItem->data(Qt::UserRole + 2).toBool();
        acquisition.petName = this->MultiTimepointAcquisitionTable->item(row, 2)
            ? this->MultiTimepointAcquisitionTable->item(row, 2)->text()
            : QObject::tr("PET acquisition");
        acquisition.petNodeID = QString::fromUtf8(petNode->GetID());
        acquisition.segmentationNodeID = QString::fromUtf8(segmentationNode->GetID());

        const QString metadataNodeID = useItem->data(Qt::UserRole + 4).toString();
        vtkMRMLNode* metadataNode = !metadataNodeID.isEmpty()
            ? scene->GetNodeByID(metadataNodeID.toUtf8().constData())
            : static_cast<vtkMRMLNode*>(petNode);
        if (!metadataNode)
        {
            metadataNode = petNode;
        }
        acquisition.metadataNodeID = QString::fromUtf8(metadataNode->GetID());

        if (!readPersistedKineticTiming(
                metadataNode, acquisition.start, acquisition.end,
                acquisition.durationSec))
        {
            finishUI();
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "Could not recover persisted timing while preparing '%1'.")
                    .arg(acquisition.petName);
            }
            return false;
        }

        acquisition.valueType = nodeAttributeText(petNode, "dPET.ValueType");
        if (acquisition.valueType.isEmpty())
        {
            acquisition.valueType = nodeAttributeText(metadataNode, "dPET.ValueType");
        }
        acquisition.decayCorrection = normalizedDecayCorrection(metadataNode);
        QString factorError;
        if (!multiAcquisitionSUVbwFactor(
                petNode, metadataNode, acquisition.decayCorrection,
                acquisition.sourceSUVbwFactor, &factorError))
        {
            finishUI();
            if (errorMessage)
            {
                *errorMessage = QObject::tr("Could not prepare quantitative normalization for '%1': %2")
                    .arg(acquisition.petName, factorError);
            }
            return false;
        }
        acquisitions.push_back(std::move(acquisition));
    }

    if (acquisitions.size() < 2)
    {
        finishUI();
        if (errorMessage) *errorMessage = QObject::tr("Select at least two acquisitions.");
        return false;
    }

    std::stable_sort(acquisitions.begin(), acquisitions.end(),
        [](const PreparedMultiTimepointAcquisition& a,
           const PreparedMultiTimepointAcquisition& b)
        {
            return a.start < b.start;
        });

    this->multiTimepointReferenceSUVbwFactor = acquisitions.front().sourceSUVbwFactor;
    this->multiTimepointReferenceDecayCorrection = acquisitions.front().decayCorrection;
    this->suvbwFactorValidated =
        std::isfinite(this->multiTimepointReferenceSUVbwFactor) &&
        this->multiTimepointReferenceSUVbwFactor > 0.0;
    if (!this->suvbwFactorValidated)
    {
        finishUI();
        if (errorMessage) *errorMessage = QObject::tr(
            "The earliest dynamic acquisition has no valid reference SUVbw factor.");
        return false;
    }

    QStringList quantitativePreparationSummary;
    for (const PreparedMultiTimepointAcquisition& acquisition : acquisitions)
    {
        quantitativePreparationSummary << QObject::tr(
            "  '%1': stored=%2 DecayCorrection=%3 sourceSUVfactor=%4")
            .arg(acquisition.petName)
            .arg(acquisition.valueType)
            .arg(acquisition.decayCorrection)
            .arg(acquisition.sourceSUVbwFactor, 0, 'g', 12);
    }
    this->logToPythonConsole(QObject::tr(
        "[SlicerDynamicPET multi-timepoint PREP][QUANT] Canonical scalar domain=SUVbw; "
        "reference dynamic SUV factor=%1 (%2). PET voxel data are unchanged; no cross-acquisition PET sequence is created.\n%3")
        .arg(this->multiTimepointReferenceSUVbwFactor, 0, 'g', 12)
        .arg(this->multiTimepointReferenceDecayCorrection)
        .arg(quantitativePreparationSummary.join(QStringLiteral("\n"))));

    const auto totalStart = Clock::now();
    std::vector<PreparedMultiTimepointObservation> observations;

    for (size_t acquisitionIndex = 0;
         acquisitionIndex < acquisitions.size(); ++acquisitionIndex)
    {
        if (q->stopRequested)
        {
            finishUI();
            if (errorMessage) *errorMessage = QObject::tr("Multi-timepoint preparation was cancelled.");
            return false;
        }

        PreparedMultiTimepointAcquisition& acquisition = acquisitions[acquisitionIndex];
        const size_t firstObservation = observations.size();
        QString localError;
        const bool ok = acquisition.dynamic
            ? this->prepareDynamicMultiTimepointAcquisition(
                acquisition, observations, &localError)
            : this->prepareStaticMultiTimepointAcquisition(
                acquisition, observations, &localError);
        if (!ok)
        {
            finishUI();
            if (errorMessage) *errorMessage = localError;
            return false;
        }

        for (size_t observationIndex = firstObservation;
             observationIndex < observations.size(); ++observationIndex)
        {
            observations[observationIndex].acquisitionIndex =
                static_cast<int>(acquisitionIndex);
        }
    }

    std::stable_sort(observations.begin(), observations.end(),
        [](const PreparedMultiTimepointObservation& a,
           const PreparedMultiTimepointObservation& b)
        {
            if (a.start != b.start)
            {
                return a.start < b.start;
            }
            return a.end < b.end;
        });

    this->preparedMultiTimepointAcquisitions = std::move(acquisitions);
    this->preparedMultiTimepointObservations = std::move(observations);
    this->multiTimepointPreparationValid =
        !this->preparedMultiTimepointAcquisitions.empty() &&
        !this->preparedMultiTimepointObservations.empty();

    if (this->multiTimepointPreparationValid)
    {
        this->setupMultiTimepointSegmentationWatchers();
    }

    const double totalMs = std::chrono::duration<double, std::milli>(
        Clock::now() - totalStart).count();

    finishUI();
    if (this->multiTimepointPreparationValid)
    {
        this->MultiTimepointStatusLabel->setText(
            this->MultiTimepointStatusLabel->text() +
            QObject::tr(" Prepared %1 observation(s); TAC can now run without representation conversion.")
                .arg(static_cast<int>(this->preparedMultiTimepointObservations.size())));
    }
    q->enableTACbutton();
    this->updateSegmentationAdvancedUI();
    this->updateSegmentationFrameUI(true);

    this->logToPythonConsole(QObject::tr(
        "[SlicerDynamicPET PERF][MULTI PREP] SUMMARY: total=%1 ms; acquisitions=%2; observations=%3. "
        "Preparation is complete before TAC; Binary is treated only as a derived/cache representation and no Binary-to-Binary conversion is requested.")
        .arg(totalMs, 0, 'f', 3)
        .arg(static_cast<int>(this->preparedMultiTimepointAcquisitions.size()))
        .arg(static_cast<int>(this->preparedMultiTimepointObservations.size())));

    if (errorMessage)
    {
        errorMessage->clear();
    }
    return this->multiTimepointPreparationValid;
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
computeMultiTimepointTAC(QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    using Clock = std::chrono::steady_clock;

    bool extractionStarted = false;
    this->multiTimepointExtractionRunning = true;

    auto resetProgress = [&]()
    {
        if (q->ProgressBar)
        {
            q->ProgressBar->setValue(0);
            q->ProgressBar->setVisible(false);
            q->ProgressBar->setMinimum(0);
            q->ProgressBar->setMaximum(100);
            q->ProgressBar->setFormat("%p%");
        }
        if (q->stopButton)
        {
            q->stopButton->setVisible(false);
        }
    };

    auto fail = [&](const QString& message)
    {
        this->multiTimepointExtractionRunning = false;
        resetProgress();
        if (extractionStarted)
        {
            q->segmentTACs.clear();
            q->segmentTACsnames.clear();
            q->timePoints.clear();
            q->durations.clear();
            q->suvFactors.clear();
            q->numberOfTimepoints = 0;
            this->setPostTACEnabled(false);
        }
        if (errorMessage)
        {
            *errorMessage = message;
        }
        return false;
    };

    if (!this->multiTimepointMode || !this->multiTimepointSelectionValidated)
    {
        return fail(QObject::tr(
            "Multi-timepoint acquisition validation has not passed."));
    }
    if (!this->multiTimepointPreparationValid ||
        this->preparedMultiTimepointAcquisitions.empty() ||
        this->preparedMultiTimepointObservations.empty())
    {
        return fail(QObject::tr(
            "Multi-timepoint acquisitions have not been prepared. Re-open Edit acquisitions and close the validated selection to run preparation."));
    }
    if (q->segmentIDs.empty())
    {
        return fail(QObject::tr("Select at least one common segment."));
    }

    vtkMRMLScene* scene = q->mrmlScene();
    vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(q->logic());
    if (!scene || !logic)
    {
        return fail(QObject::tr("The MRML scene or DynamicPET logic is unavailable."));
    }

    // Freeze the computational problem.  MRML/Qt events are still serviced for
    // progress and cancellation, but they cannot silently change which ROIs or
    // observations are being extracted halfway through the TAC.
    const std::vector<QString> selectedROIs = q->segmentIDs;
    const std::vector<PreparedMultiTimepointAcquisition> acquisitions =
        this->preparedMultiTimepointAcquisitions;
    const std::vector<PreparedMultiTimepointObservation> observations =
        this->preparedMultiTimepointObservations;

    if (selectedROIs.empty() || observations.empty())
    {
        return fail(QObject::tr("The prepared Multi-timepoint problem is empty."));
    }

    const QDateTime referenceStart = acquisitions.front().start;
    if (!referenceStart.isValid())
    {
        return fail(QObject::tr("The first prepared acquisition has no valid start time."));
    }

    q->clearTACdata();
    extractionStarted = true;
    q->timePoints.clear();
    q->durations.clear();
    q->suvFactors.clear();
    q->numberOfTimepoints = 0;
    q->PET_flatten_values.clear();
    this->suvbwFactorValidated = true;

    for (const QString& segmentName : selectedROIs)
    {
        const std::string key = segmentName.toStdString();
        q->segmentTACs[key] = std::vector<VoxelStatistics>();
        q->segmentTACsnames[key] = key;
    }

    // Prepared Multi observations are canonicalized to SUVbw acquisition by
    // acquisition.  Bq-family display/export uses the earliest dynamic
    // acquisition's SUVbw factor as the common scalar reference.
    q->dPETvalueType = "SUVbw";
    if (!std::isfinite(this->multiTimepointReferenceSUVbwFactor) ||
        this->multiTimepointReferenceSUVbwFactor <= 0.0)
    {
        return fail(QObject::tr("The prepared Multi-timepoint reference SUVbw factor is invalid."));
    }

    if (q->ProgressBar)
    {
        q->ProgressBar->setFormat(QObject::tr("Computing prepared Multi-timepoint TAC (%p%)"));
        q->ProgressBar->setMinimum(0);
        q->ProgressBar->setMaximum(static_cast<int>(observations.size()));
        q->ProgressBar->setValue(0);
        q->ProgressBar->setVisible(true);
        q->ProgressBar->show();
    }
    if (q->stopButton)
    {
        q->stopButton->setVisible(true);
        q->stopButton->show();
    }
    q->stopRequested = false;

    const auto totalStart = Clock::now();
    double nodeResolutionMs = 0.0;
    double mergedLabelmapMs = 0.0;
    double voxelStatisticsMs = 0.0;
    double uiEventMs = 0.0;
    size_t mergedLabelmapCalls = 0;
    size_t voxelStatisticsCalls = 0;
    size_t missingSegmentCalls = 0;

    auto elapsedMs = [](const Clock::time_point& start,
                        const Clock::time_point& end)
    {
        return std::chrono::duration<double, std::milli>(end - start).count();
    };

    double previousEndSec = -std::numeric_limits<double>::infinity();
    int gapCount = 0;
    QStringList extractionWarnings;
    QSet<QString> geometryWarningLoggedForAcquisition;

    for (size_t observationIndex = 0;
         observationIndex < observations.size(); ++observationIndex)
    {
        if (q->stopRequested)
        {
            break;
        }

        const PreparedMultiTimepointObservation& observation =
            observations[observationIndex];

        const auto resolveStart = Clock::now();
        vtkMRMLScalarVolumeNode* petVolume = nullptr;
        vtkMRMLSegmentationNode* segmentationNode = nullptr;
        vtkMRMLNode* metadataNode = !observation.metadataNodeID.isEmpty()
            ? scene->GetNodeByID(observation.metadataNodeID.toUtf8().constData())
            : nullptr;

        if (observation.dynamic)
        {
            vtkMRMLSequenceNode* petSequence = vtkMRMLSequenceNode::SafeDownCast(
                scene->GetNodeByID(observation.petSequenceNodeID.toUtf8().constData()));
            vtkMRMLSequenceNode* segmentationSequence = vtkMRMLSequenceNode::SafeDownCast(
                scene->GetNodeByID(observation.segmentationSequenceNodeID.toUtf8().constData()));
            if (petSequence && segmentationSequence)
            {
                const std::string indexValue = observation.sequenceIndex.toStdString();
                petVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
                    petSequence->GetDataNodeAtValue(indexValue));
                segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
                    segmentationSequence->GetDataNodeAtValue(indexValue));
            }
        }
        else
        {
            petVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
                scene->GetNodeByID(observation.petNodeID.toUtf8().constData()));
            segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
                scene->GetNodeByID(observation.segmentationNodeID.toUtf8().constData()));
        }
        nodeResolutionMs += elapsedMs(resolveStart, Clock::now());

        if (!petVolume || !segmentationNode || !segmentationNode->GetSegmentation())
        {
            return fail(QObject::tr(
                "A prepared PET/segmentation observation is no longer available for '%1'. Re-run acquisition preparation.")
                .arg(observation.acquisitionName));
        }

        const double startSec =
            static_cast<double>(referenceStart.msecsTo(observation.start)) / 1000.0;
        const double endSec =
            static_cast<double>(referenceStart.msecsTo(observation.end)) / 1000.0;
        if (!(endSec > startSec))
        {
            return fail(QObject::tr(
                "Prepared observation '%1' has an invalid time interval.")
                .arg(observation.acquisitionName));
        }
        if (startSec < previousEndSec - 1e-3)
        {
            return fail(QObject::tr(
                "Prepared observation intervals overlap near '%1'.")
                .arg(observation.acquisitionName));
        }
        if (std::isfinite(previousEndSec) && startSec > previousEndSec + 1e-3)
        {
            ++gapCount;
        }

        if (q->ProgressBar)
        {
            if (observation.dynamic)
            {
                q->ProgressBar->setFormat(QObject::tr(
                    "TAC frame %1 (%p%)")
                    .arg(observation.frameIndex + 1));
            }
            else
            {
                q->ProgressBar->setFormat(QObject::tr(
                    "TAC static (%p%)"));
            }
        }

        int nonEmptySegmentCount = 0;
        QStringList emptySegments;

        for (size_t roiIndex = 0; roiIndex < selectedROIs.size(); ++roiIndex)
        {
            const QString& commonName = selectedROIs[roiIndex];
            const std::string commonKey = commonName.toStdString();

            VoxelStatistics stats;
            stats.keep = false;
            stats.empty = true;

            const auto idIt = observation.segmentIDsByName.find(commonKey);
            if (idIt == observation.segmentIDsByName.end() || idIt->second.empty())
            {
                ++missingSegmentCalls;
                extractionWarnings << QObject::tr(
                    "Segment '%1' is absent from one prepared observation of '%2'; the observation was marked unavailable.")
                    .arg(commonName, observation.acquisitionName);
            }
            else
            {
                vtkNew<vtkStringArray> segmentArray;
                segmentArray->InsertNextValue(idIt->second);
                vtkSmartPointer<vtkOrientedImageData> labelmap =
                    vtkSmartPointer<vtkOrientedImageData>::New();

                const auto labelmapStart = Clock::now();
                vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(
                    segmentationNode,
                    petVolume,
                    segmentArray,
                    vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
                    labelmap);
                mergedLabelmapMs += elapsedMs(labelmapStart, Clock::now());
                ++mergedLabelmapCalls;

                vtkDataArray* labelScalars = labelmap && labelmap->GetPointData()
                    ? labelmap->GetPointData()->GetScalars() : nullptr;
                vtkImageData* petImage = petVolume->GetImageData();
                vtkDataArray* petScalars = petImage && petImage->GetPointData()
                    ? petImage->GetPointData()->GetScalars() : nullptr;

                if (labelScalars && petScalars &&
                    labelScalars->GetNumberOfTuples() == petScalars->GetNumberOfTuples())
                {
                    const auto statsStart = Clock::now();
                    stats = logic->ComputeVoxelStatistics(petVolume, labelmap, 1);
                    voxelStatisticsMs += elapsedMs(statsStart, Clock::now());
                    ++voxelStatisticsCalls;

                    if (!stats.empty)
                    {
                        const PreparedMultiTimepointAcquisition& sourceAcquisition =
                            acquisitions.at(static_cast<size_t>(observation.acquisitionIndex));
                        const QString sourceType = sourceAcquisition.valueType.trimmed().toUpper();
                        if (sourceType == QStringLiteral("BQML"))
                        {
                            scaleVoxelStatistics(stats, sourceAcquisition.sourceSUVbwFactor);
                        }
                    }
                }

                bool representativeObservation = !observation.dynamic;
                if (observation.dynamic)
                {
                    vtkMRMLSequenceNode* petSequence = vtkMRMLSequenceNode::SafeDownCast(
                        scene->GetNodeByID(observation.petSequenceNodeID.toUtf8().constData()));
                    const int frameCount = petSequence ? petSequence->GetNumberOfDataNodes() : 0;
                    representativeObservation =
                        observation.frameIndex == 0 ||
                        (frameCount > 0 && observation.frameIndex == frameCount / 2) ||
                        (frameCount > 0 && observation.frameIndex + 1 == frameCount);
                }

                if (roiIndex == 0 && representativeObservation)
                {
                    vtkSegmentation* segmentation = segmentationNode->GetSegmentation();
                    const std::string binaryRep =
                        vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();
                    double labelRange[2] = {0.0, 0.0};
                    if (labelmap && labelmap->GetPointData() && labelmap->GetPointData()->GetScalars())
                    {
                        labelmap->GetScalarRange(labelRange);
                    }
                    this->logToPythonConsole(QObject::tr(
                        "[SlicerDynamicPET DIAG][MULTI TAC] acquisition='%1' frame=%2 ROI='%3' localID='%4' "
                        "SEGsource='%5' BinaryAvailable=%6 mergedRange=[%7,%8] stats.count=%9 empty=%10.")
                        .arg(observation.acquisitionName)
                        .arg(observation.dynamic
                            ? QString::number(observation.frameIndex + 1)
                            : QStringLiteral("static"))
                        .arg(commonName)
                        .arg(QString::fromStdString(idIt->second))
                        .arg(QString::fromStdString(segmentation->GetSourceRepresentationName()))
                        .arg(segmentation->ContainsRepresentation(binaryRep) ? 1 : 0)
                        .arg(labelRange[0], 0, 'g', 8)
                        .arg(labelRange[1], 0, 'g', 8)
                        .arg(static_cast<qlonglong>(stats.count))
                        .arg(stats.empty ? 1 : 0));
                }
            }

            if (!stats.empty && stats.count > 0)
            {
                ++nonEmptySegmentCount;
            }
            else
            {
                emptySegments << commonName;
            }
            q->segmentTACs[commonKey].push_back(stats);
        }

        if (nonEmptySegmentCount == 0)
        {
            const QString warning = QObject::tr(
                "'%1' at %2-%3 s produced no voxels for any selected ROI in the prepared PET/segmentation geometry.")
                .arg(observation.acquisitionName)
                .arg(startSec, 0, 'g', 12)
                .arg(endSec, 0, 'g', 12);
            extractionWarnings << warning;
            if (!geometryWarningLoggedForAcquisition.contains(observation.acquisitionName))
            {
                this->logToPythonConsole(
                    QObject::tr("[SlicerDynamicPET multi-timepoint TAC] WARNING: %1")
                    .arg(warning));
                geometryWarningLoggedForAcquisition.insert(observation.acquisitionName);
            }
        }
        else if (!emptySegments.isEmpty())
        {
            extractionWarnings << QObject::tr(
                "'%1' at %2-%3 s: unavailable ROI(s): %4")
                .arg(observation.acquisitionName)
                .arg(startSec, 0, 'g', 12)
                .arg(endSec, 0, 'g', 12)
                .arg(emptySegments.join(QStringLiteral(", ")));
        }

        q->timePoints.push_back(endSec);
        q->durations.push_back(endSec - startSec);

        q->suvFactors.push_back(this->multiTimepointReferenceSUVbwFactor);

        previousEndSec = endSec;
        if (q->ProgressBar)
        {
            q->ProgressBar->setValue(static_cast<int>(observationIndex + 1));
        }

        const auto uiStart = Clock::now();
        qApp->processEvents(QEventLoop::ExcludeUserInputEvents);
        uiEventMs += elapsedMs(uiStart, Clock::now());
    }

    resetProgress();
    this->multiTimepointExtractionRunning = false;

    if (q->stopRequested)
    {
        q->segmentTACs.clear();
        q->segmentTACsnames.clear();
        q->timePoints.clear();
        q->durations.clear();
        q->suvFactors.clear();
        return fail(QObject::tr("Multi-timepoint TAC extraction was cancelled."));
    }

    q->numberOfTimepoints = static_cast<int>(q->timePoints.size());
    if (q->numberOfTimepoints == 0)
    {
        return fail(QObject::tr("No TAC observations were extracted."));
    }

    const double totalMs = elapsedMs(totalStart, Clock::now());
    const double accountedMs =
        nodeResolutionMs + mergedLabelmapMs + voxelStatisticsMs + uiEventMs;
    this->logToPythonConsole(QObject::tr(
        "[SlicerDynamicPET PERF][MULTI TAC] SUMMARY: total=%1 ms; observations=%2; frozenROIs=%3; "
        "merged-labelmap generation=%4 ms across %5 call(s); voxel statistics=%6 ms across %7 call(s); "
        "node resolution=%8 ms; UI/processEvents=%9 ms; missing-segment observations=%10; other/control=%11 ms.")
        .arg(totalMs, 0, 'f', 3)
        .arg(q->numberOfTimepoints)
        .arg(static_cast<int>(selectedROIs.size()))
        .arg(mergedLabelmapMs, 0, 'f', 3)
        .arg(static_cast<qulonglong>(mergedLabelmapCalls))
        .arg(voxelStatisticsMs, 0, 'f', 3)
        .arg(static_cast<qulonglong>(voxelStatisticsCalls))
        .arg(nodeResolutionMs, 0, 'f', 3)
        .arg(uiEventMs, 0, 'f', 3)
        .arg(static_cast<qulonglong>(missingSegmentCalls))
        .arg(std::max(0.0, totalMs - accountedMs), 0, 'f', 3));

    QStringList completelyEmptySegments;
    for (const QString& segmentName : selectedROIs)
    {
        const auto it = q->segmentTACs.find(segmentName.toStdString());
        bool hasValidVoxelObservation = false;
        if (it != q->segmentTACs.end())
        {
            for (const VoxelStatistics& stats : it->second)
            {
                if (!stats.empty && stats.count > 0 && std::isfinite(stats.mean))
                {
                    hasValidVoxelObservation = true;
                    break;
                }
            }
        }
        if (!hasValidVoxelObservation)
        {
            completelyEmptySegments << segmentName;
        }
    }
    if (!completelyEmptySegments.isEmpty())
    {
        return fail(QObject::tr(
            "Prepared Multi-timepoint ROI extraction produced no PET voxels for: %1. "
            "No TAC has been accepted.")
            .arg(completelyEmptySegments.join(QStringLiteral(", "))));
    }

    this->updateAcquisitionTimingContext(false);
    this->updateQuantitativeUnitUI();

    QString message = QObject::tr(
        "[SlicerDynamicPET multi-timepoint] Merged ROI TAC created from prepared provenance: "
        "%1 observations from %2 acquisitions, %3 preserved temporal gap(s), time span 0 -> %4 s. "
        "TAC extraction performed no representation conversion or browser adoption.")
        .arg(q->numberOfTimepoints)
        .arg(static_cast<int>(acquisitions.size()))
        .arg(gapCount)
        .arg(q->timePoints.back(), 0, 'g', 12);
    if (!extractionWarnings.isEmpty())
    {
        extractionWarnings.removeDuplicates();
        message += QStringLiteral("\n") + extractionWarnings.join(QStringLiteral("\n"));
    }
    this->logToPythonConsole(message);

    if (errorMessage)
    {
        errorMessage->clear();
    }
    return true;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
setImageSetupVisible(bool visible)
{
    const bool singleVisible = visible && !this->multiTimepointMode;
    const bool multiVisible = visible && this->multiTimepointMode;

    this->label1->setVisible(visible);
    this->PatSelector->setVisible(visible);

    this->label2->setVisible(singleVisible);
    this->StuSelector->setVisible(singleVisible);
    this->label3->setVisible(singleVisible);
    this->CTSelector->setVisible(singleVisible);
    this->label4->setVisible(singleVisible);
    this->PETSelector->setVisible(singleVisible);
    if (!singleVisible)
    {
        this->AcquisitionTimingStatusLabel->setVisible(false);
    }
    this->labseg->setVisible(singleVisible);
    this->SegSelector->setVisible(singleVisible);

    const bool segmentSelectionVisible =
        singleVisible || (multiVisible && this->multiTimepointSelectionValidated);
    this->labsegments->setVisible(segmentSelectionVisible);
    this->segmentSelectAll->setVisible(segmentSelectionVisible);
    this->SegmentCheckScrollArea->setVisible(segmentSelectionVisible);
    if (multiVisible)
    {
        this->labsegments->setText(QObject::tr("Common segments:"));
    }
    else
    {
        this->labsegments->setText(QObject::tr("Segments:"));
    }

    this->MultiTimepointSelectionLabel->setVisible(multiVisible);
    this->MultiTimepointSelectionRow->setVisible(multiVisible);
    this->TACbutton->setVisible(visible);
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateTableUnitUI()
{
    const bool isSUV =
        this->selectedTableActivityUnit() == ActivityUnit::SUVbw;

    if (this->TableSUVbwFactorLabel)
    {
        this->TableSUVbwFactorLabel->setVisible(isSUV);
    }
    if (this->TableSUVbwFactorEdit)
    {
        this->TableSUVbwFactorEdit->setVisible(isSUV);
    }

    this->updateQuantitativeUnitUI();
}

//-----------------------------------------------------------------------------
double
qSlicerDynamicPETModuleWidgetPrivate::
tissueSigmaForWeighting(
    const std::string& segmentID,
    size_t frameIndex,
    const std::string& statistic,
    const VoxelStatistics& stats) const
{
    if (this->tableBasedMode)
    {
        const auto segmentIt = this->tableSigma.find(segmentID);
        if (segmentIt != this->tableSigma.end())
        {
            const auto statIt = segmentIt->second.find(statistic);
            if (statIt != segmentIt->second.end() &&
                frameIndex < statIt->second.size())
            {
                const double sigma = statIt->second[frameIndex];
                if (std::isfinite(sigma) && sigma > 0.0)
                {
                    return sigma;
                }
            }
        }

        return std::numeric_limits<double>::quiet_NaN();
    }

    return statisticDispersionSigma(stats, statistic);
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateTableWeightingAvailability()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->tableBasedMode)
    {
        auto imageStatisticSupportsWeighting = [](QComboBox* combo)
        {
            if (!combo || combo->currentIndex() < 0)
            {
                return false;
            }
            const std::string stat =
                combo->currentData().toString().toStdString();
            return stat == "Mean" || stat == "Median" || stat == "Peak";
        };

        const bool tcmWeightedAvailable =
            imageStatisticSupportsWeighting(this->StatSelector);
        const bool mtgaWeightedAvailable =
            imageStatisticSupportsWeighting(this->StatSelectorMTGA);

        if (!tcmWeightedAvailable && this->weightFitCheckBox->isChecked())
        {
            this->weightFitCheckBox->setChecked(false);
            this->standardFitCheckBox->setChecked(true);
        }
        if (!mtgaWeightedAvailable && this->weightedFitCheckBox->isChecked())
        {
            this->weightedFitCheckBox->setChecked(false);
        }

        this->weightFitCheckBox->setEnabled(tcmWeightedAvailable);
        this->weightedFitCheckBox->setEnabled(mtgaWeightedAvailable);
        this->weightFitCheckBox->setToolTip(
            tcmWeightedAvailable
            ? QObject::tr("Weighted least squares using normalized inverse-variance proxy weights from spatial dispersion: SD for Mean, IQR/1.349 for Median, and local SUVpeak-region SD for Peak.")
            : QObject::tr("WLS is not enabled for Max because no measurement-uncertainty estimate is derived from a single maximum voxel value."));
        this->weightedFitCheckBox->setToolTip(
            mtgaWeightedAvailable
            ? QObject::tr("Weighted least squares using inverse-variance proxy weights: ROI SD for Mean, IQR/1.349 for Median, and local SUVpeak-region SD for Peak.")
            : QObject::tr("WLS is not enabled for Max because no measurement-uncertainty estimate is derived from a single maximum voxel value."));
        return;
    }

    auto hasUsableSigma =
        [this](const std::string& statistic)
        {
            if (statistic.empty() || this->tableTACState.segmentTACs.empty())
            {
                return false;
            }

            // Match the actual weighting implementation: individual zero or
            // undefined sigma values are treated as neutral weights, not as a
            // reason to disable WLS for the entire workbook.  Require at least
            // one finite positive uncertainty for each ROI so WLS remains
            // meaningful while DynamicPET round-trip workbooks containing an
            // occasional zero-dispersion frame stay usable.
            for (const auto& pair : this->tableTACState.segmentTACs)
            {
                const auto segIt = this->tableSigma.find(pair.first);
                if (segIt == this->tableSigma.end())
                {
                    return false;
                }
                const auto statIt = segIt->second.find(statistic);
                if (statIt == segIt->second.end() ||
                    statIt->second.size() != pair.second.size())
                {
                    return false;
                }

                bool roiHasPositiveSigma = false;
                for (double sigma : statIt->second)
                {
                    if (std::isfinite(sigma) && sigma > 1e-12)
                    {
                        roiHasPositiveSigma = true;
                        break;
                    }
                }
                if (!roiHasPositiveSigma)
                {
                    return false;
                }
            }
            return true;
        };

    const int tcmIndex = this->StatSelector->currentIndex();
    const std::string tcmStat =
        tcmIndex >= 0
        ? this->StatSelector->itemData(tcmIndex).toString().toStdString()
        : std::string();

    const int mtgaIndex = this->StatSelectorMTGA->currentIndex();
    const std::string mtgaStat =
        mtgaIndex >= 0
        ? this->StatSelectorMTGA->itemData(mtgaIndex).toString().toStdString()
        : std::string();

    const bool tcmHasSigma = hasUsableSigma(tcmStat);
    const bool mtgaHasSigma = hasUsableSigma(mtgaStat);

    if (!tcmHasSigma && this->weightFitCheckBox->isChecked())
    {
        this->weightFitCheckBox->setChecked(false);
        this->standardFitCheckBox->setChecked(true);
    }
    if (!mtgaHasSigma && this->weightedFitCheckBox->isChecked())
    {
        this->weightedFitCheckBox->setChecked(false);
    }

    this->weightFitCheckBox->setEnabled(tcmHasSigma);
    this->weightedFitCheckBox->setEnabled(mtgaHasSigma);

    this->weightFitCheckBox->setToolTip(
        tcmHasSigma
        ? QObject::tr("Weighted least squares using imported one-sigma TAC uncertainty. Zero or undefined uncertainty entries are assigned a neutral weight instead of disabling the whole fit.")
        : QObject::tr("Weighted least squares is unavailable because at least one ROI does not provide any positive uncertainty values for the selected TAC statistic."));

    this->weightedFitCheckBox->setToolTip(
        mtgaHasSigma
        ? QObject::tr("Weighted least squares using imported one-sigma TAC uncertainty. Zero or undefined uncertainty entries are assigned a neutral weight instead of disabling the whole fit.")
        : QObject::tr("Weighted least squares is unavailable because at least one ROI does not provide any positive uncertainty values for the selected TAC statistic."));
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
rebuildTACStatisticUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    std::vector<TACStatisticOption> fitOptions;
    std::vector<TACStatisticOption> plotOptions;

    if (this->tableBasedMode)
    {
        fitOptions = this->tableFitStatistics;
        plotOptions = this->tablePlotStatistics;
    }
    else
    {
        fitOptions = {
            {QObject::tr("Mean"), "Mean"},
            {QObject::tr("Median"), "Median"},
            {QObject::tr("Peak"), "Peak"},
            {QObject::tr("Max"), "Max"}
        };
        plotOptions = {
            {QObject::tr("Mean"), "Mean"},
            {QObject::tr("Median"), "Median"},
            {QObject::tr("Peak"), "Peak"},
            {QObject::tr("Min"), "Min"},
            {QObject::tr("Max"), "Max"},
            {QObject::tr("Volume [PET] (cm3)"), "VolumePET"}
        };
        // Distribution is evaluated on one concrete prepared observation.
        // In Single this maps to a PET/segmentation sequence frame; in Multi
        // it maps through prepared provenance to the acquisition-specific PET
        // and segmentation geometry. No cross-acquisition resampling occurs.
        plotOptions.push_back({QObject::tr("Distribution"), "Distribution"});

    }

    auto rebuildCombo =
        [](QComboBox* combo, const std::vector<TACStatisticOption>& options)
        {
            QSignalBlocker blocker(combo);
            combo->clear();
            for (const TACStatisticOption& option : options)
            {
                combo->addItem(option.label, QString::fromStdString(option.id));
            }
            if (combo->count() > 0)
            {
                combo->setCurrentIndex(0);
            }
        };

    rebuildCombo(this->StatSelector, fitOptions);
    rebuildCombo(this->StatSelectorMTGA, fitOptions);

    while (QLayoutItem* item = this->PlotStatsCheckLayout->takeAt(0))
    {
        delete item->widget();
        delete item;
    }

    for (const TACStatisticOption& option : plotOptions)
    {
        QCheckBox* cb = new QCheckBox(option.label, this->PlotStatsCheckContents);
        cb->setProperty("StatID", QString::fromStdString(option.id));
        if (option.id == "VolumePET")
        {
            cb->setToolTip(QObject::tr(
                "Physical ROI volume after rasterizing the dynamic segmentation on the PET reference grid."));
        }
        else if (option.id == "Max")
        {
            cb->setToolTip(QObject::tr(
                "Maximum voxel value inside the ROI for each frame."));
        }
        else if (option.id == "Distribution")
        {
            cb->setToolTip(QObject::tr(
                "Image mode only. Plot the voxel-value histogram for one ROI at the frame/observation selected by the Distribution slider. In Multi, the selected observation is resolved through acquisition provenance and remains in its native PET geometry. Selecting Distribution unchecks all other plot metrics and keeps only one ROI."));
        }

        QObject::connect(
            cb, &QCheckBox::toggled, q,
            [this, cb, option](bool checked)
            {
                if (option.id == "Distribution")
                {
                    if (this->distributionFrameLabel)
                    {
                        this->distributionFrameLabel->setVisible(checked);
                    }
                    if (this->distributionFrameWidget)
                    {
                        this->distributionFrameWidget->setVisible(checked);
                    }
                    if (!checked)
                    {
                        this->updatePlotDispersionAvailability();
                        return;
                    }
                    for (int i = 0; i < this->PlotStatsCheckLayout->count(); ++i)
                    {
                        QCheckBox* other = qobject_cast<QCheckBox*>(
                            this->PlotStatsCheckLayout->itemAt(i)->widget());
                        if (other && other != cb)
                        {
                            QSignalBlocker blocker(other);
                            other->setChecked(false);
                        }
                    }
                    this->enforceDistributionSelection();
                    this->updateDistributionFrameUI(false);
                    this->updatePlotDispersionAvailability();
                    this->refreshDistributionPlotIfActive();
                    return;
                }

                if (!checked)
                {
                    this->updatePlotDispersionAvailability();
                    return;
                }

                for (int i = 0; i < this->PlotStatsCheckLayout->count(); ++i)
                {
                    QCheckBox* other = qobject_cast<QCheckBox*>(
                        this->PlotStatsCheckLayout->itemAt(i)->widget());
                    if (other && other != cb &&
                        other->property("StatID").toString() == "Distribution")
                    {
                        QSignalBlocker blocker(other);
                        other->setChecked(false);
                    }
                }
                if (this->distributionFrameLabel)
                {
                    this->distributionFrameLabel->setVisible(false);
                }
                if (this->distributionFrameWidget)
                {
                    this->distributionFrameWidget->setVisible(false);
                }
                this->updatePlotDispersionAvailability();
            });

        this->PlotStatsCheckLayout->addWidget(cb);
    }

    this->updateDistributionFrameUI(false);
    const bool showDistributionFrame = this->plotDistributionSelected() && !this->tableBasedMode;
    if (this->distributionFrameLabel)
    {
        this->distributionFrameLabel->setVisible(showDistributionFrame);
    }
    if (this->distributionFrameWidget)
    {
        this->distributionFrameWidget->setVisible(showDistributionFrame);
    }
    this->updatePlotDispersionAvailability();
    this->updateTableWeightingAvailability();
    q->clearFITdata();
    q->clearFITMTGAdata();
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
plotDistributionSelected() const
{
    for (int i = 0; i < this->PlotStatsCheckLayout->count(); ++i)
    {
        QCheckBox* cb = qobject_cast<QCheckBox*>(
            this->PlotStatsCheckLayout->itemAt(i)->widget());
        if (cb && cb->isChecked() &&
            cb->property("StatID").toString() == "Distribution")
        {
            return true;
        }
    }
    return false;
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
plotStatisticSelected(const QString& statisticID) const
{
    for (int i = 0; i < this->PlotStatsCheckLayout->count(); ++i)
    {
        QCheckBox* cb = qobject_cast<QCheckBox*>(
            this->PlotStatsCheckLayout->itemAt(i)->widget());
        if (cb && cb->isChecked() &&
            cb->property("StatID").toString() == statisticID)
        {
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updatePlotDispersionAvailability()
{
    const bool incompatible =
        this->plotStatisticSelected(QStringLiteral("VolumePET")) ||
        this->plotStatisticSelected(QStringLiteral("Distribution"));

    if (incompatible && this->PlotErrorCheckbox->isChecked())
    {
        QSignalBlocker blocker(this->PlotErrorCheckbox);
        this->PlotErrorCheckbox->setChecked(false);
    }
    this->PlotErrorCheckbox->setEnabled(!incompatible);
    this->PlotErrorCheckbox->setToolTip(
        incompatible
            ? QObject::tr("Dispersion is not defined for PET volume or voxel-distribution plots.")
            : QObject::tr("Show vertical ROI-dispersion bars when available."));
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
enforceDistributionSelection()
{
    if (!this->plotDistributionSelected())
    {
        return;
    }

    QCheckBox* keep = nullptr;
    for (int i = 0; i < this->PlotsegmentCheckLayout->count(); ++i)
    {
        QCheckBox* cb = qobject_cast<QCheckBox*>(
            this->PlotsegmentCheckLayout->itemAt(i)->widget());
        if (!cb)
        {
            continue;
        }
        const std::string id = cb->property("SegmentID").toString().toStdString();
        if (!this->lastPlotSegmentID.empty() && id == this->lastPlotSegmentID)
        {
            keep = cb;
            break;
        }
        if (!keep && cb->isChecked())
        {
            keep = cb;
        }
        if (!keep)
        {
            keep = cb;
        }
    }

    if (!keep)
    {
        return;
    }

    this->lastPlotSegmentID =
        keep->property("SegmentID").toString().toStdString();

    for (int i = 0; i < this->PlotsegmentCheckLayout->count(); ++i)
    {
        QCheckBox* cb = qobject_cast<QCheckBox*>(
            this->PlotsegmentCheckLayout->itemAt(i)->widget());
        if (!cb)
        {
            continue;
        }
        QSignalBlocker blocker(cb);
        cb->setChecked(cb == keep);
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateDistributionFrameInfo()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->distributionFrameSlider || !this->distributionFrameInfoEdit)
    {
        return;
    }

    const int observationIndex = this->distributionFrameSlider->value() - 1;
    if (observationIndex < 0 ||
        observationIndex >= static_cast<int>(q->timePoints.size()))
    {
        this->distributionFrameInfoEdit->clear();
        return;
    }

    const double endSec = this->frameEndForInputSec(
        static_cast<size_t>(observationIndex));
    this->distributionFrameInfoEdit->setText(
        QObject::tr("Frame %1 | end %2 s (%3 min)")
            .arg(observationIndex + 1)
            .arg(endSec, 0, 'f', 2)
            .arg(endSec / 60.0, 0, 'f', 2));
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateDistributionFrameUI(bool resetRange)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!this->distributionFrameSlider)
    {
        return;
    }

    const int count = static_cast<int>(q->timePoints.size());
    const bool available =
        !this->tableBasedMode && count > 0 &&
        (!this->multiTimepointMode ||
         static_cast<int>(this->preparedMultiTimepointObservations.size()) == count);

    this->distributionFrameSlider->setEnabled(available);
    if (!available)
    {
        QSignalBlocker blocker(this->distributionFrameSlider);
        this->distributionFrameSlider->setMinimum(1);
        this->distributionFrameSlider->setMaximum(1);
        this->distributionFrameSlider->setValue(1);
        if (this->distributionFrameInfoEdit)
        {
            this->distributionFrameInfoEdit->clear();
        }
        return;
    }

    const bool uninitializedRange =
        count > 1 && this->distributionFrameSlider->maximum() <= 1;
    int preferred = this->distributionFrameSlider->value();
    if (this->resetDistributionFrameToFirstPending)
    {
        // A mode switch must always start Plot/Distribution navigation from
        // the first observation, even if the previous mode had a selected
        // plot point or the old slider range has not yet been rebuilt.
        preferred = 1;
    }
    else if (resetRange || uninitializedRange || preferred < 1 || preferred > count)
    {
        preferred = count;
        if (!this->multiTimepointMode && q->sequenceBrowserPETNode)
        {
            const int browserFrame = q->sequenceBrowserPETNode->GetSelectedItemNumber();
            if (browserFrame >= 0 && browserFrame < count)
            {
                preferred = browserFrame + 1;
            }
        }
        else if (q->PlotSelectedFrame >= 0 && q->PlotSelectedFrame < count)
        {
            preferred = q->PlotSelectedFrame + 1;
        }
    }

    {
        QSignalBlocker blocker(this->distributionFrameSlider);
        this->distributionFrameSlider->setMinimum(1);
        this->distributionFrameSlider->setMaximum(count);
        this->distributionFrameSlider->setValue(std::clamp(preferred, 1, count));
    }
    this->resetDistributionFrameToFirstPending = false;
    this->updateDistributionFrameInfo();
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
refreshDistributionPlotIfActive()
{
    if (!this->plotDistributionSelected())
    {
        return;
    }

    this->enforceDistributionSelection();
    if (this->lastPlotSegmentID.empty())
    {
        return;
    }

    QString error;
    if (!this->plotROIDistribution(this->lastPlotSegmentID, &error) && !error.isEmpty())
    {
        this->logToPythonConsole(
            QObject::tr("[SlicerDynamicPET Distribution] %1").arg(error));
    }
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
plotROIDistribution(
    const std::string& segmentID,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (this->tableBasedMode)
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr(
                "Voxel distributions require image data and are not available in Table mode.");
        }
        return false;
    }

    const int observationIndex = this->distributionFrameSlider
        ? this->distributionFrameSlider->value() - 1
        : q->PlotSelectedFrame;
    if (observationIndex < 0 || observationIndex >= static_cast<int>(q->timePoints.size()))
    {
        if (errorMessage) *errorMessage = QObject::tr("No PET frame/observation is available.");
        return false;
    }

    vtkMRMLScalarVolumeNode* petVolume = nullptr;
    vtkMRMLSegmentationNode* segmentationNode = nullptr;
    std::string localSegmentID = segmentID;
    ActivityUnit sourceUnit = this->petStoredActivityUnit();
    double sourceSUVbwFactor = std::numeric_limits<double>::quiet_NaN();
    QString sourceDescription;

    if (this->multiTimepointMode)
    {
        if (!this->multiTimepointPreparationValid ||
            observationIndex >= static_cast<int>(this->preparedMultiTimepointObservations.size()))
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "Multi-timepoint provenance is unavailable. Re-open and validate the acquisition selection.");
            }
            return false;
        }

        vtkMRMLScene* scene = q->mrmlScene();
        const PreparedMultiTimepointObservation& observation =
            this->preparedMultiTimepointObservations[static_cast<size_t>(observationIndex)];

        if (observation.dynamic)
        {
            vtkMRMLSequenceNode* petSequence = scene
                ? vtkMRMLSequenceNode::SafeDownCast(
                    scene->GetNodeByID(observation.petSequenceNodeID.toUtf8().constData()))
                : nullptr;
            vtkMRMLSequenceNode* segmentationSequence = scene
                ? vtkMRMLSequenceNode::SafeDownCast(
                    scene->GetNodeByID(observation.segmentationSequenceNodeID.toUtf8().constData()))
                : nullptr;
            if (petSequence && segmentationSequence)
            {
                const std::string indexValue = observation.sequenceIndex.toStdString();
                petVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
                    petSequence->GetDataNodeAtValue(indexValue));
                segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
                    segmentationSequence->GetDataNodeAtValue(indexValue));
            }
            sourceDescription = QObject::tr("observation %1, frame %2")
                .arg(observationIndex + 1)
                .arg(observation.frameIndex + 1);
        }
        else
        {
            petVolume = scene
                ? vtkMRMLScalarVolumeNode::SafeDownCast(
                    scene->GetNodeByID(observation.petNodeID.toUtf8().constData()))
                : nullptr;
            segmentationNode = scene
                ? vtkMRMLSegmentationNode::SafeDownCast(
                    scene->GetNodeByID(observation.segmentationNodeID.toUtf8().constData()))
                : nullptr;
            sourceDescription = QObject::tr("observation %1, static")
                .arg(observationIndex + 1);
        }

        const auto localIDIt = observation.segmentIDsByName.find(segmentID);
        if (localIDIt == observation.segmentIDsByName.end() || localIDIt->second.empty())
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The selected ROI is unavailable at the requested Multi-timepoint observation.");
            }
            return false;
        }
        localSegmentID = localIDIt->second;

        if (observation.acquisitionIndex < 0 ||
            observation.acquisitionIndex >= static_cast<int>(this->preparedMultiTimepointAcquisitions.size()))
        {
            if (errorMessage) *errorMessage = QObject::tr("The prepared acquisition provenance is invalid.");
            return false;
        }
        const PreparedMultiTimepointAcquisition& acquisition =
            this->preparedMultiTimepointAcquisitions[
                static_cast<size_t>(observation.acquisitionIndex)];
        const QString valueType = acquisition.valueType.trimmed().toUpper();
        if (valueType == QStringLiteral("BQML"))
        {
            sourceUnit = ActivityUnit::BqPerMl;
            sourceSUVbwFactor = acquisition.sourceSUVbwFactor;
            if (!std::isfinite(sourceSUVbwFactor) || sourceSUVbwFactor <= 0.0)
            {
                if (errorMessage) *errorMessage = QObject::tr("The source acquisition SUVbw factor is invalid.");
                return false;
            }
        }
        else if (valueType == QStringLiteral("SUVBW"))
        {
            sourceUnit = ActivityUnit::SUVbw;
        }
        else
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The prepared acquisition has unsupported quantitative type '%1'.")
                    .arg(acquisition.valueType);
            }
            return false;
        }
    }
    else
    {
        if (!q->sequencePETNode || !q->segSequenceNode)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr("Dynamic PET or segmentation sequence is unavailable.");
            }
            return false;
        }

        int frameIndex = observationIndex;
        if (frameIndex >= q->sequencePETNode->GetNumberOfDataNodes())
        {
            frameIndex = q->sequencePETNode->GetNumberOfDataNodes() - 1;
        }
        if (frameIndex < 0)
        {
            if (errorMessage) *errorMessage = QObject::tr("No PET frame is available.");
            return false;
        }

        const std::string indexValue = q->sequencePETNode->GetNthIndexValue(frameIndex);
        petVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
            q->sequencePETNode->GetDataNodeAtValue(indexValue));
        segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
            q->segSequenceNode->GetDataNodeAtValue(indexValue));
        sourceDescription = QObject::tr("frame %1").arg(frameIndex + 1);
    }

    if (!petVolume || !segmentationNode || !segmentationNode->GetSegmentation() ||
        !segmentationNode->GetSegmentation()->GetSegment(localSegmentID))
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr("The selected ROI is unavailable at the requested frame/observation.");
        }
        return false;
    }

    vtkNew<vtkStringArray> segmentArray;
    segmentArray->InsertNextValue(localSegmentID);
    vtkSmartPointer<vtkOrientedImageData> labelmap =
        vtkSmartPointer<vtkOrientedImageData>::New();
    vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(
        segmentationNode,
        petVolume,
        segmentArray,
        vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
        labelmap);

    vtkImageData* petImage = petVolume->GetImageData();
    vtkDataArray* petScalars = petImage && petImage->GetPointData()
        ? petImage->GetPointData()->GetScalars() : nullptr;
    vtkDataArray* labelScalars = labelmap && labelmap->GetPointData()
        ? labelmap->GetPointData()->GetScalars() : nullptr;
    if (!petScalars || !labelScalars ||
        petScalars->GetNumberOfTuples() != labelScalars->GetNumberOfTuples())
    {
        if (errorMessage) *errorMessage = QObject::tr("Could not extract ROI voxel values.");
        return false;
    }

    std::vector<double> values;
    values.reserve(static_cast<size_t>(petScalars->GetNumberOfTuples() / 8));
    const ActivityUnit displayUnit = this->selectedDisplayActivityUnit();
    for (vtkIdType i = 0; i < petScalars->GetNumberOfTuples(); ++i)
    {
        if (static_cast<int>(std::llround(labelScalars->GetComponent(i, 0))) != 1)
        {
            continue;
        }
        const double nativeValue = petScalars->GetComponent(i, 0);
        if (!std::isfinite(nativeValue))
        {
            continue;
        }

        double canonicalValue = nativeValue;
        ActivityUnit canonicalUnit = sourceUnit;
        if (this->multiTimepointMode && sourceUnit == ActivityUnit::BqPerMl)
        {
            // Match Multi TAC assembly: source-acquisition BQML is first
            // normalized to canonical SUVbw using that acquisition's own
            // validated factor. Display conversion then uses the common
            // reference-dynamic factor without modifying PET image voxels.
            canonicalValue = nativeValue * sourceSUVbwFactor;
            canonicalUnit = ActivityUnit::SUVbw;
        }

        double converted = nativeValue;
        if (!this->convertActivityValue(
                canonicalValue, canonicalUnit, displayUnit, converted, nullptr) ||
            !std::isfinite(converted))
        {
            continue;
        }
        values.push_back(converted);
    }
    if (values.empty())
    {
        if (errorMessage) *errorMessage = QObject::tr("The ROI contains no finite PET voxel values.");
        return false;
    }

    std::sort(values.begin(), values.end());
    const double minimum = values.front();
    const double maximum = values.back();
    int bins = 1;
    if (maximum > minimum && values.size() > 1)
    {
        const size_t n = values.size();
        const double q1 = values[n / 4];
        const double q3 = values[(3 * n) / 4];
        const double iqr = q3 - q1;
        double binWidth = 2.0 * iqr / std::cbrt(static_cast<double>(n));
        if (!(binWidth > 0.0) || !std::isfinite(binWidth))
        {
            binWidth = (maximum - minimum) / std::sqrt(static_cast<double>(n));
        }
        if (binWidth > 0.0 && std::isfinite(binWidth))
        {
            bins = static_cast<int>(std::ceil((maximum - minimum) / binWidth));
        }
        bins = std::max(10, std::min(100, bins));
    }

    const double width = bins > 1 ? (maximum - minimum) / bins : 1.0;
    std::vector<double> counts(static_cast<size_t>(bins), 0.0);
    for (double value : values)
    {
        int bin = 0;
        if (bins > 1 && width > 0.0)
        {
            bin = static_cast<int>((value - minimum) / width);
            bin = std::max(0, std::min(bins - 1, bin));
        }
        counts[static_cast<size_t>(bin)] += 1.0;
    }

    q->RemoveExistingPlotChartAndTable();
    vtkMRMLTableNode* tableNode = q->GetOrCreatePlotTable();
    vtkNew<vtkDoubleArray> centers;
    vtkNew<vtkDoubleArray> countArray;
    centers->SetName("Value");
    countArray->SetName("Voxel count");
    for (int i = 0; i < bins; ++i)
    {
        const double center = bins == 1
            ? minimum
            : minimum + (static_cast<double>(i) + 0.5) * width;
        centers->InsertNextValue(center);
        countArray->InsertNextValue(counts[static_cast<size_t>(i)]);
    }
    tableNode->AddColumn(centers);
    tableNode->AddColumn(countArray);

    vtkMRMLPlotChartNode* chartNode = q->GetOrCreatePlotChart();
    const auto nameIt = q->segmentTACsnames.find(segmentID);
    const std::string roiName = nameIt != q->segmentTACsnames.end()
        ? nameIt->second : segmentID;
    const double frameEnd = this->frameEndForInputSec(
        static_cast<size_t>(observationIndex));
    chartNode->SetTitle(
        QObject::tr("ROI distribution - %1 - %2 (end %3 s)")
            .arg(QString::fromStdString(roiName))
            .arg(sourceDescription)
            .arg(frameEnd, 0, 'g', 8)
            .toStdString().c_str());
    chartNode->SetXAxisTitle(this->activityUnitLabel(displayUnit).toStdString().c_str());
    chartNode->SetYAxisTitle("Voxel count");

    vtkSmartPointer<vtkMRMLPlotSeriesNode> series =
        vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
    q->mrmlScene()->AddNode(series);
    series->SetName(roiName.c_str());
    series->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatterBar);
    series->SetAndObserveTableNodeID(tableNode->GetID());
    series->SetXColumnName("Value");
    series->SetYColumnName("Voxel count");
    series->SetUniqueColor();
    chartNode->AddAndObservePlotSeriesNodeID(series->GetID());

    vtkMRMLLayoutNode* layoutNode = vtkMRMLLayoutNode::SafeDownCast(
        q->mrmlScene()->GetFirstNodeByClass("vtkMRMLLayoutNode"));
    if (layoutNode)
    {
        layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);
    }
    vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
        q->mrmlScene()->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
    if (plotViewNode)
    {
        plotViewNode->SetPlotChartNodeID(chartNode->GetID());
    }
    return true;
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
loadTableWorkbook(
    const QString& filePath,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (filePath.trimmed().isEmpty() || !QFileInfo::exists(filePath))
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr("The selected TAC workbook does not exist.");
        }
        return false;
    }

    PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
    const QVariant resultVariant = mainContext.call(
        "DPE_load_tac_workbook",
        QVariantList{filePath, this->selectedTableTimeMode()});

    const QVariantMap result = resultVariant.toMap();
    if (!result.value("ok").toBool())
    {
        if (errorMessage)
        {
            *errorMessage = result.value("error").toString();
        }
        return false;
    }

    TACModeState newState;
    newState.valid = true;
    this->tableSigma.clear();
    this->tableFitStatistics.clear();
    this->tablePlotStatistics.clear();
    this->tablePlotTimesSec.clear();

    const QVariantList fitStats = result.value("fit_stats").toList();
    for (const QVariant& statVariant : fitStats)
    {
        const QVariantMap stat = statVariant.toMap();
        TACStatisticOption option;
        option.id = stat.value("id").toString().toStdString();
        option.label = stat.value("label").toString();
        if (!option.id.empty())
        {
            this->tableFitStatistics.push_back(option);
        }
    }

    const QVariantList plotStats = result.value("plot_stats").toList();
    for (const QVariant& statVariant : plotStats)
    {
        const QVariantMap stat = statVariant.toMap();
        TACStatisticOption option;
        option.id = stat.value("id").toString().toStdString();
        option.label = stat.value("label").toString();
        if (!option.id.empty())
        {
            this->tablePlotStatistics.push_back(option);
        }
    }

    const QVariantList plotTimes = result.value("plot_times_sec").toList();
    for (const QVariant& value : plotTimes)
    {
        this->tablePlotTimesSec.push_back(value.toDouble());
    }

    const QVariantList rois = result.value("rois").toList();
    if (rois.isEmpty())
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr("The workbook contains no usable ROI sheets.");
        }
        return false;
    }

    auto numericOrNaN = [](const QVariantMap& row, const QString& key)
    {
        const QVariant value = row.value(key);
        if (!value.isValid() || value.isNull())
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        bool ok = false;
        const double d = value.toDouble(&ok);
        return ok ? d : std::numeric_limits<double>::quiet_NaN();
    };

    bool timingInitialized = false;
    int roiCounter = 0;

    for (const QVariant& roiVariant : rois)
    {
        const QVariantMap roi = roiVariant.toMap();
        const QString roiName = roi.value("name").toString().trimmed();
        const QVariantList rows = roi.value("rows").toList();

        if (roiName.isEmpty() || rows.isEmpty())
        {
            continue;
        }

        const std::string segmentID =
            std::string("table::") + std::to_string(roiCounter++);

        std::vector<VoxelStatistics> statsVector;
        statsVector.reserve(rows.size());

        std::map<std::string, std::vector<double>> sigmaByStat;
        sigmaByStat["Mean"].reserve(rows.size());
        sigmaByStat["Median"].reserve(rows.size());
        sigmaByStat["Peak"].reserve(rows.size());
        sigmaByStat["Max"].reserve(rows.size());

        std::vector<double> localEnds;
        std::vector<double> localDurations;
        localEnds.reserve(rows.size());
        localDurations.reserve(rows.size());

        for (const QVariant& rowVariant : rows)
        {
            const QVariantMap row = rowVariant.toMap();

            VoxelStatistics vs;
            vs.mean = numericOrNaN(row, "Mean");
            vs.median = numericOrNaN(row, "Median");
            vs.min = numericOrNaN(row, "Min");
            vs.max = numericOrNaN(row, "Max");
            vs.stddev = numericOrNaN(row, "StDev");
            vs.q1 = numericOrNaN(row, "Q1");
            vs.q3 = numericOrNaN(row, "Q3");
            vs.iqr = numericOrNaN(row, "IQR");
            vs.volume_mm3 = numericOrNaN(row, "Volume(mm3)");
            vs.volume_cm3 = numericOrNaN(row, "Volume(cm3)");
            vs.peak = numericOrNaN(row, "Peak");
            vs.peakStddev = numericOrNaN(row, "PeakStDev");
            vs.count = row.value("VoxelCount").toInt();
            vs.peakCount = row.value("PeakVoxelCount").toInt();
            vs.keep = true;
            vs.empty = false;
            statsVector.push_back(vs);

            sigmaByStat["Mean"].push_back(numericOrNaN(row, "MeanSigma"));
            sigmaByStat["Median"].push_back(numericOrNaN(row, "MedianSigma"));
            sigmaByStat["Peak"].push_back(numericOrNaN(row, "PeakSigma"));
            sigmaByStat["Max"].push_back(numericOrNaN(row, "MaxSigma"));

            localEnds.push_back(row.value("FrameEnd_s").toDouble());
            localDurations.push_back(row.value("Duration_s").toDouble());
        }

        if (!timingInitialized)
        {
            newState.timePoints = localEnds;
            newState.durations = localDurations;
            newState.numberOfTimepoints = static_cast<int>(rows.size());
            timingInitialized = true;
        }

        newState.segmentTACs[segmentID] = std::move(statsVector);
        newState.segmentTACsnames[segmentID] = roiName.toStdString();
        newState.segmentDisplayOrder.push_back(segmentID);
        this->tableSigma[segmentID] = std::move(sigmaByStat);
    }

    // Keep Table mode consistent with Single and Multi: all user-facing ROI
    // lists use one case-insensitive A->Z display order. Synthetic table::N
    // IDs and the imported TAC data remain unchanged.
    newState.segmentDisplayOrder =
        sortedSegmentIDs(newState.segmentTACsnames);

    if (newState.segmentTACs.empty() ||
        newState.timePoints.empty() ||
        newState.durations.size() != newState.timePoints.size())
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr("The workbook did not produce a valid shared TAC time grid.");
        }
        return false;
    }

    this->tableTACState = std::move(newState);
    this->tableDataLoaded = true;
    this->tableWorkbookPath = filePath;
    this->tableFramingExact = result.value("framing_exact").toBool();
    this->tableTimingSummary = result.value("timing_summary").toString();

    // A newly loaded workbook must never inherit a conversion factor from
    // the previously loaded table. Metadata below may repopulate it.
    this->TableSUVbwFactorEdit->clear();

    const QVariantMap metadata = result.value("metadata").toMap();
    this->tableWorkbookMetadata = metadata;
    QString activityUnit = metadata.value("ActivityUnit").toString().trimmed();
    if (activityUnit.isEmpty())
    {
        activityUnit = result.value("suggested_activity_unit").toString().trimmed();
    }

    if (!activityUnit.isEmpty())
    {
        QSignalBlocker blocker(this->TableActivityUnitSelector);
        if (activityUnit.compare("SUVbw", Qt::CaseInsensitive) == 0 ||
            activityUnit.compare("SUV", Qt::CaseInsensitive) == 0)
        {
            this->TableActivityUnitSelector->setCurrentIndex(2);
        }
        else if (activityUnit.compare("kBq/mL", Qt::CaseInsensitive) == 0)
        {
            this->TableActivityUnitSelector->setCurrentIndex(1);
        }
        else if (activityUnit.compare("Bq/mL", Qt::CaseInsensitive) == 0)
        {
            this->TableActivityUnitSelector->setCurrentIndex(0);
        }
    }

    const QVariant factor = metadata.value("SUVbwFactor");
    if (factor.isValid() && !factor.isNull())
    {
        bool ok = false;
        const double value = factor.toDouble(&ok);
        if (ok && std::isfinite(value) && value > 0.0)
        {
            this->TableSUVbwFactorEdit->setText(QString::number(value, 'g', 12));
        }
    }

    this->TableWorkbookPathEdit->setText(filePath);
    this->updateTableUnitUI();

    if (this->tableBasedMode)
    {
        this->restoreActiveTACState(this->tableTACState);
        this->rebuildTACStatisticUI();
        this->populatePlotSegmentCheckboxes();

        // A new workbook is a new table dataset. Do not silently carry a
        // table-ROI IF identity across files, but keep the chosen source type.
        this->tableIFID.clear();
        q->IFID.clear();
        this->populateIF();
        this->populateVOI(std::string());
        this->populateVOIMTGA(std::string());
        this->populateTimeBarMTGA(true);
        this->updateAcquisitionTimingContext(true);
        this->invalidateInputFunctionResults();
        this->updateInputFunctionStatus();
        this->setPostTACEnabled(true);
    }

    const QString framingLabel = this->tableFramingExact
        ? QObject::tr("exact framing")
        : QObject::tr("inferred framing");

    const double acquisitionEndSec =
        this->tableTACState.timePoints.empty()
        ? 0.0
        : this->tableTACState.timePoints.back();

    this->TableStatusLabel->setText(
        QObject::tr("Loaded %1 ROIs, %2 frames | end=%3 s | %4 | %5")
            .arg(static_cast<qulonglong>(this->tableTACState.segmentTACs.size()))
            .arg(this->tableTACState.numberOfTimepoints)
            .arg(acquisitionEndSec, 0, 'g', 10)
            .arg(framingLabel)
            .arg(this->tableTimingSummary));

    return true;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
clearTableWorkbook()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    this->tableDataLoaded = false;
    this->tableTACState = TACModeState{};
    this->tableSigma.clear();
    this->tableFitStatistics.clear();
    this->tablePlotStatistics.clear();
    this->tablePlotTimesSec.clear();
    this->tableWorkbookPath.clear();
    this->tableFramingExact = false;
    this->tableTimingSummary.clear();
    this->tableWorkbookMetadata.clear();
    this->acquisitionTiming = AcquisitionTimingContext{};

    if (this->TableWorkbookPathEdit)
    {
        this->TableWorkbookPathEdit->clear();
    }
    if (this->TableSUVbwFactorEdit)
    {
        this->TableSUVbwFactorEdit->clear();
    }
    if (this->TableStatusLabel)
    {
        this->TableStatusLabel->setText(QObject::tr("No table TAC workbook loaded."));
    }

    if (this->tableBasedMode)
    {
        this->tableIFID.clear();
        q->IFID.clear();
        this->clearActiveTACState();
        this->rebuildTACStatisticUI();
        this->populateIF();
        this->populatePlotSegmentCheckboxes();
        this->populateVOI(std::string());
        this->populateVOIMTGA(std::string());
        this->invalidateInputFunctionResults();
        this->updateInputFunctionStatus();
        this->setPostTACEnabled(false);
        this->updateQuantitativeUnitUI();
        q->RemoveExistingPlotChartAndTable();
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
setTableBasedMode(bool enabled)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    // Multi-timepoint is an Image-mode choice only. Enforce exclusivity here
    // as well as in the checkbox handlers so programmatic mode changes cannot
    // leave both states active.
    if (enabled && this->multiTimepointMode)
    {
        {
            QSignalBlocker blocker(this->MultiTimepointAnalysisCheckBox);
            this->MultiTimepointAnalysisCheckBox->setChecked(false);
        }
        this->setMultiTimepointMode(false);
    }

    if (enabled == this->tableBasedMode)
    {
        this->MultiTimepointAnalysisCheckBox->setVisible(!enabled);
        return;
    }

    q->RemoveExistingPlotChartAndTable();
    q->clearFITdata();
    q->clearFITMTGAdata();

    if (enabled)
    {
        // Save image-mode state before switching sources.
        this->captureActiveTACState(this->imageTACState);
        this->imageIFSourceIndex = this->IFSourceSelector->currentIndex();
        this->imageIFID = q->IFID;

        this->tableBasedMode = true;
        if (this->AcquisitionTimingStatusLabel)
        {
            this->AcquisitionTimingStatusLabel->setVisible(false);
        }
        this->setImageSetupVisible(false);
        this->TableSetupGroupBox->setVisible(true);
        this->IFSourceSelector->setItemText(0, QObject::tr("Table ROI (IDIF)"));
        this->IFLabel->setText(QObject::tr("IDIF table ROI:"));

        if (this->tableDataLoaded)
        {
            this->restoreActiveTACState(this->tableTACState);
        }
        else
        {
            this->clearActiveTACState();
        }

        q->IFID = this->tableIFID;
        this->rebuildTACStatisticUI();
        this->populateIF();

        {
            QSignalBlocker blocker(this->IFSourceSelector);
            this->IFSourceSelector->setCurrentIndex(
                std::clamp(this->tableIFSourceIndex, 0, 1));
        }

        if (this->tableIFSourceIndex == 0)
        {
            const int index = this->IFSelector->findData(
                QString::fromStdString(this->tableIFID));
            QSignalBlocker blocker(this->IFSelector);
            this->IFSelector->setCurrentIndex(index >= 0 ? index : 0);
            q->IFID = index >= 0 ? this->tableIFID : std::string();
        }
        else
        {
            q->IFID.clear();
        }

        this->populatePlotSegmentCheckboxes();
        this->populateVOI(this->tableIFSourceIndex == 0 ? q->IFID : std::string());
        this->populateVOIMTGA(this->tableIFSourceIndex == 0 ? q->IFID : std::string());
        if (this->tableDataLoaded)
        {
            this->populateTimeBarMTGA();
        }
    }
    else
    {
        if (this->tableDataLoaded)
        {
            this->captureActiveTACState(this->tableTACState);
        }
        this->tableIFSourceIndex = this->IFSourceSelector->currentIndex();
        this->tableIFID = q->IFID;

        this->tableBasedMode = false;
        this->setImageSetupVisible(true);
        this->TableSetupGroupBox->setVisible(false);
        this->IFSourceSelector->setItemText(0, QObject::tr("Segment (IDIF)"));
        this->IFLabel->setText(QObject::tr("IDIF segment:"));

        this->restoreActiveTACState(this->imageTACState);
        q->IFID = this->imageIFID;
        this->rebuildTACStatisticUI();
        this->populateIF();

        {
            QSignalBlocker blocker(this->IFSourceSelector);
            this->IFSourceSelector->setCurrentIndex(
                std::clamp(this->imageIFSourceIndex, 0, 1));
        }

        if (this->imageIFSourceIndex == 0)
        {
            const int index = this->IFSelector->findData(
                QString::fromStdString(this->imageIFID));
            QSignalBlocker blocker(this->IFSelector);
            this->IFSelector->setCurrentIndex(index >= 0 ? index : 0);
            q->IFID = index >= 0 ? this->imageIFID : std::string();
        }
        else
        {
            q->IFID.clear();
        }

        this->populatePlotSegmentCheckboxes();
        this->populateVOI(this->imageIFSourceIndex == 0 ? q->IFID : std::string());
        this->populateVOIMTGA(this->imageIFSourceIndex == 0 ? q->IFID : std::string());
        if (!q->timePoints.empty())
        {
            this->populateTimeBarMTGA();
        }

        // The subject hierarchy may have changed while table mode was active.
        q->onSubjectHierarchyChanged();

        // A Multi -> Table transition intentionally cancels the Single refresh
        // scheduled while Multi is being left. Now that Table is also being
        // left, perform the same PET/segmentation display and watcher rebuild as
        // the direct Multi -> Single path.
        this->scheduleSingleModeAcquisitionRefresh();
    }

    // Both source choices are valid in either mode: source 0 is an image
    // Segment IDIF in Image mode and a table ROI IDIF in Table mode; source 1
    // is the shared external CSV workflow. Re-enable explicitly so no stale
    // item state survives a mode switch.
    if (QStandardItemModel* sourceModel =
            qobject_cast<QStandardItemModel*>(this->IFSourceSelector->model()))
    {
        if (sourceModel->item(0)) sourceModel->item(0)->setEnabled(true);
        if (sourceModel->item(1)) sourceModel->item(1)->setEnabled(true);
    }

    // Image-only controls are hidden, not merely disabled, in table mode.
    this->SegmentationAdvancedCollapsibleButton->setVisible(
        !this->tableBasedMode && !this->multiTimepointMode);
    this->Plottacsave->setVisible(!this->tableBasedMode);
    this->direxcel->setVisible(!this->tableBasedMode);
    this->fileexcel->setVisible(!this->tableBasedMode);
    this->saveExcelButton->setVisible(!this->tableBasedMode);

    const int imagingIndex = this->PlotsTabWidget->indexOf(this->ImagingWidget);
    if (imagingIndex >= 0)
    {
#if QT_VERSION >= QT_VERSION_CHECK(5, 15, 0)
        this->PlotsTabWidget->setTabVisible(
            imagingIndex,
            !this->tableBasedMode && !this->multiTimepointMode);
#else
        this->PlotsTabWidget->setTabEnabled(
            imagingIndex,
            !this->tableBasedMode && !this->multiTimepointMode);
#endif
    }

    this->externalIFPreviewSelectedIndex = -1;
    this->updateAcquisitionTimingContext(true);
    this->invalidateInputFunctionResults();
    this->updateTableUnitUI();
    this->updateInputFunctionStatus();
    this->updateSegmentationAdvancedUI();
    this->setPostTACEnabled(!q->segmentTACs.empty());
    this->updateParametricImagingAvailability();

    // Multi-timepoint belongs exclusively to Image mode; do not show an
    // irrelevant second mode switch while Table mode is active.
    this->MultiTimepointAnalysisCheckBox->setVisible(!this->tableBasedMode);

    // Force layouts to recompute after a mode switch; this minimizes dead
    // vertical space without adding another nested/scrolling UI hierarchy.
    this->verticalLayout->invalidate();
    this->TACwidget->adjustSize();
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
initializeTableBasedUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    // All widgets in this section are defined in the .ui file so the layout
    // remains editable in Qt Designer. Do not create duplicate controls here.
    this->TableSetupGroupBox->setVisible(false);
    this->verticalLayout->setAlignment(Qt::AlignTop);

    QObject::connect(
        this->TableBasedAnalysisCheckBox,
        &QCheckBox::toggled,
        q,
        [this](bool checked)
        {
            if (checked && this->multiTimepointMode)
            {
                {
                    QSignalBlocker blocker(this->MultiTimepointAnalysisCheckBox);
                    this->MultiTimepointAnalysisCheckBox->setChecked(false);
                }
                this->setMultiTimepointMode(false);
            }
            this->setTableBasedMode(checked);
        });

    QObject::connect(
        this->TableWorkbookBrowseButton,
        &QPushButton::clicked,
        q,
        [this, q]()
        {
            const QString path = QFileDialog::getOpenFileName(
                q,
                QObject::tr("Select ROI TAC workbook"),
                this->tableWorkbookPath,
                QObject::tr("Excel workbooks (*.xlsx)"));

            if (path.isEmpty())
            {
                return;
            }

            QString error;
            if (!this->loadTableWorkbook(path, &error))
            {
                QMessageBox::warning(q, QObject::tr("Table-based TACs"), error);
            }
        });

    QObject::connect(
        this->TableWorkbookClearButton,
        &QPushButton::clicked,
        q,
        [this]()
        {
            this->clearTableWorkbook();
        });

    QObject::connect(
        this->TableTimeModeSelector,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        q,
        [this, q](int)
        {
            if (this->tableWorkbookPath.isEmpty())
            {
                return;
            }
            QString error;
            if (!this->loadTableWorkbook(this->tableWorkbookPath, &error))
            {
                QMessageBox::warning(q, QObject::tr("Table-based TACs"), error);
            }
        });

    QObject::connect(
        this->TableActivityUnitSelector,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        q,
        [this](int)
        {
            this->updateTableUnitUI();
            const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
            this->invalidateInputFunctionResults();
            this->updateInputFunctionStatus();
            if (previewOpen)
            {
                this->previewInputFunction();
            }
        });

    QObject::connect(
        this->TableSUVbwFactorEdit,
        &QLineEdit::editingFinished,
        q,
        [this]()
        {
            const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
            this->updateTableUnitUI();
            this->invalidateInputFunctionResults();
            this->updateInputFunctionStatus();
            if (previewOpen)
            {
                this->previewInputFunction();
            }
        });

    QObject::connect(
        this->StatSelector,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        q,
        [this](int)
        {
            this->updateTableWeightingAvailability();
        });

    QObject::connect(
        this->StatSelectorMTGA,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        q,
        [this](int)
        {
            this->updateTableWeightingAvailability();
        });

    this->updateTableUnitUI();
}

//-----------------------------------------------------------------------------
void qSlicerDynamicPETModuleWidgetPrivate::init()
{
  Q_Q(qSlicerDynamicPETModuleWidget);
  this->setupUi(q);

  // Distribution is tied to one concrete image observation. Keep its frame
  // selector immediately below the Plot controls and hidden unless the
  // Distribution statistic is active. It is created programmatically so the
  // feature does not depend on regenerated Ui_* members.
  this->distributionFrameLabel = new QLabel(
      QObject::tr("Distribution frame:"), this->TACCollapsibleButton);
  this->distributionFrameWidget = new QWidget(this->TACCollapsibleButton);
  QHBoxLayout* distributionFrameLayout = new QHBoxLayout(this->distributionFrameWidget);
  distributionFrameLayout->setContentsMargins(0, 0, 0, 0);
  this->distributionFrameSlider = new QSlider(Qt::Horizontal, this->distributionFrameWidget);
  this->distributionFrameSlider->setMinimum(1);
  this->distributionFrameSlider->setMaximum(1);
  this->distributionFrameSlider->setValue(1);
  this->distributionFrameSlider->setToolTip(QObject::tr(
      "Select the PET frame/observation used for the ROI voxel distribution."));
  this->distributionFrameInfoEdit = new QLineEdit(this->distributionFrameWidget);
  this->distributionFrameInfoEdit->setReadOnly(true);
  this->distributionFrameInfoEdit->setMinimumWidth(235);
  distributionFrameLayout->addWidget(this->distributionFrameSlider, 1);
  distributionFrameLayout->addWidget(this->distributionFrameInfoEdit);
  if (this->formLayout_2)
  {
      this->formLayout_2->insertRow(3, this->distributionFrameLabel, this->distributionFrameWidget);
  }
  this->distributionFrameLabel->setVisible(false);
  this->distributionFrameWidget->setVisible(false);
  QObject::connect(
      this->distributionFrameSlider, &QSlider::valueChanged, q,
      [this](int)
      {
          this->updateDistributionFrameInfo();
          this->refreshDistributionPlotIfActive();
      });

  this->PlotErrorCheckbox->setText(QObject::tr("Dispersion"));

  // Keep the segmentation tools adjacent to the image/TAC visualization they
  // operate on. Move Advanced [Segmentation] directly below Plot without
  // depending on Designer item ordering.
  if (this->SegmentationAdvancedCollapsibleButton && this->TACCollapsibleButton)
  {
      QBoxLayout* parentBox = qobject_cast<QBoxLayout*>(
          this->TACCollapsibleButton->parentWidget()
              ? this->TACCollapsibleButton->parentWidget()->layout()
              : nullptr);
      if (parentBox)
      {
          const int plotIndex = parentBox->indexOf(this->TACCollapsibleButton);
          if (plotIndex >= 0)
          {
              parentBox->removeWidget(this->SegmentationAdvancedCollapsibleButton);
              parentBox->insertWidget(plotIndex + 1, this->SegmentationAdvancedCollapsibleButton);
          }
      }
  }

  // Frame/observation navigator used only to change what PET/segmentation is
  // displayed for inspection. Editing detection is independent and comes from
  // acquisition-specific segmentation watchers. In Multi this traverses
  // prepared provenance observations, including late static acquisitions,
  // without creating a synthetic combined image sequence.
  this->segmentationFrameLabel = new QLabel(
      QObject::tr("Displayed frame:"), this->SegmentationAdvancedCollapsibleButton);
  this->segmentationFrameWidget = new QWidget(this->SegmentationAdvancedCollapsibleButton);
  QHBoxLayout* segmentationFrameLayout = new QHBoxLayout(this->segmentationFrameWidget);
  segmentationFrameLayout->setContentsMargins(0, 0, 0, 0);
  this->segmentationFrameSlider = new QSlider(Qt::Horizontal, this->segmentationFrameWidget);
  this->segmentationFrameSlider->setMinimum(1);
  this->segmentationFrameSlider->setMaximum(1);
  this->segmentationFrameSlider->setValue(1);
  this->segmentationFrameSlider->setToolTip(QObject::tr(
      "Select the PET/segmentation frame displayed for inspection. Segment Editor changes are detected independently by acquisition-specific watchers."));
  this->segmentationFrameInfoEdit = new QLineEdit(this->segmentationFrameWidget);
  this->segmentationFrameInfoEdit->setReadOnly(true);
  this->segmentationFrameInfoEdit->setMinimumWidth(215);
  segmentationFrameLayout->addWidget(this->segmentationFrameSlider, 1);
  segmentationFrameLayout->addWidget(this->segmentationFrameInfoEdit);
  if (this->SegmentationAdvancedLayout)
  {
      this->SegmentationAdvancedLayout->insertRow(
          0, this->segmentationFrameLabel, this->segmentationFrameWidget);
  }
  QObject::connect(
      this->segmentationFrameSlider, &QSlider::valueChanged, q,
      [this](int)
      {
          if (this->updatingSegmentationFrameSlider)
          {
              return;
          }
          this->updateSegmentationFrameInfo();
          this->displaySelectedSegmentationFrame();
      });

  // qSlicerDynamicPETModuleWidget.ui contains this checkbox in current
  // sources.  Resolve it by object name instead of relying on the generated
  // Ui_* class member so incremental builds using an older uic header still
  // compile.  If the Designer widget is absent, create the same control
  // programmatically at the end of the Advanced IF layout.
  this->fengExtrapolationCheckBox =
      q->findChild<QCheckBox*>(QStringLiteral("IFFengExtrapolationCheckBox"));
  if (!this->fengExtrapolationCheckBox)
  {
      this->fengExtrapolationCheckBox =
          new QCheckBox(
              QObject::tr("Allow Feng extrapolation to PET support"),
              this->IFAdvancedCollapsibleButton);
      this->fengExtrapolationCheckBox->setObjectName(
          QStringLiteral("IFFengExtrapolationCheckBox"));
      this->fengExtrapolationCheckBox->setToolTip(
          QObject::tr(
              "When Feng modeling is selected, extend the fitted analytic "
              "input function beyond the last measured blood sample only as "
              "far as required by the available PET observations. Off by "
              "default. Interpolation/smoothing and PBIF templates are never "
              "extrapolated. Disabled while PBIF calibration is active so "
              "modeled tail values cannot enter PBIF AUC calibration."));
      this->fengExtrapolationCheckBox->setChecked(false);
      if (this->IFAdvancedLayout)
      {
          this->IFAdvancedLayout->addWidget(
              this->fengExtrapolationCheckBox,
              this->IFAdvancedLayout->rowCount(),
              0, 1, 4);
      }
  }

  this->initializeTableBasedUI();
  this->initializeMultiTimepointUI();

  this->PlotLiveSegEdit->setToolTip(
      QObject::tr(
          "Track Segment Editor corrections and refresh image-derived TAC/plots. In Multi-timepoint mode each prepared PET acquisition has its own segmentation watcher; the displayed-frame slider is visualization-only."));
  this->OpenSegmentEditorButton->setEnabled(false);
  this->SaveDynamicRTStructButton->setEnabled(false);
  this->SegmentationAdvancedCollapsibleButton->setCollapsed(true);
  this->RESETbutton->setToolTip(
      QObject::tr(
          "Restore points removed from the active mode's tissue TACs and external input function."));

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
            if (this->tableBasedMode)
                this->tableIFSourceIndex = source;
            else
                this->imageIFSourceIndex = source;

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

            const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
            this->invalidateInputFunctionResults();
            this->updateInputFunctionStatus();
            if (previewOpen)
            {
                this->previewInputFunction();
            }
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

            const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
            this->invalidateInputFunctionResults();
            this->updateInputFunctionStatus();
            if (previewOpen)
            {
                this->previewInputFunction();
            }
        });


    QObject::connect(
        this->IFCurveTypeSelector,
        QOverload<int>::of(
            &QComboBox::currentIndexChanged),
        q,
        [this](int)
        {
            const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
            this->invalidateInputFunctionResults();
            this->updateInputFunctionUI();
            if (previewOpen)
            {
                this->previewInputFunction();
            }
        });

    auto onPBRChanged =
        [this]()
        {
            const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
            this->invalidateInputFunctionResults();
            this->updateInputFunctionUI();
            if (previewOpen)
            {
                this->previewInputFunction();
            }
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

          const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
          if (previewOpen)
          {
              this->previewInputFunction();
          }
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
              this->logToPythonConsole(
                  QObject::tr(
                      "[SlicerDynamicPET PBIF] Template begins after 0 s; "
                      "an internal (0 s, 0) interpolation anchor was added. "
                      "The CSV file was not modified."));
          }

          this->updateAcquisitionTimingContext(false);

          double startTime = std::max(
              this->pbifTimesSec.empty() ? 0.0 : this->pbifTimesSec.front(),
              this->currentObservedInputStartSec());

          double endTime =
              this->pbifTimesSec.empty()
              ? 0.0
              : this->pbifTimesSec.back();

          if (!q->timePoints.empty())
          {
              endTime = std::min(
                  endTime,
                  this->frameEndForInputSec(q->timePoints.size() - 1));
          }

          if (endTime > startTime)
          {
              this->PBIFCalibrationStartSpinBox->setValue(startTime);
              this->PBIFCalibrationEndSpinBox->setValue(endTime);
          }

          const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
          if (previewOpen)
          {
              this->previewInputFunction();
          }
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

          const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
          if (previewOpen)
          {
              this->previewInputFunction();
          }
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
          const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
          // Clear the cached/precomputed IF before rebuilding status.
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
          this->updateROIModelingAvailability();
          this->updateParametricImagingAvailability();
          if (previewOpen)
          {
              const bool pbifIncomplete =
                  this->PBIFOptionCheckBox->isChecked() &&
                  (this->pbifTimesSec.size() < 2 ||
                   this->pbifTimesSec.size() != this->pbifTemplateValues.size());
              if (pbifIncomplete)
              {
                  this->removePreviewGroup("InputFunctionPreview");
              }
              else
              {
                  this->previewInputFunction();
              }
          }
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
          // Display conversion never changes fitting data/results. If the IF
          // preview is already open, replace it in place so unit changes are
          // immediately visible without forcing the user to press Preview.
          const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
          this->updateInputFunctionStatus();
          if (previewOpen)
          {
              this->previewInputFunction();
          }
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
      this->fengExtrapolationCheckBox,
      &QCheckBox::toggled,
      q,
      [invalidateIFPipeline](bool)
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
      this->IFGaussianSigmaSpinBox,
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
          const bool previewOpen = this->previewGroupExists("InputFunctionPreview");
          this->invalidateInputFunctionResults();
          this->updateInputFunctionUI();
          if (previewOpen)
          {
              this->previewInputFunction();
          }
      });

  QObject::connect(
      this->timeStepEdit,
      &QLineEdit::editingFinished,
      q,
      [invalidateIFPipeline]()
      {
          invalidateIFPipeline();
      });

  QObject::connect(
      this->IFPreviewButton,
      &QPushButton::clicked,
      q,
      [this]()
      {
          this->previewInputFunction();
      });

  QObject::connect(
      this->IFExportButton,
      &QPushButton::clicked,
      q,
      [this, q]()
      {
          QString exportBaseName;
          if (q->SubjectHierarchyNode
              && q->patID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
          {
              exportBaseName = QString::fromStdString(
                  q->SubjectHierarchyNode->GetItemName(q->patID));
          }
          const QString suggestedFileName = exportBaseName.isEmpty()
              ? QStringLiteral("IF.csv")
              : exportBaseName + QStringLiteral("_IF.csv");
          const QString suggestedPath =
              this->sharedOutputDirectory.trimmed().isEmpty()
              ? suggestedFileName
              : QDir(this->sharedOutputDirectory)
                    .filePath(suggestedFileName);

          QString filePath = QFileDialog::getSaveFileName(
              q,
              QObject::tr("Export final input function"),
              suggestedPath,
              QObject::tr("CSV files (*.csv)"));
          if (filePath.isEmpty())
              return;
          if (!filePath.endsWith(".csv", Qt::CaseInsensitive))
              filePath += ".csv";

          this->propagateOutputDirectory(QFileInfo(filePath).absolutePath());

          QString error;
          if (!this->exportFinalInputFunctionCSV(filePath, &error))
          {
              QMessageBox::warning(q, QObject::tr("Export input function"), error);
          }
      });


  this->IFExportButton->setEnabled(false);

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
  QObject::connect( this->OpenSegmentEditorButton, SIGNAL(clicked(bool)),
    q, SLOT(onOpenSegmentEditor()));
  QObject::connect( this->SaveDynamicRTStructButton, SIGNAL(clicked(bool)),
    q, SLOT(onSaveDynamicRTStruct()));
  QObject::connect(
      this->DynamicRTStructDirectory,
      &ctkPathLineEdit::currentPathChanged,
      q,
      [this](const QString&) { this->updateSegmentationAdvancedUI(); });
  QObject::connect(
      this->DynamicRTStructFilename,
      &QLineEdit::textChanged,
      q,
      [this](const QString&) { this->updateSegmentationAdvancedUI(); });
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
  QObject::connect( this->PlotErrorCheckbox, SIGNAL(toggled(bool)),
    q, SLOT(onPlotbutton()));
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

  // One output directory selection is shared across all output destinations.
  // Filenames remain independent; only directory browsing is synchronized.
  for (ctkPathLineEdit* outputDirectory : {
           this->direxcel,
           this->direxceltcm,
           this->direxcelmtga,
           this->direxceltcmfitted,
           this->direxcelmtgafitted,
           this->DynamicRTStructDirectory})
  {
    QObject::connect(
        outputDirectory,
        &ctkPathLineEdit::currentPathChanged,
        q,
        [this](const QString& path)
        {
          this->propagateOutputDirectory(path);
        });
  }
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
  QObject::connect(this->standardFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onStdFitclicked()));
  QObject::connect(this->standardFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onStdFitImgclicked()));
  QObject::connect(this->weightFitCheckBox, SIGNAL(toggled(bool)),
    q, SLOT(onWFitclicked()));
  QObject::connect(this->weightFitCheckBoxImg, SIGNAL(toggled(bool)),
    q, SLOT(onWFitImgclicked()));

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


  this->SegmentationAdvancedCollapsibleButton->setCollapsed(true);
  this->TACCollapsibleButton->setCollapsed(true);
  this->TCMCollapsibleButton->setCollapsed(true);
  this->MTGACollapsibleButton->setCollapsed(true);

  // Keep Parametric Imaging compact when optional collapsible sections are
  // closed. The module's outer scroll area can still grow when sections open.
  this->verticalLayoutImg->setAlignment(Qt::AlignTop);
  this->mtgaLayoutImg->setAlignment(Qt::AlignTop);
  this->tcmLayoutImg->setAlignment(Qt::AlignTop);
  this->verticalLayoutImg->setSpacing(4);
  this->mtgaLayoutImg->setSpacing(4);
  this->tcmLayoutImg->setSpacing(4);
  this->ModelsTabWidgetImg->setSizePolicy(
      QSizePolicy::Preferred, QSizePolicy::Maximum);
  this->MTGAStatTestButton->setCollapsed(true);
  this->MTGAStatTestButton->setEnabled(false);
  this->TCMStatTestButton->setCollapsed(true);
  this->TCMStatTestButton->setEnabled(false);
  for (const QString& name : q->checkboxNames)
  {
    QCheckBox* cb = new QCheckBox(name, this->PlotStatsCheckContents);
    cb->setProperty("StatID", name);
    this->PlotStatsCheckLayout->addWidget(cb);
  }

  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  QString dpePythonScript;
  dpePythonScript += QString::fromUtf8(R"DPEPY1(
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

try:
    import openpyxl
except ImportError:
    import slicer
    slicer.util.pip_install("openpyxl")

import importlib
import importlib.metadata
import sys
import slicer

def DPE_console_message(message):
    print(message)

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

def DPE_save_multisheet_excel(filepath, sheet_data_dict, metadata=None):
    """
    Save the detailed DynamicPET TAC workbook.  ROI sheet names remain ROI
    names.  The optional _DynamicPET sheet is deliberately small and optional
    so that externally produced workbooks do not need to reproduce it.
    """
    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        for sheet, data in sheet_data_dict.items():
            df = pd.DataFrame(data)[["Time(s)", "FrameStart_s", "FrameMid_s", "FrameEnd_s", "Duration_s", "Mean", "Median", "StDev","IQR","Min", "Max", "Q1", "Q3", "Peak", "PeakStDev", "PeakVoxelCount", "VoxelCount","Volume(mm3)","Volume(cm3)"]]
            df.to_excel(writer, sheet_name=sheet, index=False)

        if metadata:
            md = pd.DataFrame(
                [{"Key": str(k), "Value": v} for k, v in metadata.items()]
            )
            md.to_excel(writer, sheet_name="_DynamicPET", index=False)


def DPE_load_tac_workbook(filepath, time_mode="auto"):
    """
    Normalize a sheet-per-ROI workbook for the C++ widget.

    Minimum per ROI sheet:
      * one numeric time column
      * one numeric TAC value column

    Exact frame start/end/duration columns are used when available.  If only a
    representative time exists, a contiguous frame schedule is inferred.  A
    generic time is interpreted as frame end in Auto mode.
    """
    import math
    import re

    def norm(name):
        return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())

    def col_lookup(df):
        return {norm(c): c for c in df.columns}

    def find_col(lookup, aliases):
        for alias in aliases:
            key = norm(alias)
            if key in lookup:
                return lookup[key]
        return None

    def numeric_list(df, column):
        if column is None:
            return None
        values = pd.to_numeric(df[column], errors="coerce")
        if values.isna().any():
            return None
        return [float(v) for v in values.tolist()]

    def optional_value(row, column):
        if column is None:
            return None
        value = row[column]
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        try:
            return float(value)
        except Exception:
            return None

    try:
        xls = pd.ExcelFile(filepath)
    except Exception as exc:
        return {"ok": False, "error": "Could not open Excel workbook: " + str(exc)}

    metadata = {}
    if "_DynamicPET" in xls.sheet_names:
        try:
            md = pd.read_excel(filepath, sheet_name="_DynamicPET")
            if "Key" in md.columns and "Value" in md.columns:
                for _, row in md.iterrows():
                    key = str(row["Key"]).strip()
                    if key and key.lower() != "nan":
                        value = row["Value"]
                        if not pd.isna(value):
                            metadata[key] = value.item() if hasattr(value, "item") else value
        except Exception:
            metadata = {}

    roi_sheet_names = [s for s in xls.sheet_names if s != "_DynamicPET"]
    if not roi_sheet_names:
        return {"ok": False, "error": "No ROI worksheets were found."}

    normalized_rois = []
    reference_ends = None
    reference_durations = None
    reference_plot_times = None
    framing_exact_all = True
    timing_summaries = []
    common_fit_ids = None
    common_plot_ids = None
    labels_by_id = {}
    suggested_units = []

    for sheet_name in roi_sheet_names:
        try:
            df = pd.read_excel(filepath, sheet_name=sheet_name)
        except Exception as exc:
            return {"ok": False, "error": f"Could not read ROI sheet '{sheet_name}': {exc}"}

        df = df.dropna(how="all").reset_index(drop=True)
        if df.empty:
            continue

        df.columns = [str(c).strip() for c in df.columns]
        lookup = col_lookup(df)

        start_col = find_col(lookup, ["FrameStart_s", "FrameStart", "Start_s", "Start"])
        mid_col = find_col(lookup, ["FrameMid_s", "FrameMid", "MidTime_s", "MidTime", "Midpoint_s", "Midpoint"])
        end_col = find_col(lookup, ["FrameEnd_s", "FrameEnd", "End_s", "End"])
        duration_col = find_col(lookup, ["Duration_s", "FrameDuration_s", "Duration", "FrameDuration"])
        generic_time_col = find_col(lookup, ["Time(s)", "Time_s", "Time", "time_sec", "Seconds"])

        starts = numeric_list(df, start_col)
        mids = numeric_list(df, mid_col)
        ends = numeric_list(df, end_col)
        durations = numeric_list(df, duration_col)
        generic_times = numeric_list(df, generic_time_col)

        exact_framing = starts is not None and ends is not None

        if exact_framing:
            if durations is None:
                durations = [b - a for a, b in zip(starts, ends)]
            if mids is None:
                mids = [a + 0.5 * d for a, d in zip(starts, durations)]
            plot_times = generic_times if generic_times is not None else mids
            timing_summary = "exact frame start/end"
        else:
            framing_exact_all = False
            rep = None
            convention = str(time_mode or "auto").lower()

            if convention == "auto":
                if mids is not None:
                    rep = mids
                    convention = "midpoint"
                elif ends is not None:
                    rep = ends
                    convention = "end"
                elif generic_times is not None:
                    rep = generic_times
                    convention = "end"
                else:
                    return {"ok": False, "error": f"ROI sheet '{sheet_name}' has no usable time column."}
            elif convention == "midpoint":
                rep = mids if mids is not None else generic_times
                if rep is None:
                    rep = ends
            elif convention == "end":
                rep = ends if ends is not None else generic_times
                if rep is None:
                    rep = mids
            else:
                return {"ok": False, "error": "Unknown table time convention."}

            if rep is None or len(rep) < 2:
                return {"ok": False, "error": f"ROI sheet '{sheet_name}' needs at least two time points when framing must be inferred."}

            if any(rep[i] <= rep[i - 1] for i in range(1, len(rep))):
                return {"ok": False, "error": f"Time values in ROI sheet '{sheet_name}' must be strictly increasing."}

            if durations is not None:
                if convention == "midpoint":
                    mids = rep
                    starts = [t - 0.5 * d for t, d in zip(rep, durations)]
                    ends = [t + 0.5 * d for t, d in zip(rep, durations)]
                else:
                    ends = rep
                    starts = [t - d for t, d in zip(rep, durations)]
                    mids = [a + 0.5 * d for a, d in zip(starts, durations)]
            elif convention == "end":
                ends = rep
                starts = [0.0] + list(rep[:-1])
                durations = [b - a for a, b in zip(starts, ends)]
                mids = [a + 0.5 * d for a, d in zip(starts, durations)]
            else:
                mids = rep
                boundaries = [max(0.0, rep[0] - 0.5 * (rep[1] - rep[0]))]
                boundaries += [0.5 * (rep[i - 1] + rep[i]) for i in range(1, len(rep))]
                boundaries += [rep[-1] + 0.5 * (rep[-1] - rep[-2])]
                starts = boundaries[:-1]
                ends = boundaries[1:]
                durations = [b - a for a, b in zip(starts, ends)]

            plot_times = rep
            timing_summary = "inferred from " + ("frame midpoint" if convention == "midpoint" else "frame end")

        if any((not math.isfinite(d)) or d <= 0.0 for d in durations):
            return {"ok": False, "error": f"ROI sheet '{sheet_name}' produced non-positive frame durations."}

        if any(s < -1e-9 for s in starts):
            return {"ok": False, "error": f"ROI sheet '{sheet_name}' produced a frame start before 0 s. Check the selected time convention."}

        # TAC columns. DynamicPET columns are recognized first; a simple foreign
        # workbook can instead supply Value/TAC/Activity/Concentration/SUV.
        mean_col = find_col(lookup, ["Mean"])
        median_col = find_col(lookup, ["Median"])
        peak_col = find_col(lookup, ["Peak", "SUVpeak"])
        max_col = find_col(lookup, ["Max", "Maximum", "SUVmax"])
  )DPEPY1");

  dpePythonScript += QString::fromUtf8(R"DPEPY2(
        generic_value_col = find_col(lookup, ["Value", "TAC", "Activity", "Concentration", "Radioactivity", "SUVbw", "SUV"])

        stat_columns = {}
        stat_labels = {}
        if mean_col is not None:
            stat_columns["Mean"] = mean_col
            stat_labels["Mean"] = "Mean"
        if median_col is not None:
            stat_columns["Median"] = median_col
            stat_labels["Median"] = "Median"
        if peak_col is not None:
            stat_columns["Peak"] = peak_col
            stat_labels["Peak"] = "Peak"
        if max_col is not None:
            stat_columns["Max"] = max_col
            stat_labels["Max"] = "Max"

        if not stat_columns and generic_value_col is not None:
            stat_columns["Mean"] = generic_value_col
            stat_labels["Mean"] = str(generic_value_col)
            n = norm(generic_value_col)
            if "suv" in n:
                suggested_units.append("SUVbw")

        if not stat_columns:
            return {"ok": False, "error": f"ROI sheet '{sheet_name}' has no recognized TAC value column (for example Mean, Value, TAC, Activity, Concentration, or SUV)."}

        # Require the chosen fit statistic to exist for every ROI. The common
        # set is exposed to the GUI after all sheets have been inspected.
        fit_ids = set(stat_columns.keys())

        min_col = find_col(lookup, ["Min", "Minimum"])
        plot_ids = set(fit_ids)
        if min_col is not None:
            plot_ids.add("Min")
            stat_labels["Min"] = "Min"
        if max_col is not None:
            plot_ids.add("Max")
            stat_labels["Max"] = "Max"

        common_fit_ids = fit_ids if common_fit_ids is None else common_fit_ids.intersection(fit_ids)
        common_plot_ids = plot_ids if common_plot_ids is None else common_plot_ids.intersection(plot_ids)
        for key, value in stat_labels.items():
            labels_by_id.setdefault(key, []).append(value)

        stdev_col = find_col(lookup, ["StDev", "StdDev", "StandardDeviation", "SD"])
        iqr_col = find_col(lookup, ["IQR", "InterquartileRange"])
        peak_stdev_col = find_col(lookup, ["PeakStDev", "PeakStdDev", "SUVpeakStDev"])
        generic_sigma_col = find_col(lookup, ["Sigma", "Uncertainty", "Error", "SEM", "StandardError"])
        max_sigma_col = find_col(lookup, ["MaxSigma", "SUVmaxSigma", "MaxUncertainty", "MaxSEM"])

        q1_col = find_col(lookup, ["Q1"])
        q3_col = find_col(lookup, ["Q3"])
        voxel_count_col = find_col(lookup, ["VoxelCount", "Count"])
        peak_count_col = find_col(lookup, ["PeakVoxelCount"])
        vol_mm3_col = find_col(lookup, ["PETVolume(mm3)", "Volume(mm3)", "Volume_mm3", "VolumeMM3"])
        vol_cm3_col = find_col(lookup, ["PETVolume(cm3)", "Volume(cm3)", "Volume_cm3", "VolumeCC", "Volume(cc)"])
        if vol_cm3_col is not None or vol_mm3_col is not None:
            plot_ids.add("VolumePET")
            stat_labels["VolumePET"] = "Volume [PET] (cm3)"

        rows = []
        for i, (_, source_row) in enumerate(df.iterrows()):
            row = {
                "FrameStart_s": float(starts[i]),
                "FrameMid_s": float(mids[i]),
                "FrameEnd_s": float(ends[i]),
                "Duration_s": float(durations[i]),
            }

            for stat_id, column in stat_columns.items():
                row[stat_id] = optional_value(source_row, column)

            row["Min"] = optional_value(source_row, min_col)
            row["Max"] = optional_value(source_row, max_col)
            row["StDev"] = optional_value(source_row, stdev_col)
            row["IQR"] = optional_value(source_row, iqr_col)
            row["Q1"] = optional_value(source_row, q1_col)
            row["Q3"] = optional_value(source_row, q3_col)
            row["PeakStDev"] = optional_value(source_row, peak_stdev_col)
            row["VoxelCount"] = optional_value(source_row, voxel_count_col)
            row["PeakVoxelCount"] = optional_value(source_row, peak_count_col)
            row["Volume(mm3)"] = optional_value(source_row, vol_mm3_col)
            row["Volume(cm3)"] = optional_value(source_row, vol_cm3_col)
            # Normalize uncertainty to one-sigma values for WLS, irrespective
            # of whether the workbook stored SD, IQR or an explicit sigma/SEM.
            mean_sigma = optional_value(source_row, stdev_col)
            if mean_sigma is None and "Mean" in stat_columns:
                mean_sigma = optional_value(source_row, generic_sigma_col)

            median_sigma = None
            iqr = optional_value(source_row, iqr_col)
            if iqr is not None:
                median_sigma = iqr / 1.3489795003921634
            elif "Median" in stat_columns:
                median_sigma = optional_value(source_row, generic_sigma_col)

            peak_sigma = optional_value(source_row, peak_stdev_col)
            if peak_sigma is None and "Peak" in stat_columns:
                peak_sigma = optional_value(source_row, generic_sigma_col)

            row["MeanSigma"] = mean_sigma
            row["MedianSigma"] = median_sigma
            max_sigma = optional_value(source_row, max_sigma_col)
            if max_sigma is None and "Max" in stat_columns:
                # Only use a generic sigma if the workbook explicitly has one;
                # image-derived SUVmax itself has no automatic uncertainty proxy.
                max_sigma = optional_value(source_row, generic_sigma_col)

            row["PeakSigma"] = peak_sigma
            row["MaxSigma"] = max_sigma
            rows.append(row)

        if reference_ends is None:
            reference_ends = list(ends)
            reference_durations = list(durations)
            reference_plot_times = list(plot_times)
        else:
            if len(ends) != len(reference_ends):
                return {"ok": False, "error": "All ROI sheets must contain the same number of observations."}
            for a, b, da, db in zip(ends, reference_ends, durations, reference_durations):
                tol = 1e-6 * max(1.0, abs(a), abs(b), abs(da), abs(db))
                if abs(a - b) > tol or abs(da - db) > tol:
                    return {"ok": False, "error": "All ROI sheets must describe the same temporal grid."}

        timing_summaries.append(timing_summary)
        normalized_rois.append({"name": str(sheet_name), "rows": rows})

    if not normalized_rois:
        return {"ok": False, "error": "No non-empty ROI worksheets were found."}

    if not common_fit_ids:
        return {"ok": False, "error": "The ROI sheets do not share a common usable TAC value type."}

    order = ["Mean", "Median", "Peak", "Max", "Min", "VolumePET"]

    def display_label(stat_id):
        labels = labels_by_id.get(stat_id, [stat_id])
        first = labels[0] if labels else stat_id
        if all(str(x) == str(first) for x in labels):
            return str(first)
        return stat_id

    fit_stats = [
        {"id": stat_id, "label": display_label(stat_id)}
        for stat_id in order
        if stat_id in common_fit_ids
    ]
    plot_stats = [
        {"id": stat_id, "label": display_label(stat_id)}
        for stat_id in order
        if stat_id in common_plot_ids
    ]

    timing_summary = timing_summaries[0] if len(set(timing_summaries)) == 1 else "mixed source columns, common normalized grid"
    suggested_activity_unit = "SUVbw" if suggested_units and all(x == "SUVbw" for x in suggested_units) else ""

    # Convert any numpy scalar metadata to plain Python/Qt-friendly values.
    clean_metadata = {}
    for key, value in metadata.items():
        if hasattr(value, "item"):
            value = value.item()
        clean_metadata[str(key)] = value

    return {
        "ok": True,
        "rois": normalized_rois,
        "fit_stats": fit_stats,
        "plot_stats": plot_stats,
        "plot_times_sec": [float(x) for x in reference_plot_times],
        "framing_exact": bool(framing_exact_all),
        "timing_summary": timing_summary,
        "metadata": clean_metadata,
        "suggested_activity_unit": suggested_activity_unit,
    }

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
  )DPEPY2");

  dpePythonScript += QString::fromUtf8(R"DPEPY3(

        # 2. Separate spatial construction sources from provenance.

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
  )DPEPY3");

  dpePythonScript += QString::fromUtf8(R"DPEPY4(
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
  )DPEPY4");

  mainContext.evalScript(dpePythonScript);

  // Small Python bridge for RTSTRUCT import/export. The DICOM work remains in
  // the scripted dRTImporter/dRTExporter helpers; C++ only passes
  // MRML node IDs and displays the result.
  mainContext.evalScript(R"PYTHON(
import importlib
import importlib.util
import os
import traceback
import slicer

def _DPE_load_dynamic_rt_module(module_name, file_name):
    errors = []
    # Support both historical flat installs and the current dRTImporterLib
    # package layout used by dRTImporter.py.
    for import_name in (module_name, 'dRTImporterLib.' + module_name):
        try:
            return importlib.import_module(import_name)
        except Exception as exc:
            errors.append(f'{import_name}: {exc}')

    candidates = []

    # dRTExporter.py is installed next to dRTImporterPlugin.py. Resolve the
    # helper from either supported plugin import layout first.
    for plugin_name in ('dRTImporterPlugin', 'dRTImporterLib.dRTImporterPlugin'):
        try:
            importer_plugin = importlib.import_module(plugin_name)
            importer_plugin_path = getattr(importer_plugin, '__file__', '')
            if importer_plugin_path:
                candidates.append(
                    os.path.join(os.path.dirname(importer_plugin_path), file_name))
        except Exception as exc:
            errors.append(f'{plugin_name}: {exc}')

    # Fallback for installations where the scripted module is known to Slicer
    # but its Python directory is not currently on sys.path.
    try:
        module_path = slicer.util.modulePath('dRTImporter')
        if module_path:
            module_dir = os.path.dirname(module_path)
            candidates.append(os.path.join(module_dir, file_name))
            candidates.append(os.path.join(module_dir, 'dRTImporterLib', file_name))
    except Exception as exc:
        errors.append(str(exc))

    module_object = getattr(slicer.modules, 'drtimporter', None)
    module_path = getattr(module_object, 'path', '') if module_object else ''
    if module_path:
        module_dir = os.path.dirname(module_path)
        candidates.append(os.path.join(module_dir, file_name))
        candidates.append(os.path.join(module_dir, 'dRTImporterLib', file_name))

    seen = set()
    for candidate in candidates:
        candidate = os.path.abspath(candidate)
        if candidate in seen:
            continue
        seen.add(candidate)
        if not os.path.isfile(candidate):
            continue
        spec = importlib.util.spec_from_file_location(module_name, candidate)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    raise ImportError(
        f"Could not locate {file_name}. Ensure it is installed beside "
        "dRTImporterPlugin.py (add dRTExporter.py to MODULE_PYTHON_SCRIPTS "
        "in dRTImporter/CMakeLists.txt, then rebuild/reinstall dRTImporter). "
        + ("Previous import errors: " + " | ".join(errors) if errors else ""))

def DPE_adopt_dynamic_rtstruct(segmentation_node_id, pet_sequence_node_id, pet_browser_node_id):
    try:
        module = _DPE_load_dynamic_rt_module('dRTImporterPlugin', 'dRTImporterPlugin.py')
        return module.adopt_dynamic_rtstruct_to_pet(
            segmentation_node_id, pet_sequence_node_id, pet_browser_node_id)
    except Exception as exc:
        return {'ok': False, 'error': str(exc) + '\n\n' + traceback.format_exc()}

def DPE_export_dynamic_rtstruct(segmentation_sequence_node_id, pet_sequence_node_id, reference_volume_node_id, output_path, overwrite=False):
    try:
        module = _DPE_load_dynamic_rt_module('dRTExporter', 'dRTExporter.py')
        path = module.export_dynamic_rtstruct_from_node_ids(
            segmentation_sequence_node_id,
            pet_sequence_node_id,
            reference_volume_node_id,
            output_path,
            overwrite=bool(overwrite),
            show_progress=True)
        return {'ok': True, 'path': path}
    except Exception as exc:
        return {'ok': False, 'error': str(exc) + '\n\n' + traceback.format_exc()}

def DPE_export_static_rtstruct(segmentation_node_id, reference_volume_node_id, output_path, overwrite=False):
    try:
        module = _DPE_load_dynamic_rt_module('dRTExporter', 'dRTExporter.py')
        path = module.export_static_rtstruct_from_node_ids(
            segmentation_node_id,
            reference_volume_node_id,
            output_path,
            overwrite=bool(overwrite),
            show_progress=True)
        return {'ok': True, 'path': path}
    except Exception as exc:
        return {'ok': False, 'error': str(exc) + '\n\n' + traceback.format_exc()}

def DPE_open_segment_editor(segmentation_node_id, source_volume_node_id=''):
    try:
        segmentation_node = slicer.mrmlScene.GetNodeByID(segmentation_node_id)
        if segmentation_node is None:
            raise ValueError('Selected segmentation node is no longer available.')
        source_volume = (
            slicer.mrmlScene.GetNodeByID(source_volume_node_id)
            if source_volume_node_id else None)
        slicer.util.selectModule('SegmentEditor')
        module_widget = slicer.modules.segmenteditor.widgetRepresentation()
        if module_widget is None:
            raise RuntimeError('Segment Editor module widget is unavailable.')
        editor = module_widget.self().editor
        editor.setSegmentationNode(segmentation_node)
        if source_volume is not None:
            editor.setSourceVolumeNode(source_volume)
        return {'ok': True}
    except Exception as exc:
        return {'ok': False, 'error': str(exc) + '\n\n' + traceback.format_exc()}
)PYTHON");

  for (const QString& name : q->ModelsNamesMTGA)
  {
    auto mtgaTooltip = [](const QString& modelName) -> QString
    {
      if (modelName == "Patlak")
        return QObject::tr("Standard Patlak. Requires plasma input from injection. A delayed tissue acquisition is allowed when a complete measured or PBIF-reconstructed input is available.");
      if (modelName == "Relative Patlak")
        return QObject::tr("Relative Patlak for late/partial acquisitions. The plasma integral restarts at the selected late-time start t*, so early input is not required; the slope is relative Ki' rather than absolute Ki. Zuo, Qi & Wang, Phys Med Biol 2018, 63:165004.");
      if (modelName == "Logan")
        return QObject::tr("Standard Logan. Requires both tissue and plasma histories from the early acquisition.");
      if (modelName == "RE")
        return QObject::tr("Reversible-equilibrium Logan (RE). Requires tissue and plasma integrals from injection and a user-selected equilibrium start for the regression.");
      if (modelName == "Relative RE")
        return QObject::tr("Relative reversible-equilibrium analysis for late acquisitions. Both tissue and plasma integrals restart at the selected equilibrium start t*. The slope is relative DV_T' rather than absolute DV_T. Tian et al., Phys Med Biol 2024, 69:165005.");
      return QString();
    };

    QCheckBox* cb = new QCheckBox(name, this->ModelsMTGACheckContents);
    cb->setToolTip(mtgaTooltip(name));
    this->ModelsMTGACheckLayout->addWidget(cb);
    QObject::connect(cb, SIGNAL(stateChanged(int)),
                q, SLOT(onModelsMTGAChanged()));
    QCheckBox* cb2 = new QCheckBox(name, this->ModelsMTGACheckContents);
    cb2->setToolTip(mtgaTooltip(name));
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
      [this](const QString& path)
      {
        this->propagateOutputDirectory(path);
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
      [this](const QString& path)
      {
        this->propagateOutputDirectory(path);
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
  this->rebuildTACStatisticUI();

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

            std::string displayName =
                shNode->GetItemName(itemID);

            if (requiredModality == "PT")
            {
                // A dynamic PET entry is the proxy of a dPET master sequence.
                // Do not use Sequences.BaseName or the proxy name as a
                // discovery criterion: both are naming/UI details and depend
                // on the Sequence Browser Rename option.
                bool isDynamicPETProxy = false;

                for (int browserIndex = 0;
                     browserIndex < scene->GetNumberOfNodesByClass(
                         "vtkMRMLSequenceBrowserNode");
                     ++browserIndex)
                {
                    vtkMRMLSequenceBrowserNode* browser =
                        vtkMRMLSequenceBrowserNode::SafeDownCast(
                            scene->GetNthNodeByClass(
                                browserIndex,
                                "vtkMRMLSequenceBrowserNode"));

                    if (!browser)
                        continue;

                    vtkMRMLSequenceNode* masterSequence =
                        browser->GetMasterSequenceNode();
                    if (!masterSequence)
                        continue;

                    vtkMRMLNode* proxyNode =
                        browser->GetProxyNode(masterSequence);
                    if (proxyNode != dataNode)
                        continue;

                    const char* proxyLoadedBy =
                        dataNode->GetAttribute(
                            "dPETImporter.LoadedBy");
                    const char* sequenceLoadedBy =
                        masterSequence->GetAttribute(
                            "dPETImporter.LoadedBy");

                    const bool loadedByDPET =
                        (proxyLoadedBy &&
                         std::string(proxyLoadedBy) ==
                             "dPETImporterPlugin") ||
                        (sequenceLoadedBy &&
                         std::string(sequenceLoadedBy) ==
                             "dPETImporterPlugin");

                    if (!loadedByDPET)
                        continue;

                    isDynamicPETProxy = true;

                    // Keep the selector label stable across frame changes.
                    // The sequence name represents the dynamic study; the
                    // proxy name may intentionally stay fixed when Rename is
                    // disabled.
                    if (masterSequence->GetName() &&
                        std::string(masterSequence->GetName()).size() > 0)
                    {
                        displayName = masterSequence->GetName();
                    }
                    break;
                }

                modalityMatches =
                    modalityMatches &&
                    isDynamicPETProxy;
            }

            if (modalityMatches)
            {
                comboBox->addItem(
                    QString::fromStdString(displayName),
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
    QObject::connect(
        checkbox, &QCheckBox::toggled, q,
        [this, checkbox](bool checked)
        {
            if (!checked)
            {
                return;
            }
            this->lastPlotSegmentID =
                checkbox->property("SegmentID").toString().toStdString();
            if (this->plotDistributionSelected())
            {
                this->enforceDistributionSelection();
                this->refreshDistributionPlotIfActive();
            }
        });
  }

  this->PlotsegmentCheckLayout->addStretch();
  this->PlotSegmentCheckContents->blockSignals(false);
}

void qSlicerDynamicPETModuleWidgetPrivate::populateTimeBarMTGA(bool resetRange)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (q->numberOfTimepoints <= 0 ||
      q->timePoints.empty() ||
      q->durations.size() != q->timePoints.size())
  {
    return;
  }

  const int lastFrame = q->numberOfTimepoints;
  const int lastValidStartFrame = std::max(1, lastFrame - 1);

  this->timeOffsetSlider->setMinimum(1);
  this->timeOffsetSlider->setMaximum(lastValidStartFrame);
  if (resetRange ||
      this->timeOffsetSlider->value() < 1 ||
      this->timeOffsetSlider->value() > lastValidStartFrame)
  {
    this->timeOffsetSlider->setValue(1);
  }

  this->timeEndSlider->setMinimum(1);
  this->timeEndSlider->setMaximum(lastFrame);
  if (resetRange ||
      this->timeEndSlider->value() < 1 ||
      this->timeEndSlider->value() > lastFrame)
  {
    this->timeEndSlider->setValue(lastFrame);
  }

  this->TCMEndSlider->setMinimum(1);
  this->TCMEndSlider->setMaximum(lastFrame);
  if (resetRange ||
      this->TCMEndSlider->value() < 1 ||
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
    const double endSec = this->frameEndForInputSec(index - 1);
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
    const double endSec = this->frameEndForInputSec(index - 1);
    this->TCMEndInfoEdit->setText(
        QObject::tr("Frame %1 | end %2 s (%3 min)")
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

void qSlicerDynamicPETModuleWidgetPrivate::populateTimeBarMTGAImg(bool resetRange)
{
  Q_Q(qSlicerDynamicPETModuleWidget);

  if (q->numberOfTimepoints <= 0 ||
      q->timePoints.empty() ||
      q->durations.size() != q->timePoints.size())
  {
    return;
  }

  const int lastFrame = q->numberOfTimepoints;
  const int lastValidStartFrame = std::max(1, lastFrame - 1);

  this->timeOffsetSliderImg->setMinimum(1);
  this->timeOffsetSliderImg->setMaximum(lastValidStartFrame);
  if (resetRange ||
      this->timeOffsetSliderImg->value() < 1 ||
      this->timeOffsetSliderImg->value() > lastValidStartFrame)
  {
    this->timeOffsetSliderImg->setValue(1);
  }

  this->timeEndSliderImg->setMinimum(1);
  this->timeEndSliderImg->setMaximum(lastFrame);
  if (resetRange ||
      this->timeEndSliderImg->value() < 1 ||
      this->timeEndSliderImg->value() > lastFrame)
  {
    this->timeEndSliderImg->setValue(lastFrame);
  }

  this->TCMEndSliderImg->setMinimum(1);
  this->TCMEndSliderImg->setMaximum(lastFrame);
  if (resetRange ||
      this->TCMEndSliderImg->value() < 1 ||
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
    const double endSec = this->frameEndForInputSec(index - 1);
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
    const double endSec = this->frameEndForInputSec(index - 1);
    this->TCMEndInfoEditImg->setText(
        QObject::tr("Frame %1 | end %2 s (%3 min)")
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

void
qSlicerDynamicPETModuleWidgetPrivate::
refreshFitRangeSliderLabels()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (q->timePoints.empty())
    {
        return;
    }

    if (this->timeOffsetSlider && this->timeOffsetSlider->value() >= 1)
    {
        q->onSliderChanged(this->timeOffsetSlider->value());
    }
    if (this->timeOffsetSliderImg && this->timeOffsetSliderImg->value() >= 1)
    {
        q->onSliderImgChanged(this->timeOffsetSliderImg->value());
    }

    auto setEndText =
        [this](QSlider* slider, QLineEdit* edit)
        {
            if (!slider || !edit) return;
            const int frame = slider->value();
            if (frame < 1) return;
            const size_t index = static_cast<size_t>(frame - 1);
            const double endSec = this->frameEndForInputSec(index);
            edit->setText(
                QObject::tr("Frame %1 | end %2 s (%3 min)")
                    .arg(frame)
                    .arg(endSec, 0, 'f', 2)
                    .arg(endSec / 60.0, 0, 'f', 2));
        };

    setEndText(this->timeEndSlider, this->timeEndInfoEdit);
    setEndText(this->TCMEndSlider, this->TCMEndInfoEdit);
    setEndText(this->timeEndSliderImg, this->timeEndInfoEditImg);
    setEndText(this->TCMEndSliderImg, this->TCMEndInfoEditImg);
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateFitRangeSliders(const InputFunctionResult& result)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const int totalFrames = static_cast<int>(q->timePoints.size());
    if (totalFrames < 2)
    {
        return;
    }

    const int firstUsableFrame =
        static_cast<int>(
            std::min(result.supportFrameStartIndex, q->timePoints.size() - 1)) + 1;

    const int lastUsableFrame =
        static_cast<int>(
            std::min(
                result.supportFrameCount > 0
                    ? result.supportFrameCount
                    : q->timePoints.size(),
                q->timePoints.size()));

    if (lastUsableFrame - firstUsableFrame + 1 < 2)
    {
        return;
    }

    const int lastStartFrame =
        std::max(firstUsableFrame, lastUsableFrame - 1);
    const int firstEndFrame =
        std::min(lastUsableFrame, firstUsableFrame + 1);

    auto updateStartSlider =
        [](QSlider* slider, int minimum, int maximum)
        {
            if (!slider) return;
            const int oldValue = slider->value();
            QSignalBlocker blocker(slider);
            slider->setMinimum(minimum);
            slider->setMaximum(maximum);
            if (oldValue < minimum || oldValue > maximum)
            {
                slider->setValue(minimum);
            }
            else
            {
                slider->setValue(oldValue);
            }
        };

    auto updateEndSlider =
        [](QSlider* slider, int minimum, int maximum)
        {
            if (!slider) return;
            const int oldValue = slider->value();
            const int oldMaximum = slider->maximum();
            const bool pinnedToEnd = oldValue >= oldMaximum;

            QSignalBlocker blocker(slider);
            slider->setMinimum(minimum);
            slider->setMaximum(maximum);

            if (pinnedToEnd || oldValue > maximum)
            {
                slider->setValue(maximum);
            }
            else if (oldValue < minimum)
            {
                slider->setValue(minimum);
            }
            else
            {
                slider->setValue(oldValue);
            }
        };

    updateStartSlider(this->timeOffsetSlider, firstUsableFrame, lastStartFrame);
    updateStartSlider(this->timeOffsetSliderImg, firstUsableFrame, lastStartFrame);

    updateEndSlider(this->timeEndSlider, firstEndFrame, lastUsableFrame);
    updateEndSlider(this->TCMEndSlider, firstEndFrame, lastUsableFrame);
    updateEndSlider(this->timeEndSliderImg, firstEndFrame, lastUsableFrame);
    updateEndSlider(this->TCMEndSliderImg, firstEndFrame, lastUsableFrame);

    this->refreshFitRangeSliderLabels();
}

void
qSlicerDynamicPETModuleWidgetPrivate::
resetAcquisitionTimingDisplay()
{
    this->acquisitionTiming = AcquisitionTimingContext{};
    if (this->AcquisitionTimingStatusLabel)
    {
        this->AcquisitionTimingStatusLabel->setText(
            QObject::tr("Acquisition timing: —"));
        this->AcquisitionTimingStatusLabel->setVisible(false);
    }
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
      *errorMessage = this->multiTimepointMode
          ? QObject::tr(
              "Input-function observation count does not match the merged Multi-timepoint TAC.")
          : QObject::tr(
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

    const std::vector<bool>& externalKeep = this->activeExternalIFKeep();
    const size_t retainedExternalSamples =
        externalKeep.size() == this->externalIFTimesSec.size()
        ? static_cast<size_t>(std::count(externalKeep.begin(), externalKeep.end(), true))
        : this->externalIFTimesSec.size();

    const bool primarySourceSelected =
        segmentSource
        ? !q->IFID.empty()
        : (this->externalIFTimesSec.size() >= 2 &&
           this->externalIFTimesSec.size() ==
               this->externalIFConcentrations.size() &&
           retainedExternalSamples >= 2);

    // Keep source selection available, but do not expose an active
    // processing pipeline when there is no actual IF selected.
    this->IFAdvancedCollapsibleButton->setEnabled(
        primarySourceSelected);

    if (!primarySourceSelected)
    {
        this->IFAdvancedCollapsibleButton->setCollapsed(true);
    }

    const int sourceProcessingIndex =
        this->IFSourceProcessingSelector->currentIndex();
    const bool lowessSelected = sourceProcessingIndex == 1;
    const bool gaussianSelected = sourceProcessingIndex == 2;
    const bool fengSelected = sourceProcessingIndex == 3;
    this->IFLowessSpanLabel->setVisible(lowessSelected);
    this->IFLowessSpanSpinBox->setVisible(lowessSelected);
    this->IFGaussianSigmaLabel->setVisible(gaussianSelected);
    this->IFGaussianSigmaSpinBox->setVisible(gaussianSelected);
    this->fengExtrapolationCheckBox->setVisible(fengSelected);
    this->fengExtrapolationCheckBox->setEnabled(
        fengSelected && !this->PBIFOptionCheckBox->isChecked());

    this->IFLabel->setText(
        this->tableBasedMode
        ? QObject::tr("IDIF table ROI:")
        : QObject::tr("IDIF segment:"));
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
    this->IFExportButton->setEnabled(
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
            if (this->IFSourceProcessingSelector->currentIndex() == 3)
            {
                status += QObject::tr(" | Feng tau=%1 min, lambdas=%2/%3/%4 min^-1")
                    .arg(result.fengParameters.tau, 0, 'g', 4)
                    .arg(result.fengParameters.lambda1, 0, 'g', 4)
                    .arg(result.fengParameters.lambda2, 0, 'g', 4)
                    .arg(result.fengParameters.lambda3, 0, 'g', 4);
                if (result.fengExtrapolationApplied)
                {
                    status += QObject::tr(" | Feng extrapolated %1 -> %2 s")
                        .arg(result.sourceMeasuredEndTimeSec, 0, 'g', 7)
                        .arg(result.sourceModeledEndTimeSec, 0, 'g', 7);
                }
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
            if (result.parentFractionExtrapolationApplied)
            {
                status += QObject::tr(" (%1 -> %2 s modeled)")
                    .arg(result.parentFractionMeasuredEndTimeSec, 0, 'g', 7)
                    .arg(result.parentFractionModeledEndTimeSec, 0, 'g', 7);
            }
        }

        if (parentPlasmaSource &&
            !result.hasWholeBlood)
        {
            status +=
                QObject::tr(
                    " | TCM requires companion whole blood");
        }

        // One authoritative usable-frame range drives all ROI/imaging
        // TCM/MTGA sliders. End sliders remain pinned to the last supported
        // frame if the user had previously left them at the end.
        this->updateFitRangeSliders(result);

        this->logToPythonConsole(
            QObject::tr("[SlicerDynamicPET IF] %1").arg(status));
        this->IFStatusLabel->clear();
        this->IFStatusLabel->setToolTip(status);
        this->IFStatusLabel->setVisible(false);
    }
    else
    {
        const bool pbifConfigurationIncomplete =
            this->PBIFOptionCheckBox->isChecked() &&
            (this->pbifTimesSec.size() < 2 ||
             this->pbifTimesSec.size() != this->pbifTemplateValues.size());

        if (pbifConfigurationIncomplete)
        {
            // Enabling PBIF before choosing a template is an incomplete setup,
            // not an error worth occupying the GUI with. A real error is still
            // reported if the user later tries to fit/preview with invalid data.
            this->IFStatusLabel->clear();
            this->IFStatusLabel->setToolTip(
                QObject::tr("Select a PBIF template CSV to complete PBIF setup."));
            this->IFStatusLabel->setVisible(false);
        }
        else
        {
            this->IFStatusLabel->setText(error);
            this->IFStatusLabel->setToolTip(error);
            this->IFStatusLabel->setVisible(!error.isEmpty());
        }
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
        this->standardFitCheckBoxImg->setChecked(true);

        this->weightedFitCheckBoxImg->setEnabled(false);
        this->weightFitCheckBoxImg->setEnabled(false);
    }
    else
    {
        this->weightedFitCheckBoxImg->setEnabled(true);
        this->weightFitCheckBoxImg->setEnabled(true);
    }

    this->updateKineticModelAvailability();
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

    if (this->tableBasedMode || this->multiTimepointMode)
    {
        this->PlotsTabWidget->setTabEnabled(imagingIndex, false);
        return;
    }

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
updateKineticModelAvailability()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    // Checking PBIF before choosing a template is only an incomplete setup
    // state. Preserve the current model selections until a real PBIF exists.
    if (this->PBIFOptionCheckBox->isChecked() &&
        (this->pbifTimesSec.size() < 2 ||
         this->pbifTimesSec.size() != this->pbifTemplateValues.size()))
    {
        return;
    }

    InputFunctionResult result;
    QString ignoredError;
    const bool inputValid =
        this->buildCurrentInputFunction(result, false, &ignoredError);

    const bool fullInput = inputValid && result.inputCoversFromInjection;
    const bool delayedTissue = this->acquisitionTiming.delayedAcquisition;
    const bool lateOrPartial =
        inputValid &&
        (delayedTissue || !fullInput || result.supportFrameStartIndex > 0);

    auto updateLayout =
        [&](QLayout* layout, bool tcmLayout)
        {
            if (!layout) return;
            for (int i = 0; i < layout->count(); ++i)
            {
                QCheckBox* cb = qobject_cast<QCheckBox*>(layout->itemAt(i)->widget());
                if (!cb) continue;
                const QString name = cb->text();
                bool enabled = inputValid;
                QString availability;

                if (tcmLayout)
                {
                    enabled = inputValid && fullInput && result.hasWholeBlood;

                    if (enabled && delayedTissue)
                    {
                        availability = result.inputCoverageReconstructedByPBIF
                            ? QObject::tr("Delayed tissue acquisition: PBIF reconstructs plasma/whole-blood input from injection, so compartment modeling is available. A warning is shown before fitting because the unobserved early tissue response can weaken microparameter identifiability.")
                            : QObject::tr("Delayed tissue acquisition: complete plasma/whole-blood input from injection is available, so compartment modeling is available. A warning is shown before fitting because the unobserved early tissue response can weaken microparameter identifiability.");
                    }
                    else
                    {
                        availability = enabled
                            ? QObject::tr("Complete plasma and whole-blood input from injection is available.")
                            : (!result.hasWholeBlood
                               ? QObject::tr("Compartment modeling requires total whole blood in addition to plasma. Provide/derive a compatible whole-blood input.")
                               : QObject::tr("Compartment modeling requires a complete input function from injection. Use a full external IF or a PBIF that reconstructs the missing early input."));
                    }
                }
                else if (name == "Patlak")
                {
                    enabled = inputValid && fullInput;
                    availability = enabled
                        ? QObject::tr("Standard Patlak is available because the plasma input covers injection onward.")
                        : QObject::tr("Standard Patlak requires the plasma integral from injection. Use Relative Patlak when early input is unavailable.");
                }
                else if (name == "Logan" || name == "RE")
                {
                    enabled = inputValid && fullInput && !delayedTissue && result.supportFrameStartIndex == 0;
                    availability = enabled
                        ? QObject::tr("Standard reversible graphical analysis is available because early tissue and plasma histories are present.")
                        : QObject::tr("Standard Logan/RE require early tissue and plasma integrals. Use Relative RE for a delayed or partial acquisition.");
                }
                else if (name == "Relative Patlak")
                {
                    enabled = lateOrPartial;
                    availability = enabled
                        ? QObject::tr("Relative Patlak is enabled for the current delayed/partial acquisition. The selected MTGA start is t* and output is Ki'. Zuo, Qi & Wang, Phys Med Biol 2018, 63:165004.")
                        : QObject::tr("Relative Patlak is intended for delayed or partial acquisitions; standard Patlak is available for this complete early-time dataset.");
                }
                else if (name == "Relative RE")
                {
                    enabled = lateOrPartial;
                    availability = enabled
                        ? QObject::tr("Relative RE is enabled for the current delayed/partial acquisition. Both integrals restart at the selected equilibrium start and output is DV_T'. Tian et al., Phys Med Biol 2024.")
                        : QObject::tr("Relative RE is intended for delayed or partial acquisitions; standard RE is available for this complete early-time dataset.");
                }

                cb->setEnabled(enabled);
                if (!cb->property("BaseScientificTooltipInitialized").toBool())
                {
                    cb->setProperty("BaseScientificTooltip", cb->toolTip());
                    cb->setProperty("BaseScientificTooltipInitialized", true);
                }
                const QString scientificTip = cb->property("BaseScientificTooltip").toString();
                cb->setToolTip(scientificTip + (scientificTip.isEmpty() ? QString() : QString("\n\n")) + availability);
                if (!enabled && cb->isChecked())
                {
                    cb->setChecked(false);
                }
            }
        };

    updateLayout(this->ModelsMTGACheckLayout, false);
    updateLayout(this->ModelsCheckLayoutMTGAImg, false);
    updateLayout(this->ModelsCheckLayout, true);
    updateLayout(this->ModelsCheckLayoutTCMImg, true);
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

std::vector<bool>&
qSlicerDynamicPETModuleWidgetPrivate::
activeExternalIFKeep()
{
    return this->tableBasedMode
        ? this->tableExternalIFKeep
        : this->imageExternalIFKeep;
}

const std::vector<bool>&
qSlicerDynamicPETModuleWidgetPrivate::
activeExternalIFKeep() const
{
    return this->tableBasedMode
        ? this->tableExternalIFKeep
        : this->imageExternalIFKeep;
}

void
qSlicerDynamicPETModuleWidgetPrivate::
updateAcquisitionTimingContext(bool logMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (this->AcquisitionTimingStatusLabel)
    {
        this->AcquisitionTimingStatusLabel->setVisible(false);
    }

    AcquisitionTimingContext ctx;

    if (this->tableBasedMode)
    {
        if (q->timePoints.empty() || q->durations.empty())
        {
            this->acquisitionTiming = ctx;
            if (this->AcquisitionTimingStatusLabel)
            {
                this->AcquisitionTimingStatusLabel->setText(
                    QObject::tr("Acquisition timing: —"));
            }
            return;
        }

        const double firstFrameStartSec =
            q->timePoints.front() - q->durations.front();

        bool metadataOffsetOk = false;
        const double metadataOffset =
            this->tableWorkbookMetadata
                .value("AcquisitionStartPostInjection_s")
                .toDouble(&metadataOffsetOk);

        if (metadataOffsetOk &&
            std::isfinite(metadataOffset) &&
            metadataOffset > 60.0)
        {
            const QString timeConvention =
                this->tableWorkbookMetadata.value("TimeConvention").toString();
            const bool metadataTimesPostInjection =
                timeConvention.contains("PostInjection", Qt::CaseInsensitive);

            // DynamicPET exports delayed studies directly on the post-injection
            // time axis. Older workbooks may still be scan-relative; support
            // both without ever applying the injection offset twice.
            ctx.timingAvailable = true;
            ctx.delayedAcquisition = true;
            ctx.tableTimesAlreadyPostInjection =
                metadataTimesPostInjection || firstFrameStartSec > 60.0;
            ctx.rawInjectionToAcquisitionOffsetSec = metadataOffset;
            ctx.acquisitionStartPostInjectionSec =
                ctx.tableTimesAlreadyPostInjection && firstFrameStartSec > 60.0
                ? firstFrameStartSec
                : metadataOffset;
            ctx.source = QObject::tr("_DynamicPET timing metadata");
        }
        else
        {
            ctx.timingAvailable = true;
            ctx.tableTimesAlreadyPostInjection = firstFrameStartSec > 60.0;
            ctx.rawInjectionToAcquisitionOffsetSec = firstFrameStartSec;
            ctx.delayedAcquisition = firstFrameStartSec > 60.0;
            ctx.acquisitionStartPostInjectionSec =
                ctx.delayedAcquisition ? firstFrameStartSec : 0.0;
            ctx.source = QObject::tr("table frame timing");
        }

        this->acquisitionTiming = ctx;
    }
    else
    {
        vtkMRMLNode* timingNode = nullptr;
        if (this->multiTimepointMode)
        {
            vtkMRMLScene* scene = q->mrmlScene();
            if (scene && !this->multiTimepointReferenceMetadataNodeID.isEmpty())
            {
                timingNode = scene->GetNodeByID(
                    this->multiTimepointReferenceMetadataNodeID.toUtf8().constData());
            }
        }
        else
        {
            timingNode = q->sequencePETNode;
        }

        if (!timingNode)
        {
            this->acquisitionTiming = ctx;
            if (this->AcquisitionTimingStatusLabel)
            {
                this->AcquisitionTimingStatusLabel->setText(
                    QObject::tr("Acquisition timing: —"));
            }
            return;
        }

        QString injectionText = nodeOrKineticMetadataText(
            timingNode,
            "RadiopharmaceuticalStartDateTime",
            QStringLiteral("RadiopharmaceuticalStartDateTime"));
        if (injectionText.trimmed().isEmpty())
        {
            injectionText = nodeOrKineticMetadataText(
                timingNode,
                "RadionuclideStartDateTime",
                QStringLiteral("RadionuclideStartDateTime"));
        }

        QString firstFrameText = nodeAttributeText(
            timingNode, "dPET.FirstFrameAcquisitionDateTime");
        if (firstFrameText.isEmpty())
        {
            firstFrameText = nodeAttributeText(
                timingNode, "dPET.AcquisitionStartDateTime");
        }

        const QDateTime injectionDT = parseDICOMDateTimeText(injectionText);
        const QDateTime firstFrameDT = parseDICOMDateTimeText(firstFrameText);

        bool rawOffsetOk = false;
        double rawOffsetSec = nodeOrKineticMetadataText(
            timingNode,
            "dPET.InjectionToAcquisitionOffsetSec",
            QStringLiteral("InjectionToAcquisitionOffsetSec"))
            .toDouble(&rawOffsetOk);

        if ((!rawOffsetOk || !std::isfinite(rawOffsetSec)) &&
            injectionDT.isValid() && firstFrameDT.isValid())
        {
            rawOffsetSec = static_cast<double>(
                injectionDT.msecsTo(firstFrameDT)) / 1000.0;
            rawOffsetOk = true;
        }

        if (rawOffsetOk && std::isfinite(rawOffsetSec))
        {
            ctx.timingAvailable = true;
            ctx.rawInjectionToAcquisitionOffsetSec = rawOffsetSec;
            ctx.injectionDateTime = injectionText;
            ctx.firstFrameDateTime = firstFrameText;
            ctx.source = nodeOrKineticMetadataText(
                timingNode,
                "dPET.InjectionDateTimeSource",
                QStringLiteral("InjectionDateTimeSource"));
            if (ctx.source.isEmpty())
            {
                ctx.source = QObject::tr("DICOM");
            }

            // Keep the same conservative rule used in Single Image mode. A
            // clear >=5 min delay shifts the merged scan-relative timeline onto
            // the post-injection input-function clock; smaller differences do
            // not redefine the kinetic origin.
            constexpr double ClearDelayThresholdSec = 300.0;
            if (rawOffsetSec >= ClearDelayThresholdSec)
            {
                ctx.delayedAcquisition = true;
                ctx.acquisitionStartPostInjectionSec = rawOffsetSec;
            }
        }

        this->acquisitionTiming = ctx;
    }

    if (this->AcquisitionTimingStatusLabel)
    {
        const bool showTimingLabel =
            !this->tableBasedMode &&
            !this->multiTimepointMode &&
            q->sequencePETNode != nullptr;
        this->AcquisitionTimingStatusLabel->setVisible(showTimingLabel);
        QString text = QObject::tr("Acquisition timing: —");
        QString tooltip = QObject::tr(
            "Timing classification used by IF/PBIF processing. "
            "For delayed DICOM studies, the displayed delay is the interval from "
            "the recorded radiopharmaceutical start time to the START of the first "
            "acquired PET frame (not its midpoint or end). Early-time means the "
            "analysis starts near t=0.");

        if (this->acquisitionTiming.timingAvailable ||
            (this->tableBasedMode && !q->timePoints.empty()))
        {
            if (this->acquisitionTiming.delayedAcquisition)
            {
                text = QObject::tr("Acquisition timing: Delayed start (+%1 min)")
                    .arg(this->acquisitionTiming.acquisitionStartPostInjectionSec / 60.0, 0, 'f', 1);
            }
            else
            {
                text = QObject::tr("Acquisition timing: Early-time (t≈0)");
                tooltip += QObject::tr(
                    " For DICOM image mode, a near-start injection timestamp is not "
                    "treated as quantitatively trusted; scan start simply remains analysis t=0.");
            }
        }

        this->AcquisitionTimingStatusLabel->setText(text);
        this->AcquisitionTimingStatusLabel->setToolTip(tooltip);
    }

    if (logMessage &&
        !this->pbifTimesSec.empty() &&
        !q->timePoints.empty())
    {
        const double startTime = std::max(
            this->pbifTimesSec.front(),
            this->currentObservedInputStartSec());
        const double endTime = std::min(
            this->pbifTimesSec.back(),
            this->frameEndForInputSec(q->timePoints.size() - 1));
        if (endTime > startTime)
        {
            QSignalBlocker startBlocker(this->PBIFCalibrationStartSpinBox);
            QSignalBlocker endBlocker(this->PBIFCalibrationEndSpinBox);
            this->PBIFCalibrationStartSpinBox->setValue(startTime);
            this->PBIFCalibrationEndSpinBox->setValue(endTime);
        }
    }

    if (!logMessage || !this->acquisitionTiming.timingAvailable)
    {
        return;
    }

    if (this->acquisitionTiming.delayedAcquisition)
    {
        this->logToPythonConsole(
            QObject::tr(
                "[SlicerDynamicPET timing] Delayed tissue acquisition detected: "
                "first acquired tissue frame starts about %1 s post injection (%2). "
                "IF/PBIF operations use the post-injection clock; standard/relative "
                "kinetic-model availability is updated from the current input support.")
                .arg(this->acquisitionTiming.acquisitionStartPostInjectionSec, 0, 'g', 10)
                .arg(this->acquisitionTiming.source));
    }
    else if (!this->tableBasedMode)
    {
        this->logToPythonConsole(
            QObject::tr(
                "[SlicerDynamicPET timing] Injection timing metadata found, but it is "
                "not used to shift kinetic time because no clear >=5 min post-injection "
                "acquisition delay was established."));
    }
}


double
qSlicerDynamicPETModuleWidgetPrivate::
frameTimeShiftForInputSec() const
{
    if (!this->acquisitionTiming.delayedAcquisition ||
        this->acquisitionTiming.tableTimesAlreadyPostInjection)
    {
        return 0.0;
    }
    return this->acquisitionTiming.acquisitionStartPostInjectionSec;
}


double
qSlicerDynamicPETModuleWidgetPrivate::
frameEndForInputSec(size_t frameIndex) const
{
    Q_Q(const qSlicerDynamicPETModuleWidget);
    if (frameIndex >= q->timePoints.size())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return q->timePoints[frameIndex] + this->frameTimeShiftForInputSec();
}


double
qSlicerDynamicPETModuleWidgetPrivate::
frameStartForInputSec(size_t frameIndex) const
{
    Q_Q(const qSlicerDynamicPETModuleWidget);
    if (frameIndex >= q->timePoints.size() || frameIndex >= q->durations.size())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return this->frameEndForInputSec(frameIndex) - q->durations[frameIndex];
}


double
qSlicerDynamicPETModuleWidgetPrivate::
frameMidForInputSec(size_t frameIndex) const
{
    Q_Q(const qSlicerDynamicPETModuleWidget);
    if (frameIndex >= q->timePoints.size() || frameIndex >= q->durations.size())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return this->frameEndForInputSec(frameIndex) - 0.5 * q->durations[frameIndex];
}


double
qSlicerDynamicPETModuleWidgetPrivate::
currentObservedInputStartSec() const
{
    Q_Q(const qSlicerDynamicPETModuleWidget);
    const int source = this->IFSourceSelector ? this->IFSourceSelector->currentIndex() : 0;

    if (source == 1 && !this->externalIFTimesSec.empty())
    {
        if (this->externalIFZeroAnchorAdded && this->externalIFTimesSec.size() > 1)
        {
            return this->externalIFTimesSec[1];
        }
        return this->externalIFTimesSec.front();
    }

    if (!q->timePoints.empty() && !q->durations.empty())
    {
        return this->frameStartForInputSec(0);
    }
    return 0.0;
}


void
qSlicerDynamicPETModuleWidgetPrivate::
updateSegmentationAdvancedUI()
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    const bool imageMode = !this->tableBasedMode;
    vtkMRMLSegmentationNode* segmentationNode = nullptr;
    if (imageMode
        && !this->multiTimepointMode
        && q->mrmlScene()
        && q->segID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    {
        vtkMRMLSubjectHierarchyNode* shNode =
            vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(q->mrmlScene());
        if (shNode)
        {
            segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
                shNode->GetItemDataNode(q->segID));
        }
    }

    const bool singleHasSegmentation = segmentationNode != nullptr;
    const bool multiHasPreparedSegmentation =
        imageMode
        && this->multiTimepointMode
        && this->multiTimepointPreparationValid
        && !this->preparedMultiTimepointObservations.empty();
    const bool hasDynamicSegmentation =
        singleHasSegmentation
        && q->segSequenceNode
        && q->segSequenceNode->GetNumberOfDataNodes() > 0
        && q->sequencePETNode
        && q->sequenceBrowserPETNode;
    const bool hasPreparedMultiExport =
        imageMode
        && this->multiTimepointMode
        && this->multiTimepointPreparationValid
        && !this->preparedMultiTimepointAcquisitions.empty();
    const bool hasRTStructExport =
        hasDynamicSegmentation || hasPreparedMultiExport;

    this->SegmentationAdvancedCollapsibleButton->setVisible(imageMode);
    const bool editableSegmentationAvailable =
        singleHasSegmentation || multiHasPreparedSegmentation;
    this->PlotLiveSegEdit->setEnabled(editableSegmentationAvailable);
    this->OpenSegmentEditorButton->setEnabled(editableSegmentationAvailable);
    this->DynamicRTStructDirectory->setEnabled(hasRTStructExport);
    this->DynamicRTStructFilename->setEnabled(hasRTStructExport);
    this->SaveDynamicRTStructButton->setText(
        hasPreparedMultiExport
            ? QObject::tr("Save Multi RTSTRUCTs")
            : QObject::tr("Save Dynamic RTSTRUCT"));
    this->SaveDynamicRTStructButton->setToolTip(
        hasPreparedMultiExport
            ? QObject::tr(
                "Save one RT Structure Set per prepared acquisition. Dynamic PET acquisitions use the temporal RTSTRUCT convention; static PET acquisitions use a conventional static RTSTRUCT. The filename field is used as the batch base name.")
            : QObject::tr(
                "Save the selected dynamic segmentation sequence as one temporal RT Structure Set."));
    const bool hasOutputPath =
        !this->DynamicRTStructDirectory->currentPath().trimmed().isEmpty()
        && !this->DynamicRTStructFilename->text().trimmed().isEmpty();
    this->SaveDynamicRTStructButton->setEnabled(
        hasRTStructExport && hasOutputPath);

    if (this->segmentationFrameLabel)
    {
        this->segmentationFrameLabel->setVisible(imageMode);
    }
    if (this->segmentationFrameWidget)
    {
        this->segmentationFrameWidget->setVisible(imageMode);
    }
    this->updateSegmentationFrameUI(false);
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateSegmentationFrameUI(bool resetRange)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!this->segmentationFrameSlider || !this->segmentationFrameInfoEdit)
    {
        return;
    }

    int count = 0;
    int desired = 1;
    if (this->multiTimepointMode)
    {
        count = static_cast<int>(this->preparedMultiTimepointObservations.size());
    }
    else if (q->sequencePETNode && q->sequenceBrowserPETNode)
    {
        count = q->sequencePETNode->GetNumberOfDataNodes();
        const int selected = q->sequenceBrowserPETNode->GetSelectedItemNumber();
        if (selected >= 0)
        {
            desired = selected + 1;
        }
    }

    this->updatingSegmentationFrameSlider = true;
    this->segmentationFrameSlider->setEnabled(count > 0);
    this->segmentationFrameSlider->setMinimum(1);
    this->segmentationFrameSlider->setMaximum(std::max(1, count));
    if (count <= 0)
    {
        this->segmentationFrameSlider->setValue(1);
        this->segmentationFrameInfoEdit->clear();
    }
    else
    {
        if (!resetRange)
        {
            desired = std::max(1, std::min(this->segmentationFrameSlider->value(), count));
        }
        else
        {
            desired = std::max(1, std::min(desired, count));
        }
        this->segmentationFrameSlider->setValue(desired);
    }
    this->updatingSegmentationFrameSlider = false;
    this->updateSegmentationFrameInfo();
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
updateSegmentationFrameInfo()
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!this->segmentationFrameSlider || !this->segmentationFrameInfoEdit)
    {
        return;
    }

    const int index = this->segmentationFrameSlider->value() - 1;
    if (index < 0)
    {
        this->segmentationFrameInfoEdit->clear();
        return;
    }

    double endSec = std::numeric_limits<double>::quiet_NaN();
    if (this->multiTimepointMode &&
        index < static_cast<int>(this->preparedMultiTimepointObservations.size()) &&
        !this->preparedMultiTimepointAcquisitions.empty() &&
        this->preparedMultiTimepointAcquisitions.front().start.isValid())
    {
        endSec = this->preparedMultiTimepointAcquisitions.front().start.msecsTo(
            this->preparedMultiTimepointObservations[static_cast<size_t>(index)].end) / 1000.0;
    }
    else if (index < static_cast<int>(q->timePoints.size()))
    {
        endSec = this->frameEndForInputSec(static_cast<size_t>(index));
    }

    if (std::isfinite(endSec))
    {
        this->segmentationFrameInfoEdit->setText(
            QObject::tr("Frame %1 | end %2 s (%3 min)")
                .arg(index + 1)
                .arg(endSec, 0, 'f', 2)
                .arg(endSec / 60.0, 0, 'f', 2));
    }
    else
    {
        this->segmentationFrameInfoEdit->setText(
            QObject::tr("Frame %1").arg(index + 1));
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
displaySelectedSegmentationFrame()
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!q->mrmlScene() || !this->segmentationFrameSlider)
    {
        return;
    }

    vtkMRMLScalarVolumeNode* displayedPET = nullptr;

    if (!this->multiTimepointMode)
    {
        if (!q->sequenceBrowserPETNode || !q->sequencePETNode)
        {
            return;
        }
        const int frameIndex = this->segmentationFrameSlider->value() - 1;
        if (frameIndex < 0 || frameIndex >= q->sequencePETNode->GetNumberOfDataNodes())
        {
            return;
        }
        q->sequenceBrowserPETNode->SetSelectedItemNumber(frameIndex);
        displayedPET = vtkMRMLScalarVolumeNode::SafeDownCast(
            q->sequenceBrowserPETNode->GetProxyNode(q->sequencePETNode));
    }
    else
    {
        const int observationIndex = this->segmentationFrameSlider->value() - 1;
        if (observationIndex < 0 ||
            observationIndex >= static_cast<int>(this->preparedMultiTimepointObservations.size()))
        {
            return;
        }

        const PreparedMultiTimepointObservation& observation =
            this->preparedMultiTimepointObservations[static_cast<size_t>(observationIndex)];
        if (observation.dynamic)
        {
            if (observation.acquisitionIndex < 0 ||
                observation.acquisitionIndex >= static_cast<int>(this->preparedMultiTimepointAcquisitions.size()))
            {
                return;
            }
            vtkMRMLSequenceBrowserNode* browser = vtkMRMLSequenceBrowserNode::SafeDownCast(
                q->mrmlScene()->GetNodeByID(
                    this->preparedMultiTimepointAcquisitions[
                        static_cast<size_t>(observation.acquisitionIndex)].petBrowserNodeID.toUtf8().constData()));
            vtkMRMLSequenceNode* petSequence = vtkMRMLSequenceNode::SafeDownCast(
                q->mrmlScene()->GetNodeByID(observation.petSequenceNodeID.toUtf8().constData()));
            if (browser && petSequence && observation.frameIndex >= 0)
            {
                browser->SetSelectedItemNumber(observation.frameIndex);
                displayedPET = vtkMRMLScalarVolumeNode::SafeDownCast(
                    browser->GetProxyNode(petSequence));
            }
        }
        else
        {
            displayedPET = vtkMRMLScalarVolumeNode::SafeDownCast(
                q->mrmlScene()->GetNodeByID(observation.petNodeID.toUtf8().constData()));
        }
    }

    if (!displayedPET)
    {
        return;
    }

    vtkSlicerApplicationLogic* appLogic = qSlicerApplication::application()->applicationLogic();
    if (!appLogic || !appLogic->GetMRMLScene())
    {
        return;
    }
    for (int i = 0;
         i < appLogic->GetMRMLScene()->GetNumberOfNodesByClass("vtkMRMLSliceCompositeNode");
         ++i)
    {
        vtkMRMLSliceCompositeNode* compositeNode = vtkMRMLSliceCompositeNode::SafeDownCast(
            appLogic->GetMRMLScene()->GetNthNodeByClass(i, "vtkMRMLSliceCompositeNode"));
        if (compositeNode)
        {
            compositeNode->SetBackgroundVolumeID(displayedPET->GetID());
        }
    }
}



//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
clearMultiTimepointSegmentationWatchers()
{
    for (auto& watcher : this->multiSegWatchers)
    {
        if (watcher)
        {
            watcher->Clear();
        }
    }
    this->multiSegWatchers.clear();
    this->multiDirtySegEdits.clear();
    this->pendingMultiSegStructureChanges.clear();
    this->multiSegStructureUpdateQueued = false;
    this->processingMultiSegmentationChanges = false;
    if (this->multiSegEditTimer)
    {
        this->multiSegEditTimer->stop();
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
setupMultiTimepointSegmentationWatchers()
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    this->clearMultiTimepointSegmentationWatchers();

    vtkMRMLScene* scene = q->mrmlScene();
    if (!scene || !this->multiTimepointMode || !this->multiTimepointPreparationValid)
    {
        return;
    }

    if (!this->multiSegEditTimer)
    {
        this->multiSegEditTimer = new QTimer(q);
        this->multiSegEditTimer->setSingleShot(true);
        this->multiSegEditTimer->setInterval(120);
        QObject::connect(
            this->multiSegEditTimer, &QTimer::timeout, q,
            [this]() { this->processQueuedMultiTimepointSegmentEdits(); });
    }

    this->multiSegWatchers.reserve(this->preparedMultiTimepointAcquisitions.size());
    for (size_t acquisitionIndex = 0;
         acquisitionIndex < this->preparedMultiTimepointAcquisitions.size();
         ++acquisitionIndex)
    {
        const PreparedMultiTimepointAcquisition& acquisition =
            this->preparedMultiTimepointAcquisitions[acquisitionIndex];
        vtkSmartPointer<SegmentationChangeWatcher> watcher =
            vtkSmartPointer<SegmentationChangeWatcher>::New();
        watcher->ContextIndex = static_cast<int>(acquisitionIndex);
        watcher->DetectionOnly = true;
        watcher->GetLogic = [q]()
        {
            return vtkSlicerDynamicPETLogic::SafeDownCast(q->logic());
        };
        watcher->OnSegmentContentChanged =
            [this](int contextIndex, int frameIndex, const std::string& segmentID)
            {
                this->queueMultiTimepointSegmentEdit(
                    contextIndex, frameIndex, segmentID);
            };
        watcher->OnSegmentStructureChangedDetailed =
            [this](int contextIndex,
                   int frameIndex,
                   vtkMRMLSegmentationNode* sourceNode,
                   const std::string& segmentID,
                   SegmentationChangeWatcher::StructureChangeType changeType)
            {
                this->queueMultiTimepointStructureChange(
                    contextIndex, frameIndex, sourceNode, segmentID, changeType);
            };

        vtkMRMLSegmentationNode* proxySegmentation =
            vtkMRMLSegmentationNode::SafeDownCast(
                scene->GetNodeByID(acquisition.segmentationNodeID.toUtf8().constData()));

        if (acquisition.dynamic)
        {
            watcher->Mode = SegmentationChangeWatcher::AcquisitionMode::Dynamic;
            watcher->browser = vtkMRMLSequenceBrowserNode::SafeDownCast(
                scene->GetNodeByID(acquisition.petBrowserNodeID.toUtf8().constData()));
            vtkMRMLSequenceNode* segmentationSequence = vtkMRMLSequenceNode::SafeDownCast(
                scene->GetNodeByID(acquisition.segmentationSequenceNodeID.toUtf8().constData()));
            watcher->SegmentationSequence = segmentationSequence;
            if (proxySegmentation)
            {
                watcher->ObserveSegmentationNode(
                    proxySegmentation,
                    SegmentationChangeWatcher::UseBrowserFrameIndex);
            }
            if (segmentationSequence)
            {
                const int frameCount = segmentationSequence->GetNumberOfDataNodes();
                for (int frame = 0; frame < frameCount; ++frame)
                {
                    vtkMRMLSegmentationNode* frameNode = vtkMRMLSegmentationNode::SafeDownCast(
                        segmentationSequence->GetNthDataNode(frame));
                    if (frameNode)
                    {
                        watcher->ObserveSegmentationNode(frameNode, frame);
                    }
                }
            }
        }
        else
        {
            watcher->Mode = SegmentationChangeWatcher::AcquisitionMode::Static;
            watcher->browser = nullptr;
            if (proxySegmentation)
            {
                watcher->ObserveSegmentationNode(
                    proxySegmentation,
                    SegmentationChangeWatcher::StaticFrameIndex);
            }
        }

        this->multiSegWatchers.push_back(watcher);
    }

    this->logToPythonConsole(QObject::tr(
        "[SlicerDynamicPET Segment tracking] Prepared %1 acquisition-specific watcher(s): dynamic acquisitions observe their proxy + temporal segmentation frames; static acquisitions observe one segmentation node.")
        .arg(static_cast<int>(this->multiSegWatchers.size())));
}

//-----------------------------------------------------------------------------
int
qSlicerDynamicPETModuleWidgetPrivate::
preparedObservationIndexForAcquisitionFrame(
    int acquisitionIndex,
    int frameIndex) const
{
    for (size_t i = 0; i < this->preparedMultiTimepointObservations.size(); ++i)
    {
        const PreparedMultiTimepointObservation& observation =
            this->preparedMultiTimepointObservations[i];
        if (observation.acquisitionIndex != acquisitionIndex)
        {
            continue;
        }
        if (observation.dynamic)
        {
            if (observation.frameIndex == frameIndex)
            {
                return static_cast<int>(i);
            }
        }
        else if (frameIndex == SegmentationChangeWatcher::StaticFrameIndex)
        {
            return static_cast<int>(i);
        }
    }
    return -1;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
queueMultiTimepointSegmentEdit(
    int acquisitionIndex,
    int frameIndex,
    const std::string& segmentID)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!this->multiTimepointMode || !this->multiTimepointPreparationValid ||
        this->multiTimepointPreparationRunning || this->multiTimepointExtractionRunning ||
        this->processingMultiSegmentationChanges || segmentID.empty())
    {
        return;
    }

    this->multiDirtySegEdits.insert(
        std::make_tuple(acquisitionIndex, frameIndex, segmentID));

    // Coalesce an event burst from one brush stroke into one live refresh.
    // Do not restart an active timer: continuous painting therefore refreshes
    // at most about every 120 ms instead of once per VTK modification event.
    if (this->multiSegEditTimer && !this->multiSegEditTimer->isActive())
    {
        this->multiSegEditTimer->start();
    }
    else if (!this->multiSegEditTimer)
    {
        QTimer::singleShot(0, q, [this]()
        {
            this->processQueuedMultiTimepointSegmentEdits();
        });
    }
}

//-----------------------------------------------------------------------------
bool
qSlicerDynamicPETModuleWidgetPrivate::
recomputePreparedMultiTimepointSegmentObservation(
    int observationIndex,
    const std::string& commonSegmentName,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!q->mrmlScene() ||
        observationIndex < 0 ||
        observationIndex >= static_cast<int>(this->preparedMultiTimepointObservations.size()))
    {
        if (errorMessage) *errorMessage = QObject::tr("Prepared observation is unavailable.");
        return false;
    }

    const PreparedMultiTimepointObservation& observation =
        this->preparedMultiTimepointObservations[static_cast<size_t>(observationIndex)];
    if (observation.acquisitionIndex < 0 ||
        observation.acquisitionIndex >= static_cast<int>(this->preparedMultiTimepointAcquisitions.size()))
    {
        if (errorMessage) *errorMessage = QObject::tr("Prepared acquisition provenance is invalid.");
        return false;
    }

    const auto localIt = observation.segmentIDsByName.find(commonSegmentName);
    if (localIt == observation.segmentIDsByName.end() || localIt->second.empty())
    {
        if (errorMessage) *errorMessage = QObject::tr("The edited segment is no longer mapped at this observation.");
        return false;
    }

    auto tacIt = q->segmentTACs.find(commonSegmentName);
    if (tacIt == q->segmentTACs.end() ||
        observationIndex >= static_cast<int>(tacIt->second.size()))
    {
        // No computed TAC exists yet; the watcher stays active but there is
        // nothing quantitative to refresh until the user computes TAC.
        return true;
    }

    vtkMRMLScene* scene = q->mrmlScene();
    vtkMRMLScalarVolumeNode* petVolume = nullptr;
    vtkMRMLSegmentationNode* segmentationNode = nullptr;
    if (observation.dynamic)
    {
        vtkMRMLSequenceNode* petSequence = vtkMRMLSequenceNode::SafeDownCast(
            scene->GetNodeByID(observation.petSequenceNodeID.toUtf8().constData()));
        vtkMRMLSequenceNode* segmentationSequence = vtkMRMLSequenceNode::SafeDownCast(
            scene->GetNodeByID(observation.segmentationSequenceNodeID.toUtf8().constData()));
        if (petSequence && segmentationSequence)
        {
            const std::string indexValue = observation.sequenceIndex.toStdString();
            petVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
                petSequence->GetDataNodeAtValue(indexValue));
            segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
                segmentationSequence->GetDataNodeAtValue(indexValue));
        }
    }
    else
    {
        petVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
            scene->GetNodeByID(observation.petNodeID.toUtf8().constData()));
        segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
            scene->GetNodeByID(observation.segmentationNodeID.toUtf8().constData()));
    }

    vtkSlicerDynamicPETLogic* logic = vtkSlicerDynamicPETLogic::SafeDownCast(q->logic());
    if (!petVolume || !segmentationNode || !segmentationNode->GetSegmentation() || !logic)
    {
        if (errorMessage) *errorMessage = QObject::tr("The edited PET/segmentation observation is unavailable.");
        return false;
    }

    vtkNew<vtkStringArray> segmentArray;
    segmentArray->InsertNextValue(localIt->second);
    vtkSmartPointer<vtkOrientedImageData> labelmap =
        vtkSmartPointer<vtkOrientedImageData>::New();
    vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(
        segmentationNode,
        petVolume,
        segmentArray,
        vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
        labelmap);

    VoxelStatistics stats;
    vtkDataArray* labelScalars = labelmap && labelmap->GetPointData()
        ? labelmap->GetPointData()->GetScalars() : nullptr;
    vtkImageData* petImage = petVolume->GetImageData();
    vtkDataArray* petScalars = petImage && petImage->GetPointData()
        ? petImage->GetPointData()->GetScalars() : nullptr;
    if (!labelScalars || !petScalars ||
        labelScalars->GetNumberOfTuples() != petScalars->GetNumberOfTuples())
    {
        stats.keep = false;
        stats.empty = true;
    }
    else
    {
        stats = logic->ComputeVoxelStatistics(petVolume, labelmap, 1);
        if (!stats.empty)
        {
            const PreparedMultiTimepointAcquisition& acquisition =
                this->preparedMultiTimepointAcquisitions[
                    static_cast<size_t>(observation.acquisitionIndex)];
            if (acquisition.valueType.trimmed().toUpper() == QStringLiteral("BQML"))
            {
                scaleVoxelStatistics(stats, acquisition.sourceSUVbwFactor);
            }
        }
    }

    const VoxelStatistics previousStats =
        tacIt->second[static_cast<size_t>(observationIndex)];
    stats.keep = stats.empty
        ? false
        : (previousStats.empty ? true : previousStats.keep);
    tacIt->second[static_cast<size_t>(observationIndex)] = stats;
    if (errorMessage) errorMessage->clear();
    return true;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
processQueuedMultiTimepointSegmentEdits()
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (this->multiDirtySegEdits.empty())
    {
        return;
    }

    const auto dirty = this->multiDirtySegEdits;
    this->multiDirtySegEdits.clear();

    if (!this->multiTimepointMode || !this->multiTimepointPreparationValid ||
        this->multiTimepointPreparationRunning || this->multiTimepointExtractionRunning ||
        this->processingMultiSegmentationChanges)
    {
        return;
    }

    this->processingMultiSegmentationChanges = true;
    bool anyTACChanged = false;
    bool inputSegmentChanged = false;
    QSet<int> changedObservationIndices;
    QStringList warnings;

    for (const auto& key : dirty)
    {
        const int acquisitionIndex = std::get<0>(key);
        const int frameIndex = std::get<1>(key);
        const std::string localSegmentID = std::get<2>(key);
        const int observationIndex = this->preparedObservationIndexForAcquisitionFrame(
            acquisitionIndex, frameIndex);
        if (observationIndex < 0)
        {
            continue;
        }

        const PreparedMultiTimepointObservation& observation =
            this->preparedMultiTimepointObservations[
                static_cast<size_t>(observationIndex)];
        std::string commonName;
        for (const auto& nameAndID : observation.segmentIDsByName)
        {
            if (nameAndID.second == localSegmentID)
            {
                commonName = nameAndID.first;
                break;
            }
        }
        if (commonName.empty())
        {
            continue;
        }

        QString error;
        if (!this->recomputePreparedMultiTimepointSegmentObservation(
                observationIndex, commonName, &error))
        {
            if (!error.isEmpty()) warnings << error;
            continue;
        }

        auto tacIt = q->segmentTACs.find(commonName);
        if (tacIt != q->segmentTACs.end() &&
            observationIndex < static_cast<int>(tacIt->second.size()))
        {
            anyTACChanged = true;
            changedObservationIndices.insert(observationIndex);
            if (commonName == q->IFID)
            {
                inputSegmentChanged = true;
            }
        }
    }

    if (anyTACChanged)
    {
        q->clearFITdata();
        q->clearFITMTGAdata();
        if (inputSegmentChanged)
        {
            this->invalidateInputFunctionResults();
            this->updateInputFunctionStatus();
        }
        else
        {
            q->enableFITbutton();
            q->enableFITMTGAbutton();
        }

        if (this->PlotLiveSegEdit && this->PlotLiveSegEdit->isChecked())
        {
            if (this->plotDistributionSelected())
            {
                const int displayedObservation = this->distributionFrameSlider
                    ? this->distributionFrameSlider->value() - 1 : -1;
                if (changedObservationIndices.contains(displayedObservation))
                {
                    this->refreshDistributionPlotIfActive();
                }
            }
            else
            {
                q->onPlotbutton();
            }
        }
    }

    for (const QString& warning : warnings)
    {
        this->logToPythonConsole(
            QObject::tr("[SlicerDynamicPET Segment tracking] %1").arg(warning));
    }
    this->processingMultiSegmentationChanges = false;
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
queueMultiTimepointStructureChange(
    int acquisitionIndex,
    int frameIndex,
    vtkMRMLSegmentationNode* sourceNode,
    const std::string& segmentID,
    SegmentationChangeWatcher::StructureChangeType changeType)
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    if (!this->multiTimepointMode || !this->multiTimepointPreparationValid ||
        this->processingMultiSegmentationChanges || !sourceNode || segmentID.empty())
    {
        return;
    }

    PendingMultiSegStructureChange pending;
    pending.acquisitionIndex = acquisitionIndex;
    pending.frameIndex = frameIndex;
    pending.sourceNode = sourceNode;
    pending.segmentID = segmentID;
    pending.changeType = changeType;
    this->pendingMultiSegStructureChanges.push_back(std::move(pending));

    if (!this->multiSegStructureUpdateQueued)
    {
        this->multiSegStructureUpdateQueued = true;
        QTimer::singleShot(0, q, [this]()
        {
            this->processQueuedMultiTimepointStructureChanges();
        });
    }
}

//-----------------------------------------------------------------------------
void
qSlicerDynamicPETModuleWidgetPrivate::
processQueuedMultiTimepointStructureChanges()
{
    Q_Q(qSlicerDynamicPETModuleWidget);
    this->multiSegStructureUpdateQueued = false;
    if (this->pendingMultiSegStructureChanges.empty())
    {
        return;
    }

    const std::vector<PendingMultiSegStructureChange> changes =
        this->pendingMultiSegStructureChanges;
    this->pendingMultiSegStructureChanges.clear();

    if (!this->multiTimepointMode || !this->multiTimepointPreparationValid ||
        this->multiTimepointPreparationRunning || this->multiTimepointExtractionRunning)
    {
        return;
    }

    this->processingMultiSegmentationChanges = true;
    for (const PendingMultiSegStructureChange& change : changes)
    {
        if (change.acquisitionIndex < 0 ||
            change.acquisitionIndex >= static_cast<int>(this->multiSegWatchers.size()))
        {
            continue;
        }
        vtkMRMLSegmentationNode* sourceNode = change.sourceNode.GetPointer();
        SegmentationChangeWatcher* watcher =
            this->multiSegWatchers[static_cast<size_t>(change.acquisitionIndex)];
        if (watcher)
        {
            watcher->ApplyDeferredStructureChange(
                sourceNode, change.segmentID, change.changeType);
        }
    }

    // Structure changes can alter the common ROI intersection. Do not rebuild
    // selections inside the VTK callback. Once Slicer has returned here, clear
    // stale quantitative results, revalidate the acquisition table, and
    // reprepare the same acquisition set. This is intentionally conservative
    // and avoids dangling segment IDs after add/remove/rename.
    q->clearTACdata();
    q->timePoints.clear();
    q->durations.clear();
    q->suvFactors.clear();
    q->numberOfTimepoints = 0;
    this->multiTimepointPreparationValid = false;
    this->processingMultiSegmentationChanges = false;

    this->updateMultiTimepointSelectionStatus();
    if (this->multiTimepointSelectionValidated)
    {
        QString error;
        if (!this->prepareMultiTimepointAcquisitions(&error))
        {
            this->MultiTimepointStatusLabel->setText(
                QObject::tr("Segmentation structure changed; preparation failed: %1").arg(error));
        }
        else
        {
            this->populateMultiTimepointCommonSegmentCheckboxes();
            this->MultiTimepointStatusLabel->setText(
                this->MultiTimepointStatusLabel->text() +
                QObject::tr(" Segmentation structure refreshed; recompute TAC."));
        }
    }
    else
    {
        this->populateMultiTimepointCommonSegmentCheckboxes();
    }
    q->enableTACbutton();
    this->updateSegmentationAdvancedUI();
}

void
qSlicerDynamicPETModuleWidgetPrivate::
propagateOutputDirectory(const QString& path)
{
    if (this->propagatingOutputDirectory)
    {
        return;
    }

    const QString trimmed = path.trimmed();
    if (trimmed.isEmpty())
    {
        return;
    }

    const QString cleanPath = QDir::cleanPath(trimmed);
    this->sharedOutputDirectory = cleanPath;
    this->propagatingOutputDirectory = true;

    for (ctkPathLineEdit* editor : {
             this->direxcel,
             this->direxceltcm,
             this->direxcelmtga,
             this->direxceltcmfitted,
             this->direxcelmtgafitted,
             this->DynamicRTStructDirectory,
             this->MTGADICOMDirectoryImg,
             this->TCMDICOMDirectoryImg})
    {
        if (editor && editor->currentPath() != cleanPath)
        {
            editor->setCurrentPath(cleanPath);
        }
    }

    this->propagatingOutputDirectory = false;
}


void
qSlicerDynamicPETModuleWidgetPrivate::
logToPythonConsole(const QString& message) const
{
    PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
    mainContext.call(
        "DPE_console_message",
        QVariantList{message});
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
exportFinalInputFunctionCSV(
    const QString& filePath,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    InputFunctionResult result;
    QString error;
    if (!this->buildCurrentInputFunction(result, false, &error))
    {
        if (errorMessage)
            *errorMessage = error;
        return false;
    }

    if (result.supportFrameCount < 1 ||
        result.supportFrameCount > q->timePoints.size())
    {
        if (errorMessage)
            *errorMessage = QObject::tr("The current input function has no supported TAC frames to export.");
        return false;
    }

    bool stepOk = false;
    double sampleStepSec = this->timeStepEdit->text().toDouble(&stepOk);
    if (!stepOk || !std::isfinite(sampleStepSec) || sampleStepSec <= 0.0)
        sampleStepSec = 1.0;
    sampleStepSec = std::max(0.01, sampleStepSec);

    const double endTimeSec =
        this->frameEndForInputSec(result.supportFrameCount - 1);
    if (!std::isfinite(endTimeSec) || endTimeSec <= 0.0)
    {
        if (errorMessage)
            *errorMessage = QObject::tr("The supported input-function end time is invalid.");
        return false;
    }

    const std::string plasmaInterpolation =
        result.pbifApplied ? std::string("linear") : this->selectedIFInterpolation();

    const double startTimeSec =
        result.inputCoversFromInjection
        ? 0.0
        : (std::isfinite(result.earliestAvailableInputTimeSec)
           ? std::max(0.0, result.earliestAvailableInputTimeSec)
           : 0.0);

    std::vector<double> exportTimes;
    std::vector<double> exportValuesNative;
    const size_t reserveCount =
        static_cast<size_t>(std::ceil((endTimeSec - startTimeSec) / sampleStepSec)) + 2;
    exportTimes.reserve(reserveCount);
    exportValuesNative.reserve(reserveCount);

    auto finalPlasmaAt =
        [&](double timeSec) -> double
        {
            double plasma = std::numeric_limits<double>::quiet_NaN();
            if (!result.nativePlasmaTimesSec.empty() &&
                result.nativePlasmaTimesSec.size() == result.nativePlasmaValues.size())
            {
                plasma = this->interpolateInputFunction(
                    result.nativePlasmaTimesSec,
                    result.nativePlasmaValues,
                    timeSec,
                    plasmaInterpolation);
            }
            else
            {
                plasma = this->evaluateFrameCurve(
                    result.framePlasma,
                    timeSec,
                    plasmaInterpolation);
            }

            if (!std::isfinite(plasma))
                return plasma;

            if (result.applyParentFraction)
            {
                const double fraction = this->interpolateInputFunction(
                    result.parentFractionTimesSec,
                    result.parentFractionValues,
                    timeSec,
                    "linear");
                if (!std::isfinite(fraction))
                    return std::numeric_limits<double>::quiet_NaN();
                plasma *= fraction;
            }
            return plasma;
        };

    for (double t = startTimeSec; t < endTimeSec; t += sampleStepSec)
    {
        const double value = finalPlasmaAt(t);
        if (!std::isfinite(value))
        {
            if (errorMessage)
                *errorMessage = QObject::tr("The final input function produced a non-finite value at %1 s.")
                    .arg(t, 0, 'g', 10);
            return false;
        }
        exportTimes.push_back(t);
        exportValuesNative.push_back(value);
    }

    if (exportTimes.empty() ||
        std::fabs(exportTimes.back() - endTimeSec) > 1e-9)
    {
        const double value = finalPlasmaAt(endTimeSec);
        if (!std::isfinite(value))
        {
            if (errorMessage)
                *errorMessage = QObject::tr("The final input function produced a non-finite value at the supported end time.");
            return false;
        }
        exportTimes.push_back(endTimeSec);
        exportValuesNative.push_back(value);
    }

    std::vector<double> exportValues;
    if (!this->convertActivityVector(
            exportValuesNative,
            this->petStoredActivityUnit(),
            this->selectedDisplayActivityUnit(),
            exportValues,
            errorMessage))
    {
        return false;
    }

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate))
    {
        if (errorMessage)
            *errorMessage = QObject::tr("Could not open the selected CSV for writing: %1")
                .arg(file.errorString());
        return false;
    }

    QTextStream stream(&file);
    stream.setRealNumberNotation(QTextStream::SmartNotation);
    stream.setRealNumberPrecision(12);
    stream << "time_s,concentration\n";
    for (size_t i = 0; i < exportTimes.size(); ++i)
    {
        stream << exportTimes[i] << ',' << exportValues[i] << '\n';
    }
    file.close();

    const QString domain =
        (result.plasmaIsParent || result.applyParentFraction)
        ? QObject::tr("parent plasma")
        : QObject::tr("total plasma");

    this->logToPythonConsole(
        QObject::tr("[SlicerDynamicPET IF] Exported final %1 to %2 | %3 samples | step=%4 s | unit=%5")
            .arg(domain)
            .arg(QDir::toNativeSeparators(filePath))
            .arg(exportTimes.size())
            .arg(sampleStepSec, 0, 'g', 8)
            .arg(this->activityUnitLabel(this->selectedDisplayActivityUnit())));

    return true;
}

ParentFractionModel
qSlicerDynamicPETModuleWidgetPrivate::
selectedParentFractionModel() const
{
    switch (this->ParentFractionProcessingSelector->currentIndex())
    {
      case 1:
        return ParentFractionModel::Hill;
      case 2:
        return ParentFractionModel::ExtendedHill;
      case 3:
        return ParentFractionModel::ExponentialPlateau;
      default:
        return ParentFractionModel::Linear;
    }
}

bool
qSlicerDynamicPETModuleWidgetPrivate::
buildProcessedParentFraction(
    double requiredEndTimeSec,
    std::vector<double>& processedTimesSec,
    std::vector<double>& processedValues,
    ParentFractionFitParameters* fitParameters,
    QString* processingLabel,
    QString* fitSummary,
    QString* errorMessage)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    processedTimesSec.clear();
    processedValues.clear();
    if (fitParameters)
        *fitParameters = ParentFractionFitParameters{};
    if (processingLabel)
        processingLabel->clear();
    if (fitSummary)
        fitSummary->clear();

    if (this->parentFractionTimesSec.size() < 2 ||
        this->parentFractionTimesSec.size() != this->parentFractionValues.size())
    {
        if (errorMessage)
            *errorMessage = QObject::tr("Load a valid time_s,parent_fraction CSV first.");
        return false;
    }

    const std::vector<double>& retainedTimes = this->parentFractionTimesSec;
    const std::vector<double>& retainedValues = this->parentFractionValues;

    const ParentFractionModel model = this->selectedParentFractionModel();
    if (model == ParentFractionModel::Linear)
    {
        processedTimesSec = retainedTimes;
        processedValues = retainedValues;
        if (processingLabel)
            *processingLabel = QObject::tr("Linear interpolation");
        return true;
    }

    // The automatically inserted (0,1) point defines interpolation support,
    // but it is not a measured metabolite sample and must not bias a nonlinear fit.
    std::vector<double> fitTimes = retainedTimes;
    std::vector<double> fitValues = retainedValues;
    if (this->parentFractionZeroAnchorAdded &&
        !fitTimes.empty() && fitTimes.front() == 0.0)
    {
        fitTimes.erase(fitTimes.begin());
        fitValues.erase(fitValues.begin());
    }

    vtkSlicerDynamicPETLogic* logic =
        vtkSlicerDynamicPETLogic::SafeDownCast(q->logic());
    if (!logic)
    {
        if (errorMessage)
            *errorMessage = QObject::tr("DynamicPET logic is unavailable.");
        return false;
    }

    ParentFractionFitParameters params;
    std::vector<double> fittedObservations;
    std::string fitError;
    if (!logic->FitParentFraction(
            fitTimes,
            fitValues,
            model,
            params,
            fittedObservations,
            &fitError))
    {
        if (errorMessage)
        {
            *errorMessage = QObject::tr("Parent-fraction fit failed: %1")
                .arg(QString::fromStdString(fitError));
        }
        return false;
    }

    const double endTimeSec = std::max(
        requiredEndTimeSec,
        retainedTimes.back());
    const double denseStepSec = 1.0;
    for (double t = 0.0; t < endTimeSec; t += denseStepSec)
    {
        const double value = logic->EvaluateParentFraction(t, model, params);
        if (!std::isfinite(value))
        {
            if (errorMessage)
                *errorMessage = QObject::tr("Parent-fraction model produced a non-finite value.");
            return false;
        }
        processedTimesSec.push_back(t);
        processedValues.push_back(value);
    }
    processedTimesSec.push_back(endTimeSec);
    processedValues.push_back(logic->EvaluateParentFraction(endTimeSec, model, params));

    QString label;
    QString summary;
    if (model == ParentFractionModel::Hill)
    {
        label = QObject::tr("Hill fit");
        summary = QObject::tr("Hill: A=%1, B=%2, C=%3, SSE=%4")
            .arg(params.A, 0, 'g', 5)
            .arg(params.B, 0, 'g', 5)
            .arg(params.C, 0, 'g', 5)
            .arg(params.SSE, 0, 'g', 5);
    }
    else if (model == ParentFractionModel::ExtendedHill)
    {
        label = QObject::tr("Extended Hill fit");
        summary = QObject::tr("Extended Hill: A=%1, B=%2, C=%3, D=%4, E=%5 min, SSE=%6")
            .arg(params.A, 0, 'g', 5)
            .arg(params.B, 0, 'g', 5)
            .arg(params.C, 0, 'g', 5)
            .arg(params.D, 0, 'g', 5)
            .arg(params.E, 0, 'g', 5)
            .arg(params.SSE, 0, 'g', 5);
    }
    else
    {
        label = QObject::tr("Exponential-to-plateau fit");
        summary = QObject::tr("Exponential: plateau=%1, k=%2 min^-1, SSE=%3")
            .arg(params.plateau, 0, 'g', 5)
            .arg(params.rate, 0, 'g', 5)
            .arg(params.SSE, 0, 'g', 5);
    }

    this->logToPythonConsole(
        QObject::tr("[SlicerDynamicPET parent fraction] %1").arg(summary));
    if (fitParameters)
        *fitParameters = params;
    if (processingLabel)
        *processingLabel = label;
    if (fitSummary)
        *fitSummary = summary;
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
    this->imageExternalIFKeep.assign(this->externalIFTimesSec.size(), true);
    this->tableExternalIFKeep.assign(this->externalIFTimesSec.size(), true);
    this->externalIFPreviewIndexMap.clear();
    this->externalIFPreviewTimesSec.clear();
    this->externalIFPreviewDisplayValues.clear();
    this->externalIFPreviewSelectedIndex = -1;
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
    if (interpolationType == "pchip")
    {
        return pchipInterpolate(times, values, targetTime);
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
                    times, values, 0.5 * (a + b), interpolationType);
            if (!std::isfinite(y))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
            integral += y * (b - a);
        }
        else if (interpolationType == "pchip")
        {
            const double ya = this->interpolateInputFunction(times, values, a, interpolationType);
            const double ym = this->interpolateInputFunction(times, values, 0.5 * (a + b), interpolationType);
            const double yb = this->interpolateInputFunction(times, values, b, interpolationType);
            if (!std::isfinite(ya) || !std::isfinite(ym) || !std::isfinite(yb))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
            integral += (b - a) * (ya + 4.0 * ym + yb) / 6.0;
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
        q->timePoints.empty() ||
        endTime <= startTime ||
        startTime < q->timePoints.front() - q->durations.front() - 1e-9 ||
        endTime > q->timePoints.back() + 1e-9)
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    double integral = 0.0;

    for (size_t i = 0;
         i < q->durations.size();
         ++i)
    {
        const double frameEnd = q->timePoints[i];
        const double frameStart = frameEnd - q->durations[i];

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
        targetTime < this->frameStartForInputSec(0) - 1e-9 ||
        targetTime > this->frameEndForInputSec(frameValues.size() - 1) + 1e-9)
    {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    if (interpolationType == "const")
    {
        for (size_t i = 0;
             i < q->durations.size();
             ++i)
        {
            const double frameEnd = this->frameEndForInputSec(i);

            if (targetTime < frameEnd ||
                i + 1 == q->durations.size())
            {
                return frameValues[i];
            }
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
            this->frameMidForInputSec(i));
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

    if (interpolationType == "pchip")
    {
        std::vector<double> splineTimes;
        std::vector<double> splineValues;
        splineTimes.reserve(midpoints.size() + 2);
        splineValues.reserve(frameValues.size() + 2);
        splineTimes.push_back(this->frameStartForInputSec(0));
        splineValues.push_back(frameValues.front());
        for (size_t i = 0; i < midpoints.size(); ++i)
        {
            splineTimes.push_back(midpoints[i]);
            splineValues.push_back(frameValues[i]);
        }
        splineTimes.push_back(this->frameEndForInputSec(frameValues.size() - 1));
        splineValues.push_back(frameValues.back());
        return pchipInterpolate(splineTimes, splineValues, targetTime);
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

double
qSlicerDynamicPETModuleWidgetPrivate::
integrateModelPlasmaOverInterval(
    const InputFunctionResult& result,
    double startTimeSec,
    double endTimeSec)
{
    Q_Q(qSlicerDynamicPETModuleWidget);

    if (!(endTimeSec > startTimeSec))
    {
        return 0.0;
    }

    const std::string interpolation = this->selectedIFInterpolation();

    // Preferred path: integrate the actual processed/native plasma curve.
    // This retains external/PBIF temporal detail across periods where no PET
    // tissue image was acquired.
    if (result.nativePlasmaTimesSec.size() >= 2 &&
        result.nativePlasmaTimesSec.size() == result.nativePlasmaValues.size() &&
        startTimeSec >= result.nativePlasmaTimesSec.front() - 1e-6 &&
        endTimeSec <= result.nativePlasmaTimesSec.back() + 1e-6)
    {
        if (result.plasmaIsParent || !result.applyParentFraction)
        {
            return this->integrateInputFunctionOverInterval(
                result.nativePlasmaTimesSec,
                result.nativePlasmaValues,
                startTimeSec,
                endTimeSec,
                interpolation);
        }

        if (result.parentFractionTimesSec.size() >= 2 &&
            result.parentFractionTimesSec.size() == result.parentFractionValues.size() &&
            startTimeSec >= result.parentFractionTimesSec.front() - 1e-6 &&
            endTimeSec <= result.parentFractionTimesSec.back() + 1e-6)
        {
            const double average =
                this->averagePlasmaTimesParentFractionOverInterval(
                    result.nativePlasmaTimesSec,
                    result.nativePlasmaValues,
                    false,
                    interpolation,
                    result.parentFractionTimesSec,
                    result.parentFractionValues,
                    startTimeSec,
                    endTimeSec);
            return std::isfinite(average)
                ? average * (endTimeSec - startTimeSec)
                : std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Frame-derived fallback. Acquired frame averages contribute their exact
    // rectangular area. A genuine acquisition gap is the only interval that
    // is inferred, using a linear bridge between the neighbouring plasma
    // frame averages. This deliberately does not collapse the gap.
    const std::vector<double>& values = result.frameModelPlasma;
    if (values.size() != q->durations.size() ||
        values.empty() || q->timePoints.size() != values.size())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double firstStart = this->frameStartForInputSec(0);
    const double lastEnd = this->frameEndForInputSec(values.size() - 1);
    if (startTimeSec < firstStart - 1e-6 || endTimeSec > lastEnd + 1e-6)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto integrateLinearBridge = [](
        double gapStart, double gapEnd,
        double valueStart, double valueEnd,
        double a, double b) -> double
    {
        if (!(gapEnd > gapStart) || !(b > a))
        {
            return 0.0;
        }
        const double span = gapEnd - gapStart;
        auto primitive = [&](double t)
        {
            const double x = t - gapStart;
            return valueStart * x +
                0.5 * (valueEnd - valueStart) * x * x / span;
        };
        return primitive(b) - primitive(a);
    };

    double integral = 0.0;
    for (size_t i = 0; i < values.size(); ++i)
    {
        const double frameStart = this->frameStartForInputSec(i);
        const double frameEnd = this->frameEndForInputSec(i);

        if (i > 0)
        {
            const double previousEnd = this->frameEndForInputSec(i - 1);
            if (frameStart > previousEnd + 1e-9)
            {
                const double a = std::max(startTimeSec, previousEnd);
                const double b = std::min(endTimeSec, frameStart);
                if (b > a)
                {
                    if (!std::isfinite(values[i - 1]) || !std::isfinite(values[i]))
                    {
                        return std::numeric_limits<double>::quiet_NaN();
                    }
                    integral += integrateLinearBridge(
                        previousEnd, frameStart,
                        values[i - 1], values[i], a, b);
                }
            }
        }

        const double a = std::max(startTimeSec, frameStart);
        const double b = std::min(endTimeSec, frameEnd);
        if (b > a)
        {
            if (!std::isfinite(values[i]))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
            integral += values[i] * (b - a);
        }

        if (frameEnd >= endTimeSec)
        {
            break;
        }
    }

    return integral;
}

double
qSlicerDynamicPETModuleWidgetPrivate::
initialModelPlasmaIntegralSec(
    const InputFunctionResult& result,
    double endTimeSec)
{
    if (!(endTimeSec > 0.0))
    {
        return 0.0;
    }

    if (result.nativePlasmaTimesSec.size() < 2 ||
        result.nativePlasmaTimesSec.size() != result.nativePlasmaValues.size())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (result.nativePlasmaTimesSec.front() > 1e-6 ||
        result.nativePlasmaTimesSec.back() + 1e-6 < endTimeSec)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const std::string interpolation = this->selectedIFInterpolation();

    // If plasma is already parent plasma, or no metabolite correction is
    // active, integrate the native model input directly.
    if (result.plasmaIsParent || !result.applyParentFraction)
    {
        return this->integrateInputFunctionOverInterval(
            result.nativePlasmaTimesSec,
            result.nativePlasmaValues,
            0.0,
            endTimeSec,
            interpolation);
    }

    if (result.parentFractionTimesSec.size() < 2 ||
        result.parentFractionTimesSec.size() != result.parentFractionValues.size() ||
        result.parentFractionTimesSec.front() > 1e-6 ||
        result.parentFractionTimesSec.back() + 1e-6 < endTimeSec)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double average = this->averagePlasmaTimesParentFractionOverInterval(
            result.nativePlasmaTimesSec,
            result.nativePlasmaValues,
            false,
            interpolation,
            result.parentFractionTimesSec,
            result.parentFractionValues,
            0.0,
            endTimeSec);

    return std::isfinite(average) ? average * endTimeSec
                                  : std::numeric_limits<double>::quiet_NaN();
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

    this->updateAcquisitionTimingContext(false);

    const double petEndTimeSec =
        this->frameEndForInputSec(q->timePoints.size() - 1);
    double modelEndTimeSec = petEndTimeSec;

    const std::string inputInterpolation =
        this->selectedIFInterpolation();

    const IFCurveDomain sourceDomain =
        this->selectedIFCurveDomain();

    auto fillFrameAverages =
        [&](const std::vector<double>& times,
            const std::vector<double>& values,
            const std::string& interpolation,
            std::vector<double>& output,
            size_t frameCount,
            size_t firstFrameIndex) -> bool
        {
            frameCount = std::min(frameCount, q->durations.size());
            firstFrameIndex = std::min(firstFrameIndex, frameCount);
            output.assign(
                q->durations.size(),
                std::numeric_limits<double>::quiet_NaN());

            for (size_t i = firstFrameIndex; i < frameCount; ++i)
            {
                const double frameEnd = this->frameEndForInputSec(i);
                const double frameStart = frameEnd - q->durations[i];

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

                output[i] = value;
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
        modelEndTimeSec = this->frameEndForInputSec(result.supportFrameCount - 1);
        result.earliestAvailableInputTimeSec = this->frameStartForInputSec(0);
        result.inputCoversFromInjection =
            result.earliestAvailableInputTimeSec <= 60.0 + 1e-6;

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
                    this->frameMidForInputSec(i));
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
                    this->frameMidForInputSec(i);

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
                    this->frameMidForInputSec(i);
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
                result.processedSourcePreviewTimesSec.push_back(this->frameEndForInputSec(i));
                result.processedSourcePreviewValues.push_back(smoothed[i]);
            }
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel =
                QObject::tr("LOWESS (span %1)")
                    .arg(this->IFLowessSpanSpinBox->value(), 0, 'f', 2);
        }
        else if (sourceProcessing == 2)
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
                    this->frameMidForInputSec(i);
                evaluationTimes.push_back(midpoint);
                if (result.frameKeep[i])
                {
                    retainedTimes.push_back(midpoint);
                    retainedValues.push_back(observedWholeBlood[i]);
                }
            }

            const std::vector<double> smoothed =
                gaussianKernelPredict(
                    retainedTimes,
                    retainedValues,
                    evaluationTimes,
                    this->IFGaussianSigmaSpinBox->value());
            if (smoothed.size() != result.supportFrameCount)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "Gaussian input-function smoothing failed.");
                }
                return false;
            }
            for (size_t i = 0; i < result.supportFrameCount; ++i)
            {
                result.frameWholeBlood[i] = smoothed[i];
                result.processedSourcePreviewTimesSec.push_back(this->frameEndForInputSec(i));
                result.processedSourcePreviewValues.push_back(smoothed[i]);
            }
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel =
                QObject::tr("Gaussian (sigma %1 s)")
                    .arg(this->IFGaussianSigmaSpinBox->value(), 0, 'g', 5);
        }
        else if (sourceProcessing == 3)
        {
            const double measuredSourceEndSec = modelEndTimeSec;
            const bool extendFengToPETSupport =
                this->fengExtrapolationCheckBox->isChecked() &&
                !this->PBIFOptionCheckBox->isChecked() &&
                measuredSourceEndSec + 1e-6 < petEndTimeSec;

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
                    this->frameMidForInputSec(i));
                retainedValues.push_back(observedWholeBlood[i]);
                retainedFrameStarts.push_back(this->frameStartForInputSec(i));
                retainedFrameEnds.push_back(this->frameEndForInputSec(i));
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

            if (extendFengToPETSupport)
            {
                result.supportFrameCount = q->durations.size();
                modelEndTimeSec = petEndTimeSec;
                result.fengExtrapolationApplied = true;
                result.sourceMeasuredEndTimeSec = measuredSourceEndSec;
                result.sourceModeledEndTimeSec = modelEndTimeSec;
                std::fill(result.frameKeep.begin(), result.frameKeep.end(), true);
                this->logToPythonConsole(
                    QObject::tr("[SlicerDynamicPET IF] Feng extrapolation enabled: measured source support ends at %1 s; analytic model is extended to PET support at %2 s.")
                        .arg(measuredSourceEndSec, 0, 'g', 8)
                        .arg(modelEndTimeSec, 0, 'g', 8));
            }

            for (size_t i = 0; i < result.supportFrameCount; ++i)
            {
                result.frameWholeBlood[i] =
                    logic->AverageFengInputFunction(
                        this->frameStartForInputSec(i),
                        this->frameEndForInputSec(i),
                        result.fengParameters);
            }

            const double supportEndSec = modelEndTimeSec;
            // Dense representation of the analytic Feng curve is independent
            // of the user-selected TCM integration step. The TCM backend may
            // subsequently sample this continuous representation more finely.
            const double denseStepSec = 0.25;
            const double denseStartSec =
                result.inputCoversFromInjection
                ? 0.0
                : std::max(0.0, result.earliestAvailableInputTimeSec);

            for (double t = denseStartSec; t < supportEndSec; t += denseStepSec)
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
                this->frameMidForInputSec(i);

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

        const std::vector<bool>& externalKeep = this->activeExternalIFKeep();
        std::vector<double> retainedExternalTimes;
        std::vector<double> retainedExternalValues;
        retainedExternalTimes.reserve(this->externalIFTimesSec.size());
        retainedExternalValues.reserve(this->externalIFConcentrations.size());

        for (size_t i = 0; i < this->externalIFTimesSec.size(); ++i)
        {
            const bool retained =
                externalKeep.size() == this->externalIFTimesSec.size()
                ? externalKeep[i]
                : true;
            if (!retained)
            {
                continue;
            }
            retainedExternalTimes.push_back(this->externalIFTimesSec[i]);
            retainedExternalValues.push_back(this->externalIFConcentrations[i]);
        }

        if (retainedExternalTimes.size() < 2)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "At least two external input-function observations must remain.");
            }
            return false;
        }

        double retainedIFObservedStartSec = retainedExternalTimes.front();
        if (this->externalIFZeroAnchorAdded)
        {
            const auto firstRealIt = std::find_if(
                retainedExternalTimes.begin(),
                retainedExternalTimes.end(),
                [](double t) { return t > 1e-6; });
            if (firstRealIt != retainedExternalTimes.end())
            {
                retainedIFObservedStartSec = *firstRealIt;
            }
        }
        const double retainedIFEndSec = retainedExternalTimes.back();
        result.earliestAvailableInputTimeSec = retainedIFObservedStartSec;
        result.inputCoversFromInjection =
            retainedIFObservedStartSec <= 60.0 + 1e-6;

        // A real first sample within the first minute is accepted as an
        // early-time input and the internal (0,0) anchor may support the very
        // first frames.  A much later first sample remains genuinely partial;
        // the artificial anchor must never make it look complete.
        const double retainedIFCoverageStartSec =
            result.inputCoversFromInjection ? 0.0 : retainedIFObservedStartSec;

        result.supportFrameStartIndex = q->timePoints.size();
        result.supportFrameCount = 0;
        for (size_t i = 0; i < q->timePoints.size(); ++i)
        {
            const double frameStart = this->frameStartForInputSec(i);
            const double frameEnd = this->frameEndForInputSec(i);
            if (frameStart + 1e-6 >= retainedIFCoverageStartSec &&
                frameEnd <= retainedIFEndSec + 1e-6)
            {
                if (result.supportFrameStartIndex == q->timePoints.size())
                {
                    result.supportFrameStartIndex = i;
                }
                result.supportFrameCount = i + 1;
            }
        }

        if (result.supportFrameStartIndex >= result.supportFrameCount ||
            result.supportFrameCount - result.supportFrameStartIndex < 2)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The retained input function does not cover at least two TAC frames.");
            }
            return false;
        }

        modelEndTimeSec = this->frameEndForInputSec(result.supportFrameCount - 1);
        if (modelEndTimeSec + 1e-6 < petEndTimeSec)
        {
            this->logToPythonConsole(
                QObject::tr("[SlicerDynamicPET IF] Retained IF ends at %1 s; fitting is shortened to TAC frame %2 ending at %3 s.")
                    .arg(retainedIFEndSec, 0, 'g', 8)
                    .arg(result.supportFrameCount)
                    .arg(modelEndTimeSec, 0, 'g', 8));
        }

        // External activity curves are converted to the active model's native
        // quantitative scale before any biological-domain processing.
        std::vector<double> externalPatientValues;
        if (!this->convertActivityVector(
                retainedExternalValues,
                this->selectedExternalIFActivityUnit(),
                this->petStoredActivityUnit(),
                externalPatientValues,
                errorMessage))
        {
            return false;
        }

        std::vector<double> processedPatientTimes =
            retainedExternalTimes;
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
            processedPatientValues =
                gaussianKernelPredict(
                    processedPatientTimes,
                    processedPatientValues,
                    processedPatientTimes,
                    this->IFGaussianSigmaSpinBox->value());
            if (processedPatientValues.size() != processedPatientTimes.size())
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "Gaussian input-function smoothing failed.");
                }
                return false;
            }
            result.sourceProcessingApplied = true;
            result.sourceProcessingLabel =
                QObject::tr("Gaussian (sigma %1 s)")
                    .arg(this->IFGaussianSigmaSpinBox->value(), 0, 'g', 5);
        }
        else if (sourceProcessing == 3)
        {
            const double measuredSourceEndSec = retainedIFEndSec;
            const bool extendFengToPETSupport =
                this->fengExtrapolationCheckBox->isChecked() &&
                !this->PBIFOptionCheckBox->isChecked() &&
                measuredSourceEndSec + 1e-6 < petEndTimeSec;

            if (sourceDomain == IFCurveDomain::ParentPlasma)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "Feng source modeling is intended for whole-blood or total-plasma input functions, not an already metabolite-corrected parent-plasma curve. Use None, LOWESS, or Gaussian smoothing for this source.");
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

            if (extendFengToPETSupport)
            {
                result.supportFrameCount = q->timePoints.size();
                modelEndTimeSec = petEndTimeSec;
                result.fengExtrapolationApplied = true;
                result.sourceMeasuredEndTimeSec = measuredSourceEndSec;
                result.sourceModeledEndTimeSec = modelEndTimeSec;
                this->logToPythonConsole(
                    QObject::tr("[SlicerDynamicPET IF] Feng extrapolation enabled: measured source support ends at %1 s; analytic model is extended to PET support at %2 s.")
                        .arg(measuredSourceEndSec, 0, 'g', 8)
                        .arg(modelEndTimeSec, 0, 'g', 8));
            }

            processedPatientTimes.clear();
            processedPatientValues.clear();
            // Dense representation of the analytic Feng curve is independent
            // of the user-selected TCM integration step. The TCM backend may
            // subsequently sample this continuous representation more finely.
            const double denseStepSec = 0.25;
            const double denseStartSec =
                result.inputCoversFromInjection
                ? 0.0
                : std::max(0.0, result.earliestAvailableInputTimeSec);
            for (double t = denseStartSec; t < modelEndTimeSec; t += denseStepSec)
            {
                processedPatientTimes.push_back(t);
                processedPatientValues.push_back(
                    logic->EvaluateFengInputFunction(t, result.fengParameters));
            }
            processedPatientTimes.push_back(modelEndTimeSec);
            processedPatientValues.push_back(
                logic->EvaluateFengInputFunction(
                    modelEndTimeSec,
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
                    modelEndTimeSec;

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
                            "the retained fit range.");
                }
                return false;
            }
        }

        const bool pbifWillReplacePatientCurve =
            this->PBIFOptionCheckBox->isChecked() && !result.plasmaIsParent;

        const bool plasmaFramesAvailable = fillFrameAverages(
            result.nativePlasmaTimesSec,
            result.nativePlasmaValues,
            inputInterpolation,
            result.framePlasma,
            result.supportFrameCount,
            result.supportFrameStartIndex);

        if (!plasmaFramesAvailable && !pbifWillReplacePatientCurve)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The external input could not be sampled over its supported tissue-frame interval.");
            }
            return false;
        }

        if (result.hasWholeBlood)
        {
            const bool wholeBloodFramesAvailable = fillFrameAverages(
                result.nativeWholeBloodTimesSec,
                result.nativeWholeBloodValues,
                inputInterpolation,
                result.frameWholeBlood,
                result.supportFrameCount,
                result.supportFrameStartIndex);

            if (!wholeBloodFramesAvailable && !pbifWillReplacePatientCurve)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "The external whole-blood input could not be sampled over its supported tissue-frame interval.");
                }
                return false;
            }
        }

        result.frameKeep.assign(q->durations.size(), false);
        std::fill(
            result.frameKeep.begin() + static_cast<std::ptrdiff_t>(result.supportFrameStartIndex),
            result.frameKeep.begin() + static_cast<std::ptrdiff_t>(result.supportFrameCount),
            true);
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
                        "Select a valid PBIF time_s,template_value CSV.");
            }
            return false;
        }

        const double patientSourceSupportStart =
            this->frameStartForInputSec(
                std::min(result.supportFrameStartIndex, q->timePoints.size() - 1));
        const double patientSourceSupportEnd =
            this->frameEndForInputSec(
                std::min(
                    result.supportFrameCount > 0 ? result.supportFrameCount - 1 : size_t{0},
                    q->timePoints.size() - 1));

        // PBIF is never extrapolated beyond its supplied support. Determine
        // the tissue frames that are fully covered by the real template
        // samples. An artificial (0,0) anchor does not extend real support.
        const double pbifObservedStartSec =
            this->pbifZeroAnchorAdded && this->pbifTimesSec.size() > 1
            ? this->pbifTimesSec[1]
            : this->pbifTimesSec.front();
        const bool pbifCoversFromInjection =
            pbifObservedStartSec <= 60.0 + 1e-6;
        const double pbifCoverageStartSec =
            pbifCoversFromInjection ? 0.0 : pbifObservedStartSec;

        size_t pbifSupportStartIndex = q->timePoints.size();
        size_t pbifSupportFrameCount = 0;
        for (size_t i = 0; i < q->timePoints.size(); ++i)
        {
            const double frameStart = this->frameStartForInputSec(i);
            const double frameEnd = this->frameEndForInputSec(i);
            if (frameStart + 1e-6 >= pbifCoverageStartSec &&
                frameEnd <= this->pbifTimesSec.back() + 1e-6)
            {
                if (pbifSupportStartIndex == q->timePoints.size())
                {
                    pbifSupportStartIndex = i;
                }
                pbifSupportFrameCount = i + 1;
            }
        }

        if (pbifSupportStartIndex >= pbifSupportFrameCount ||
            pbifSupportFrameCount - pbifSupportStartIndex < 2)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The PBIF template does not cover at least two tissue frames.");
            }
            return false;
        }

        // Once calibrated, the PBIF template is the authoritative model
        // input. Its own temporal support - not the shorter patient
        // calibration curve - determines which tissue frames are usable.
        result.supportFrameStartIndex = pbifSupportStartIndex;
        result.supportFrameCount = pbifSupportFrameCount;
        if (result.supportFrameCount - result.supportFrameStartIndex < 2)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The common PBIF/tissue support contains fewer than two frames.");
            }
            return false;
        }

        modelEndTimeSec =
            this->frameEndForInputSec(result.supportFrameCount - 1);

        if (result.frameKeep.size() == q->durations.size())
        {
            std::fill(result.frameKeep.begin(), result.frameKeep.end(), false);
            std::fill(result.frameKeep.begin() + static_cast<std::ptrdiff_t>(result.supportFrameStartIndex),
                      result.frameKeep.begin() + static_cast<std::ptrdiff_t>(result.supportFrameCount),
                      true);
        }

        if (pbifSupportFrameCount < q->timePoints.size())
        {
            this->logToPythonConsole(
                QObject::tr(
                    "[SlicerDynamicPET PBIF] Template ends at %1 s; "
                    "the usable acquisition is shortened to frame %2 ending at %3 s. "
                    "PBIF extrapolation is not performed.")
                    .arg(this->pbifTimesSec.back(), 0, 'g', 8)
                    .arg(result.supportFrameCount)
                    .arg(modelEndTimeSec, 0, 'g', 8));
        }

        // Calibrate on the actual overlap of patient and template coverage.
        // If a previously stored GUI interval falls outside that overlap,
        // clamp it silently instead of presenting a warning on mode changes.
        double patientCoverageStart = patientSourceSupportStart;
        double patientCoverageEnd = patientSourceSupportEnd;
        if (source != 0)
        {
            const std::vector<double>& patientCoverageTimes =
                result.nativePlasmaTimesSec.empty()
                ? result.nativeWholeBloodTimesSec
                : result.nativePlasmaTimesSec;
            if (!patientCoverageTimes.empty())
            {
                patientCoverageStart = std::max(patientCoverageStart, patientCoverageTimes.front());
                patientCoverageEnd = std::min(patientCoverageEnd, patientCoverageTimes.back());
            }
        }

        const double overlapStart =
            std::max(patientCoverageStart, pbifObservedStartSec);
        const double overlapEnd =
            std::min(patientCoverageEnd, this->pbifTimesSec.back());

        if (!(overlapEnd > overlapStart + 1e-6))
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "Patient input and PBIF template have no usable common calibration interval.");
            }
            return false;
        }

        double calibrationStart =
            std::max(this->PBIFCalibrationStartSpinBox->value(), overlapStart);
        double calibrationEnd =
            std::min(this->PBIFCalibrationEndSpinBox->value(), overlapEnd);
        if (!(calibrationEnd > calibrationStart + 1e-6))
        {
            calibrationStart = overlapStart;
            calibrationEnd = overlapEnd;
        }

        if (std::abs(calibrationStart - this->PBIFCalibrationStartSpinBox->value()) > 1e-6 ||
            std::abs(calibrationEnd - this->PBIFCalibrationEndSpinBox->value()) > 1e-6)
        {
            QSignalBlocker startBlocker(this->PBIFCalibrationStartSpinBox);
            QSignalBlocker endBlocker(this->PBIFCalibrationEndSpinBox);
            this->PBIFCalibrationStartSpinBox->setValue(calibrationStart);
            this->PBIFCalibrationEndSpinBox->setValue(calibrationEnd);
            this->logToPythonConsole(
                QObject::tr(
                    "[SlicerDynamicPET PBIF] Calibration interval adjusted to the common "
                    "patient/template support: %1-%2 s.")
                    .arg(calibrationStart, 0, 'g', 8)
                    .arg(calibrationEnd, 0, 'g', 8));
        }

        const PBIFTemplateDomain templateChoice =
            this->selectedPBIFTemplateDomain();

        // The template domain is an intrinsic property of the supplied PBIF.
        // Always convert/use the patient calibration curve in that same domain.
        const IFCurveDomain calibrationDomain =
            templateChoice == PBIFTemplateDomain::TotalPlasma
            ? IFCurveDomain::TotalPlasma
            : IFCurveDomain::WholeBlood;

        double patientAUC =
            std::numeric_limits<double>::quiet_NaN();

        if (source == 0)
        {
            const std::vector<double>& patientFrames =
                calibrationDomain == IFCurveDomain::WholeBlood
                ? result.frameWholeBlood
                : result.framePlasma;

            const double frameTimeShift =
                this->frameTimeShiftForInputSec();
            patientAUC =
                this->integrateFrameAverageCurveOverInterval(
                    patientFrames,
                    calibrationStart - frameTimeShift,
                    calibrationEnd - frameTimeShift);

            result.pbifPatientCalibrationTimesSec.reserve(
                patientFrames.size());
            result.pbifPatientCalibrationValues =
                patientFrames;

            for (size_t i = 0;
                 i < patientFrames.size();
                 ++i)
            {
                result.pbifPatientCalibrationTimesSec.push_back(
                    this->frameEndForInputSec(i));
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

        if (pbifCoversFromInjection)
        {
            result.inputCoversFromInjection = true;
            result.inputCoverageReconstructedByPBIF = true;
            result.earliestAvailableInputTimeSec = pbifObservedStartSec;
        }
        else
        {
            result.inputCoversFromInjection = false;
            result.inputCoverageReconstructedByPBIF = false;
            result.earliestAvailableInputTimeSec = pbifObservedStartSec;
        }
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
        if (!fillFrameAverages(
                result.nativePlasmaTimesSec,
                result.nativePlasmaValues,
                "linear",
                result.framePlasma,
                result.supportFrameCount,
                result.supportFrameStartIndex) ||
            !fillFrameAverages(
                result.nativeWholeBloodTimesSec,
                result.nativeWholeBloodValues,
                "linear",
                result.frameWholeBlood,
                result.supportFrameCount,
                result.supportFrameStartIndex))
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
    //    Linear interpolation uses the supplied measurements directly; if that
    //    curve ends before the current IF support, the supported acquisition is
    //    shortened. Parametric models may extrapolate smoothly to the IF support.
    // ------------------------------------------------------------

    if (this->MetaboliteCorrectionCheckBox->isChecked())
    {
        if (result.plasmaIsParent)
        {
            if (errorMessage)
            {
                *errorMessage = QObject::tr(
                    "The selected external curve already represents parent plasma; metabolite correction must be off.");
            }
            return false;
        }

        std::vector<double> processedParentTimes;
        std::vector<double> processedParentValues;
        ParentFractionFitParameters parentFit;
        QString parentLabel;
        QString parentSummary;
        if (!this->buildProcessedParentFraction(
                modelEndTimeSec,
                processedParentTimes,
                processedParentValues,
                &parentFit,
                &parentLabel,
                &parentSummary,
                errorMessage))
        {
            return false;
        }

        result.parentFractionMeasuredEndTimeSec = this->parentFractionTimesSec.back();
        result.parentFractionModeledEndTimeSec = processedParentTimes.back();
        if (this->selectedParentFractionModel() != ParentFractionModel::Linear &&
            result.parentFractionModeledEndTimeSec >
                result.parentFractionMeasuredEndTimeSec + 1e-6)
        {
            result.parentFractionExtrapolationApplied = true;
            this->logToPythonConsole(
                QObject::tr("[SlicerDynamicPET parent fraction] Measured support ends at %1 s; %2 extrapolates the parent fraction to %3 s.")
                    .arg(result.parentFractionMeasuredEndTimeSec, 0, 'g', 8)
                    .arg(parentLabel)
                    .arg(result.parentFractionModeledEndTimeSec, 0, 'g', 8));
        }

        if (this->selectedParentFractionModel() == ParentFractionModel::Linear &&
            processedParentTimes.back() + 1e-6 < modelEndTimeSec)
        {
            size_t parentSupportFrameCount = result.supportFrameStartIndex;
            for (size_t i = result.supportFrameStartIndex; i < result.supportFrameCount; ++i)
            {
                if (this->frameEndForInputSec(i) <= processedParentTimes.back() + 1e-6)
                    parentSupportFrameCount = i + 1;
                else
                    break;
            }

            if (parentSupportFrameCount - result.supportFrameStartIndex < 2)
            {
                if (errorMessage)
                {
                    *errorMessage = QObject::tr(
                        "The parent-fraction curve does not cover at least two TAC frames.");
                }
                return false;
            }

            result.supportFrameCount = parentSupportFrameCount;
            modelEndTimeSec = this->frameEndForInputSec(result.supportFrameCount - 1);
            this->logToPythonConsole(
                QObject::tr("[SlicerDynamicPET parent fraction] Linear parent fraction ends at %1 s; fitting is shortened to TAC frame %2 ending at %3 s.")
                    .arg(processedParentTimes.back(), 0, 'g', 8)
                    .arg(result.supportFrameCount)
                    .arg(modelEndTimeSec, 0, 'g', 8));
        }

        result.applyParentFraction = true;
        result.parentFractionTimesSec = std::move(processedParentTimes);
        result.parentFractionValues = std::move(processedParentValues);

        result.frameModelPlasma.assign(
            q->durations.size(),
            std::numeric_limits<double>::quiet_NaN());

        const bool plasmaIsFrameCurve = result.nativePlasmaTimesSec.empty();
        const std::string plasmaInterpolation = result.pbifApplied ? "linear" : inputInterpolation;

        for (size_t i = result.supportFrameStartIndex; i < result.supportFrameCount; ++i)
        {
            const double frameEnd = this->frameEndForInputSec(i);
            const double frameStart = frameEnd - q->durations[i];
            const std::vector<double>& plasmaValues =
                plasmaIsFrameCurve ? result.framePlasma : result.nativePlasmaValues;

            const double value = this->averagePlasmaTimesParentFractionOverInterval(
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
                    *errorMessage = QObject::tr("Could not calculate parent-plasma frame averages.");
                return false;
            }
            result.frameModelPlasma[i] = value;
        }
    }
    else
    {
        result.frameModelPlasma =
            result.framePlasma;
    }

    if (result.supportFrameCount == 0)
    {
        result.supportFrameStartIndex = 0;
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

    std::vector<double> sigmas;
    double fallbackSigma = std::numeric_limits<double>::quiet_NaN();
    if (weighted)
    {
        sigmas.reserve(it->second.size());
        for (const VoxelStatistics& vs : it->second)
        {
            sigmas.push_back(statisticDispersionSigma(vs, statistic));
        }
        fallbackSigma = medianValidSigma(sigmas);
    }

    for (size_t i = 0; i < it->second.size(); ++i)
    {
        const VoxelStatistics& vs = it->second[i];
        double weight = 1.0;

        if (weighted)
        {
            weight = inverseVarianceWeightFromSigma(sigmas[i], fallbackSigma);
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
  if (currentSelectedID.empty() && !q->plotMTGAVOI.empty())
  {
    currentSelectedID = q->plotMTGAVOI;
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
  if (currentSelectedID.empty() && !q->plotTCMVOI.empty())
  {
    currentSelectedID = q->plotTCMVOI;
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
  QStringList headers = { "", "Ki / Ki'", "DV / DV'", "Intercept", "R2", "AIC", "MASE"};
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

  const int comparableModelCount =
      static_cast<int>(modelsForSegment.size());

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

  InputFunctionResult currentInput;
  QString ignoredInputError;
  const bool inputValid =
      this->buildCurrentInputFunction(
          currentInput,
          false,
          &ignoredInputError);

  // Delayed tissue acquisition remains a relative-MTGA comparison even when
  // PBIF reconstructs a complete input from injection: PBIF can restore the
  // missing INPUT history, but it cannot restore the unobserved early TISSUE
  // response. The same relative pair is also used for any other partial-input
  // configuration.
  const bool relativeMode =
      this->acquisitionTiming.delayedAcquisition ||
      (inputValid &&
       (!currentInput.inputCoversFromInjection ||
        currentInput.supportFrameStartIndex > 0));

  const std::string primaryModel =
      relativeMode ? "Relative Patlak" : "Patlak";

  const bool hasPrimary =
      hasResult(primaryModel);

  const bool hasLogan =
      !relativeMode && hasResult("Logan");

  const bool hasRE =
      !relativeMode && hasResult("RE");

  const bool hasRelativeRE =
      relativeMode && hasResult("Relative RE");

  const bool available =
      hasPrimary &&
      (relativeMode ? hasRelativeRE : (hasLogan || hasRE));

  if (relativeMode)
  {
    this->MTGAOptimizationInfoLabelImg->setText(
        QObject::tr(
            "Compare Relative Patlak with Relative RE voxel by voxel. "
            "For delayed/partial acquisitions the red channel contains relative Ki' "
            "and the blue channel contains relative DV_T'."));
    this->MTGAReversibleModelLabelImg->setText(
        QObject::tr("Relative reversible model"));
    this->MTGAChannelInfoLabelImg->setText(
        QObject::tr(
            "Ki' and DV_T' are relative quantitative maps for delayed/partial acquisitions. "
            "Their initial display window/level is determined automatically by Slicer. "
            "Adjust the display ranges in the Volumes module if desired, then refresh the RGB visualization."));
    this->GenerateMTGAOptimizedImgButton->setText(
        QObject::tr("Generate optimized Ki'/DV_T' + RGB"));
    this->RefreshMTGARGBButtonImg->setText(
        QObject::tr("Refresh RGB from Ki'/DV_T' display ranges"));
  }
  else
  {
    this->MTGAOptimizationInfoLabelImg->setText(
        QObject::tr(
            "Compare Patlak with a reversible model voxel by voxel. The selected "
            "Patlak Ki contributes to the red channel and the selected reversible DV to the blue channel."));
    this->MTGAReversibleModelLabelImg->setText(
        QObject::tr("Reversible model"));
    this->MTGAChannelInfoLabelImg->setText(
        QObject::tr(
            "Ki and DV are quantitative scalar maps. Their initial display window/level is determined automatically by Slicer. "
            "Adjust the display ranges in the Volumes module if desired, then refresh the RGB visualization to apply the current Ki and DV display ranges."));
    this->GenerateMTGAOptimizedImgButton->setText(
        QObject::tr("Generate optimized Ki/DV + RGB"));
    this->RefreshMTGARGBButtonImg->setText(
        QObject::tr("Refresh RGB from Ki/DV display ranges"));
  }

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

  if (relativeMode)
  {
    if (hasRelativeRE)
    {
      this->MTGAReversibleModelComboImg
          ->addItem(
              "Relative RE",
              "Relative RE");
    }
  }
  else
  {
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

    if (rgbName.startsWith("MTGA Optimized KiPrime"))
    {
      rgbName.replace(
          "MTGA Optimized KiPrime",
          "MTGA Relative Selection RGB");
    }
    else
    {
      rgbName.replace(
          "MTGA Optimized Ki",
          "MTGA Selection RGB");
    }

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
  // Required fitted models. Delayed/partial acquisition uses the relative
  // graphical pair; early complete acquisition keeps the standard pair.
  // ------------------------------------------------------------------------

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

  const bool relativeMode =
      reversibleModel == "Relative RE";

  const std::string primaryModel =
      relativeMode
      ? "Relative Patlak"
      : "Patlak";

  auto primaryIt =
      q->MTGAImgOutcomes.find(
          primaryModel);

  if (primaryIt ==
          q->MTGAImgOutcomes.end() ||
      primaryIt->second.empty())
  {
    return;
  }


  auto reversibleIt =
      q->MTGAImgOutcomes.find(
          reversibleModel);

  if (reversibleIt ==
          q->MTGAImgOutcomes.end() ||
      reversibleIt->second.empty())
  {
    return;
  }


  const auto& primary =
      primaryIt->second;

  const auto& reversible =
      reversibleIt->second;


  if (primary.size() !=
          reversible.size() ||
      primary.size() !=
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
      << primaryModel << " voxels=" << primary.size()
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

  const QString primaryQString =
      QString::fromStdString(primaryModel);

  const QString kiField =
      relativeMode ? "KiPrime" : "Ki";

  const QString dvField =
      relativeMode ? "DVPrime" : "DV";

  QString optimizationSuffix =
      primaryQString +
      " vs " +
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
      primary.size();

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


  size_t primarySelectedCount = 0;
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
        primary[v];

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

      ++primarySelectedCount;
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
      primarySelectedCount +
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
      << primaryModel << "=" << primarySelectedCount
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
      << "| primary =" << primaryQString
      << "=" << primarySelectedCount
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
            "MTGA Optimized " + kiField + " - " + optimizationSuffix,
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
            "MTGA Optimized " + dvField + " - " + optimizationSuffix,
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
              "SlicerDynamicPET.MTGA.PrimaryModel",
              primaryQString
                  .toUtf8()
                  .constData());

          node->SetAttribute(
              "SlicerDynamicPET.MTGA.ReversibleModel",
              reversibleQString
                  .toUtf8()
                  .constData());

          node->SetAttribute(
              "SlicerDynamicPET.MTGA.RelativeMode",
              relativeMode ? "1" : "0");

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
          "MTGA Optimized - " + kiField + "...");

      q->ProgressBar->setVisible(true);

      QApplication::processEvents();

      this->exportParametricMapDICOM(
          refPETNode,
          this->MTGAOptimizedKiValues,
          "MTGA",
          "MTGAOptimized",
          kiField.toStdString(),
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
              "MTGA Optimized - %1 was not exported "
              "because Framing Norm is %2 s.\n\n"
              "Its physical unit cannot be represented "
              "honestly as seconds or minutes without "
              "rescaling the numerical values.")
              .arg(kiField)
              .arg(framingNorm));
    }


    q->ProgressBar->setMinimum(0);
    q->ProgressBar->setMaximum(0);

    q->ProgressBar->setFormat(
        "Saving DICOM PMAP: "
        "MTGA Optimized - " + dvField + "...");

    q->ProgressBar->setVisible(true);

    QApplication::processEvents();

    this->exportParametricMapDICOM(
        refPETNode,
        this->MTGAOptimizedDVValues,
        "MTGA",
        "MTGAOptimized",
        dvField.toStdString(),
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

  if (modelID == "Patlak" || modelID == "Relative Patlak")
  {
    fields =
    {
      modelID == "Relative Patlak" ? "KiPrime" : "Ki",
      "Intercept",
      "AIC",
      "MASE",
      "R2",
      "chi2"
    };
  }
  else if (modelID == "Logan" ||
           modelID == "RE" ||
           modelID == "Relative RE")
  {
    fields =
    {
      modelID == "Relative RE" ? "DVPrime" : "DV",
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

    if (modelID == "Patlak")               modelIndex = 0;
    else if (modelID == "Relative Patlak") modelIndex = 1;
    else if (modelID == "Logan")           modelIndex = 2;
    else if (modelID == "RE")              modelIndex = 3;
    else if (modelID == "Relative RE")     modelIndex = 4;

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
      if ((modelID == "Patlak" || modelID == "Relative Patlak") &&
          (field == "Ki" || field == "KiPrime"))
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
           modelID == "RE" ||
           modelID == "Relative RE") &&
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
  else if (field == "KiPrime")
  {
    quantityCode = "SDP_KIP";
    quantityMeaning = "Relative Patlak slope Ki prime";
  }
  else if (field == "DV")
  {
    quantityCode = "SDP_DV";
    quantityMeaning = "Distribution volume";
  }
  else if (field == "DVPrime")
  {
    quantityCode = "SDP_DVP";
    quantityMeaning = "Relative equilibrium slope DV T prime";
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
  else if (modelID == "Relative Patlak")
    methodCode = "SDP_RPATLAK";
  else if (modelID == "Logan")
    methodCode = "SDP_LOGAN";
  else if (modelID == "RE")
    methodCode = "SDP_RE";
  else if (modelID == "Relative RE")
    methodCode = "SDP_RRE";
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

  QString exportBaseName;
  if (q->SubjectHierarchyNode
      && q->patID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    exportBaseName = QString::fromStdString(
        q->SubjectHierarchyNode->GetItemName(q->patID));
  }
  if (exportBaseName.trimmed().isEmpty())
  {
    exportBaseName = QStringLiteral("SlicerDynamicPET");
  }

  const QString outputPath =
      outputDir.filePath(
          QString(
              "%1_PMAP_%2_%3_%4.dcm")
              .arg(
                  exportBaseName,
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

    const BodySupportSource effectiveSource = requestedSource;

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
    q->SegWatcher->Clear();
    q->SegWatcher->browser = nullptr;
  }

  q->durations.clear();
  q->timePoints.clear();
  q->suvFactors.clear();

  q->numberOfTimepoints = 0;

  q->petID = newPetID;

  this->updateBodySupportUI();
  this->updateSegmentationAdvancedUI();
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
    "Patlak", "Relative Patlak", "Logan", "RE", "Relative RE"
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
  if (d->isTableBasedMode() || d->multiTimepointModeTransitionRunning)
  {
    // Table mode is independent of the MRML tree.  Mode transitions must also
    // be atomic because clearTACdata() removes plot nodes and can synchronously
    // emit Subject Hierarchy events.
    return;
  }

  if (d->isMultiTimepointMode())
  {
    // Multi preparation and plotting legitimately create/remove many MRML
    // nodes. Rebuilding the patient/acquisition/common-ROI UI for each of
    // those events caused lost checkbox state, slow plotting, invalidated
    // preparation caches, and re-entrant mode-switch crashes.
    //
    // The acquisition table already has explicit refresh points: entering
    // Multi and opening "Select/Edit acquisitions".  If no patient is
    // selected yet then refreshing the patient list is still useful (e.g.
    // after a DICOM import), but an established Multi analysis is left alone.
    if (this->patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    {
      d->populatePatientComboBox();
    }
    return;
  }

  // A segmentation (or another study child) can be announced by Subject
  // Hierarchy before Slicer has finished assigning/reparenting its hierarchy
  // item. Rebuilding synchronously can therefore miss a node that has just
  // been added while this module is open. Coalesce the event burst and refresh
  // once after Slicer has completed the current MRML/Subject Hierarchy update.
  if (d->subjectHierarchyRefreshQueued)
  {
    return;
  }
  d->subjectHierarchyRefreshQueued = true;

  QTimer::singleShot(0, this, [this]()
  {
    Q_D(qSlicerDynamicPETModuleWidget);
    d->subjectHierarchyRefreshQueued = false;

    if (!this->IsActive ||
        d->isTableBasedMode() ||
        d->isMultiTimepointMode() ||
        d->multiTimepointModeTransitionRunning)
    {
      return;
    }

    d->populatePatientComboBox();
  });
}

void qSlicerDynamicPETModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  this->Superclass::setMRMLScene(scene);
  Q_D(qSlicerDynamicPETModuleWidget);

  d->clearMultiTimepointSegmentationWatchers();
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
      if (d->isTableBasedMode() || d->isMultiTimepointMode())
      {
        return;
      }
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
    if (d->isTableBasedMode() || d->isMultiTimepointMode())
    {
      return;
    }
    // Any edited dynamic segmentation invalidates ROI fits.
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
  this->SegWatcher->OnSegmentContentChanged =
      [this, d](int, int frameIndex, const std::string& segmentID)
  {
    // Segment Editor undo/redo may still be restoring its segmentation history
    // while vtkSegmentation emits the content event. Never rasterize/recompute
    // TAC from that VTK callback. Let Slicer finish the restoration first, then
    // update the exact frame that generated the event on the Qt event loop.
    const std::string deferredSegmentID = segmentID;
    QTimer::singleShot(0, this, [this, d, frameIndex, deferredSegmentID]()
    {
      if (d->isTableBasedMode() || d->isMultiTimepointMode() ||
          !this->SegWatcher || !this->SubjectHierarchyNode ||
          this->segID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
      {
        return;
      }

      vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(
          this->SubjectHierarchyNode->GetItemDataNode(this->segID));
      if (!segNode)
      {
        return;
      }

      this->SegWatcher->ApplyDeferredLegacyContentChange(
          segNode, deferredSegmentID, frameIndex);
    });
  };

}




void qSlicerDynamicPETModuleWidget::onPatChanged (int index) {
  Q_D(qSlicerDynamicPETModuleWidget);
  d->resetAcquisitionTimingDisplay();
  this->resetPETSelection();
  this->patID = d->PatSelector->itemData(index).value<vtkIdType>();
  if (d->isMultiTimepointMode())
  {
    this->stuID = vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID;
    d->populateMultiTimepointAcquisitionTable();
    if (this->patID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    {
      QTimer::singleShot(0, this, [d]()
      {
        d->showMultiTimepointSelectionDialog();
      });
    }
  }
  else
  {
    d->populateStudyComboBox(this->patID);
  }
  if (this->patID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    // No patient name is available: keep neutral, usable fallback names for
    // every user-editable export field instead of resetting TAC only.
    d->fileexcel->setText(QStringLiteral("TAC.xlsx"));
    d->fileexceltcm->setText(QStringLiteral("TCMparameters.xlsx"));
    d->fileexceltcmfitted->setText(QStringLiteral("TCMfitted.xlsx"));
    d->fileexcelmtga->setText(QStringLiteral("MTGAparameters.xlsx"));
    d->fileexcelmtgafitted->setText(QStringLiteral("MTGAfitted.xlsx"));
    d->DynamicRTStructFilename->setText(QStringLiteral("RTSTRUCT.dcm"));
    return;
  }

  const QString name = QString::fromStdString(
      this->SubjectHierarchyNode->GetItemName(this->patID));

  // Keep all editable export filenames on the same convention:
  // <patient name>_<export type>.<extension>
  d->fileexcel->setText(name + QStringLiteral("_TAC.xlsx"));
  d->fileexceltcm->setText(name + QStringLiteral("_TCMparameters.xlsx"));
  d->fileexceltcmfitted->setText(name + QStringLiteral("_TCMfitted.xlsx"));
  d->fileexcelmtga->setText(name + QStringLiteral("_MTGAparameters.xlsx"));
  d->fileexcelmtgafitted->setText(name + QStringLiteral("_MTGAfitted.xlsx"));
  d->DynamicRTStructFilename->setText(name + QStringLiteral("_RTSTRUCT.dcm"));
}


void qSlicerDynamicPETModuleWidget::onStuChanged (int index) {
  Q_D(qSlicerDynamicPETModuleWidget);
  d->resetAcquisitionTimingDisplay();
  this->stuID = d->StuSelector->itemData(index).value<vtkIdType>();
  this->resetPETSelection();
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

    // PET TACs are independent of the selected CT. CT remains available
    // only for the optional parametric-imaging body-support mask.

    this->enableTACbutton();
}

void qSlicerDynamicPETModuleWidget::resetPETSelection()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->resetInputFunctionData(true);
  this->sequencePETNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  if (this->SegWatcher)
  {
    this->SegWatcher->Clear();
    this->SegWatcher->browser = nullptr;
  }

  this->durations.clear();
  this->timePoints.clear();
  this->suvFactors.clear();
  this->numberOfTimepoints = 0;
  d->suvbwFactorValidated = false;
  d->resetAcquisitionTimingDisplay();

  // A PET change is a new temporal problem. Do not carry slider values or
  // labels from the previous acquisition while no PET is selected.
  for (QSlider* slider : {d->timeOffsetSlider, d->timeEndSlider, d->TCMEndSlider,
                          d->timeOffsetSliderImg, d->timeEndSliderImg, d->TCMEndSliderImg})
  {
    if (!slider) continue;
    QSignalBlocker blocker(slider);
    slider->setRange(1, 1);
    slider->setValue(1);
  }
  for (QLineEdit* edit : {d->frameEdit, d->timeSecEdit, d->timeMinEdit,
                          d->timeEndInfoEdit, d->TCMEndInfoEdit,
                          d->frameEditImg, d->timeSecEditImg, d->timeMinEditImg,
                          d->timeEndInfoEditImg, d->TCMEndInfoEditImg})
  {
    if (edit) edit->clear();
  }

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
  d->resetAcquisitionTimingDisplay();
  d->resetInputFunctionData(true);
  const vtkIdType newPetID =
      d->PETSelector
          ->itemData(index)
          .value<vtkIdType>();

  d->setPETItemID(newPetID);

  this->sequencePETNode = nullptr;
  this->sequenceBrowserPETNode = nullptr;
  if (this->SegWatcher)
  {
    this->SegWatcher->Clear();
    this->SegWatcher->browser = nullptr;
  }

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

    // Resolve the selected PET from the actual master-sequence/proxy
    // relationship. This is independent of the Sequence Browser Rename
    // option and of the current proxy node name.
    vtkMRMLSequenceNode* seqNode = browser->GetMasterSequenceNode();
    if (!seqNode)
      continue;

    vtkMRMLNode* proxyNode = browser->GetProxyNode(seqNode);
    if (proxyNode != petNode)
      continue;

    foundSeqNode = seqNode;
    foundBrowser = browser;

    // If proxy was not valid, check sequence provenance.
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
  d->updateAcquisitionTimingContext(true);
  d->populateTimeBarMTGA(true);
  d->populateTimeBarMTGAImg(true);

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

  // A new Single segmentation context must start with exactly one observer
  // set. Clear any previous segmentation/proxy observations before rebuilding
  // the selected segmentation sequence below.
  if (this->SegWatcher)
  {
    this->SegWatcher->Clear();
    this->SegWatcher->browser = this->sequenceBrowserPETNode;
  }

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
  d->updateSegmentationAdvancedUI();

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

  // A dynamic RTSTRUCT may have been imported before its PET.  In that case
  // the dRTImporter deliberately creates a standalone browser.  Adopt that
  // existing temporal sequence into the selected PET browser instead of
  // treating the proxy as a static segmentation and cloning one frame.
  const char* importedDynamicAttribute =
      segNode->GetAttribute("dRTImporter.DynamicRTStruct");
  const bool importedDynamicSegmentation =
      importedDynamicAttribute
      && QString::fromUtf8(importedDynamicAttribute) == QStringLiteral("1");

  if (!seqNode && importedDynamicSegmentation)
  {
    PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
    const QVariant resultVariant = mainContext.call(
        "DPE_adopt_dynamic_rtstruct",
        QVariantList{
            QString::fromUtf8(segNode->GetID()),
            QString::fromUtf8(this->sequencePETNode->GetID()),
            QString::fromUtf8(this->sequenceBrowserPETNode->GetID())});
    const QVariantMap result = resultVariant.toMap();

    if (!result.value("ok").toBool())
    {
      this->ProgressBar->hide();
      this->ProgressBar->setValue(0);
      this->ProgressBar->setFormat("%p%");
      d->populateSegmentCheckboxes(this->segID);
      d->updateSegmentationAdvancedUI();
      QMessageBox::warning(
          this,
          tr("Dynamic RTSTRUCT"),
          tr("The imported dynamic segmentation could not be matched to the selected PET.\n\n%1")
              .arg(result.value("error").toString()));
      return;
    }

    const QByteArray sequenceNodeID =
        result.value("sequence_node_id").toString().toUtf8();
    seqNode = vtkMRMLSequenceNode::SafeDownCast(
        scene->GetNodeByID(sequenceNodeID.constData()));
    if (!seqNode)
    {
      this->ProgressBar->hide();
      this->ProgressBar->setValue(0);
      this->ProgressBar->setFormat("%p%");
      d->updateSegmentationAdvancedUI();
      QMessageBox::warning(
          this,
          tr("Dynamic RTSTRUCT"),
          tr("The imported dynamic segmentation was matched to the PET, but its aligned sequence could not be recovered."));
      return;
    }
  }

  if (!seqNode)
  {
    const std::string stableSegmentationName =
        shNode->GetItemName(segID);

    // Preserve the original, known-good Sequence Browser setup order.
    // Attach the existing segmentation proxy before populating the temporal
    // sequence.  This lets Sequence Browser establish its proxy bookkeeping
    // without forcing index/name updates on the already-loaded PET proxy.
    vtkSmartPointer<vtkMRMLSequenceNode> newSeqNode =
      vtkSmartPointer<vtkMRMLSequenceNode>::New();

    newSeqNode->SetName(
      stableSegmentationName.c_str());

    scene->AddNode(newSeqNode);

    this->sequenceBrowserPETNode->AddProxyNode(
      segNode,
      newSeqNode,
      false);

    this->sequenceBrowserPETNode->SetSaveChanges(
      newSeqNode,
      true);

    std::string indexValue;
    for (int i = 0; i < this->numberOfTimepoints; ++i)
    {
      indexValue = this->sequencePETNode->GetNthIndexValue(i);

      if (!newSeqNode->GetDataNodeAtValue(indexValue))
      {
        newSeqNode->SetDataNodeAtValue(segNode, indexValue);
      }

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
      qApp->processEvents();
    }

    // The synchronized segmentation proxy is now fully established and the
    // temporal sequence has been populated. Disable automatic proxy renaming
    // only now, so browsing changes content but not the user-visible name.
    this->sequenceBrowserPETNode->SetOverwriteProxyName(
      newSeqNode,
      false);
    segNode->SetName(stableSegmentationName.c_str());
    shNode->SetItemName(segID, stableSegmentationName);

    this->SegWatcher->ObserveSegmentationNode(segNode);
    seqNode = newSeqNode;
  }

  this->segSequenceNode = seqNode;
  this->SegWatcher->ObserveSegmentationNode(segNode);

  d->populateSegmentCheckboxes(this->segID);
  this->enableTACbutton();
  d->updateSegmentationAdvancedUI();
  d->updateSegmentationFrameUI(true);

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
  if (d->isMultiTimepointMode())
  {
    d->syncMultiTimepointSelectedSegments();
    return;
  }
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

void qSlicerDynamicPETModuleWidget::onOpenSegmentEditor()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  vtkMRMLScene* scene = this->mrmlScene();
  vtkMRMLSegmentationNode* segNode = nullptr;
  vtkMRMLScalarVolumeNode* sourceVolume = nullptr;

  if (d->isMultiTimepointMode())
  {
    if (!scene || !d->multiTimepointPreparationValid ||
        d->preparedMultiTimepointObservations.empty())
    {
      QMessageBox::warning(
          this, tr("Segment Editor"),
          tr("Prepare the Multi-timepoint acquisitions first."));
      return;
    }

    const int observationIndex = d->segmentationFrameSlider
        ? d->segmentationFrameSlider->value() - 1 : 0;
    if (observationIndex < 0 ||
        observationIndex >= static_cast<int>(d->preparedMultiTimepointObservations.size()))
    {
      QMessageBox::warning(this, tr("Segment Editor"), tr("No displayed observation is available."));
      return;
    }

    // The slider only selects what is displayed/opened. Once Segment Editor is
    // active, acquisition-specific watchers detect edits independently of this
    // slider and resolve the actual changed acquisition/frame.
    const PreparedMultiTimepointObservation& observation =
        d->preparedMultiTimepointObservations[static_cast<size_t>(observationIndex)];
    const PreparedMultiTimepointAcquisition& acquisition =
        d->preparedMultiTimepointAcquisitions[
            static_cast<size_t>(observation.acquisitionIndex)];

    if (observation.dynamic)
    {
      vtkMRMLSequenceBrowserNode* browser = vtkMRMLSequenceBrowserNode::SafeDownCast(
          scene->GetNodeByID(acquisition.petBrowserNodeID.toUtf8().constData()));
      vtkMRMLSequenceNode* petSequence = vtkMRMLSequenceNode::SafeDownCast(
          scene->GetNodeByID(observation.petSequenceNodeID.toUtf8().constData()));
      segNode = vtkMRMLSegmentationNode::SafeDownCast(
          scene->GetNodeByID(acquisition.segmentationNodeID.toUtf8().constData()));
      if (browser && observation.frameIndex >= 0)
      {
        browser->SetSelectedItemNumber(observation.frameIndex);
      }
      if (browser && petSequence)
      {
        sourceVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
            browser->GetProxyNode(petSequence));
      }
    }
    else
    {
      segNode = vtkMRMLSegmentationNode::SafeDownCast(
          scene->GetNodeByID(observation.segmentationNodeID.toUtf8().constData()));
      sourceVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
          scene->GetNodeByID(observation.petNodeID.toUtf8().constData()));
    }
  }
  else
  {
    vtkMRMLSubjectHierarchyNode* shNode = scene
        ? vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene) : nullptr;
    segNode = (shNode && this->segID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        ? vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(this->segID))
        : nullptr;
    if (this->sequenceBrowserPETNode && this->sequencePETNode)
    {
      sourceVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
          this->sequenceBrowserPETNode->GetProxyNode(this->sequencePETNode));
    }
  }

  if (!segNode)
  {
    QMessageBox::warning(
        this, tr("Segment Editor"), tr("Select a segmentation first."));
    d->updateSegmentationAdvancedUI();
    return;
  }

  // Segment Editor can modify only Binary labelmap source representations.
  // Multi preparation intentionally keeps Planar contour/other source data
  // untouched for read-only TAC extraction. Convert only when the user
  // explicitly opens Segment Editor, so editing never triggers Slicer's
  // confirmation dialog while non-editing workflows preserve the original
  // source representation for as long as possible.
  vtkSegmentation* segmentation = segNode->GetSegmentation();
  const std::string binaryRep =
      vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();
  if (segmentation && segmentation->GetSourceRepresentationName() != binaryRep)
  {
    if (sourceVolume)
    {
      segNode->SetReferenceImageGeometryParameterFromVolumeNode(sourceVolume);
    }
    if (!segNode->SetSourceRepresentationToBinaryLabelmap())
    {
      QMessageBox::warning(
          this,
          tr("Segment Editor"),
          tr("Could not convert '%1' to an editable Binary labelmap source representation.")
              .arg(segNode->GetName() ? QString::fromUtf8(segNode->GetName()) : tr("Segmentation")));
      return;
    }
    d->logToPythonConsole(
        tr("[SlicerDynamicPET Segment tracking] Converted '%1' to Binary labelmap source representation for editing.")
            .arg(segNode->GetName() ? QString::fromUtf8(segNode->GetName()) : tr("Segmentation")));
  }

  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  const QVariant resultVariant = mainContext.call(
      "DPE_open_segment_editor",
      QVariantList{
          QString::fromUtf8(segNode->GetID()),
          sourceVolume ? QString::fromUtf8(sourceVolume->GetID()) : QString()});
  const QVariantMap result = resultVariant.toMap();
  if (!result.value("ok").toBool())
  {
    QMessageBox::warning(
        this,
        tr("Segment Editor"),
        tr("Could not open Segment Editor.\n\n%1")
            .arg(result.value("error").toString()));
  }
}

void qSlicerDynamicPETModuleWidget::onSaveDynamicRTStruct()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  QString directory = d->DynamicRTStructDirectory->currentPath().trimmed();
  QString fileName = d->DynamicRTStructFilename->text().trimmed();
  if (directory.isEmpty() || fileName.isEmpty())
  {
    QMessageBox::warning(
        this,
        tr("Save RTSTRUCT"),
        tr("Choose an output directory and filename."));
    return;
  }

  QDir outputDirectory(directory);
  if (!outputDirectory.exists() && !QDir().mkpath(directory))
  {
    QMessageBox::warning(
        this,
        tr("Save RTSTRUCT"),
        tr("The output directory could not be created:\n%1").arg(directory));
    return;
  }

  // Multi export deliberately follows prepared acquisition provenance, not the
  // merged observation timeline. Therefore N prepared acquisitions produce N
  // RTSTRUCT files. Dynamic acquisitions retain the temporal dRT convention;
  // static acquisitions are exported as ordinary non-temporal RTSTRUCTs.
  if (d->isMultiTimepointMode())
  {
    if (!d->multiTimepointPreparationValid
        || d->preparedMultiTimepointAcquisitions.empty())
    {
      QMessageBox::warning(
          this,
          tr("Save Multi RTSTRUCTs"),
          tr("Prepare the selected Multi-timepoint acquisitions before exporting segmentations."));
      d->updateSegmentationAdvancedUI();
      return;
    }

    QString baseName = fileName;
    if (baseName.endsWith(QStringLiteral(".dcm"), Qt::CaseInsensitive))
    {
      baseName.chop(4);
    }
    baseName = baseName.trimmed();
    if (baseName.isEmpty())
    {
      baseName = QStringLiteral("RTSTRUCT");
    }

    struct PendingRTStructExport
    {
      int acquisitionIndex{-1};
      bool dynamic{false};
      QString acquisitionName;
      QString segmentationNodeID;
      QString segmentationSequenceNodeID;
      QString petNodeID;
      QString petSequenceNodeID;
      QString outputPath;
    };

    std::vector<PendingRTStructExport> exports;
    exports.reserve(d->preparedMultiTimepointAcquisitions.size());
    QStringList existingPaths;

    for (size_t acquisitionIndex = 0;
         acquisitionIndex < d->preparedMultiTimepointAcquisitions.size();
         ++acquisitionIndex)
    {
      const PreparedMultiTimepointAcquisition& acquisition =
          d->preparedMultiTimepointAcquisitions[acquisitionIndex];

      PendingRTStructExport pending;
      pending.acquisitionIndex = static_cast<int>(acquisitionIndex);
      pending.dynamic = acquisition.dynamic;
      pending.acquisitionName = acquisition.petName;
      pending.segmentationNodeID = acquisition.segmentationNodeID;
      pending.segmentationSequenceNodeID = acquisition.segmentationSequenceNodeID;
      pending.petNodeID = acquisition.petNodeID;
      pending.petSequenceNodeID = acquisition.petSequenceNodeID;

      const QString acquisitionTag = QStringLiteral("Acq%1_%2")
          .arg(static_cast<int>(acquisitionIndex) + 1, 2, 10, QLatin1Char('0'))
          .arg(acquisition.dynamic
              ? QStringLiteral("Dynamic")
              : QStringLiteral("Static"));
      pending.outputPath = outputDirectory.filePath(
          baseName + QStringLiteral("_") + acquisitionTag + QStringLiteral(".dcm"));

      vtkMRMLScene* scene = this->mrmlScene();
      vtkMRMLNode* petNode = scene
          ? scene->GetNodeByID(pending.petNodeID.toUtf8().constData())
          : nullptr;
      vtkMRMLNode* segmentationNode = scene
          ? scene->GetNodeByID(pending.segmentationNodeID.toUtf8().constData())
          : nullptr;
      if (!petNode || !segmentationNode)
      {
        QMessageBox::warning(
            this,
            tr("Save Multi RTSTRUCTs"),
            tr("Prepared acquisition %1 ('%2') lost its PET or segmentation node. Re-prepare Multi-timepoint acquisitions before exporting.")
                .arg(static_cast<int>(acquisitionIndex) + 1)
                .arg(acquisition.petName));
        return;
      }
      if (pending.dynamic)
      {
        vtkMRMLNode* petSequence = scene->GetNodeByID(
            pending.petSequenceNodeID.toUtf8().constData());
        vtkMRMLNode* segmentationSequence = scene->GetNodeByID(
            pending.segmentationSequenceNodeID.toUtf8().constData());
        if (!petSequence || !segmentationSequence)
        {
          QMessageBox::warning(
              this,
              tr("Save Multi RTSTRUCTs"),
              tr("Prepared dynamic acquisition %1 ('%2') lost its PET or segmentation sequence. Re-prepare Multi-timepoint acquisitions before exporting.")
                  .arg(static_cast<int>(acquisitionIndex) + 1)
                  .arg(acquisition.petName));
          return;
        }
      }

      if (QFileInfo::exists(pending.outputPath))
      {
        existingPaths << pending.outputPath;
      }
      exports.push_back(std::move(pending));
    }

    if (!existingPaths.isEmpty())
    {
      QString existingText = existingPaths.mid(0, 6).join(QStringLiteral("\n"));
      if (existingPaths.size() > 6)
      {
        existingText += tr("\n... and %1 more file(s)")
            .arg(existingPaths.size() - 6);
      }
      const QMessageBox::StandardButton answer = QMessageBox::question(
          this,
          tr("Save Multi RTSTRUCTs"),
          tr("The following batch output file(s) already exist:\n%1\n\nReplace the existing files and continue?")
              .arg(existingText),
          QMessageBox::Yes | QMessageBox::No,
          QMessageBox::No);
      if (answer != QMessageBox::Yes)
      {
        return;
      }
    }

    PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
    QStringList savedPaths;
    for (const PendingRTStructExport& pending : exports)
    {
      QVariant resultVariant;
      if (pending.dynamic)
      {
        resultVariant = mainContext.call(
            "DPE_export_dynamic_rtstruct",
            QVariantList{
                pending.segmentationSequenceNodeID,
                pending.petSequenceNodeID,
                pending.petNodeID,
                pending.outputPath,
                true});
      }
      else
      {
        resultVariant = mainContext.call(
            "DPE_export_static_rtstruct",
            QVariantList{
                pending.segmentationNodeID,
                pending.petNodeID,
                pending.outputPath,
                true});
      }

      const QVariantMap result = resultVariant.toMap();
      if (!result.value("ok").toBool())
      {
        QString partialText;
        if (!savedPaths.isEmpty())
        {
          partialText = tr("\n\nAlready saved before the failure:\n%1")
              .arg(savedPaths.join(QStringLiteral("\n")));
        }
        QMessageBox::warning(
            this,
            tr("Save Multi RTSTRUCTs"),
            tr("RTSTRUCT export failed for acquisition %1 ('%2').\n\n%3%4")
                .arg(pending.acquisitionIndex + 1)
                .arg(pending.acquisitionName)
                .arg(result.value("error").toString())
                .arg(partialText));
        return;
      }
      savedPaths << result.value("path").toString();
    }

    d->propagateOutputDirectory(directory);
    QMessageBox::information(
        this,
        tr("Save Multi RTSTRUCTs"),
        tr("Saved %1 acquisition-specific RTSTRUCT file(s):\n\n%2")
            .arg(savedPaths.size())
            .arg(savedPaths.join(QStringLiteral("\n"))));
    return;
  }

  if (!this->segSequenceNode
      || !this->sequencePETNode
      || !this->sequenceBrowserPETNode)
  {
    QMessageBox::warning(
        this,
        tr("Save Dynamic RTSTRUCT"),
        tr("Select a dynamic PET and segmentation first."));
    d->updateSegmentationAdvancedUI();
    return;
  }

  if (!fileName.endsWith(QStringLiteral(".dcm"), Qt::CaseInsensitive))
  {
    fileName += QStringLiteral(".dcm");
    d->DynamicRTStructFilename->setText(fileName);
  }

  const QString outputPath = outputDirectory.filePath(fileName);
  if (QFileInfo::exists(outputPath))
  {
    const QMessageBox::StandardButton answer = QMessageBox::question(
        this,
        tr("Save Dynamic RTSTRUCT"),
        tr("The file already exists:\n%1\n\nReplace it?").arg(outputPath),
        QMessageBox::Yes | QMessageBox::No,
        QMessageBox::No);
    if (answer != QMessageBox::Yes)
    {
      return;
    }
  }

  vtkMRMLScene* scene = this->mrmlScene();
  vtkMRMLSubjectHierarchyNode* shNode =
      scene ? vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene) : nullptr;
  vtkMRMLScalarVolumeNode* referenceVolume = nullptr;
  if (shNode && this->ctID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    referenceVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
        shNode->GetItemDataNode(this->ctID));
  }
  if (!referenceVolume)
  {
    referenceVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
        this->sequenceBrowserPETNode->GetProxyNode(this->sequencePETNode));
  }
  if (!referenceVolume)
  {
    QMessageBox::warning(
        this,
        tr("Save Dynamic RTSTRUCT"),
        tr("No CT or current PET volume is available as export geometry reference."));
    return;
  }

  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  const QVariant resultVariant = mainContext.call(
      "DPE_export_dynamic_rtstruct",
      QVariantList{
          QString::fromUtf8(this->segSequenceNode->GetID()),
          QString::fromUtf8(this->sequencePETNode->GetID()),
          QString::fromUtf8(referenceVolume->GetID()),
          outputPath,
          true});
  const QVariantMap result = resultVariant.toMap();
  if (!result.value("ok").toBool())
  {
    QMessageBox::warning(
        this,
        tr("Save Dynamic RTSTRUCT"),
        tr("Dynamic RTSTRUCT export failed.\n\n%1")
            .arg(result.value("error").toString()));
    return;
  }

  d->propagateOutputDirectory(QFileInfo(outputPath).absolutePath());
  QMessageBox::information(
      this,
      tr("Save Dynamic RTSTRUCT"),
      tr("Dynamic RTSTRUCT saved successfully:\n%1")
          .arg(result.value("path").toString()));
}

void qSlicerDynamicPETModuleWidget::clearTACdata()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  // Image hierarchy/segmentation changes must never erase a table dataset.
  // Table data are cleared explicitly from the Table Setup controls.
  if (d->isTableBasedMode())
  {
    return;
  }

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
  d->rebuildTACStatisticUI();
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
  if (d->isTableBasedMode())
  {
    d->TACbutton->setEnabled(false);
    return;
  }
  if (d->isMultiTimepointMode())
  {
    d->TACbutton->setEnabled(
        d->multiTimepointSelectionValidated &&
        d->multiTimepointPreparationValid &&
        !this->segmentIDs.empty());
    return;
  }
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

  // Error bars use auxiliary table nodes because vtkMRMLPlotSeriesNode does
  // not expose a native error-column property. Remove them with the chart.
  vtkCollection* tableNodes =
      this->mrmlScene()->GetNodesByClass("vtkMRMLTableNode");
  if (tableNodes)
  {
    std::vector<vtkMRMLNode*> removeNodes;
    for (int i = 0; i < tableNodes->GetNumberOfItems(); ++i)
    {
      vtkMRMLNode* node = vtkMRMLNode::SafeDownCast(
          tableNodes->GetItemAsObject(i));
      if (node && node->GetName() &&
          QString::fromUtf8(node->GetName()).startsWith("DynamicPET.ErrorTable."))
      {
        removeNodes.push_back(node);
      }
    }
    tableNodes->Delete();
    for (vtkMRMLNode* node : removeNodes)
    {
      this->mrmlScene()->RemoveNode(node);
    }
  }
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

  if (d->isMultiTimepointMode())
  {
    QString error;
    if (!d->computeMultiTimepointTAC(&error))
    {
      if (!error.isEmpty() && error != tr("Multi-timepoint TAC extraction was cancelled."))
      {
        QMessageBox::warning(this, tr("Multi-timepoint TAC"), error);
        d->logToPythonConsole(
            tr("[SlicerDynamicPET multi-timepoint] TAC extraction failed: %1").arg(error));
      }
      return;
    }

    d->populatePlotSegmentCheckboxes();
    d->rebuildTACStatisticUI();
    d->populateIF();
    d->populateTimeBarMTGA(true);
    d->populateTimeBarMTGAImg(true);
    d->updateInputFunctionStatus();

    const bool tacReady =
        !this->segmentTACs.empty() &&
        !this->segmentTACsnames.empty();
    d->setPostTACEnabled(tacReady);
    d->updateParametricImagingAvailability();
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
    bool needsComputation =
        this->segmentTACsnames.find(segmentID) == this->segmentTACsnames.end();

    if (needsComputation)
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
  using SinglePerfClock = std::chrono::steady_clock;
  double singleFlattenMs = 0.0;
  bool singleFlattenExecuted = false;
  if (this->PET_flatten_values.empty()) {
    const auto flattenStart = SinglePerfClock::now();
    logic->Image2Flatten(this->petID, this->PET_flatten_values, this->PETdims, this->numberOfTimepoints, this->ProgressBar, this->stopButton, this->stopRequested);
    singleFlattenMs = std::chrono::duration<double, std::milli>(
        SinglePerfClock::now() - flattenStart).count();
    singleFlattenExecuted = true;
  }
  if (this->stopRequested) {
    return;
  }
  if (!segmentsToCompute.empty())
  {
    // Any new/updated segmentation TAC also invalidates lazily cached CT
    // volumes for that dynamic segmentation. Keep normal TAC computation PET-only.
  }

  double singleTacCallMs = 0.0;
  if (!segmentsToCompute.empty())
  {
    const auto tacStart = SinglePerfClock::now();
    logic->TAC(this->sequencePETNode, this->segSequenceNode, segmentsToCompute, this->segmentTACs, this->segmentTACsnames, this->ProgressBar, this->stopButton, this->stopRequested);
    singleTacCallMs = std::chrono::duration<double, std::milli>(
        SinglePerfClock::now() - tacStart).count();
  }

  if (this->sequencePETNode)
  {
    const char* logicSummary = this->sequencePETNode->GetAttribute(
        "SlicerDynamicPET.TACPerfSummary");
    const char* frameSummary = this->sequencePETNode->GetAttribute(
        "SlicerDynamicPET.TACPerfFrames");

    d->logToPythonConsole(tr(
        "[SlicerDynamicPET PERF][SINGLE] Widget timing: Image2Flatten=%1 ms (%2); logic->TAC wall=%3 ms; requestedROIs=%4.")
        .arg(singleFlattenMs, 0, 'f', 3)
        .arg(singleFlattenExecuted ? tr("executed") : tr("cached/skipped"))
        .arg(singleTacCallMs, 0, 'f', 3)
        .arg(static_cast<int>(segmentsToCompute.size())));

    if (logicSummary && *logicSummary)
    {
      d->logToPythonConsole(QString::fromUtf8(logicSummary));
    }
    if (frameSummary && *frameSummary)
    {
      d->logToPythonConsole(QString::fromUtf8(frameSummary));
    }
    this->sequencePETNode->RemoveAttribute("SlicerDynamicPET.TACPerfSummary");
    this->sequencePETNode->RemoveAttribute("SlicerDynamicPET.TACPerfFrames");
  }
  else
  {
    d->logToPythonConsole(tr(
        "[SlicerDynamicPET PERF][SINGLE] Widget timing: Image2Flatten=%1 ms (%2); logic->TAC wall=%3 ms; requestedROIs=%4; no PET sequence node available for stage summary.")
        .arg(singleFlattenMs, 0, 'f', 3)
        .arg(singleFlattenExecuted ? tr("executed") : tr("cached/skipped"))
        .arg(singleTacCallMs, 0, 'f', 3)
        .arg(static_cast<int>(segmentsToCompute.size())));
  }

  this->ProgressBar->setValue(0);
  this->ProgressBar->setVisible(false);
  if (this->stopRequested) {
    return;
  }
  // CT-referenced ROI volumes are intentionally not computed here.

  d->updateAcquisitionTimingContext(false);
  d->populatePlotSegmentCheckboxes();
  d->rebuildTACStatisticUI();
  d->populateIF();
  d->populateTimeBarMTGA(true);
  d->populateTimeBarMTGAImg(true);

  // Rebuild timing-dependent IF support after TAC extraction so slider ranges
  // immediately reflect PBIF/external-IF coverage. No checkbox retoggle should
  // ever be required to refresh the usable frame range.
  d->updateInputFunctionStatus();

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
  if (d->isMultiTimepointMode())
  {
    d->syncMultiTimepointSelectedSegments();
  }
  else
  {
    d->populateSegmentCheckboxes(this->segID);
  }
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
        const double frameEndSec = d->frameEndForInputSec(i);
        const double frameStartSec = d->frameStartForInputSec(i);
        const double frameMidSec = d->frameMidForInputSec(i);

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
      if (modelName == "Relative Patlak")
      {
        row["KiPrime"] = params.Ki;
        row["Ki"] = QVariant();
      }
      else
      {
        row["Ki"] = params.Ki;
        row["KiPrime"] = QVariant();
      }
      if (modelName == "Relative RE")
      {
        row["DVPrime"] = params.DV;
        row["DV"] = QVariant();
      }
      else
      {
        row["DV"] = params.DV;
        row["DVPrime"] = QVariant();
      }
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
      const double frameEndSec = d->frameEndForInputSec(i);
      const double frameStartSec = d->frameStartForInputSec(i);
      const double frameMidSec = d->frameMidForInputSec(i);
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

  QVariantMap metadata;
  metadata["FormatVersion"] = 1;
  metadata["ActivityUnit"] = d->activityUnitLabel(d->petStoredActivityUnit());
  d->updateAcquisitionTimingContext(false);
  metadata["TimeConvention"] =
      d->acquisitionTiming.delayedAcquisition
      ? "FrameEndPostInjection"
      : "FrameEnd";
  metadata["SourceMode"] = d->isTableBasedMode() ? "TableBased" : "ImageBased";

  double suvbwFactor = 0.0;
  if (d->getCommonSUVbwFactor(suvbwFactor, nullptr))
  {
    metadata["SUVbwFactor"] = suvbwFactor;
  }

  if (!d->isTableBasedMode() && this->SubjectHierarchyNode)
  {
    auto addSHAttribute =
        [this, &metadata](vtkIdType itemID, const char* attributeName)
        {
          if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID ||
              !this->SubjectHierarchyNode->HasItemAttribute(itemID, attributeName))
          {
            return;
          }
          const std::string value =
              this->SubjectHierarchyNode->GetItemAttribute(itemID, attributeName);
          if (!value.empty())
          {
            metadata[QString::fromUtf8(attributeName)] =
                QString::fromStdString(value);
          }
        };

    for (const char* attr : {
             "DICOM.PatientSex",
             "DICOM.PatientAge",
             "DICOM.PatientWeight",
             "DICOM.PatientSize"})
    {
      addSHAttribute(this->patID, attr);
    }
    for (const char* attr : {
             "DICOM.StudyDate",
             "DICOM.StudyTime",
             "DICOM.StudyDescription",
             "DICOM.StudyID",
             "DICOM.StudyInstanceUID"})
    {
      addSHAttribute(this->stuID, attr);
    }

    if (this->sequencePETNode)
    {
      auto addSequenceAttribute =
          [this, &metadata](const char* outputKey, const char* attributeName)
          {
            const char* value = this->sequencePETNode->GetAttribute(attributeName);
            if (value && *value)
            {
              metadata[QString::fromUtf8(outputKey)] = QString::fromUtf8(value);
            }
          };

      addSequenceAttribute(
          "RadiopharmaceuticalStartDateTime",
          "RadiopharmaceuticalStartDateTime");
      addSequenceAttribute(
          "RadionuclideStartDateTime",
          "RadionuclideStartDateTime");
      addSequenceAttribute(
          "RadionuclideTotalDose",
          "RadionuclideTotalDose");
      addSequenceAttribute(
          "FirstFrameAcquisitionDateTime",
          "dPET.FirstFrameAcquisitionDateTime");
      addSequenceAttribute(
          "InjectionDateTimeSource",
          "dPET.InjectionDateTimeSource");
      addSequenceAttribute(
          "InjectionToAcquisitionOffset_s",
          "dPET.InjectionToAcquisitionOffsetSec");
    }

    d->updateAcquisitionTimingContext(false);
    if (d->acquisitionTiming.delayedAcquisition)
    {
      metadata["AcquisitionStartPostInjection_s"] =
          d->acquisitionTiming.acquisitionStartPostInjectionSec;
    }
  }

  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();
  const QVariant result = mainContext.call(
      "DPE_save_multisheet_excel",
      QVariantList{ fullPath, segmentDict, metadata });
  Q_UNUSED(result);
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
  Q_D(qSlicerDynamicPETModuleWidget);
  vtkMRMLScene* scene = this->mrmlScene();
  if (!scene || !mrmlPlotSeriesIDs || !selectionCol)
    return;

  // IF Preview selection is deliberately tolerant of overlapping curves.
  // If VTK hit-tests a derived line instead of the visible Original-source
  // marker, map the selected x coordinate to the nearest displayed raw CSV
  // observation. This makes Delete reliable without hiding processed curves.
  for (vtkIdType i = 0; i < mrmlPlotSeriesIDs->GetNumberOfValues(); ++i)
  {
    const std::string seriesID = mrmlPlotSeriesIDs->GetValue(i);
    vtkMRMLPlotSeriesNode* seriesNode = vtkMRMLPlotSeriesNode::SafeDownCast(
        scene->GetNodeByID(seriesID));
    if (!seriesNode)
      continue;

    const char* group = seriesNode->GetAttribute("SlicerDynamicPET.IFPreviewGroup");
    if (!group || std::string(group) != "InputFunctionPreview")
      continue;

    if (i >= selectionCol->GetNumberOfItems())
      return;

    vtkIdTypeArray* selectedPoints = vtkIdTypeArray::SafeDownCast(
        selectionCol->GetItemAsObject(i));
    if (!selectedPoints || selectedPoints->GetNumberOfTuples() < 1)
      return;

    const vtkIdType pointIndex =
        selectedPoints->GetValue(selectedPoints->GetNumberOfTuples() - 1);

    // Only external CSV observations are removable. Segment/Table ROI IDIF
    // points continue to use the normal TAC plot keep-mask workflow.
    if (d->IFSourceSelector->currentIndex() != 1 ||
        d->externalIFPreviewIndexMap.empty())
    {
      d->externalIFPreviewSelectedIndex = -1;
      return;
    }

    const char* sourceObservations =
        seriesNode->GetAttribute("SlicerDynamicPET.SourceObservations");
    const char* observationRole =
        seriesNode->GetAttribute("SlicerDynamicPET.SourceObservationRole");
    const bool directSourceSelection =
        sourceObservations && std::string(sourceObservations) == "1" &&
        observationRole && std::string(observationRole) == "ExternalIF";

    if (directSourceSelection && pointIndex >= 0 &&
        static_cast<size_t>(pointIndex) < d->externalIFPreviewIndexMap.size())
    {
      const size_t previewIndex = static_cast<size_t>(pointIndex);
      d->externalIFPreviewSelectedIndex =
          static_cast<int>(d->externalIFPreviewIndexMap[previewIndex]);
      if (previewIndex < d->externalIFPreviewTimesSec.size())
      {
        d->logToPythonConsole(
            QObject::tr("[SlicerDynamicPET IF] Selected external IF observation at %1 s; press Delete to exclude it.")
                .arg(d->externalIFPreviewTimesSec[previewIndex], 0, 'g', 10));
      }
      return;
    }

    vtkMRMLTableNode* tableNode = seriesNode->GetTableNode();
    vtkTable* table = tableNode ? tableNode->GetTable() : nullptr;
    const std::string xColumnName = seriesNode->GetXColumnName();
    if (!table || xColumnName.empty() || pointIndex < 0)
    {
      d->externalIFPreviewSelectedIndex = -1;
      return;
    }

    vtkAbstractArray* xArray = table->GetColumnByName(xColumnName.c_str());
    if (!xArray || pointIndex >= xArray->GetNumberOfTuples())
    {
      d->externalIFPreviewSelectedIndex = -1;
      return;
    }

    const double selectedX = xArray->GetVariantValue(pointIndex).ToDouble();
    const bool xInMinutes = xColumnName.find("min") != std::string::npos;
    const double xScale = xInMinutes ? (1.0 / 60.0) : 1.0;

    size_t nearestPreviewIndex = 0;
    double nearestDx = std::numeric_limits<double>::infinity();
    for (size_t j = 0; j < d->externalIFPreviewTimesSec.size(); ++j)
    {
      const double rawX = d->externalIFPreviewTimesSec[j] * xScale;
      const double dx = std::fabs(selectedX - rawX);
      if (dx < nearestDx)
      {
        nearestDx = dx;
        nearestPreviewIndex = j;
      }
    }

    // Use a local-spacing tolerance so dense early IF samples are not confused
    // with their neighbours while late sparse samples remain easy to select.
    const double nearestX =
        d->externalIFPreviewTimesSec[nearestPreviewIndex] * xScale;
    double localSpacing = std::numeric_limits<double>::infinity();
    if (nearestPreviewIndex > 0)
    {
      localSpacing = std::min(
          localSpacing,
          nearestX - d->externalIFPreviewTimesSec[nearestPreviewIndex - 1] * xScale);
    }
    if (nearestPreviewIndex + 1 < d->externalIFPreviewTimesSec.size())
    {
      localSpacing = std::min(
          localSpacing,
          d->externalIFPreviewTimesSec[nearestPreviewIndex + 1] * xScale - nearestX);
    }
    if (!std::isfinite(localSpacing) || localSpacing <= 0.0)
      localSpacing = 1.0 * xScale;

    const double tolerance = std::max(1e-9, 0.40 * localSpacing);
    if (nearestDx <= tolerance &&
        nearestPreviewIndex < d->externalIFPreviewIndexMap.size())
    {
      d->externalIFPreviewSelectedIndex =
          static_cast<int>(d->externalIFPreviewIndexMap[nearestPreviewIndex]);
      d->logToPythonConsole(
          QObject::tr("[SlicerDynamicPET IF] Selected external IF observation at %1 s; press Delete to exclude it.")
              .arg(d->externalIFPreviewTimesSec[nearestPreviewIndex], 0, 'g', 10));
    }
    else
    {
      d->externalIFPreviewSelectedIndex = -1;
    }
    return;
  }

  if (!this->checkdisplayedDynamicPET() || this->MapPlotSeriesNodeIDToPlot.empty())
    return;

  QSet<QPair<QString, vtkIdType>> newSelection;

  int psf_value = -1;
  std::string psv_value;
  std::string lastseriesID;
  for (vtkIdType i = 0; i < mrmlPlotSeriesIDs->GetNumberOfValues(); ++i)
  {
    QString seriesID = QString::fromStdString(mrmlPlotSeriesIDs->GetValue(i));
    vtkIdTypeArray* selectedPoints = vtkIdTypeArray::SafeDownCast(
        selectionCol->GetItemAsObject(i));

    if (!selectedPoints || selectedPoints->GetNumberOfTuples() < 1)
      continue;

    vtkIdType pointIndex = selectedPoints->GetValue(selectedPoints->GetNumberOfTuples()-1);
    QPair<QString, vtkIdType> candidate(seriesID, pointIndex);
    vtkMRMLPlotSeriesNode* seriesNode = vtkMRMLPlotSeriesNode::SafeDownCast(
        scene->GetNodeByID(seriesID.toStdString()));
    if (!seriesNode)
        continue;
    vtkMRMLTableNode* tableNode = seriesNode->GetTableNode();
    if (!tableNode)
        continue;
    vtkTable* table = tableNode->GetTable();
    if (!table)
        continue;
    const std::string labelName = seriesNode->GetLabelColumnName();
    if (labelName.empty())
        continue;
    vtkAbstractArray* labelArray = table->GetColumnByName(labelName.c_str());
    vtkStringArray* strArray = vtkStringArray::SafeDownCast(labelArray);
    if (!strArray || pointIndex >= strArray->GetNumberOfValues())
        continue;

    QString labelValue = QString::fromStdString(strArray->GetValue(pointIndex));
    QStringList parts = labelValue.split(',');
    if (!parts.isEmpty())
    {
      QString framePart = parts[0].trimmed();
      framePart.remove("Frame:");
      psf_value = framePart.trimmed().toInt();
      const auto idIt = this->ColNameToSegmentID.find(seriesNode->GetName());
      if (idIt != this->ColNameToSegmentID.end())
      {
        psv_value = idIt->second;
      }
    }
    if (!newSelection.contains(candidate))
      newSelection.insert(candidate);
    if (!this->lastSelection.contains(candidate))
    {
      this->PlotSelectedFrame = psf_value;
      this->PlotSelectedVOI = psv_value;
    }
    else
    {
      if (this->PlotSelectedFrame == psf_value && this->PlotSelectedVOI == psv_value)
      {
        lastseriesID = seriesID.toStdString();
        continue;
      }
      if (this->MapPlotSeriesNodeIDToPlot.contains(seriesID))
      {
        vtkPlot* vtkplot = this->MapPlotSeriesNodeIDToPlot.value(seriesID);
        vtkSmartPointer<vtkIdTypeArray> emptySelection = vtkSmartPointer<vtkIdTypeArray>::New();
        vtkplot->SetSelection(emptySelection);
      }
    }
  }

  if (newSelection.isEmpty())
  {
    this->PlotSelectedFrame = -1;
    this->PlotSelectedVOI.clear();
  }
  else if (!(newSelection.size() == 1 && !lastseriesID.empty()) && !lastseriesID.empty())
  {
    const QString lastID = QString::fromStdString(lastseriesID);
    if (this->MapPlotSeriesNodeIDToPlot.contains(lastID))
    {
      vtkPlot* vtkplot = this->MapPlotSeriesNodeIDToPlot.value(lastID);
      vtkSmartPointer<vtkIdTypeArray> emptySelection = vtkSmartPointer<vtkIdTypeArray>::New();
      vtkplot->SetSelection(emptySelection);
    }
  }

  this->lastSelection = newSelection;
}

void qSlicerDynamicPETModuleWidget::onPlotbutton()
{
  using PlotClock = std::chrono::steady_clock;
  const auto plotBuildStart = PlotClock::now();
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

  // Get selected statistics. The visible label may come from an external
  // workbook (for example "Activity"), while StatID remains the canonical
  // internal statistic used by the existing fitting/plotting code.
  std::vector<std::string> PlotSelectedStats;
  std::map<std::string, std::string> PlotStatLabels;
  for (int i = 0; i < d->PlotStatsCheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->PlotStatsCheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      QString statID = checkbox->property("StatID").toString();
      if (statID.isEmpty())
      {
        statID = checkbox->text();
      }
      const std::string id = statID.toStdString();
      PlotSelectedStats.push_back(id);
      PlotStatLabels[id] = checkbox->text().toStdString();
    }
  }

  if (PlotSelectedIDs.empty() || PlotSelectedStats.empty())
    return;

  if (std::find(PlotSelectedStats.begin(), PlotSelectedStats.end(), "Distribution") !=
      PlotSelectedStats.end())
  {
    if (d->isTableBasedMode())
    {
      QMessageBox::information(
          this, tr("ROI distribution"),
          tr("Voxel distributions require image data and are not available in Table mode."));
      return;
    }
    d->enforceDistributionSelection();
    d->updateDistributionFrameUI(false);

    std::string segmentID = d->lastPlotSegmentID;
    if (segmentID.empty() || this->segmentTACs.find(segmentID) == this->segmentTACs.end())
    {
      segmentID = PlotSelectedIDs.front();
    }
    QString distributionError;
    if (!d->plotROIDistribution(segmentID, &distributionError))
    {
      QMessageBox::warning(this, tr("ROI distribution"), distributionError);
    }
    return;
  }

  // Clear previous plot/chart/table
  this->RemoveExistingPlotChartAndTable();

  // Create or get table
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  // Add time column
  vtkNew<vtkDoubleArray> timeArray;
  timeArray->SetName("Time (min)");
  vtkNew<vtkStringArray> labelArray;
  labelArray->SetName("ToolTipLabelTAC");

  auto plotTimeSec = [this, d](size_t index)
  {
    double t = this->timePoints[index];
    if (d->isTableBasedMode() &&
        d->tablePlotTimesSec.size() == this->timePoints.size() &&
        index < d->tablePlotTimesSec.size())
    {
      t = d->tablePlotTimesSec[index];
    }

    if (d->acquisitionTiming.delayedAcquisition &&
        !d->acquisitionTiming.tableTimesAlreadyPostInjection)
    {
      t += d->acquisitionTiming.acquisitionStartPostInjectionSec;
    }
    return t;
  };

  for (int i = 0; i < this->timePoints.size(); ++i)
  {
    const double frameEndSec = d->frameEndForInputSec(static_cast<size_t>(i));
    const double frameStartSec =
        frameEndSec - this->durations[i];
    const double frameMidSec =
        frameStartSec + 0.5 * this->durations[i];

    timeArray->InsertNextValue(plotTimeSec(static_cast<size_t>(i)) / 60.0);

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

  const auto isVolumeStatistic = [](const std::string& name)
  {
    return name == "VolumePET";
  };
  bool allActivity = true;
  bool allVolume = true;
  for (const std::string& stat : PlotSelectedStats)
  {
    const bool activity =
        stat == "Mean" || stat == "Median" || stat == "Peak" ||
        stat == "Min" || stat == "Max";
    allActivity = allActivity && activity;
    allVolume = allVolume && isVolumeStatistic(stat);
  }
  if (allActivity)
  {
    chartNode->SetYAxisTitle(
        d->activityUnitLabel(displayUnit).toStdString().c_str());
  }
  else if (allVolume)
  {
    chartNode->SetYAxisTitle("ROI volume (cm3)");
  }
  else
  {
    chartNode->SetYAxisTitle("Value (mixed units)");
  }

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
      const std::string statDisplayName =
          PlotStatLabels.count(statName) ? PlotStatLabels[statName] : statName;
      std::string colName = segmentName + " - " + statDisplayName;
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
          else if (statName == "VolumePET")  value = vs.volume_cm3;
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
                x1 = plotTimeSec(static_cast<size_t>(next_ivs));
                if (statName == "Mean") nextValue = vsNext.mean;
                else if (statName == "Median") nextValue = vsNext.median;
                else if (statName == "Peak") nextValue = vsNext.peak;
                else if (statName == "VoxelCount") nextValue = vsNext.count;
                else if (statName == "Min") nextValue = vsNext.min;
                else if (statName == "Max") nextValue = vsNext.max;
                else if (statName == "VolumePET") nextValue = vsNext.volume_cm3;
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
                x0 = plotTimeSec(static_cast<size_t>(prev_ivs));
                if (statName == "Mean") prevValue = vsPrev.mean;
                else if (statName == "Median") prevValue = vsPrev.median;
                else if (statName == "Peak") prevValue = vsPrev.peak;
                else if (statName == "VoxelCount") prevValue = vsPrev.count;
                else if (statName == "Min") prevValue = vsPrev.min;
                else if (statName == "Max") prevValue = vsPrev.max;
                else if (statName == "VolumePET") prevValue = vsPrev.volume_cm3;
                prevValue = convertPlotActivity(prevValue, statName);
                break;
            }
          }
          if (std::isnan(prevValue)) {
            statArrayLine->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
            continue;
          }

          double x  = plotTimeSec(static_cast<size_t>(ivs));
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

      if (d->PlotErrorCheckbox && d->PlotErrorCheckbox->isChecked())
      {
        vtkSmartPointer<vtkMRMLTableNode> errorTable =
            vtkSmartPointer<vtkMRMLTableNode>::New();
        errorTable->SetName(
            (std::string("DynamicPET.ErrorTable.") + segmentID + "." + statName).c_str());
        scene->AddNode(errorTable);

        vtkNew<vtkDoubleArray> errorX;
        vtkNew<vtkDoubleArray> errorY;
        errorX->SetName("X");
        errorY->SetName("Y");

        for (size_t ivs = 0; ivs < this->segmentTACs[segmentID].size(); ++ivs)
        {
          const VoxelStatistics& vs = this->segmentTACs[segmentID][ivs];
          if (!vs.keep)
          {
            continue;
          }

          double center = std::numeric_limits<double>::quiet_NaN();
          if (statName == "Mean") center = vs.mean;
          else if (statName == "Median") center = vs.median;
          else if (statName == "Peak") center = vs.peak;
          else if (statName == "Max") center = vs.max;

          double sigma = d->tissueSigmaForWeighting(
              segmentID, ivs, statName, vs);
          if (!std::isfinite(center) || !std::isfinite(sigma) || sigma <= 0.0)
          {
            continue;
          }

          center = convertPlotActivity(center, statName);
          sigma = convertPlotActivity(sigma, statName);
          if (!std::isfinite(center) || !std::isfinite(sigma))
          {
            continue;
          }

          const double xMin = plotTimeSec(ivs) / 60.0;
          errorX->InsertNextValue(xMin);
          errorY->InsertNextValue(center - sigma);
          errorX->InsertNextValue(xMin);
          errorY->InsertNextValue(center + sigma);
          errorX->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
          errorY->InsertNextValue(std::numeric_limits<double>::quiet_NaN());
        }

        if (errorX->GetNumberOfTuples() > 0)
        {
          errorTable->AddColumn(errorX);
          errorTable->AddColumn(errorY);

          vtkSmartPointer<vtkMRMLPlotSeriesNode> errorSeries =
              vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
          scene->AddNode(errorSeries);
          errorSeries->SetName("");
          errorSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
          errorSeries->SetAndObserveTableNodeID(errorTable->GetID());
          errorSeries->SetXColumnName("X");
          errorSeries->SetYColumnName("Y");
          errorSeries->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
          errorSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleSolid);
          errorSeries->SetLineWidth(1.0);
          errorSeries->SetColor(lineSeries->GetColor());
          chartNode->AddAndObservePlotSeriesNodeID(errorSeries->GetID());
        }
        else
        {
          scene->RemoveNode(errorTable);
        }
      }
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
                         this, SLOT(onSelectedPoint(vtkStringArray*, vtkCollection*)),
                         Qt::UniqueConnection);

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

  const double plotBuildMs =
      std::chrono::duration<double, std::milli>(PlotClock::now() - plotBuildStart).count();
  d->logToPythonConsole(
      QObject::tr("[SlicerDynamicPET PERF][PLOT] mode=%1; curves=%2; statistics=%3; observations=%4; total=%5 ms.")
          .arg(d->isMultiTimepointMode() ? QStringLiteral("MULTI") :
               (d->isTableBasedMode() ? QStringLiteral("TABLE") : QStringLiteral("SINGLE")))
          .arg(static_cast<int>(PlotSelectedIDs.size()))
          .arg(static_cast<int>(PlotSelectedStats.size()))
          .arg(static_cast<int>(this->timePoints.size()))
          .arg(plotBuildMs, 0, 'f', 3));
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

  if (d->isTableBasedMode())
    d->tableIFID = this->IFID;
  else
    d->imageIFID = this->IFID;

  const bool previewOpen = d->previewGroupExists("InputFunctionPreview");

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

  if (previewOpen && index >= 0)
  {
    d->previewInputFunction();
  }
}

void qSlicerDynamicPETModuleWidget::onVOISelectionChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std :: string segmentID = d->VOISelector->itemData(index).toString().toStdString();
  this->plotTCMVOI = segmentID;
  d->populateResultsTable(segmentID);

  if (!d->syncingResultVOISelection)
  {
    d->syncingResultVOISelection = true;
    const int peerIndex = d->VOISelectorMTGA->findData(QString::fromStdString(segmentID));
    if (peerIndex >= 0)
    {
      QSignalBlocker blocker(d->VOISelectorMTGA);
      d->VOISelectorMTGA->setCurrentIndex(peerIndex);
      this->plotMTGAVOI = segmentID;
      d->populateResultsMTGATable(segmentID);
    }
    d->syncingResultVOISelection = false;
  }
}

void qSlicerDynamicPETModuleWidget::onVOIMTGASelectionChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std :: string segmentID = d->VOISelectorMTGA->itemData(index).toString().toStdString();
  this->plotMTGAVOI = segmentID;
  d->populateResultsMTGATable(segmentID);

  if (!d->syncingResultVOISelection)
  {
    d->syncingResultVOISelection = true;
    const int peerIndex = d->VOISelector->findData(QString::fromStdString(segmentID));
    if (peerIndex >= 0)
    {
      QSignalBlocker blocker(d->VOISelector);
      d->VOISelector->setCurrentIndex(peerIndex);
      this->plotTCMVOI = segmentID;
      d->populateResultsTable(segmentID);
    }
    d->syncingResultVOISelection = false;
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
  this->onVOISegmentsChanged();
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
  this->onVOIMTGASegmentsChanged();
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
  if (!d->standardFitCheckBox->isChecked())
  {
    d->standardFitCheckBox->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onStdFitImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->weightFitCheckBoxImg->blockSignals(true);
  d->weightFitCheckBoxImg->setChecked(false);
  d->weightFitCheckBoxImg->blockSignals(false);
  if (!d->standardFitCheckBoxImg->isChecked())
  {
    d->standardFitCheckBoxImg->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onWFitclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->standardFitCheckBox->blockSignals(true);
  d->standardFitCheckBox->setChecked(false);
  d->standardFitCheckBox->blockSignals(false);
  if (!d->weightFitCheckBox->isChecked())
  {
    d->standardFitCheckBox->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onWFitImgclicked()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  d->standardFitCheckBoxImg->blockSignals(true);
  d->standardFitCheckBoxImg->setChecked(false);
  d->standardFitCheckBoxImg->blockSignals(false);
  if (!d->weightFitCheckBoxImg->isChecked())
  {
    d->standardFitCheckBoxImg->setChecked(true);
  }
}

void qSlicerDynamicPETModuleWidget::onSliderChanged(int index)
{
  Q_D(qSlicerDynamicPETModuleWidget);
  const double timeSec =
      d->frameEndForInputSec(static_cast<size_t>(index - 1));
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
      d->frameEndForInputSec(static_cast<size_t>(index - 1));
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

  InputFunctionResult ifResult;
  QString ifError;

  if (!d->buildCurrentInputFunction(
          ifResult,
          true,
          &ifError) ||
      !ifResult.inputCoversFromInjection)
  {
      d->FITbuttonTCMImg->setEnabled(false);
      d->FITbuttonTCMImg->setToolTip(
          tr("Voxelwise TCM requires plasma and whole-blood input that covers injection onward. For a delayed acquisition, provide a full external IF or enable a PBIF template that reconstructs the missing early input."));
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

  if (d->acquisitionTiming.delayedAcquisition)
  {
    d->FITbuttonTCMImg->setToolTip(
        ifResult.inputCoverageReconstructedByPBIF
        ? tr("PBIF reconstructs the missing early input, so voxelwise TCM is available. A confirmation warning will be shown because the early tissue response was not observed and microparameter identifiability may be reduced.")
        : tr("Complete input from injection is available, so voxelwise TCM is available for this delayed tissue acquisition. A confirmation warning will be shown because the early tissue response was not observed and microparameter identifiability may be reduced."));
  }
  else
  {
    d->FITbuttonTCMImg->setToolTip(
        tr("Fit selected TCM models voxelwise and generate their parametric maps using the selected output options."));
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
  std::vector<std::string> selectedIDs;
  QSet<QString> selectedSet;

  for (int i = 0; i < d->VOICheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->VOICheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      const QString id = checkbox->property("SegmentID").toString();
      selectedSet.insert(id);
      selectedIDs.push_back(id.toStdString());
    }
  }
  this->VOIsegmentIDs = selectedIDs;

  if (!d->syncingVOICheckSelection)
  {
    d->syncingVOICheckSelection = true;
    std::vector<std::string> peerIDs;
    for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->VOIMTGACheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (!checkbox)
      {
        continue;
      }
      const QString id = checkbox->property("SegmentID").toString();
      const bool checked = selectedSet.contains(id);
      {
        QSignalBlocker blocker(checkbox);
        checkbox->setChecked(checked);
      }
      if (checked)
      {
        peerIDs.push_back(id.toStdString());
      }
    }
    this->VOIMTGAsegmentIDs = peerIDs;
    d->syncingVOICheckSelection = false;
  }

  this->enableFITbutton();
  this->enableFITMTGAbutton();
}

void qSlicerDynamicPETModuleWidget::onVOIMTGASegmentsChanged()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  std::vector<std::string> selectedIDs;
  QSet<QString> selectedSet;

  for (int i = 0; i < d->VOIMTGACheckLayout->count(); ++i)
  {
    QLayoutItem* item = d->VOIMTGACheckLayout->itemAt(i);
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
    if (checkbox && checkbox->isChecked())
    {
      const QString id = checkbox->property("SegmentID").toString();
      selectedSet.insert(id);
      selectedIDs.push_back(id.toStdString());
    }
  }
  this->VOIMTGAsegmentIDs = selectedIDs;

  if (!d->syncingVOICheckSelection)
  {
    d->syncingVOICheckSelection = true;
    std::vector<std::string> peerIDs;
    for (int i = 0; i < d->VOICheckLayout->count(); ++i)
    {
      QLayoutItem* item = d->VOICheckLayout->itemAt(i);
      QCheckBox* checkbox = qobject_cast<QCheckBox*>(item->widget());
      if (!checkbox)
      {
        continue;
      }
      const QString id = checkbox->property("SegmentID").toString();
      const bool checked = selectedSet.contains(id);
      {
        QSignalBlocker blocker(checkbox);
        checkbox->setChecked(checked);
      }
      if (checked)
      {
        peerIDs.push_back(id.toStdString());
      }
    }
    this->VOIsegmentIDs = peerIDs;
    d->syncingVOICheckSelection = false;
  }

  this->enableFITbutton();
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

    std::vector<double> segmentSigmas;
    double fallbackSigma = std::numeric_limits<double>::quiet_NaN();
    if (d->weightFitCheckBox->isChecked())
    {
      segmentSigmas.reserve(statsVec.size());
      for (size_t i = 0; i < statsVec.size(); ++i)
      {
        segmentSigmas.push_back(
            d->tissueSigmaForWeighting(
                segmentName, i, currentSelectedStatID, statsVec[i]));
      }
      fallbackSigma = medianValidSigma(segmentSigmas);
    }

    for (int ivs=0; ivs<statsVec.size(); ++ivs)
    {
      const auto& vs = statsVec[ivs];
      double value;
      if (currentSelectedStatID == "Mean") {
        value = vs.mean;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Median") {
        value = vs.median;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Peak") {
        value = vs.peak;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Max") {
        value = vs.max;
        if (d->weightFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
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

  if (!ifResult.inputCoversFromInjection)
  {
    QMessageBox::warning(
        this,
        tr("TCM input support"),
        tr("Compartment modeling requires an input function that covers injection onward. "
           "Provide a full external input function or enable a PBIF template that reconstructs the missing early input."));
    return;
  }

  if (d->acquisitionTiming.delayedAcquisition)
  {
    const QMessageBox::StandardButton choice = QMessageBox::warning(
        this,
        tr("Delayed tissue acquisition"),
        tr("This PET acquisition starts after injection, so the early tissue kinetics were not observed.\n\n"
           "Because a complete input function from injection is available, the compartment model can be propagated "
           "through the unobserved pre-scan interval and compared with the acquired late tissue frames. However, "
           "individual kinetic rate constants may be weakly identifiable or unstable without the early tissue response.\n\n"
           "Continue with TCM fitting?"),
        QMessageBox::Yes | QMessageBox::Cancel,
        QMessageBox::Cancel);
    if (choice != QMessageBox::Yes)
    {
      return;
    }
  }

  const double acquisitionStartSec =
      d->acquisitionTiming.delayedAcquisition
      ? d->frameStartForInputSec(0)
      : 0.0;

  if (d->acquisitionTiming.delayedAcquisition)
  {
    d->logToPythonConsole(
        tr("[SlicerDynamicPET timing] Tissue acquisition starts %1 s after injection. "
           "A complete input function is available, so TCM propagation begins at injection; "
           "late-only tissue sampling may reduce microparameter identifiability.")
            .arg(acquisitionStartSec, 0, 'g', 8));
  }

  // Preserve the actual acquisition windows.  This is identical to the old
  // cumulative schedule for a normal contiguous dynamic scan, but it also
  // allows Table/Multi-timepoint observations to contain real temporal gaps.
  std::vector<double> frameStartTimesSec(static_cast<size_t>(Nframe));
  std::vector<double> frameEndTimesSec(static_cast<size_t>(Nframe));
  bool hasTCMTemporalGaps = false;
  for (size_t i = 0; i < static_cast<size_t>(Nframe); ++i)
  {
    frameStartTimesSec[i] = d->frameStartForInputSec(i);
    frameEndTimesSec[i] = d->frameEndForInputSec(i);
    if (!std::isfinite(frameStartTimesSec[i]) ||
        !std::isfinite(frameEndTimesSec[i]) ||
        !(frameEndTimesSec[i] > frameStartTimesSec[i]) ||
        (i > 0 && frameStartTimesSec[i] < frameEndTimesSec[i - 1] - 1e-6))
    {
      QMessageBox::warning(
          this,
          tr("TCM timing"),
          tr("The PET observation schedule contains invalid or overlapping frame intervals."));
      return;
    }
    if (i > 0 && frameStartTimesSec[i] > frameEndTimesSec[i - 1] + 1e-6)
    {
      hasTCMTemporalGaps = true;
    }
  }

  if (hasTCMTemporalGaps)
  {
    d->logToPythonConsole(
        tr("[SlicerDynamicPET timing] TCM fit contains separated acquisition windows. "
           "The compartment state and input function are propagated continuously through each gap; "
           "no tissue residual is created where no PET observation exists."));
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

    // Keep the original tissue-observation visibility mask for plotting.
    // The selected TCM end controls fitting only; valid observations after
    // that point remain visible as excluded/out-of-fit measurements.

    if (fitFrameCount < 2)
    {
      const auto nameIt = segmentTACsnames.find(segmentID);
      const QString roiName =
          nameIt != segmentTACsnames.end()
          ? QString::fromStdString(nameIt->second)
          : QString::fromStdString(segmentID);

      QMessageBox::warning(
          this,
          tr("TCM fit range"),
          tr("ROI '%1': At least two retained tissue/IF frames are required for TCM fitting.")
              .arg(roiName));
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
    std::vector<double> frameStartTimesFit(
        frameStartTimesSec.begin(),
        frameStartTimesSec.begin() + fitFrameCount);
    std::vector<double> frameEndTimesFit(
        frameEndTimesSec.begin(),
        frameEndTimesSec.begin() + fitFrameCount);
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                       ifResult.plasmaIsParent,
                       acquisitionStartSec,
                       &frameStartTimesFit,
                       &frameEndTimesFit
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
                : &ifResult.nativeWholeBloodValues,
            acquisitionStartSec,
            &frameStartTimesFit,
            &frameEndTimesFit);
      }
      else
      {
        std::cerr << "Unknown model ID: " << modelID << std::endl;
        return;
      }

      // Diagnostic is reported on exactly the observations used by the optimizer.
      // The plotted prediction is extended only afterwards, so fit diagnostics
      // remain completely independent of the display/prediction support.

      // Diagnostic reported on exactly the observations used by the optimizer.
      // This makes it easy to distinguish a genuinely weighted-AIC preference
      // from a plotting/data-support mismatch.
      double unweightedSSE = 0.0;
      double weightedSSE = 0.0;
      size_t diagnosticCount = 0;
      double* diagnosticFit = this->segmentTCMfits[segmentID][modelID];
      if (diagnosticFit)
      {
        for (size_t i = 0; i < fitFrameCount; ++i)
        {
          if (i >= weightFit.size() || !(weightFit[i] > 0.0))
          {
            continue;
          }
          const double residual = tacFit[i][0] - diagnosticFit[i];
          if (!std::isfinite(residual))
          {
            continue;
          }
          unweightedSSE += residual * residual;
          weightedSSE += weightFit[i] * residual * residual;
          ++diagnosticCount;
        }
      }
      const auto parameterIt = this->segmentTCM[segmentID].find(modelID);
      if (parameterIt != this->segmentTCM[segmentID].end())
      {
        d->logToPythonConsole(
            tr("[SlicerDynamicPET TCM] ROI=%1 | model=%2 | frames=%3 | weighting=%4 | AIC=%5 | SSE=%6 | weightedSSE=%7")
                .arg(QString::fromStdString(this->segmentTACsnames[segmentID]))
                .arg(QString::fromStdString(modelID))
                .arg(diagnosticCount)
                .arg(d->weightFitCheckBox->isChecked() ? "WLS" : "OLS")
                .arg(parameterIt->second.AIC, 0, 'g', 8)
                .arg(unweightedSSE, 0, 'g', 8)
                .arg(weightedSSE, 0, 'g', 8));
      }

      const size_t predictionFrameCount =
          std::min(ifSupportCount, static_cast<size_t>(Nframe));
      std::vector<double> fullPrediction;
      std::string predictionError;
      bool predictionOK = false;

      if (predictionFrameCount >= fitFrameCount && parameterIt != this->segmentTCM[segmentID].end())
      {
        // callTCM stores parameters that were fixed during optimization as NaN.
        // Restore the model-definition fixed values before display-only forward
        // prediction; fitted parameters themselves are never changed/refit.
        TCMParameters predictionParams = parameterIt->second;
        if (modelID == "1TCM" || modelID == "1TiCM" ||
            modelID == "2TCM" || modelID == "2TiCM")
        {
          predictionParams.td = 0.0;
        }
        if (modelID == "1TiCM" || modelID == "1TidCM")
        {
          predictionParams.k2 = 0.0;
        }
        if (modelID == "2TiCM" || modelID == "2TidCM")
        {
          predictionParams.k4 = 0.0;
        }
        if (modelID == "Liver DBIF")
        {
          predictionParams.td = 0.0;
        }

        std::vector<std::vector<double>> cpPrediction(
            Cp.begin(), Cp.begin() + predictionFrameCount);
        std::vector<std::vector<double>> cwbPrediction(
            Cwb.begin(), Cwb.begin() + predictionFrameCount);
        std::vector<std::vector<double>> framingPrediction(
            framing.begin(), framing.begin() + predictionFrameCount);
        std::vector<double> frameStartPrediction(
            frameStartTimesSec.begin(), frameStartTimesSec.begin() + predictionFrameCount);
        std::vector<double> frameEndPrediction(
            frameEndTimesSec.begin(), frameEndTimesSec.begin() + predictionFrameCount);

        try
        {
          if (modelID == "Liver DBIF")
          {
            predictionOK = logic->PredictLiverTCM(
                cwbPrediction, framingPrediction, predictionParams,
                dk, timestep, fullPrediction, interpolationType,
                ifResult.nativeWholeBloodTimesSec.empty() ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                ifResult.nativeWholeBloodValues.empty() ? nullptr : &ifResult.nativeWholeBloodValues,
                acquisitionStartSec, &frameStartPrediction, &frameEndPrediction,
                &predictionError);
          }
          else
          {
            const int predictionCompartments =
                (!modelID.empty() && modelID[0] == '2') ? 2 : 1;
            predictionOK = logic->PredictTCM(
                cpPrediction, cwbPrediction, framingPrediction,
                predictionParams, predictionCompartments, dk, timestep,
                fullPrediction, interpolationType,
                ifResult.nativePlasmaTimesSec.empty() ? nullptr : &ifResult.nativePlasmaTimesSec,
                ifResult.nativePlasmaValues.empty() ? nullptr : &ifResult.nativePlasmaValues,
                ifResult.nativeWholeBloodTimesSec.empty() ? nullptr : &ifResult.nativeWholeBloodTimesSec,
                ifResult.nativeWholeBloodValues.empty() ? nullptr : &ifResult.nativeWholeBloodValues,
                ifResult.parentFractionTimesSec.empty() ? nullptr : &ifResult.parentFractionTimesSec,
                ifResult.parentFractionValues.empty() ? nullptr : &ifResult.parentFractionValues,
                ifResult.plasmaIsParent, acquisitionStartSec,
                &frameStartPrediction, &frameEndPrediction, &predictionError);
          }
        }
        catch (const std::exception& e)
        {
          predictionOK = false;
          predictionError = e.what();
        }
      }

      if (!predictionOK && predictionFrameCount > fitFrameCount)
      {
        d->logToPythonConsole(
            tr("[SlicerDynamicPET TCM] Could not extend %1 prediction beyond the fit end: %2")
                .arg(QString::fromStdString(modelID))
                .arg(QString::fromStdString(predictionError)));
      }
      else if (predictionOK && predictionFrameCount > fitFrameCount)
      {
        d->logToPythonConsole(
            tr("[SlicerDynamicPET TCM] %1 fitted through frame %2; prediction propagated through IF-supported frame %3 without refitting.")
                .arg(QString::fromStdString(modelID))
                .arg(fitFrameCount)
                .arg(predictionFrameCount));
      }

      replaceFittedCurveForPlot(
          this->segmentTCMfits[segmentID][modelID],
          fitFrameCount,
          static_cast<size_t>(Nframe),
          predictionOK ? &fullPrediction : nullptr);
    }

    const auto liverFitIt = this->segmentTCMfits[segmentID].find("Liver DBIF");
    const auto tdFitIt = this->segmentTCMfits[segmentID].find("2TdCM");
    if (liverFitIt != this->segmentTCMfits[segmentID].end() &&
        tdFitIt != this->segmentTCMfits[segmentID].end() &&
        liverFitIt->second && tdFitIt->second)
    {
      const bool distinctBuffers = liverFitIt->second != tdFitIt->second;
      double maxAbsDifference = 0.0;
      double maxMagnitude = 0.0;
      for (size_t i = 0; i < fitFrameCount; ++i)
      {
        const double a = liverFitIt->second[i];
        const double b = tdFitIt->second[i];
        if (std::isfinite(a) && std::isfinite(b))
        {
          maxAbsDifference = std::max(maxAbsDifference, std::abs(a - b));
          maxMagnitude = std::max(maxMagnitude, std::max(std::abs(a), std::abs(b)));
        }
      }
      const double relativeDifference =
          maxAbsDifference / std::max(1.0, maxMagnitude);
      d->logToPythonConsole(
          tr("[SlicerDynamicPET liver] Liver DBIF vs 2TdCM | independent buffers=%1 | max curve difference=%2 | relative=%3")
              .arg(distinctBuffers ? "yes" : "NO")
              .arg(maxAbsDifference, 0, 'g', 6)
              .arg(relativeDifference, 0, 'g', 6));
      if (!distinctBuffers)
      {
        qCritical() << "Liver DBIF and 2TdCM unexpectedly share the same fitted-curve buffer.";
      }
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

  if (!ifResult.inputCoversFromInjection)
  {
      QMessageBox::warning(
          this,
          tr("TCM input support"),
          tr("Voxelwise compartment modeling requires an input function that covers injection onward. "
             "Provide a full external input function or a PBIF that reconstructs the missing early input."));
      return;
  }

  if (d->acquisitionTiming.delayedAcquisition)
  {
    const QMessageBox::StandardButton choice = QMessageBox::warning(
        this,
        tr("Delayed tissue acquisition"),
        tr("This PET acquisition starts after injection, so the early tissue kinetics were not observed.\n\n"
           "A complete input function allows voxelwise compartment propagation from injection, but the missing early "
           "tissue response can make individual kinetic rate constants weakly identifiable or unstable.\n\n"
           "Continue with voxelwise TCM fitting?"),
        QMessageBox::Yes | QMessageBox::Cancel,
        QMessageBox::Cancel);
    if (choice != QMessageBox::Yes)
    {
      return;
    }
  }

  const double acquisitionStartSec =
      d->acquisitionTiming.delayedAcquisition
      ? d->frameStartForInputSec(0)
      : 0.0;

  if (d->acquisitionTiming.delayedAcquisition)
  {
    d->logToPythonConsole(
        tr("[SlicerDynamicPET timing] Tissue acquisition starts %1 s after injection. "
           "A complete input function is available, so TCM propagation begins at injection; "
           "late-only tissue sampling may reduce microparameter identifiability.")
            .arg(acquisitionStartSec, 0, 'g', 8));
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
        appendDouble(key, acquisitionStartSec);
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
    ifResult.plasmaIsParent,
    acquisitionStartSec
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

    std::vector<double> segmentSigmas;
    double fallbackSigma = std::numeric_limits<double>::quiet_NaN();
    if (d->weightedFitCheckBox->isChecked())
    {
      segmentSigmas.reserve(statsVec.size());
      for (size_t i = 0; i < statsVec.size(); ++i)
      {
        segmentSigmas.push_back(
            d->tissueSigmaForWeighting(
                segmentName, i, currentSelectedStatID, statsVec[i]));
      }
      fallbackSigma = medianValidSigma(segmentSigmas);
    }

    for (int ivs=0; ivs<statsVec.size(); ++ivs)
    {
      const auto& vs = statsVec[ivs];
      double value;
      if (currentSelectedStatID == "Mean") {
        value = vs.mean;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Median") {
        value = vs.median;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Peak") {
        value = vs.peak;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
        } else {
          wgtVec[segmentName].push_back(1.);
        }
      }
      else if (currentSelectedStatID == "Max") {
        value = vs.max;
        if (d->weightedFitCheckBox->isChecked()) {
          wgtVec[segmentName].push_back(
              inverseVarianceWeightFromSigma(
                  segmentSigmas[static_cast<size_t>(ivs)], fallbackSigma));
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

  const double framingNorm = d->framingNormEdit->text().toDouble();
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

  const size_t startFrameIndex =
      static_cast<size_t>(d->timeOffsetSlider->value() - 1);
  const size_t firstUsableIndex =
      std::min(ifResult.supportFrameStartIndex, static_cast<size_t>(Nframe));

  // Explicit observation schedule shared by all ROI MTGA methods.  For a
  // standard contiguous dynamic scan this is numerically identical to the
  // previous cumulative-duration clock.  With separated acquisitions, the
  // real unobserved intervals are retained instead of collapsed.
  std::vector<double> frameStartTimesSec(static_cast<size_t>(Nframe));
  std::vector<double> frameEndTimesSec(static_cast<size_t>(Nframe));
  std::vector<double> frameMidTimesSec(static_cast<size_t>(Nframe));
  bool hasTemporalGaps = false;
  for (size_t i = 0; i < static_cast<size_t>(Nframe); ++i)
  {
      frameStartTimesSec[i] = d->frameStartForInputSec(i);
      frameEndTimesSec[i] = d->frameEndForInputSec(i);
      frameMidTimesSec[i] = d->frameMidForInputSec(i);

      if (!std::isfinite(frameStartTimesSec[i]) ||
          !std::isfinite(frameEndTimesSec[i]) ||
          !std::isfinite(frameMidTimesSec[i]) ||
          !(frameEndTimesSec[i] > frameStartTimesSec[i]) ||
          (i > 0 && frameStartTimesSec[i] < frameEndTimesSec[i - 1] - 1e-6))
      {
          QMessageBox::warning(
              this,
              tr("MTGA timing"),
              tr("The PET observation schedule contains invalid or overlapping frame intervals."));
          return;
      }

      if (i > 0 && frameStartTimesSec[i] > frameEndTimesSec[i - 1] + 1e-6)
      {
          hasTemporalGaps = true;
      }
  }

  if (hasTemporalGaps)
  {
    d->logToPythonConsole(
        tr("[SlicerDynamicPET timing] MTGA fit contains separated acquisition windows; temporal gaps are preserved."));

    const bool tissueIntegralModelSelected =
        std::find(this->modelsMTGAID.begin(), this->modelsMTGAID.end(), "Logan") !=
            this->modelsMTGAID.end() ||
        std::find(this->modelsMTGAID.begin(), this->modelsMTGAID.end(), "RE") !=
            this->modelsMTGAID.end() ||
        std::find(this->modelsMTGAID.begin(), this->modelsMTGAID.end(), "Relative RE") !=
            this->modelsMTGAID.end();
    if (tissueIntegralModelSelected)
    {
      d->logToPythonConsole(
          tr("[SlicerDynamicPET timing] Logan/RE tissue integrals bridge only the unobserved acquisition gaps "
             "by linear interpolation between neighbouring tissue frame averages."));
    }
  }

  const double timeOffset =
      startFrameIndex < frameMidTimesSec.size()
      ? frameMidTimesSec[startFrameIndex] / framingNorm
      : 0.0;

  double initialPlasmaIntegralNormalized = 0.0;
  const bool standardPatlakSelected =
      std::find(this->modelsMTGAID.begin(), this->modelsMTGAID.end(), "Patlak")
      != this->modelsMTGAID.end();
  if (standardPatlakSelected &&
      d->acquisitionTiming.delayedAcquisition &&
      ifResult.inputCoversFromInjection)
  {
      const double acquisitionStartSec = d->frameStartForInputSec(0);
      const double integralSec =
          d->initialModelPlasmaIntegralSec(ifResult, acquisitionStartSec);
      if (!std::isfinite(integralSec))
      {
          QMessageBox::warning(
              this,
              tr("Patlak input support"),
              tr("The complete pre-scan model plasma integral could not be reconstructed. "
                 "Standard Patlak is unavailable for this delayed acquisition."));
          return;
      }
      initialPlasmaIntegralNormalized = integralSec / framingNorm;
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

  // When a true acquisition gap exists, preferentially integrate the fully
  // processed model plasma curve on its native clock. This preserves IF area
  // through the unobserved interval. If that native representation is not
  // available, the logic layer falls back to an explicit linear bridge
  // between neighbouring plasma frame averages.
  std::vector<double> plasmaIntegralAtMidSec;
  if (hasTemporalGaps && ifSupportCount > 0)
  {
      const double integralOriginSec =
          (d->acquisitionTiming.delayedAcquisition &&
           ifResult.inputCoversFromInjection)
          ? 0.0
          : frameStartTimesSec.front();

      plasmaIntegralAtMidSec.reserve(ifSupportCount);
      bool nativeIntegralAvailable = true;
      for (size_t i = 0; i < ifSupportCount; ++i)
      {
          const double integral =
              d->integrateModelPlasmaOverInterval(
                  ifResult, integralOriginSec, frameMidTimesSec[i]);
          if (!std::isfinite(integral))
          {
              nativeIntegralAvailable = false;
              break;
          }
          plasmaIntegralAtMidSec.push_back(integral);
      }

      if (!nativeIntegralAvailable)
      {
          plasmaIntegralAtMidSec.clear();
          d->logToPythonConsole(
              tr("[SlicerDynamicPET timing] Full native plasma integration across the acquisition gap was unavailable; "
                 "MTGA will use the documented linear bridge between neighbouring plasma frame averages."));
      }
  }

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
        const size_t supportStart =
            std::min(firstUsableIndex, values.size());
        if (supportStart >= values.size())
        {
          return false;
        }

        const bool hasRemoved =
            std::any_of(
                stats.begin() + static_cast<std::ptrdiff_t>(supportStart),
                stats.begin() + static_cast<std::ptrdiff_t>(values.size()),
                [](const VoxelStatistics& vs)
                {
                  return !vs.keep;
                });

        if (!hasRemoved)
        {
          return true;
        }

        if (values.size() - supportStart < 2 || !stats[supportStart].keep)
        {
          QMessageBox::warning(
              this,
              tr("MTGA point exclusion"),
              tr("The first tissue observation inside the usable IF/tissue interval must remain available. "
                 "A removed terminal observation is handled as a shorter acquisition."));
          return false;
        }

        std::vector<double> retainedTimes;
        std::vector<double> retainedValues;

        for (size_t i = supportStart; i < values.size(); ++i)
        {
          if (!stats[i].keep)
          {
            continue;
          }

          retainedTimes.push_back(
              frameMidTimesSec[i]);
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

        for (size_t i = supportStart; i < values.size(); ++i)
        {
          if (stats[i].keep)
          {
            continue;
          }

          const double t =
              frameMidTimesSec[i];

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

    // Regression end controls the fitted subset only. Keep the full tissue/IF
    // support available to the graphical transformation so plots can show
    // observations outside the selected regression window.
    size_t tissueSupportCount = 0;
    for (size_t i = ifSupportCount; i > 0; --i)
    {
      if (statsIt->second[i - 1].keep)
      {
        tissueSupportCount = i;
        break;
      }
    }

    const size_t modelSupportCount =
        std::min(ifSupportCount, tissueSupportCount);
    const size_t regressionEndCount =
        std::min(requestedMTGAEndCount, modelSupportCount);

    if (regressionEndCount < 2 ||
        startFrameIndex < firstUsableIndex ||
        startFrameIndex >= regressionEndCount ||
        regressionEndCount - startFrameIndex < 2)
    {
      const auto nameIt = segmentTACsnames.find(segmentID);
      const QString roiName =
          nameIt != segmentTACsnames.end()
          ? QString::fromStdString(nameIt->second)
          : QString::fromStdString(segmentID);

      QMessageBox::warning(
          this,
          tr("MTGA fit range"),
          tr("ROI '%1': The selected MTGA start/end range is outside the usable IF/tissue acquisition range.")
              .arg(roiName));
      continue;
    }

    auto tac_flatten = extractColumn(tacVOI);
    tac_flatten.resize(modelSupportCount);

    std::vector<double> cpFit(
        ifResult.frameModelPlasma.begin(),
        ifResult.frameModelPlasma.begin() + modelSupportCount);
    std::vector<double> framingFit(
        framing_flatten.begin(),
        framing_flatten.begin() + modelSupportCount);
    std::vector<double> frameStartTimesFit(
        frameStartTimesSec.begin(),
        frameStartTimesSec.begin() + modelSupportCount);
    std::vector<double> frameEndTimesFit(
        frameEndTimesSec.begin(),
        frameEndTimesSec.begin() + modelSupportCount);
    std::vector<double> plasmaIntegralFit;
    const std::vector<double>* plasmaIntegralFitPtr = nullptr;
    if (plasmaIntegralAtMidSec.size() >= modelSupportCount)
    {
        plasmaIntegralFit.assign(
            plasmaIntegralAtMidSec.begin(),
            plasmaIntegralAtMidSec.begin() + modelSupportCount);
        plasmaIntegralFitPtr = &plasmaIntegralFit;
    }
    std::vector<double> weightFit(
        wgtVec[segmentID].begin(),
        wgtVec[segmentID].begin() + modelSupportCount);
    for (size_t i = regressionEndCount; i < weightFit.size(); ++i)
    {
        weightFit[i] = 0.0;
    }

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
                      max_iter,
                      initialPlasmaIntegralNormalized,
                      &frameStartTimesFit,
                      &frameEndTimesFit,
                      plasmaIntegralFitPtr
                      );
      }
      else if (modelID == "Relative Patlak") {
        logic->RelativePatlak(tac_flatten,
                              cpFit,
                              framingFit,
                              this->segmentMTGA[segmentID]["Relative Patlak"],
                              wgt,
                              timeOffset,
                              framingNorm,
                              robust,
                              std,
                              huber_tune,
                              tol,
                              max_iter,
                              firstUsableIndex,
                              &frameStartTimesFit,
                              &frameEndTimesFit,
                              plasmaIntegralFitPtr);
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
                     max_iter,
                     &frameStartTimesFit,
                     &frameEndTimesFit,
                     plasmaIntegralFitPtr
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
                  max_iter,
                  &frameStartTimesFit,
                  &frameEndTimesFit,
                  plasmaIntegralFitPtr
                 );
      }
      else if (modelID == "Relative RE") {
        logic->RelativeRE(tac_flatten,
                          cpFit,
                          framingFit,
                          this->segmentMTGA[segmentID]["Relative RE"],
                          wgt,
                          timeOffset,
                          framingNorm,
                          robust,
                          std,
                          huber_tune,
                          tol,
                          max_iter,
                          firstUsableIndex,
                          &frameStartTimesFit,
                          &frameEndTimesFit,
                          plasmaIntegralFitPtr);
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

  const size_t firstUsableIndex =
      std::min(ifResult.supportFrameStartIndex, static_cast<size_t>(Nframe));
  const size_t startFrameIndex =
      static_cast<size_t>(d->timeOffsetSliderImg->value() - 1);

  if (fitFrameCount < 2 ||
      startFrameIndex < firstUsableIndex ||
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
  double timeOffset = 0.0;
  for (size_t i = 0; i < startFrameIndex && i < durations.size(); ++i)
  {
      timeOffset += durations[i] / framingNorm;
  }
  if (startFrameIndex < durations.size())
  {
      timeOffset += 0.5 * durations[startFrameIndex] / framingNorm;
  }

  double initialPlasmaIntegralNormalized = 0.0;
  const bool standardPatlakSelected =
      std::find(this->modelsMTGAImgID.begin(), this->modelsMTGAImgID.end(), "Patlak")
      != this->modelsMTGAImgID.end();
  if (standardPatlakSelected &&
      d->acquisitionTiming.delayedAcquisition &&
      ifResult.inputCoversFromInjection)
  {
      const double acquisitionStartSec = d->frameStartForInputSec(0);
      const double integralSec =
          d->initialModelPlasmaIntegralSec(ifResult, acquisitionStartSec);
      if (!std::isfinite(integralSec))
      {
          QMessageBox::warning(
              this,
              tr("Patlak input support"),
              tr("The complete pre-scan model plasma integral could not be reconstructed. "
                 "Standard Patlak is unavailable for this delayed acquisition."));
          return;
      }
      initialPlasmaIntegralNormalized = integralSec / framingNorm;
  }

  const bool robust = d->robustFitCheckBoxImg->isChecked();
  const bool std = false;  // MTGA parametric imaging intentionally uses physical graphical variables.
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
            + "|standardize=0";

        appendDouble(key, framingNorm);
        appendDouble(key, timeOffset);
        appendDouble(key, initialPlasmaIntegralNormalized);
        key += "|dataStartIndex=" + QString::number(static_cast<qulonglong>(firstUsableIndex));
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
          numThreads,
          initialPlasmaIntegralNormalized,
          firstUsableIndex);

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
  int enabledCount = 0;
  int enabledCheckedCount = 0;
  for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsCheckLayout->itemAt(i)->widget());
    if (!checkbox || !checkbox->isEnabled()) continue;
    ++enabledCount;
    if (checkbox->isChecked()) ++enabledCheckedCount;
  }

  const bool clearEnabled = enabledCount > 0 && enabledCheckedCount == enabledCount;
  this->modelsID.clear();
  for (int i = 0; i < d->ModelsCheckLayout->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsCheckLayout->itemAt(i)->widget());
    if (!checkbox) continue;
    checkbox->blockSignals(true);
    if (checkbox->isEnabled())
    {
      checkbox->setChecked(!clearEnabled);
      if (!clearEnabled)
      {
        this->modelsID.push_back(checkbox->text().toStdString());
      }
    }
    else
    {
      checkbox->setChecked(false);
    }
    checkbox->blockSignals(false);
  }
  d->ModelsCheckContents->blockSignals(false);
  d->updateLiverParameterUI();
  this->enableFITbutton();
}


void qSlicerDynamicPETModuleWidget::onModelsMTGAAllbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  d->ModelsMTGACheckContents->blockSignals(true);
  int enabledCount = 0;
  int enabledCheckedCount = 0;
  for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsMTGACheckLayout->itemAt(i)->widget());
    if (!checkbox || !checkbox->isEnabled()) continue;
    ++enabledCount;
    if (checkbox->isChecked()) ++enabledCheckedCount;
  }

  const bool clearEnabled = enabledCount > 0 && enabledCheckedCount == enabledCount;
  this->modelsMTGAID.clear();
  for (int i = 0; i < d->ModelsMTGACheckLayout->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsMTGACheckLayout->itemAt(i)->widget());
    if (!checkbox) continue;
    checkbox->blockSignals(true);
    if (checkbox->isEnabled())
    {
      checkbox->setChecked(!clearEnabled);
      if (!clearEnabled)
      {
        this->modelsMTGAID.push_back(checkbox->text().toStdString());
      }
    }
    else
    {
      checkbox->setChecked(false);
    }
    checkbox->blockSignals(false);
  }
  d->ModelsMTGACheckContents->blockSignals(false);
  this->enableFITMTGAbutton();
}


void qSlicerDynamicPETModuleWidget::onModelsSelectAllMTGAImgbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  d->ModelsCheckContentsMTGAImg->blockSignals(true);
  int enabledCount = 0;
  int enabledCheckedCount = 0;
  for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsCheckLayoutMTGAImg->itemAt(i)->widget());
    if (!checkbox || !checkbox->isEnabled()) continue;
    ++enabledCount;
    if (checkbox->isChecked()) ++enabledCheckedCount;
  }

  const bool clearEnabled = enabledCount > 0 && enabledCheckedCount == enabledCount;
  this->modelsMTGAImgID.clear();
  for (int i = 0; i < d->ModelsCheckLayoutMTGAImg->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsCheckLayoutMTGAImg->itemAt(i)->widget());
    if (!checkbox) continue;
    checkbox->blockSignals(true);
    if (checkbox->isEnabled())
    {
      checkbox->setChecked(!clearEnabled);
      if (!clearEnabled)
      {
        this->modelsMTGAImgID.push_back(checkbox->text().toStdString());
      }
    }
    else
    {
      checkbox->setChecked(false);
    }
    checkbox->blockSignals(false);
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
    if (checkbox && checkbox->isEnabled() && checkbox->isChecked())
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
    if (checkbox && checkbox->isEnabled() && checkbox->isChecked())
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
    if (checkbox && checkbox->isEnabled() && checkbox->isChecked())
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
    if (checkbox && checkbox->isEnabled() && checkbox->isChecked())
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
  int enabledCount = 0;
  int enabledCheckedCount = 0;
  for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsCheckLayoutTCMImg->itemAt(i)->widget());
    if (!checkbox || !checkbox->isEnabled()) continue;
    ++enabledCount;
    if (checkbox->isChecked()) ++enabledCheckedCount;
  }

  const bool clearEnabled = enabledCount > 0 && enabledCheckedCount == enabledCount;
  this->modelsTCMImgID.clear();
  for (int i = 0; i < d->ModelsCheckLayoutTCMImg->count(); ++i)
  {
    QCheckBox* checkbox = qobject_cast<QCheckBox*>(d->ModelsCheckLayoutTCMImg->itemAt(i)->widget());
    if (!checkbox) continue;
    checkbox->blockSignals(true);
    if (checkbox->isEnabled())
    {
      checkbox->setChecked(!clearEnabled);
      if (!clearEnabled)
      {
        this->modelsTCMImgID.push_back(checkbox->text().toStdString());
      }
    }
    else
    {
      checkbox->setChecked(false);
    }
    checkbox->blockSignals(false);
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
      if (modelName == "Liver DBIF")
      {
        lineSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleDash);
        lineSeries->SetLineWidth(3.0);
      }
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
  std::string selectedVOI = this->plotMTGAVOI;

  QString modelName_qstr = d->MTGASelector->currentText();
  if (modelName_qstr.toStdString().empty())
    return;
  std::string modelName = modelName_qstr.toStdString();

  auto voiIt = this->segmentMTGA.find(selectedVOI);
  if (voiIt == this->segmentMTGA.end())
    return;
  auto modelIt = voiIt->second.find(modelName);
  if (modelIt == voiIt->second.end())
    return;

  const MTGAParameters& params = modelIt->second;
  const bool haveFullPlot =
      !params.plotX.empty() &&
      params.plotX.size() == params.plotY.size() &&
      params.plotX.size() == params.plotFitted.size() &&
      params.plotX.size() == params.plotFrame.size() &&
      params.plotX.size() == params.plotIncluded.size();

  const std::vector<double>& plotX = haveFullPlot ? params.plotX : params.x;
  const std::vector<double>& plotY = haveFullPlot ? params.plotY : params.y;
  const std::vector<double>& plotFit = haveFullPlot ? params.plotFitted : params.fitted;
  const std::vector<int>& plotFrame = haveFullPlot ? params.plotFrame : params.frame;

  if (plotX.empty() || plotY.size() != plotX.size() || plotFit.size() != plotX.size())
    return;

  this->RemoveExistingPlotChartAndTable();
  vtkSmartPointer<vtkMRMLTableNode> tableNode = this->GetOrCreatePlotTable();

  vtkNew<vtkDoubleArray> xArray;
  vtkNew<vtkDoubleArray> includedArray;
  vtkNew<vtkDoubleArray> excludedArray;
  vtkNew<vtkStringArray> labelArray;
  xArray->SetName("X");
  includedArray->SetName("Data");
  excludedArray->SetName("Excluded from regression");
  labelArray->SetName("ToolTipData");

  for (size_t i = 0; i < plotX.size(); ++i)
  {
    const bool included = haveFullPlot
        ? params.plotIncluded[i]
        : (i < params.keep.size() ? params.keep[i] : true);
    const int frameIndex =
        i < plotFrame.size() ? plotFrame[i] - 1 : static_cast<int>(i);
    bool originalObservationAvailable = true;
    const auto statsIt = this->segmentTACs.find(selectedVOI);
    if (statsIt != this->segmentTACs.end() &&
        frameIndex >= 0 &&
        frameIndex < static_cast<int>(statsIt->second.size()))
    {
      originalObservationAvailable = statsIt->second[static_cast<size_t>(frameIndex)].keep;
    }

    xArray->InsertNextValue(plotX[i]);
    includedArray->InsertNextValue(
        included && originalObservationAvailable
        ? plotY[i]
        : std::numeric_limits<double>::quiet_NaN());
    excludedArray->InsertNextValue(
        !included && originalObservationAvailable
        ? plotY[i]
        : std::numeric_limits<double>::quiet_NaN());

    std::ostringstream oss;
    oss << "Frame: " << frameIndex;
    if (frameIndex >= 0 && frameIndex < static_cast<int>(this->timePoints.size()))
    {
      oss << ", Time(s): " << this->timePoints[frameIndex]
          << ", Time(min): " << this->timePoints[frameIndex] / 60.0;
    }
    if (!originalObservationAvailable)
    {
      oss << ", observation removed";
    }
    else
    {
      oss << (included ? ", included in regression" : ", excluded from regression");
    }
    labelArray->InsertNextValue(oss.str());
  }
  tableNode->AddColumn(xArray);
  tableNode->AddColumn(includedArray);
  tableNode->AddColumn(excludedArray);
  tableNode->AddColumn(labelArray);

  vtkMRMLPlotChartNode* chartNode = this->GetOrCreatePlotChart();

  vtkNew<vtkDoubleArray> fitArray;
  fitArray->SetName(modelName.c_str());
  for (double value : plotFit)
  {
    fitArray->InsertNextValue(value);
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

  this->ColNameToSegmentID["Data"] = selectedVOI;
  vtkSmartPointer<vtkMRMLPlotSeriesNode> scatterSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(scatterSeries);
  scatterSeries->SetName("Included in regression");
  scatterSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
  scatterSeries->SetAndObserveTableNodeID(tableNode->GetID());
  scatterSeries->SetXColumnName("X");
  scatterSeries->SetYColumnName("Data");
  scatterSeries->SetLabelColumnName("ToolTipData");
  scatterSeries->SetUniqueColor();
  scatterSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(scatterSeries->GetID());

  this->ColNameToSegmentID["Excluded from regression"] = selectedVOI;
  vtkSmartPointer<vtkMRMLPlotSeriesNode> excludedSeries = vtkSmartPointer<vtkMRMLPlotSeriesNode>::New();
  scene->AddNode(excludedSeries);
  excludedSeries->SetName("Excluded from regression");
  excludedSeries->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
  excludedSeries->SetAndObserveTableNodeID(tableNode->GetID());
  excludedSeries->SetXColumnName("X");
  excludedSeries->SetYColumnName("Excluded from regression");
  excludedSeries->SetLabelColumnName("ToolTipData");
  excludedSeries->SetUniqueColor();
  excludedSeries->SetLineStyle(vtkMRMLPlotSeriesNode::LineStyleNone);
  chartNode->AddAndObservePlotSeriesNodeID(excludedSeries->GetID());

  chartNode->SetTitle((modelName + " - " + this->segmentTACsnames[selectedVOI]).c_str());
  if (modelName == "Patlak" || modelName == "Relative Patlak") {
    chartNode->SetXAxisTitle("intCp/Cp");
    chartNode->SetYAxisTitle("Ct/Cp");
  } else if (modelName == "Logan") {
    chartNode->SetXAxisTitle("intCp/Ct");
    chartNode->SetYAxisTitle("intCt/Ct");
  } else if (modelName == "RE" || modelName == "Relative RE") {
    chartNode->SetXAxisTitle("intCp/Cp");
    chartNode->SetYAxisTitle("intCt/Cp");
  } else {
    std::cerr << "Unknown model: " << modelName << std::endl;
  }

  auto* layoutNode = vtkMRMLLayoutNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLLayoutNode"));
  if (layoutNode)
    layoutNode->SetViewArrangement(vtkMRMLLayoutNode::SlicerLayoutConventionalPlotView);

  vtkMRMLPlotViewNode* plotViewNode = vtkMRMLPlotViewNode::SafeDownCast(
      scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"));
  if (plotViewNode) {
    plotViewNode->SetPlotChartNodeID(chartNode->GetID());
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

void qSlicerDynamicPETModuleWidget::onDeleteKeyPressed()
{
  Q_D(qSlicerDynamicPETModuleWidget);

  // External IF preview deletion is mode-specific. The raw CSV is shared,
  // but Image and Table mode have independent keep masks.
  vtkMRMLPlotViewNode* activeView = this->mrmlScene()
      ? vtkMRMLPlotViewNode::SafeDownCast(
            this->mrmlScene()->GetFirstNodeByClass("vtkMRMLPlotViewNode"))
      : nullptr;
  vtkMRMLPlotChartNode* activeChart =
      (activeView && activeView->GetPlotChartNodeID())
      ? vtkMRMLPlotChartNode::SafeDownCast(
            this->mrmlScene()->GetNodeByID(activeView->GetPlotChartNodeID()))
      : nullptr;
  const bool inputPreviewDisplayed =
      activeChart && activeChart->GetName() &&
      std::string(activeChart->GetName()) ==
          "DynamicPET.InputFunctionPreview.Chart";
  if (inputPreviewDisplayed &&
      d->externalIFPreviewSelectedIndex >= 0 &&
      d->IFSourceSelector->currentIndex() == 1)
  {
    const size_t rawIndex = static_cast<size_t>(d->externalIFPreviewSelectedIndex);
    std::vector<bool>& keep = d->activeExternalIFKeep();
    if (keep.size() != d->externalIFTimesSec.size())
    {
      keep.assign(d->externalIFTimesSec.size(), true);
    }

    if (rawIndex < keep.size())
    {
      if (d->externalIFZeroAnchorAdded && rawIndex == 0)
      {
        d->logToPythonConsole(
            QObject::tr("[SlicerDynamicPET IF] The assumed (0 s, 0) anchor is protected."));
      }
      else
      {
        const size_t retainedCount = static_cast<size_t>(
            std::count(keep.begin(), keep.end(), true));
        if (keep[rawIndex] && retainedCount > 2)
        {
          keep[rawIndex] = false;
          d->logToPythonConsole(
              QObject::tr("[SlicerDynamicPET IF] Excluded external IF observation at %1 s.")
                  .arg(d->externalIFTimesSec[rawIndex], 0, 'g', 10));
          d->invalidateInputFunctionResults();
          d->updateInputFunctionStatus();
          this->clearFITdata();
          this->clearFITMTGAdata();
          d->previewInputFunction();
        }
      }
    }

    d->externalIFPreviewSelectedIndex = -1;
    return;
  }

  if (!inputPreviewDisplayed)
  {
    d->externalIFPreviewSelectedIndex = -1;
  }
  if (!this->checkdisplayedDynamicPET())
    return;
  if (this->PlotSelectedFrame < 0 || this->PlotSelectedVOI.empty())
    return;

  const auto tacIt = this->segmentTACs.find(this->PlotSelectedVOI);
  if (tacIt == this->segmentTACs.end() ||
      static_cast<size_t>(this->PlotSelectedFrame) >= tacIt->second.size())
    return;

  vtkGenericWarningMacro(
      "Segment " << this->segmentTACsnames[this->PlotSelectedVOI]
      << " removed at frame " << this->PlotSelectedFrame);
  this->segmentTACs[this->PlotSelectedVOI][this->PlotSelectedFrame].keep = false;
  this->segmentTACs[this->PlotSelectedVOI][this->PlotSelectedFrame].empty = false;

  // Persist immediately in the active source state so deletion never leaks
  // across Image/Table mode or disappears on a mode toggle.
  if (d->isTableBasedMode())
    d->captureActiveTACState(d->tableTACState);
  else
    d->captureActiveTACState(d->imageTACState);

  if (this->PlotSelectedVOI == this->IFID)
  {
    d->invalidateInputFunctionResults();
    d->updateInputFunctionStatus();
  }

  this->clearFITdata();
  this->clearFITMTGAdata();
  this->onPlotbutton();
  this->PlotSelectedFrame = -1;
  this->PlotSelectedVOI.clear();
}

void qSlicerDynamicPETModuleWidget::onResetbutton()
{
  Q_D(qSlicerDynamicPETModuleWidget);
  vtkGenericWarningMacro("Restoring removed points");

  bool changed = false;
  bool currentInputWasModified = false;

  for (auto& segmentPair : this->segmentTACs)
  {
    const bool isCurrentInput =
        d->IFSourceSelector->currentIndex() == 0 &&
        segmentPair.first == this->IFID;

    for (VoxelStatistics& vs : segmentPair.second)
    {
      if (!vs.empty && !vs.keep)
      {
        if (isCurrentInput)
          currentInputWasModified = true;
        vs.keep = true;
        changed = true;
      }
    }
  }

  std::vector<bool>& externalKeep = d->activeExternalIFKeep();
  if (!externalKeep.empty())
  {
    const bool hadRemoved = std::any_of(
        externalKeep.begin(), externalKeep.end(), [](bool keep){ return !keep; });
    if (hadRemoved)
    {
      std::fill(externalKeep.begin(), externalKeep.end(), true);
      d->logToPythonConsole(
          QObject::tr("[SlicerDynamicPET IF] Restored excluded external IF observations for the active mode."));
      changed = true;
      if (d->IFSourceSelector->currentIndex() == 1)
        currentInputWasModified = true;
    }
  }

  if (d->isTableBasedMode())
    d->captureActiveTACState(d->tableTACState);
  else
    d->captureActiveTACState(d->imageTACState);

  if (currentInputWasModified)
  {
    d->invalidateInputFunctionResults();
    d->updateInputFunctionStatus();
  }

  if (changed)
  {
    this->clearFITdata();
    this->clearFITMTGAdata();
  }

  vtkMRMLScene* scene = this->mrmlScene();
  vtkMRMLPlotViewNode* plotViewNode = scene
      ? vtkMRMLPlotViewNode::SafeDownCast(scene->GetFirstNodeByClass("vtkMRMLPlotViewNode"))
      : nullptr;
  vtkMRMLPlotChartNode* currentPlot =
      (scene && plotViewNode && plotViewNode->GetPlotChartNodeID())
      ? vtkMRMLPlotChartNode::SafeDownCast(scene->GetNodeByID(plotViewNode->GetPlotChartNodeID()))
      : nullptr;

  if (currentPlot && currentPlot->GetName() &&
      std::string(currentPlot->GetName()) == "DynamicPET.InputFunctionPreview.Chart" &&
      d->IFSourceSelector->currentIndex() == 1)
  {
    d->previewInputFunction();
  }
  else if (this->checkdisplayedDynamicPET())
  {
    this->onPlotbutton();
  }
}
