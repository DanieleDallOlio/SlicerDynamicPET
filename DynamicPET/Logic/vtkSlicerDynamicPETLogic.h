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

// .NAME vtkSlicerDynamicPETLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerDynamicPETLogic_h
#define __vtkSlicerDynamicPETLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"
#include <QApplication>

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerDynamicPETModuleLogicExport.h"
#include "vtkMRMLSubjectHierarchyNode.h"
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLSegmentationNode.h>
#include <vtkMRMLSequenceNode.h>
#include <vtkMRMLSequenceBrowserNode.h>
#include "vtkOrientedImageDataResample.h"
#include <vtkSegmentation.h>
#include <QPointer>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMRMLTableNode.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkOrientedImageData.h>
#include "vtkSlicerSegmentationsModuleLogic.h"
#include <QProgressBar>
#include <kmodels.h>
#include <utils.hpp>
#ifdef HAVE_OPENMP
// #pragma message("HAVE_OPENMP is defined!")
#include <omp.h>
// #else
// #pragma message("HAVE_OPENMP is NOT defined!")
#endif
#include <atomic>
#include <QMetaObject>
#include <QPushButton>
#include <vtkSmartPointer.h>
#include <vtkImageLogic.h>

struct VoxelStatistics
{
  int count = 0;
  double mean = 0.0;
  double median = 0.0;
  double min = 0.;
  double max = 0.;
  double stddev = 0.0;
  double q1 = 0.0;
  double q3 = 0.0;
  double iqr = 0.0;
  double volume_mm3 = 0.0;
  double volume_cm3 = 0.0;
  double peak = 0.0;
  double peakStddev = 0.0;
  int peakCount = 0;
  bool keep = true;
  bool empty = false;
};

struct MTGAParameters
{
  double Ki = std::numeric_limits<double>::quiet_NaN();
  double Intercept = std::numeric_limits<double>::quiet_NaN();
  double DV = std::numeric_limits<double>::quiet_NaN();
  double AIC = std::numeric_limits<double>::quiet_NaN();
  double MASE = std::numeric_limits<double>::quiet_NaN();
  double R2 = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int dof = 0;
  std::vector<double> x, y, fitted, weights, r;
  std::vector<int> frame;
  std::vector<bool> keep;
};


enum TCMBoundFlag : unsigned int
{
  TCM_BOUND_NONE      = 0u,
  TCM_BOUND_VB_LOWER  = 1u << 0,
  TCM_BOUND_VB_UPPER  = 1u << 1,
  TCM_BOUND_K1_LOWER  = 1u << 2,
  TCM_BOUND_K1_UPPER  = 1u << 3,
  TCM_BOUND_K2_LOWER  = 1u << 4,
  TCM_BOUND_K2_UPPER  = 1u << 5,
  TCM_BOUND_K3_LOWER  = 1u << 6,
  TCM_BOUND_K3_UPPER  = 1u << 7,
  TCM_BOUND_K4_LOWER  = 1u << 8,
  TCM_BOUND_K4_UPPER  = 1u << 9,
  TCM_BOUND_TD_LOWER  = 1u << 10,
  TCM_BOUND_TD_UPPER  = 1u << 11,
  TCM_BOUND_KA_LOWER  = 1u << 12,
  TCM_BOUND_KA_UPPER  = 1u << 13,
  TCM_BOUND_FA_LOWER  = 1u << 14,
  TCM_BOUND_FA_UPPER  = 1u << 15
};

struct FengParameters
{
  double tau = std::numeric_limits<double>::quiet_NaN();
  double A1 = std::numeric_limits<double>::quiet_NaN();
  double A2 = std::numeric_limits<double>::quiet_NaN();
  double A3 = std::numeric_limits<double>::quiet_NaN();
  double lambda1 = std::numeric_limits<double>::quiet_NaN();
  double lambda2 = std::numeric_limits<double>::quiet_NaN();
  double lambda3 = std::numeric_limits<double>::quiet_NaN();
  double SSE = std::numeric_limits<double>::quiet_NaN();
};

struct TCMParameters
{
  double K1 = std::numeric_limits<double>::quiet_NaN();
  double k2 = std::numeric_limits<double>::quiet_NaN();
  double k3 = std::numeric_limits<double>::quiet_NaN();
  double k4 = std::numeric_limits<double>::quiet_NaN();
  double ka = std::numeric_limits<double>::quiet_NaN();
  double fa = std::numeric_limits<double>::quiet_NaN();
  double vb = std::numeric_limits<double>::quiet_NaN();
  double td = std::numeric_limits<double>::quiet_NaN();
  double Ki = std::numeric_limits<double>::quiet_NaN();
  double DV = std::numeric_limits<double>::quiet_NaN();
  double AIC = std::numeric_limits<double>::quiet_NaN();
  double MASE = std::numeric_limits<double>::quiet_NaN();
  double BIC = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double loglik = std::numeric_limits<double>::quiet_NaN();
  int dof = 0;
  unsigned int boundFlags = TCM_BOUND_NONE;
  std::vector<double> weights, r;
  std::vector<bool> keep;
};


struct ModelComparisonResult
{
  std::string type;
  // double statistic = std::numeric_limits<double>::quiet_NaN();
  double p_value = std::numeric_limits<double>::quiet_NaN();
};

enum class VuongCorrection { None, AIC, BIC };

enum class Tail { TwoSided, Model1Greater, Model2Greater };

enum class PETCompositeMode { UnweightedSum, DurationWeightedSum };

enum class BodySupportCombination { Union, Intersection };

class VTK_SLICER_DYNAMICPET_MODULE_LOGIC_EXPORT vtkSlicerDynamicPETLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerDynamicPETLogic *New();
  vtkTypeMacro(vtkSlicerDynamicPETLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  void computeTAC(vtkIdType ctNode, vtkIdType petNode, vtkIdType segNode, std::vector<QString> segments, std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs, std::map<std::string, std::string>& segmentTACsnames, QProgressBar* ProgressBar, QPushButton* stopButton, std::atomic<bool>& stopRequested);
  void setupSeg(vtkMRMLSegmentationNode* segNode);
  VoxelStatistics ComputeVoxelStatistics(vtkMRMLScalarVolumeNode* PETVolume, vtkImageData* labelmap, int labelValue = 1);
  void TAC(vtkMRMLSequenceNode* sequencePETNode, vtkMRMLSequenceNode* segSequenceNode, std::vector<QString> segmentsID, std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs, std::map<std::string, std::string>& segmentTACsnames, QProgressBar* ProgressBar, QPushButton* stopButton, std::atomic<bool>& stopRequested);
  double computeLogLik(const std::vector<double>& y,
                       const std::vector<double>& fitted,
                       const std::vector<double>* weights);
  bool FitFengInputFunction(
      const std::vector<double>& timesSec,
      const std::vector<double>& values,
      const std::vector<double>* frameStartSec,
      const std::vector<double>* frameEndSec,
      bool observationsAreFrameAverages,
      FengParameters& params,
      std::vector<double>& fittedObservationValues,
      std::string* errorMessage = nullptr);
  double EvaluateFengInputFunction(
      double timeSec,
      const FengParameters& params) const;
  double AverageFengInputFunction(
      double frameStartSec,
      double frameEndSec,
      const FengParameters& params) const;
  void callTCM(
     std::vector<std::vector<double>> tac,
     std::vector<std::vector<double>> Cp,
     std::vector<std::vector<double>> Cwb,
     std::vector<std::vector<double>> framing,
     long int Nframe,
     long int Nvox,
     double* kinit,
     double* lb,
     double* ub,
     const bool* sens,
     const double dk,
     const double timestep,
     const int maxiter,
     const int n_tc,
     TCMParameters& params,
     double*& fitted_curve,
     const std::vector<double>* wgt,
     const std::string& interpolationType = "linear",
     const std::vector<double>* nativePlasmaTimesSec = nullptr,
     const std::vector<double>* nativePlasmaValues = nullptr,
     const std::vector<double>* nativeWholeBloodTimesSec = nullptr,
     const std::vector<double>* nativeWholeBloodValues = nullptr,
     const std::vector<double>* parentFractionTimesSec = nullptr,
     const std::vector<double>* parentFractionValues = nullptr,
     bool plasmaIsParent = false);

  // ROI-only optimization-derived liver dual-blood-input model.
  // The supplied arterial input is total whole blood. The model internally
  // derives the portal-vein component and estimates ka and fa.
  void callLiverTCM(
     const std::vector<std::vector<double>>& tac,
     const std::vector<std::vector<double>>& Cwb,
     const std::vector<std::vector<double>>& framing,
     long int Nframe,
     double* kinit,
     double* lb,
     double* ub,
     const bool* sens,
     const double dk,
     const double timestep,
     const int maxiter,
     TCMParameters& params,
     double*& fitted_curve,
     const std::vector<double>* wgt,
     const std::string& interpolationType = "linear",
     const std::vector<double>* nativeWholeBloodTimesSec = nullptr,
     const std::vector<double>* nativeWholeBloodValues = nullptr);

  // void getFittedTCM(double *& fitted_curve,
  //                   std :: vector< std :: vector<double> > Cp,
  //                   std :: vector< std :: vector<double> > framing,
  //                   long int Nframe,
  //                   long int Nvox,
  //                   double* kinit,
  //                   double* lb,
  //                   double* ub,
  //                   const bool* sens,
  //                   const double dk,
  //                   const double timestep,
  //                   const double pbrp[],
  //                   const int maxiter,
  //                   const int n_tc,
  //                   TCMParameters& params
  //                   );
  double computeAIC(const std::vector<double>& obs,
                    const std::vector<double>& est,
                    int numpar,
                    const std::vector<double>* wgt = nullptr,
                    bool aicc = true);
  double computeBIC(const std::vector<double>& obs,
                    const std::vector<double>& est,
                    int numpar,
                    const std::vector<double>* wgt);
  double computeR2(const std::vector<double>& obs,
                   const std::vector<double>& est,
                   const std::vector<double>* wgt);
  double computeChi2(const std::vector<double>& y,
                     const std::vector<double>& fitted,
                     const std::vector<double>* weights);
  double computeVuongP(const std::vector<double>& r1,
                       const std::vector<double>& r2,
                       const std::vector<double>* wgt,
                       int k1, int k2,
                       VuongCorrection corr,
                       Tail tail);
  double boundaryLRTPvalue(double LR, int r_b, int r_i);
  ModelComparisonResult compareModels(
      const std::string& modelA,
      const std::string& modelB,
      const TCMParameters& m1,
      const TCMParameters& m2
  );
  double computeLRTP(const double logLik1, const double logLik2, int df2, int df1);
  double MASE(
         const std::vector<double>& Actual,
         const std::vector<double>& Predicted,
         const std::vector<double>* wgt = nullptr
         );
  void Patlak(
       const std::vector<double>& tac,
       const std::vector<double>& Cp,
       const std::vector<double>& framing,
       MTGAParameters & params,
       const std::vector<double>* wgt = nullptr, // nullptr if not using weights
       const double timeOffset = 0.,
       const double framingNorm = 60.,
       bool robust = false,
       bool std = true,
       double huber_tune = 1.345,
       double tol = 1e-6,
       int max_iter = 50
       );
  void Logan(
       const std::vector<double>& tac,
       const std::vector<double>& Cp,
       const std::vector<double>& framing,
       MTGAParameters & params,
       const std::vector<double>* wgt = nullptr, // nullptr if not using weights
       const double timeOffset = 0.,
       const double framingNorm = 60.,
       bool robust = false,
       bool std = true,
       double huber_tune = 1.345,
       double tol = 1e-6,
       int max_iter = 50
       );
  void RE(
       const std::vector<double>& tac,
       const std::vector<double>& Cp,
       const std::vector<double>& framing,
       MTGAParameters & params,
       const std::vector<double>* wgt = nullptr, // nullptr if not using weights
       const double timeOffset = 0.,
       const double framingNorm = 60.,
       bool robust = false,
       bool std = true,
       double huber_tune = 1.345,
       double tol = 1e-6,
       int max_iter = 50
       );
  void Image2Flatten(
       vtkIdType petID,
       std::vector<std::vector<double>>& flatten_voxels_values,
       int (&dims)[3],
       int& numberOfTimepoints,
       QProgressBar* ProgressBar,
       QPushButton* stopButton,
       std::atomic<bool>& stopRequested
      );
  vtkSmartPointer<vtkOrientedImageData>
  CreateFullPETSupportMask(
      vtkMRMLScalarVolumeNode* referencePETNode);
  vtkSmartPointer<vtkOrientedImageData>
  CreateCTBodySupportMask(
      vtkMRMLScalarVolumeNode* ctNode,
      vtkMRMLScalarVolumeNode* referencePETNode,
      double ctThresholdHU = -500.0,
      double bodyMarginMm = 5.0,
      bool fillHoles = true);
  vtkSmartPointer<vtkOrientedImageData>
  CreatePETBodySupportMask(
      const std::vector<std::vector<double>>& voxels,
      const int dims[3],
      vtkMRMLScalarVolumeNode* referencePETNode,
      const std::vector<double>& frameDurations,
      PETCompositeMode compositeMode,
      double bodyMarginMm,
      bool fillHoles,
      double minComponentFractionOfLargest,
      double* thresholdOut = nullptr,
      bool* usedOtsuFallbackOut = nullptr);
  vtkMRMLSegmentationNode*
  CreateOrUpdateBodySupportPreview(
      vtkOrientedImageData* mask,
      vtkMRMLScalarVolumeNode* referencePETNode,
      const std::string& nodeName =
          "SlicerDynamicPET - Fitting Support");
  vtkSmartPointer<vtkOrientedImageData>
  CombineBodySupportMasks(
      vtkOrientedImageData* firstMask,
      vtkOrientedImageData* secondMask,
      BodySupportCombination combination);
  vtkMRMLScalarVolumeNode* Flatten2Image(
      const std::vector<double>& flatten_values,
      const int dims[3],
      const std::string& name
  );
  void Patlak4Img(
      const std::vector<std::vector<double>>& voxels,
      const std::vector<double>& Cp,
      const std::vector<double>& framing,
      const std::vector<double>* wgt_global,
      double timeOffset,
      double framingNorm,
      bool robust,
      bool standardize,
      double huber_tune,
      double tol,
      int max_iter,
      std::vector<MTGAParameters>& outputParams,
      const std::vector<int>& fitVoxelIndices,
      std::atomic<bool>& stopRequested,
      int numThreads,
      std::function<void(int)> progressCallback);
  void Logan4Img(
      const std::vector<std::vector<double>>& voxels,
      const std::vector<double>& Cp,
      const std::vector<double>& framing,
      const std::vector<double>* wgt_global,
      double timeOffset,
      double framingNorm,
      bool robust,
      bool standardize,
      double huber_tune,
      double tol,
      int max_iter,
      std::vector<MTGAParameters>& outputParams,
      const std::vector<int>& fitVoxelIndices,
      std::atomic<bool>& stopRequested,
      int numThreads,
      std::function<void(int)> progressCallback);
  void RE4Img(
      const std::vector<std::vector<double>>& voxels,
      const std::vector<double>& Cp,
      const std::vector<double>& framing,
      const std::vector<double>* wgt_global,
      double timeOffset,
      double framingNorm,
      bool robust,
      bool standardize,
      double huber_tune,
      double tol,
      int max_iter,
      std::vector<MTGAParameters>& outputParams,
      const std::vector<int>& fitVoxelIndices,
      std::atomic<bool>& stopRequested,
      int numThreads,
      std::function<void(int)> progressCallback);
  std::vector<double> ExtractParameter(
      const std::vector<MTGAParameters>& outputParams,
      const std::string& field);
  void CreateMTGAParametricImages(
      const std::vector<MTGAParameters>& outputParams,
      const int dims[3],
      const std::vector<std::string>& fields,
      const std::string& modelID,
      vtkMRMLScalarVolumeNode* refNode,
      vtkMRMLSubjectHierarchyNode* refSH,
      vtkIdType refID
    );
  void callTCMImg(
      const std::vector<std::vector<double>>& voxels,   // [Nvoxels][Nframe]
      const std::vector<double>& Cp,                    // plasma [Nframe]
      const std::vector<double>& Cwb,                   // total whole blood [Nframe]
      const std::vector<double>& framing,               // [Nframe]
      double* kinit,
      double* lb,
      double* ub,
      const bool* sens,
      const double dk,
      const double timestep,
      const int maxiter,
      const int n_tc,
      std::vector<TCMParameters>& outputParams,
      const std::string& modelID,
      const std::vector<int>& fitVoxelIndices,
      std::atomic<bool>& stopRequested,
      const std::vector<double>* wgt_global = nullptr,
      int numThreads = 1,
      std::function<void(int)> progressCallback = nullptr,
      std::function<bool()> stopCallback = nullptr,
      const std::string& interpolationType = "linear",
      const std::vector<double>* nativePlasmaTimesSec = nullptr,
      const std::vector<double>* nativePlasmaValues = nullptr,
      const std::vector<double>* nativeWholeBloodTimesSec = nullptr,
      const std::vector<double>* nativeWholeBloodValues = nullptr,
      const std::vector<double>* parentFractionTimesSec = nullptr,
      const std::vector<double>* parentFractionValues = nullptr,
      bool plasmaIsParent = false
      );
  std::vector<double> ExtractParameter(
      const std::vector<TCMParameters>& outputParams,
      const std::string& field);
  void CreateTCMParametricImages(
      const std::vector<TCMParameters>& outputParams,
      const int dims[3],
      const std::vector<std::string>& fields,
      const std::string& modelID,
      vtkMRMLScalarVolumeNode* refNode,
      vtkMRMLSubjectHierarchyNode* refSH,
      vtkIdType refID
    );
protected:
  vtkSlicerDynamicPETLogic();
  ~vtkSlicerDynamicPETLogic() override;

  void SetMRMLSceneInternal(vtkMRMLScene* newScene) override;
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  void RegisterNodes() override;
  void UpdateFromMRMLScene() override;
  void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;
  void OnMRMLSceneNodeRemoved(vtkMRMLNode* node) override;
private:

  vtkSlicerDynamicPETLogic(const vtkSlicerDynamicPETLogic&); // Not implemented
  void operator=(const vtkSlicerDynamicPETLogic&); // Not implemented
};

#endif
