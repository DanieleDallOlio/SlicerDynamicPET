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

// .NAME vtkSlicerKMAPLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerKMAPLogic_h
#define __vtkSlicerKMAPLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"
#include <QApplication>

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerKMAPModuleLogicExport.h"
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

struct TCMParameters
{
  double K1 = std::numeric_limits<double>::quiet_NaN();
  double k2 = std::numeric_limits<double>::quiet_NaN();
  double k3 = std::numeric_limits<double>::quiet_NaN();
  double k4 = std::numeric_limits<double>::quiet_NaN();
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
  std::vector<double> weights, r;
  std::vector<bool> keep;
};

enum class VuongCorrection { None, AIC, BIC };

enum class Tail { TwoSided, Model1Greater, Model2Greater };

class VTK_SLICER_KMAP_MODULE_LOGIC_EXPORT vtkSlicerKMAPLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerKMAPLogic *New();
  vtkTypeMacro(vtkSlicerKMAPLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  void computeTAC(vtkIdType ctNode, vtkIdType petNode, vtkIdType segNode, std::vector<QString> segments, std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs, std::map<std::string, std::string>& segmentTACsnames, QProgressBar* ProgressBar, QPushButton* stopButton, std::atomic<bool>& stopRequested);
  void setupSeg(vtkMRMLSegmentationNode* segNode);
  VoxelStatistics ComputeVoxelStatistics(vtkMRMLScalarVolumeNode* PETVolume, vtkImageData* labelmap, int labelValue = 1);
  void TAC(vtkMRMLSequenceNode* sequencePETNode, vtkMRMLSequenceNode* segSequenceNode, std::vector<QString> segmentsID, std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs, std::map<std::string, std::string>& segmentTACsnames, QProgressBar* ProgressBar, QPushButton* stopButton, std::atomic<bool>& stopRequested);
  double computeLogLik(const std::vector<double>& y,
                       const std::vector<double>& fitted,
                       const std::vector<double>* weights);
  void callTCM(std :: vector< std :: vector<double> > tac,
               std :: vector< std :: vector<double> > Cp,
               std :: vector< std :: vector<double> > framing,
               long int Nframe,
               long int Nvox,
               double* kinit,
               double* lb,
               double* ub,
               const bool* sens,
               const double dk,
               const double timestep,
               const double pbrp[],
               const int maxiter,
               const int n_tc,
               TCMParameters& params,
               double*& fitted_curve,
               const std::vector<double>* wgt
               );
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
      std::atomic<bool>& stopRequested,
      QProgressBar* progressBar = nullptr,
      int numThreads = 1,
      QPushButton* stopButton = nullptr
  );
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
      std::atomic<bool>& stopRequested,
      QProgressBar* progressBar = nullptr,
      int numThreads = 1,
      QPushButton* stopButton = nullptr
  );
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
      std::atomic<bool>& stopRequested,
      QProgressBar* progressBar = nullptr,
      int numThreads = 1,
      QPushButton* stopButton = nullptr
  );
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
      const std::vector<double>& Cp,                    // [Nframe]
      const std::vector<double>& framing,               // [Nframe]
      double* kinit,
      double* lb,
      double* ub,
      const bool* sens,
      const double dk,
      const double timestep,
      const double pbrp[],
      const int maxiter,
      const int n_tc,
      std::vector<TCMParameters>& outputParams,
      const std::string& modelID,
      std::atomic<bool>& stopRequested,
      const std::vector<double>* wgt_global = nullptr,
      // QProgressBar* progressBar = nullptr,
      int numThreads = 1,
      // QPushButton* stopButton = nullptr,
      std::function<void(int)> progressCallback = nullptr,
      std::function<bool()> stopCallback = nullptr
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
  vtkSlicerKMAPLogic();
  ~vtkSlicerKMAPLogic() override;

  void SetMRMLSceneInternal(vtkMRMLScene* newScene) override;
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  void RegisterNodes() override;
  void UpdateFromMRMLScene() override;
  void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;
  void OnMRMLSceneNodeRemoved(vtkMRMLNode* node) override;
private:

  vtkSlicerKMAPLogic(const vtkSlicerKMAPLogic&); // Not implemented
  void operator=(const vtkSlicerKMAPLogic&); // Not implemented
};

#endif
