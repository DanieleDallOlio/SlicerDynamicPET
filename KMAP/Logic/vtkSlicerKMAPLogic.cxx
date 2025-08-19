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

// KMAP Logic includes
#include "vtkSlicerKMAPLogic.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// STD includes
#include <cassert>

// for linear regression
#include <Eigen/Dense>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerKMAPLogic);

//----------------------------------------------------------------------------
vtkSlicerKMAPLogic::vtkSlicerKMAPLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerKMAPLogic::~vtkSlicerKMAPLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerKMAPLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerKMAPLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerKMAPLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}


void vtkSlicerKMAPLogic::computeTAC(vtkIdType ctID,
                                    vtkIdType petID,
                                    vtkIdType segID,
                                    std::vector<QString> segmentsID,
                                    std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs,
                                    std::map<std::string, std::string>& segmentTACsnames,
                                    QProgressBar* ProgressBar
                                  )
{
  vtkMRMLScene* scene = this->GetMRMLScene();
  if (scene==nullptr) {
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = vtkMRMLSubjectHierarchyNode::GetSubjectHierarchyNode(scene);
  if (!shNode) {
    return;
  }
  // Fetch CT
  vtkMRMLScalarVolumeNode* ctNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(ctID));
  if (!ctNode) {
    return;
  }
  // std :: cout << ctID << std :: endl;
  // std :: string name = shNode->GetItemName(ctID);
  // std::cout << "CT: " << name
  //           << ", ID: " << ctID << std::endl;
  // Fetch PET
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(petID));
  if (!petNode) {
    return;
  }
  // name = shNode->GetItemName(petID);
  // std::cout << "PET: " << name
  //           << ", ID: " << petID << std::endl;
  // Fetch Segmentation
  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(segID));
  if (!segNode) {
    return;
  }
  // name = shNode->GetItemName(segID);
  // std::cout << "SEG: " << name
  //           << ", ID: " << segID << std::endl;
  // Make sure the source of the segmentation is a binary label map, alongside a created closed surface
  setupSeg(segNode);

  // std::cout << "Selected segments:" << std::endl;
  // for (const QString& id : segmentsID)
  // {
  //   std::string segmentName = seg->GetSegment(id.toStdString())->GetName();
  //   std::cout << segmentName << " - " << id.toStdString() << std::endl;
  // }

  // Collect the sequence for the dynamic PET
  vtkMRMLSequenceNode* sequencePETNode = nullptr;
  vtkMRMLSequenceBrowserNode* sequenceBrowserPETNode = nullptr;
  for (int i = 0; i < scene->GetNumberOfNodesByClass("vtkMRMLSequenceBrowserNode"); ++i)
  {
    vtkMRMLSequenceBrowserNode* browser = vtkMRMLSequenceBrowserNode::SafeDownCast(scene->GetNthNodeByClass(i, "vtkMRMLSequenceBrowserNode"));
    if (!browser)
      continue;

    // Check if this browser is using our PET node as a proxy node
    vtkMRMLSequenceNode* seqNode = browser->GetSequenceNode(petNode);
    if (seqNode)
    {
      sequencePETNode = seqNode;
      sequenceBrowserPETNode = browser;
      break;
    }
  }
  if (!sequencePETNode || !sequenceBrowserPETNode)
  {
    std::cerr << "Could not find sequence or browser node for PET." << std::endl;
    return;
  }
  int numberOfTimepoints = sequencePETNode->GetNumberOfDataNodes();

  // Setup the right sequence for the segmentation
  vtkMRMLSequenceNode* segSequenceNode = sequenceBrowserPETNode->GetSequenceNode(segNode);
  if (!segSequenceNode)
  {
    segSequenceNode = vtkMRMLSequenceNode::New();
    segSequenceNode->SetName(shNode->GetItemName(segID).c_str());
    scene->AddNode(segSequenceNode);
    sequenceBrowserPETNode->AddProxyNode(segNode, segSequenceNode, false);
  }
  std::string indexValue;
  sequenceBrowserPETNode->SetSaveChanges(segSequenceNode, true);
  for (int i = 0; i < numberOfTimepoints; ++i) {
    indexValue = sequencePETNode->GetNthIndexValue(i);
    if (!segSequenceNode->GetDataNodeAtValue(indexValue))
    {
      segSequenceNode->SetDataNodeAtValue(segNode, indexValue);
    }
  }

  // Get TAC
  this->TAC(sequencePETNode, segSequenceNode, segmentsID, segmentTACs, segmentTACsnames, ProgressBar);


  // Fetch segments
  // vtkSegmentation* seg = segNode->GetSegmentation();
  // if (!seg) {
  //   return;
  // }
  // if (segmentsID.empty()) {
  //   return;
  // }


}



void vtkSlicerKMAPLogic::setupSeg(vtkMRMLSegmentationNode* segNode)
{
  if (!segNode)
  {
    std::cerr << "setupSeg: Segmentation node is null!" << std::endl;
    return;
  }

  vtkSegmentation* segmentation = segNode->GetSegmentation();
  if (!segmentation)
  {
    std::cerr << "setupSeg: Failed to get vtkSegmentation from node!" << std::endl;
    return;
  }

  // 1. Make sure "Binary labelmap" is available as a representation
  const std::string labelmapRep = vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();
  segmentation->CreateRepresentation(labelmapRep);

  // 2. Set "Binary labelmap" as the master representation (i.e., source)
  segmentation->SetSourceRepresentationName(labelmapRep);

  // 3. Ensure "Closed surface" representation is present
  const std::string closedSurfRep = vtkSegmentationConverter::GetSegmentationClosedSurfaceRepresentationName();
  segmentation->CreateRepresentation(closedSurfRep);

  // std::cout << "setupSeg: Segmentation representations ensured." << std::endl;
}


VoxelStatistics vtkSlicerKMAPLogic::ComputeVoxelStatistics(vtkImageData* petImage, vtkImageData* labelmap, int labelValue)
{
  VoxelStatistics stats;
  std::vector<double> values;

  int dims[3];
  double spacing[3];
  petImage->GetDimensions(dims);
  petImage->GetSpacing(spacing);

  double voxelVolume = spacing[0] * spacing[1] * spacing[2];

  vtkDataArray* petArray = petImage->GetPointData()->GetScalars();
  vtkDataArray* labelArray = labelmap->GetPointData()->GetScalars();

  for (int z = 0; z < dims[2]; ++z)
  {
    for (int y = 0; y < dims[1]; ++y)
    {
      for (int x = 0; x < dims[0]; ++x)
      {
        int ijk[3] = { x, y, z };
        vtkIdType idx = petImage->ComputePointId(ijk);
        int label = static_cast<int>(labelArray->GetComponent(idx, 0));
        if (label == labelValue)
        {
          double val = petArray->GetComponent(idx, 0);
          stats.count++;
          stats.mean += val;
          stats.min = std::min(stats.min, val);
          stats.max = std::max(stats.max, val);
          values.push_back(val);
        }
      }
    }
  }

  if (stats.count > 0)
  {
    stats.mean /= stats.count;
    stats.volume_mm3 = voxelVolume * stats.count;
    stats.volume_cm3 = stats.volume_mm3 / 1000.0;

    double variance = 0.0;
    for (double val : values)
      variance += (val - stats.mean) * (val - stats.mean);
    stats.stddev = std::sqrt(variance / stats.count);

    std::sort(values.begin(), values.end());
    int n = stats.count;
    stats.median = (n % 2 == 0) ? (values[n / 2 - 1] + values[n / 2]) / 2.0 : values[n / 2];
    stats.q1 = values[n / 4];
    stats.q3 = values[3 * n / 4];
    stats.iqr = stats.q3 - stats.q1;
  }

  return stats;
}

void vtkSlicerKMAPLogic::TAC(vtkMRMLSequenceNode* sequencePETNode,
                             vtkMRMLSequenceNode* segSequenceNode,
                             std::vector<QString> segmentsID,
                             std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs,
                             std::map<std::string, std::string>& segmentTACsnames,
                             QProgressBar* ProgressBar
                         )
{

  if (!sequencePETNode || !segSequenceNode)
  {
    std::cerr << "Invalid input nodes!" << std::endl;
    return;
  }

  int numberOfTimepoints = sequencePETNode->GetNumberOfDataNodes();
  if (numberOfTimepoints == 0)
  {
    std::cerr << "Empty sequence node!" << std::endl;
    return;
  }

  std::string index0 = sequencePETNode->GetNthIndexValue(0);
  vtkMRMLSegmentationNode* segmentationAt0 = vtkMRMLSegmentationNode::SafeDownCast(segSequenceNode->GetDataNodeAtValue(index0));
  if (!segmentationAt0)
  {
    std::cerr << "First segmentation node is invalid!" << std::endl;
    return;
  }

  for (const QString& id : segmentsID)
  {
    const std::string& segmentID = id.toStdString();
    auto* segment = segmentationAt0->GetSegmentation()->GetSegment(segmentID);
    if (!segment)
    {
      std::cerr << "Segment not found for ID: " << id.toStdString() << std::endl;
      continue;
    }
    segmentTACsnames[segmentID] = segment->GetName();
    segmentTACs[segmentID] = std::vector<VoxelStatistics>();
    segmentTACs[segmentID].reserve(numberOfTimepoints);
  }

  for (int i = 0; i < numberOfTimepoints; ++i)
  {
    std::string indexValue = sequencePETNode->GetNthIndexValue(i);
    auto* PETVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
      sequencePETNode->GetDataNodeAtValue(indexValue));
    auto* segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(
      segSequenceNode->GetDataNodeAtValue(indexValue));
    if (!PETVolume || !segmentationNode)
    {
      std::cerr << "Missing data for timepoint " << i << std::endl;
      continue;
    }

    for (size_t s = 0; s < segmentsID.size(); ++s)
    {
      const std::string& segmentID = segmentsID[s].toStdString();
      const std::string& segmentName = segmentTACsnames[segmentID];

      vtkNew<vtkStringArray> segmentArray;
      segmentArray->InsertNextValue(segmentID);

      vtkSmartPointer<vtkOrientedImageData> labelmap = vtkSmartPointer<vtkOrientedImageData>::New();
      vtkSlicerSegmentationsModuleLogic::GenerateMergedLabelmapInReferenceGeometry(segmentationNode,
                                                                                   PETVolume,
                                                                                   segmentArray,
                                                                                   vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
                                                                                   labelmap);
      if (!labelmap)
      {
        std::cerr << "Failed to generate labelmap for segment: " << segmentID << " at timepoint " << i << std::endl;
        segmentTACs[segmentName].emplace_back(VoxelStatistics{});  // Insert empty stats
        continue;
      }
      VoxelStatistics stats = ComputeVoxelStatistics(
        PETVolume->GetImageData(), labelmap, 1);

      segmentTACs[segmentID].emplace_back(stats);
    }
    ProgressBar->setValue(static_cast<double>(i + 1) / numberOfTimepoints*100.);
    qApp->processEvents();
  }

}


void vtkSlicerKMAPLogic::callTCM(std :: vector< std :: vector<double> > tac,
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
                                 double *& fitted_curve,
                                 const std::vector<double>* wgt
                                 )
{
  const int nth = 1;

  // Basic validation
  if (containsNaN(tac)) return error_nan("TAC");
  if (containsNaN(Cp)) return error_nan("Cp");
  if (containsNaN(framing)) return error_nan("framing");

  if (tac.size() != framing.size()) return error_size("TAC", "framing", tac.size(), framing.size());
  if (tac.size() != Cp.size()) return error_size("TAC", "Cp", tac.size(), Cp.size());

  // Allocate weights
  double* wt = new double[Nframe];
  if (wgt == nullptr)
  {
    std::fill(wt, wt + Nframe, 1.0);
  } else {
    if (wgt->size() != static_cast<size_t>(Nframe))
    {
      throw std::runtime_error(
          "Weight vector length (" + std::to_string(wgt->size()) +
          ") does not match Nframe (" + std::to_string(Nframe) + ")."
      );
    }
    std::copy(wgt->begin(), wgt->end(), wt);
  }

  // Cumulative sum
  double* cumsum = new double[Nframe];
  double cum = 0.0;
  for (long int i = 0; i < Nframe; ++i)
  {
    cum += framing[i][0];
    cumsum[i] = cum;
  }

  double **scant = new double * [Nframe];
  double t, pbr;
  std :: vector< std :: vector<double> > cwb; // whole blood concentration
  scant[0L] = new double[2];
  scant[0L][0L] = 0.;
  scant[0L][1L] = cumsum[0L];
  t = (scant[0L][0L] + scant[0L][1L]) * 0.5;
  pbr = pbrp[0L] * exp(-pbrp[1L] * t / 60L) + pbrp[2L];
  std :: vector<double> cwb_r;
  cwb_r.push_back(Cp[0L][0L] / pbr);
  cwb.push_back(cwb_r);
  for (long int i = 1L; i < Nframe; ++i){
    cwb_r.clear();
    scant[i] = new double[2];
    scant[i][0L] = cumsum[i-1];
    scant[i][1L] = cumsum[i];
    t = (scant[i][0L] + scant[i][1L]) * 0.5;
    pbr = pbrp[0L] * exp(-pbrp[1L] * t / 60L) + pbrp[2L];
    cwb_r.push_back(Cp[i][0L] / pbr);
    cwb.push_back(cwb_r);
  }

  long int N_cp;
  double *Cp_new = finesample(scant, Cp, Nframe, N_cp, timestep, "linear");
  double *cwb_new = finesample(scant, cwb, Nframe, N_cp, timestep, "linear");

  double * tac_flatten = new double[Nframe*Nvox];
  for (int i=0; i<Nframe; ++i) {
    for (int j=0; j<Nvox; ++j) {
      tac_flatten[i + j * Nframe] = tac[i][j];
    }
  }
  double * scant_flatten = new double[Nframe*2];
  for (int i=0; i<Nframe; ++i) {
    for (int j=0; j<2; ++j) {
      scant_flatten[i + j * Nframe] = scant[i][j];
    }
  }

  fitted_curve   = new double[Nframe*Nvox];
  if (n_tc == 1) {
    double *fitted_params = new double[4*Nvox];
    kfit_1tcm_mex_omp(tac_flatten,
                      Nframe,
                      Nvox,
                      Nvox,
                      wt,
                      scant_flatten,
                      Cp_new,
                      cwb_new,
                      dk,
                      kinit,
                      4,
                      1,
                      lb,
                      ub,
                      sens,
                      maxiter,
                      timestep,
                      nth,
                      fitted_params,
                      fitted_curve
                      );
    params.vb = fitted_params[0];
    params.K1 = fitted_params[1];
    params.k2 = fitted_params[2];
    params.td = fitted_params[3];
    params.Ki = params.K1;
    params.DV = params.K1/(params.k2 + 1e-16);
    // std :: cout << "fitted_params" << std ::endl;
    // print_vec(fitted_params, 4*Nvox);
    // std :: cout << "fitted_curve" << std ::endl;
    // print_vec(fitted_curve, Nframe*Nvox);
  } else if (n_tc == 2) {
    double *fitted_params = new double[6*Nvox];
    kfit_2tcm_mex_omp(tac_flatten,
                      Nframe,
                      Nvox,
                      Nvox,
                      wt,
                      scant_flatten,
                      Cp_new,
                      cwb_new,
                      dk,
                      kinit,
                      6,
                      1,
                      lb,
                      ub,
                      sens,
                      maxiter,
                      timestep,
                      nth,
                      fitted_params,
                      fitted_curve
                      );
    params.vb = fitted_params[0];
    params.K1 = fitted_params[1];
    params.k2 = fitted_params[2];
    params.k3 = fitted_params[3];
    params.k4 = fitted_params[4];
    params.td = fitted_params[5];
    params.Ki = params.K1 * params.k3 / (params.k2 + params.k3);
    params.DV = params.K1/(params.k2 + 1e-16) * (1 + params.k3/(params.k4 + 1e-16));
    // std :: cout << "fitted_params" << std ::endl;
    // print_vec(fitted_params, 6*Nvox);
    // std :: cout << "fitted_curve" << std ::endl;
    // print_vec(fitted_curve, Nframe*Nvox);

  } else {
    std :: cerr << "Forbidden number of tissue compartment: " << n_tc << std :: endl;
  }

  // --- Cleanup ---
  // delete[] wt;
  // delete[] cumsum;
  // delete[] Cp_interp;
  // delete[] cwb_interp;
  // delete[] tac_flatten;
  // delete[] scant_flatten;
  // delete[] fitted_curve;
  // delete[] fitted_params;
  delete[] wt;
  return;
}



// void vtkSlicerKMAPLogic::getFittedTCM(double *& fitted_curve,
//                                       std :: vector< std :: vector<double> > Cp,
//                                       std :: vector< std :: vector<double> > framing,
//                                       long int Nframe,
//                                       long int Nvox,
//                                       double* kinit,
//                                       double* lb,
//                                       double* ub,
//                                       const bool* sens,
//                                       const double dk,
//                                       const double timestep,
//                                       const double pbrp[],
//                                       const int maxiter,
//                                       const int n_tc,
//                                       TCMParameters& params
//                                       )
// {
//   fitted_curve   = new double[Nframe*Nvox];
//   const int nth = 1;
//
//   // Basic validation
//   if (containsNaN(Cp)) return error_nan("Cp");
//   if (containsNaN(framing)) return error_nan("framing");
//
//
//   // Allocate weights
//   double* wt = new double[Nframe];
//   std::fill(wt, wt + Nframe, 1.0);
//
//   // Cumulative sum
//   double* cumsum = new double[Nframe];
//   double cum = 0.0;
//   for (long int i = 0; i < Nframe; ++i)
//   {
//     cum += framing[i][0];
//     cumsum[i] = cum;
//   }
//
//   double **scant = new double * [Nframe];
//   double t, pbr;
//   std :: vector< std :: vector<double> > cwb; // whole blood concentration
//   scant[0L] = new double[2];
//   scant[0L][0L] = 0.;
//   scant[0L][1L] = cumsum[0L];
//   t = (scant[0L][0L] + scant[0L][1L]) * 0.5;
//   pbr = pbrp[0L] * exp(-pbrp[1L] * t / 60L) + pbrp[2L];
//   std :: vector<double> cwb_r;
//   cwb_r.push_back(Cp[0L][0L] / pbr);
//   cwb.push_back(cwb_r);
//   for (long int i = 1L; i < Nframe; ++i){
//     cwb_r.clear();
//     scant[i] = new double[2];
//     scant[i][0L] = cumsum[i-1];
//     scant[i][1L] = cumsum[i];
//     t = (scant[i][0L] + scant[i][1L]) * 0.5;
//     pbr = pbrp[0L] * exp(-pbrp[1L] * t / 60L) + pbrp[2L];
//     cwb_r.push_back(Cp[i][0L] / pbr);
//     cwb.push_back(cwb_r);
//   }
//
//   long int N_cp;
//   double *Cp_new = finesample(scant, Cp, Nframe, N_cp, timestep, "linear");
//   double *cwb_new = finesample(scant, cwb, Nframe, N_cp, timestep, "linear");
//
//   double * scant_flatten = new double[Nframe*2];
//   for (int i=0; i<Nframe; ++i) {
//     for (int j=0; j<2; ++j) {
//       scant_flatten[i + j * Nframe] = scant[i][j];
//     }
//   }
//
//   if (n_tc == 1) {
//     double *fitted_params = new double[4*Nvox];
//     fitted_params[0] = params.vb;
//     fitted_params[1] = params.K1;
//     fitted_params[2] = params.k2;
//     fitted_params[3] = params.td;
//     kconv_1tcm_tac(fitted_params, dk, scant_flatten, timestep, Cp_new,
//                    cwb_new, Nframe, Nvox, fitted_curve);
//   } else if (n_tc == 2) {
//     double *fitted_params = new double[6*Nvox];
//     fitted_params[0] = params.vb;
//     fitted_params[1] = params.K1;
//     fitted_params[2] = params.k2;
//     fitted_params[3] = params.k3;
//     fitted_params[4] = params.k4;
//     fitted_params[5] = params.td;
//     kconv_2tcm_tac(fitted_params, dk, scant_flatten, timestep, Cp_new,
//                    cwb_new, Nframe, Nvox, fitted_curve);
//   } else {
//     std :: cerr << "Forbidden number of tissue compartment: " << n_tc << std :: endl;
//   }
//
//   // --- Cleanup ---
//   // delete[] wt;
//   // delete[] cumsum;
//   // delete[] Cp_interp;
//   // delete[] cwb_interp;
//   // delete[] scant_flatten;
//   // delete[] fitted_curve;
//   // delete[] fitted_params;
//   return;
// }

double vtkSlicerKMAPLogic::computeAIC(const std::vector<double>& obs,
                                      const std::vector<double>& est,
                                      int numpar,
                                      const std::vector<double>* wgt,
                                      bool aicc)
{
    size_t numfrm = obs.size();
    if (est.size() != numfrm)
    {
        throw std::invalid_argument("obs and est must have same length");
    }

    // Default weights = 1
    std::vector<double> weights;
    if (wgt == nullptr)
    {
        weights.assign(numfrm, 1.0);
    }
    else
    {
        if (wgt->size() != numfrm)
            throw std::invalid_argument("weights size must match obs");
        weights = *wgt;
    }

    // Sum of squared errors with weights
    double ss = 0.0;
    for (size_t i = 0; i < numfrm; ++i)
    {
        double diff = obs[i] - est[i];
        ss += weights[i] * diff * diff;
    }

    int numpar_new = numpar + 1; // intercept term
    double AIC = numfrm * std::log(ss / static_cast<double>(numfrm))
                 + 2.0 * numpar_new;

    // Apply AICc correction if needed
    if (aicc && (static_cast<double>(numfrm) / numpar < 40.0))
    {
        AIC += 2.0 * numpar_new * (numpar_new + 1)
               / (static_cast<double>(numfrm) - numpar_new - 1.0);
    }

    return AIC;
}

// Compute MASE
double vtkSlicerKMAPLogic::MASE(const std::vector<double>& Actual,
                                const std::vector<double>& Predicted,
                                const std::vector<double>* wgt)
{
    size_t N = Actual.size();
    if (Predicted.size() != N)
    {
        throw std::invalid_argument("Actual and Predicted must have same length");
    }
    if (N < 2)
    {
        throw std::invalid_argument("Need at least 2 observations for MASE");
    }

    // Default weights = 1
    std::vector<double> weights;
    if (wgt == nullptr)
    {
        weights.assign(N, 1.0);
    }
    else
    {
        if (wgt->size() != N)
            throw std::invalid_argument("weights size must match Actual");
        weights = *wgt;
    }

    // Numerator: weighted mean absolute error
    double num_sum = 0.0;
    double num_wsum = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        num_sum += weights[i] * std::abs(Actual[i] - Predicted[i]);
        num_wsum += weights[i];
    }
    double num = num_sum / num_wsum;

    // Denominator: weighted mean absolute difference of Actual
    double den_sum = 0.0;
    double den_wsum = 0.0;
    for (size_t i = 1; i < N; ++i)
    {
        den_sum += weights[i] * std::abs(Actual[i] - Actual[i - 1]);
        den_wsum += weights[i];
    }
    double den = den_sum / den_wsum;

    return num / den;
}


void vtkSlicerKMAPLogic::Patlak(const std::vector<double>& tac,
                                const std::vector<double>& Cp,
                                const std::vector<double>& framing,
                                MTGAParameters & params,
                                const std::vector<double>* wgt, // nullptr if not using weights
                                const double timeOffset,
                                const double framingNorm,
                                bool robust,
                                bool std,
                                double huber_tune,
                                double tol,
                                int max_iter
                              )
{
  size_t N = tac.size();
  std::vector<double> outX, outY, fittedValues;
  std::vector<int> outframe;

  // normalize framing to minutes
  std::vector<double> frameScaled(N);
  for (size_t i = 0; i < N; ++i)
      frameScaled[i] = framing[i] / framingNorm;

  // cumulative sum of framing
  std::vector<double> timeAlong(N);
  std::partial_sum(frameScaled.begin(), frameScaled.end(), timeAlong.begin());

  // cumulative trapezoid for Cp
  std::vector<double> intCp(N, 0.0);
  double acc = 0.0;
  for (size_t i = 1; i < N; ++i)
  {
      acc += 0.5 * (Cp[i] + Cp[i-1]) * frameScaled[i];
      intCp[i] = acc;
  }

  // build X, Y with time filter
  std::vector<double> wgt_adj;
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset)
      {
          double x = intCp[i] / (Cp[i] + 1e-16);
          double y = tac[i]   / (Cp[i] + 1e-16);
          outX.push_back(x);
          outY.push_back(y);
          outframe.push_back(i+1);
          if (wgt) wgt_adj.push_back((*wgt)[i]);
      }
  }

  size_t n = outX.size();
  // optional standardization
  std::vector<double> x_data_tmp, y_data_tmp;
  double meanX = 0.0, meanY = 0.0, stdX = 1.0, stdY = 1.0;
  if (std)
  {
      x_data_tmp = outX;
      y_data_tmp = outY;

      meanX = std::accumulate(outX.begin(), outX.end(), 0.0) / n;
      meanY = std::accumulate(outY.begin(), outY.end(), 0.0) / n;
      stdX = std::sqrt(std::inner_product(outX.begin(), outX.end(), outX.begin(), 0.0) / n - meanX*meanX);
      stdY = std::sqrt(std::inner_product(outY.begin(), outY.end(), outY.begin(), 0.0) / n - meanY*meanY);

      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = (outX[i] - meanX) / stdX;
          outY[i] = (outY[i] - meanY) / stdY;
      }
  }

  // regression (robust or normal OLS)
  Eigen::Map<const Eigen::VectorXd> Xv(outX.data(), n);
  Eigen::Map<const Eigen::VectorXd> Yv(outY.data(), n);
  Eigen::MatrixXd A(n, 2);
  A.col(0) = Eigen::VectorXd::Ones(n);
  A.col(1) = Xv;

  Eigen::VectorXd coeff(2);
  if (robust) {
    // simple robust regression using Iteratively Reweighted Least Squares (Huber)
    Eigen::VectorXd weights = Eigen::VectorXd::Ones(n);
    Eigen::MatrixXd W = Eigen::MatrixXd::Identity(n, n);
    Eigen::VectorXd residuals(n), prev_coeff(2);
    for (int iter = 0; iter < max_iter; ++iter) {
      // Update weight matrix
      W.diagonal() = weights;
      // Weighted least squares step
      coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
      // Check for convergence
      if (iter > 0 && (coeff - prev_coeff).norm() < tol) break;
      prev_coeff = coeff;
      // Update residuals and weights
      residuals = Yv - A * coeff;
      for (int i = 0; i < n; ++i) {
          double r = std::abs(residuals(i));
          weights(i) = (r <= huber_tune) ? 1.0 : huber_tune / std::max(r, 1e-8);  // avoid div by 0
      }
    }
    wgt_adj.assign(weights.data(), weights.data() + weights.size());
    wgt = &wgt_adj;
  } else {
    if (!wgt)
    {
      // Ordinary Least Squares
      coeff = (A.transpose() * A).ldlt().solve(A.transpose() * Yv);
    }
    else
    {
      // Weighted Least Squares
      Eigen::Map<const Eigen::VectorXd> Wv(wgt_adj.data(), n);
      Eigen::MatrixXd W = Wv.asDiagonal();
      coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
    }
  }

  double intercept = coeff(0);
  double slope = coeff(1);

  // fitted values
  fittedValues.resize(n);
  for (size_t i = 0; i < n; ++i)
      fittedValues[i] = intercept + slope * outX[i];

  // de-standardize if needed
  if (std)
  {
      double devyoverdevx = stdY / stdX;
      slope   = slope * devyoverdevx;
      intercept = meanY - slope * meanX + intercept;
      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = outX[i]*stdX + meanX;
          outY[i] = outY[i]*stdY + meanY;
          fittedValues[i] = fittedValues[i]*stdY + meanY;
      }
  }

  // Ki and Intercept
  params.Ki = slope;
  params.Intercept = intercept;

  // x, y and fitted values
  params.x = outX;
  params.y = outY;
  params.frame = outframe;
  params.fitted = fittedValues;

  // AIC and MASE
  params.AIC = computeAIC(outY, fittedValues, 2, wgt ? &wgt_adj : nullptr);
  params.MASE = MASE(outY, fittedValues, wgt ? &wgt_adj : nullptr);
}

void vtkSlicerKMAPLogic::Logan(const std::vector<double>& tac,
                               const std::vector<double>& Cp,
                               const std::vector<double>& framing,
                               MTGAParameters & params,
                               const std::vector<double>* wgt, // nullptr if not using weights
                               const double timeOffset,
                               const double framingNorm,
                               bool robust,
                               bool std,
                               double huber_tune,
                               double tol,
                               int max_iter
)
{
  size_t N = tac.size();
  std::vector<double> outX, outY, fittedValues;
  std::vector<int> outframe;

  // normalize framing to minutes
  std::vector<double> frameScaled(N);
  for (size_t i = 0; i < N; ++i)
      frameScaled[i] = framing[i] / framingNorm;

  // cumulative sum of framing
  std::vector<double> timeAlong(N);
  std::partial_sum(frameScaled.begin(), frameScaled.end(), timeAlong.begin());


  // cumulative trapezoid for Cp
  std::vector<double> intCp(N, 0.0);
  double acc = 0.0;
  for (size_t i = 1; i < N; ++i)
  {
      acc += 0.5 * (Cp[i] + Cp[i-1]) * frameScaled[i];
      intCp[i] = acc;
  }

  // cumulative trapezoid for TAC
  std::vector<double> intCt(N, 0.0);
  acc = 0.0;
  for (size_t i = 1; i < N; ++i)
  {
      acc += 0.5 * (tac[i] + tac[i-1]) * frameScaled[i];
      intCt[i] = acc;
  }

  // Build X, Y with time filter
  std::vector<double> wgt_adj;
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset)
      {
          double x = intCp[i] / (tac[i] + 1e-16);
          double y = intCt[i] / (tac[i] + 1e-16);
          outX.push_back(x);
          outY.push_back(y);
          outframe.push_back(i+1);
          if (wgt) wgt_adj.push_back((*wgt)[i]);
      }
  }

  size_t n = outX.size();
  // Optional standardization
  std::vector<double> x_data_tmp, y_data_tmp;
  double meanX = 0.0, meanY = 0.0, stdX = 1.0, stdY = 1.0;
  if (std)
  {
      x_data_tmp = outX;
      y_data_tmp = outY;

      meanX = std::accumulate(outX.begin(), outX.end(), 0.0) / n;
      meanY = std::accumulate(outY.begin(), outY.end(), 0.0) / n;
      stdX = std::sqrt(std::inner_product(outX.begin(), outX.end(), outX.begin(), 0.0) / n - meanX * meanX);
      stdY = std::sqrt(std::inner_product(outY.begin(), outY.end(), outY.begin(), 0.0) / n - meanY * meanY);

      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = (outX[i] - meanX) / stdX;
          outY[i] = (outY[i] - meanY) / stdY;
      }
  }

  // Regression (robust or normal OLS)
  Eigen::Map<const Eigen::VectorXd> Xv(outX.data(), n);
  Eigen::Map<const Eigen::VectorXd> Yv(outY.data(), n);
  Eigen::MatrixXd A(n, 2);
  A.col(0) = Eigen::VectorXd::Ones(n);
  A.col(1) = Xv;

  Eigen::VectorXd coeff(2);
  if (robust)
  {
      // Iteratively Reweighted Least Squares (Huber)
      Eigen::VectorXd weights = Eigen::VectorXd::Ones(n);
      Eigen::MatrixXd W = Eigen::MatrixXd::Identity(n, n);
      Eigen::VectorXd residuals(n), prev_coeff(2);
      for (int iter = 0; iter < max_iter; ++iter)
      {
          W.diagonal() = weights;
          coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
          if (iter > 0 && (coeff - prev_coeff).norm() < tol) break;
          prev_coeff = coeff;
          residuals = Yv - A * coeff;
          for (int i = 0; i < n; ++i)
          {
              double r = std::abs(residuals(i));
              weights(i) = (r <= huber_tune) ? 1.0 : huber_tune / std::max(r, 1e-8);
          }
      }
      wgt_adj.assign(weights.data(), weights.data() + weights.size());
      wgt = &wgt_adj;
  }
  else
  {
      if (!wgt)
      {
          // Ordinary Least Squares
          coeff = (A.transpose() * A).ldlt().solve(A.transpose() * Yv);
      }
      else
      {
          // Weighted Least Squares
          Eigen::Map<const Eigen::VectorXd> Wv(wgt_adj.data(), n);
          Eigen::MatrixXd W = Wv.asDiagonal();
          coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
      }
  }

  double intercept = coeff(0);
  double slope = coeff(1);

  // Fitted values
  fittedValues.resize(n);
  for (size_t i = 0; i < n; ++i)
      fittedValues[i] = intercept + slope * outX[i];

  // De-standardize if needed
  if (std)
  {
      double devyoverdevx = stdY / stdX;
      slope   = slope * devyoverdevx;
      intercept = meanY - slope * meanX + intercept;
      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = outX[i]*stdX + meanX;
          outY[i] = outY[i]*stdY + meanY;
          fittedValues[i] = fittedValues[i]*stdY + meanY;
      }
  }

  // DV and Intercept
  params.DV = slope;
  params.Intercept = intercept;

  // x, y and fitted values
  params.x = outX;
  params.y = outY;
  params.frame = outframe;
  params.fitted = fittedValues;

  // AIC and MASE
  params.AIC = computeAIC(outY, fittedValues, 2, wgt ? &wgt_adj : nullptr);
  params.MASE = MASE(outY, fittedValues, wgt ? &wgt_adj : nullptr);
}

void vtkSlicerKMAPLogic::RE(const std::vector<double>& tac,
                            const std::vector<double>& Cp,
                            const std::vector<double>& framing,
                            MTGAParameters & params,
                            const std::vector<double>* wgt, // nullptr if not using weights
                            const double timeOffset,
                            const double framingNorm,
                            bool robust,
                            bool std,
                            double huber_tune,
                            double tol,
                            int max_iter
)
{
  size_t N = tac.size();
  std::vector<double> outX, outY, fittedValues;
  std::vector<int> outframe;

  // normalize framing to minutes
  std::vector<double> frameScaled(N);
  for (size_t i = 0; i < N; ++i)
      frameScaled[i] = framing[i] / framingNorm;

  // cumulative sum of framing
  std::vector<double> timeAlong(N);
  std::partial_sum(frameScaled.begin(), frameScaled.end(), timeAlong.begin());

  // cumulative trapezoid for Cp
  std::vector<double> intCp(N, 0.0);
  double acc = 0.0;
  for (size_t i = 1; i < N; ++i)
  {
      acc += 0.5 * (Cp[i] + Cp[i-1]) * frameScaled[i];
      intCp[i] = acc;
  }

  // cumulative trapezoid for TAC
  std::vector<double> intCt(N, 0.0);
  acc = 0.0;
  for (size_t i = 1; i < N; ++i)
  {
      acc += 0.5 * (tac[i] + tac[i-1]) * frameScaled[i];
      intCt[i] = acc;
  }

  // Build X, Y with time filter
  std::vector<double> wgt_adj;
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset)
      {
          double x = intCp[i] / (Cp[i] + 1e-16);
          double y = intCt[i] / (Cp[i] + 1e-16);
          outX.push_back(x);
          outY.push_back(y);
          outframe.push_back(i+1);
          if (wgt) wgt_adj.push_back((*wgt)[i]);
      }
  }

  size_t n = outX.size();
  // Optional standardization
  std::vector<double> x_data_tmp, y_data_tmp;
  double meanX = 0.0, meanY = 0.0, stdX = 1.0, stdY = 1.0;
  if (std)
  {
      x_data_tmp = outX;
      y_data_tmp = outY;

      meanX = std::accumulate(outX.begin(), outX.end(), 0.0) / n;
      meanY = std::accumulate(outY.begin(), outY.end(), 0.0) / n;
      stdX = std::sqrt(std::inner_product(outX.begin(), outX.end(), outX.begin(), 0.0) / n - meanX * meanX);
      stdY = std::sqrt(std::inner_product(outY.begin(), outY.end(), outY.begin(), 0.0) / n - meanY * meanY);

      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = (outX[i] - meanX) / stdX;
          outY[i] = (outY[i] - meanY) / stdY;
      }
  }

  // Regression (robust or normal OLS)
  Eigen::Map<const Eigen::VectorXd> Xv(outX.data(), n);
  Eigen::Map<const Eigen::VectorXd> Yv(outY.data(), n);
  Eigen::MatrixXd A(n, 2);
  A.col(0) = Eigen::VectorXd::Ones(n);
  A.col(1) = Xv;

  Eigen::VectorXd coeff(2);
  if (robust)
  {
      // Iteratively Reweighted Least Squares (Huber)
      Eigen::VectorXd weights = Eigen::VectorXd::Ones(n);
      Eigen::MatrixXd W = Eigen::MatrixXd::Identity(n, n);
      Eigen::VectorXd residuals(n), prev_coeff(2);
      for (int iter = 0; iter < max_iter; ++iter)
      {
          W.diagonal() = weights;
          coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
          if (iter > 0 && (coeff - prev_coeff).norm() < tol) break;
          prev_coeff = coeff;
          residuals = Yv - A * coeff;
          for (int i = 0; i < n; ++i)
          {
              double r = std::abs(residuals(i));
              weights(i) = (r <= huber_tune) ? 1.0 : huber_tune / std::max(r, 1e-8);
          }
      }
      wgt_adj.assign(weights.data(), weights.data() + weights.size());
      wgt = &wgt_adj;
  }
  else
  {
      if (!wgt)
      {
          // Ordinary Least Squares
          coeff = (A.transpose() * A).ldlt().solve(A.transpose() * Yv);
      }
      else
      {
          // Weighted Least Squares
          Eigen::Map<const Eigen::VectorXd> Wv(wgt_adj.data(), n);
          Eigen::MatrixXd W = Wv.asDiagonal();
          coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
      }
  }

  double intercept = coeff(0);
  double slope = coeff(1);

  // Fitted values
  fittedValues.resize(n);
  for (size_t i = 0; i < n; ++i)
      fittedValues[i] = intercept + slope * outX[i];

  // De-standardize if needed
  if (std)
  {
      double devyoverdevx = stdY / stdX;
      slope   = slope * devyoverdevx;
      intercept = meanY - slope * meanX + intercept;
      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = outX[i]*stdX + meanX;
          outY[i] = outY[i]*stdY + meanY;
          fittedValues[i] = fittedValues[i]*stdY + meanY;
      }
  }

  // DV and Intercept
  params.DV = slope;
  params.Intercept = intercept;

  // x, y and fitted values
  params.x = outX;
  params.y = outY;
  params.frame = outframe;
  params.fitted = fittedValues;

  // AIC and MASE
  params.AIC = computeAIC(outY, fittedValues, 2, wgt ? &wgt_adj : nullptr);
  params.MASE = MASE(outY, fittedValues, wgt ? &wgt_adj : nullptr);
}
