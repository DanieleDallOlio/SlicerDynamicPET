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
#include <chrono>
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
#include <QThread>

// define M_PI in case of Win
#ifdef _WIN32
    #ifndef M_PI
        #define M_PI 3.14159265358979323846
    #endif
#endif


static const std::map<std::string, std::set<std::string>> MODEL_PARAMS = {
  {"1TiCM",  {"K1", "vb"}},
  {"1TCM",   {"K1", "k2", "vb"}},
  {"1TidCM", {"K1", "d", "vb"}},
  {"1TdCM",  {"K1", "k2", "d", "vb"}},
  {"2TiCM",  {"K1", "k2", "k3", "vb"}},
  {"2TCM",   {"K1", "k2", "k3", "k4", "vb"}},
  {"2TidCM", {"K1", "k2", "k3", "d", "vb"}},
  {"2TdCM",  {"K1", "k2", "k3", "k4", "d", "vb"}}
};

static bool isSubset(const std::set<std::string>& a,
                     const std::set<std::string>& b)
{
  return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

static void countConstraints(const std::set<std::string>& restricted,
                             const std::set<std::string>& full,
                             int& r_b,
                             int& r_i)
{
  r_b = 0;
  r_i = 0;

  for (const auto& p : full)
  {
    if (restricted.count(p)) continue;

    if (p == "d") r_i++;
    else if (!p.empty() && p[0] == 'k') r_b++;
  }
}


static double median(std::vector<double> v){
    size_t n=v.size();
    std::nth_element(v.begin(), v.begin()+n/2, v.end());
    double m=v[n/2];
    if(n%2==0){
        auto it=std::max_element(v.begin(), v.begin()+n/2);
        m=0.5*(m+*it);
    }
    return m;
}

static bool check_if_constant(const std::vector<double>& x, const double thres = 0.001){
    if(x.size()<2) return true;
    double c = median(x);              // robust center
    double maxdev = 0.0;
    for(double v: x) maxdev = std::max(maxdev, std::abs(v - c));
    return maxdev < thres;
}

static double norm_cdf(double x)
{
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double chi2_cdf(double x, int df)
{
    // Regularized gamma function approximation
    // Here we use std::tgamma and std::lgamma for stability
    // P(k/2, x/2) = lower_gamma(k/2, x/2) / Gamma(k/2)
    double k = df * 0.5;
    double t = x * 0.5;

    // Series expansion for lower incomplete gamma (simple for moderate df)
    double sum = 1.0 / k;
    double term = sum;
    for (int n = 1; n < 100; ++n) {
        term *= t / (k + n);
        sum += term;
        if (term < 1e-12) break;
    }

    double result = std::exp(-t + k * std::log(t) - std::lgamma(k)) * sum;
    return result;
}

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
                                    QProgressBar* ProgressBar,
                                    QPushButton* stopButton,
                                    std::atomic<bool>& stopRequested
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
  // Fetch PET
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(petID));
  if (!petNode) {
    return;
  }
  // Fetch Segmentation
  vtkMRMLSegmentationNode* segNode = vtkMRMLSegmentationNode::SafeDownCast(shNode->GetItemDataNode(segID));
  if (!segNode) {
    return;
  }
  setupSeg(segNode);


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
    std::cerr << "Could not find sequence for the segmentation PET." << std::endl;
    return;
  }

  // Get TAC
  this->TAC(sequencePETNode, segSequenceNode, segmentsID, segmentTACs, segmentTACsnames, ProgressBar, stopButton, stopRequested);


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

   using Clock = std::chrono::steady_clock;

  // 1. Make sure "Binary labelmap" is available as a representation
  const std::string labelmapRep = vtkSegmentationConverter::GetSegmentationBinaryLabelmapRepresentationName();

  // 2. Set "Binary labelmap" as the master representation (i.e., source)
  auto start = Clock::now();
  if (!segmentation->ContainsRepresentation(labelmapRep))
  {
    segmentation->CreateRepresentation(labelmapRep);
  }
  auto end = Clock::now();

  start = Clock::now();
  if (segmentation->GetSourceRepresentationName()
      != labelmapRep)
  {
    segmentation->SetSourceRepresentationName(labelmapRep);
  }
  end = Clock::now();

  // Do NOT eagerly create Closed surface here.
  // 3. Ensure "Closed surface" representation is present
  // const std::string closedSurfRep = vtkSegmentationConverter::GetSegmentationClosedSurfaceRepresentationName();
  // segmentation->CreateRepresentation(closedSurfRep);

}


VoxelStatistics vtkSlicerKMAPLogic::ComputeVoxelStatistics(vtkMRMLScalarVolumeNode* PETVolume, vtkImageData* labelmap, int labelValue)
{
  vtkImageData* petImage = PETVolume->GetImageData();
  VoxelStatistics stats;
  std::vector<double> values;

  int dims[3];
  double spacing[3];
  petImage->GetDimensions(dims);
  PETVolume->GetSpacing(spacing);

  double voxelVolume = spacing[0] * spacing[1] * spacing[2];

  vtkDataArray* petArray = petImage->GetPointData()->GetScalars();
  vtkDataArray* labelArray = labelmap->GetPointData()->GetScalars();

  int max_ijk[3] = {0,0,0};

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
          if (val > stats.max) {
            stats.max = val;
            max_ijk[0] = x;
            max_ijk[1] = y;
            max_ijk[2] = z;
          }
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

  // ================= SUVPEAK =================
  // --- 1) Sphere radius (1cc) ---
  double V_mm3 = 1000.0;
  double radius_mm = std::cbrt((3.0 * V_mm3) / (4.0 * M_PI));

  int rx = std::ceil(radius_mm / spacing[0]);
  int ry = std::ceil(radius_mm / spacing[1]);
  int rz = std::ceil(radius_mm / spacing[2]);

  // --- 2) Compute mean inside sphere ---
  double sumPeak = 0.0;
  int countPeak = 0;
  for (int z = max_ijk[2] - rz; z <= max_ijk[2] + rz; ++z)
  for (int y = max_ijk[1] - ry; y <= max_ijk[1] + ry; ++y)
  for (int x = max_ijk[0] - rx; x <= max_ijk[0] + rx; ++x)
  {
    if (x < 0 || y < 0 || z < 0 ||
        x >= dims[0] || y >= dims[1] || z >= dims[2])
      continue;
    double dx = (x - max_ijk[0]) * spacing[0];
    double dy = (y - max_ijk[1]) * spacing[1];
    double dz = (z - max_ijk[2]) * spacing[2];

    double dist2 = dx*dx + dy*dy + dz*dz;
    if (dist2 > radius_mm * radius_mm)
      continue;

    int ijk[3] = {x,y,z};
    vtkIdType idx = petImage->ComputePointId(ijk);

    int label = static_cast<int>(labelArray->GetComponent(idx, 0));
    if (label != labelValue)
      continue;

    double val = petArray->GetComponent(idx, 0);

    sumPeak += val;
    countPeak++;
  }

  if (countPeak > 0)
    stats.peak = sumPeak / countPeak;
  else
    stats.peak = std::numeric_limits<double>::quiet_NaN();

  return stats;
}

double vtkSlicerKMAPLogic::boundaryLRTPvalue(double LR, int r_b, int r_i)
{
  double p = 0.0;

  for (int j = 0; j <= r_b; ++j)
  {
    double weight = std::pow(0.5, r_b) * std::tgamma(r_b + 1) /
                    (std::tgamma(j + 1) * std::tgamma(r_b - j + 1));

    int df = r_i + j;

    double tail = (df == 0)
        ? (LR <= 0.0 ? 1.0 : 0.0)
        : (1.0 - chi2_cdf(LR, df));

    p += weight * tail;
  }

  return std::max(0.0, std::min(1.0, p));
}

ModelComparisonResult vtkSlicerKMAPLogic::compareModels(
    const std::string& modelA,
    const std::string& modelB,
    const TCMParameters& m1,
    const TCMParameters& m2
)
{
  ModelComparisonResult res;

  const auto& paramsA = MODEL_PARAMS.at(modelA);
  const auto& paramsB = MODEL_PARAMS.at(modelB);

  bool A_in_B = isSubset(paramsA, paramsB);
  bool B_in_A = isSubset(paramsB, paramsA);

  // =====================
  // CASE 1: NESTED
  // =====================
  if (A_in_B || B_in_A)
  {
    const std::set<std::string>* restricted;
    const std::set<std::string>* full;

    double LL_r, LL_f;

    if (A_in_B)
    {
      restricted = &paramsA;
      full = &paramsB;
      LL_r = m1.loglik;
      LL_f = m2.loglik;
    }
    else
    {
      restricted = &paramsB;
      full = &paramsA;
      LL_r = m2.loglik;
      LL_f = m1.loglik;
    }

    double LR = 2.0 * (LL_f - LL_r);
    if (LR < 0.0) LR = 0.0;

    int r_b = 0, r_i = 0;
    countConstraints(*restricted, *full, r_b, r_i);

    res.type = "LRT";
    // res.statistic = LR;
    res.p_value = this->boundaryLRTPvalue(LR, r_b, r_i);

    return res;
  }

  // =====================
  // CASE 2: NON-NESTED
  // =====================
  else
  {
    // average weights (your current logic)
    std::vector<double> wgt(m1.weights.size());
    for (size_t i = 0; i < wgt.size(); ++i)
      wgt[i] = 0.5 * (m1.weights[i] + m2.weights[i]);

    const std::vector<double>* wgt_ptr = &wgt;

    double p = this->computeVuongP(
        m1.r,
        m2.r,
        wgt_ptr,
        m1.dof,
        m2.dof,
        VuongCorrection::BIC,
        Tail::TwoSided
    );

    res.type = "Vuong";
    // res.statistic = std::numeric_limits<double>::quiet_NaN(); // optional
    res.p_value = p;

    return res;
  }
}

void vtkSlicerKMAPLogic::TAC(vtkMRMLSequenceNode* sequencePETNode,
                             vtkMRMLSequenceNode* segSequenceNode,
                             std::vector<QString> segmentsID,
                             std::map<std::string, std::vector<VoxelStatistics>>& segmentTACs,
                             std::map<std::string, std::string>& segmentTACsnames,
                             QProgressBar* ProgressBar,
                             QPushButton* stopButton,
                             std::atomic<bool>& stopRequested
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

    segmentTACs[segmentID].clear();
    segmentTACs[segmentID].resize(numberOfTimepoints);
  }

  if (ProgressBar) {
    ProgressBar->setFormat("Computing TAC (%p%)");
    ProgressBar->setVisible(true);
    ProgressBar->setMinimum(0);
    ProgressBar->setMaximum(100);
    ProgressBar->setValue(0);
    stopButton->setVisible(true);
    stopButton->show();
    qApp->processEvents();
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

    #ifdef HAVE_OPENMP
    int max_hw_threads = omp_get_num_procs();
    omp_set_num_threads(max_hw_threads);
    #pragma omp parallel for
    #endif
    for (int s = 0; s < segmentsID.size(); ++s)
    {
      if (stopRequested) continue;
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
      VoxelStatistics stats;
      if (!labelmap)
      {
        std::cerr << "Failed to generate labelmap for segment: " << segmentID << " at timepoint " << i << std::endl;
        stats.keep = false;
      } else {
        stats = ComputeVoxelStatistics(PETVolume, labelmap, 1);
      }

      segmentTACs[segmentID][i] = stats;
    }
    if (ProgressBar){
      ProgressBar->setValue(static_cast<double>(i + 1) / numberOfTimepoints*100.);
      // qApp->processEvents();
    }
    if (stopRequested) {
      stopButton->setVisible(false);
      segmentTACs.clear();
      qApp->processEvents();
      break;
    }
  }
  if (ProgressBar) {
    stopButton->setVisible(false);
    qApp->processEvents();
  }
}

double vtkSlicerKMAPLogic::computeLogLik(const std::vector<double>& y,
                                         const std::vector<double>& fitted,
                                         const std::vector<double>* wgt)
{
    size_t n = y.size();

    std::vector<double> weights;
    if (wgt == nullptr)
        weights.assign(n, 1.0);
    else {
        if (wgt->size() != n)
            throw std::invalid_argument("weights size must match obs");
        weights = *wgt;
    }

    if (y.size() != fitted.size() || y.size() != weights.size()) {
        throw std::invalid_argument("Input vectors must have same length");
    }

    double logLik = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        double r = y[i] - fitted[i];
        double wi = weights[i];
        if (wi <= 0) continue; // skip invalid weights

        logLik += -0.5 * (std::log(2.0 * M_PI / wi) + wi * r * r);
    }
    return logLik;
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

  params.keep.resize(Nframe);
  for (int iz = 0; iz < Nframe; ++iz) {
    params.keep[iz] = wt[iz]!=0. ? true : false;
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
    params.vb = sens[0] ? fitted_params[0] : std::numeric_limits<double>::quiet_NaN();
    params.K1 = sens[1] ? fitted_params[1] : std::numeric_limits<double>::quiet_NaN();
    params.k2 = sens[2] ? fitted_params[2] : std::numeric_limits<double>::quiet_NaN();
    params.td = sens[3] ? fitted_params[3] : std::numeric_limits<double>::quiet_NaN();
    params.Ki = params.K1;
    params.DV = params.K1/(params.k2 + 1e-16);
    params.weights.assign(wt, wt + Nframe);
    params.dof = std::count(sens, sens + 4, true);
    std::vector<double> fittedTCMvalues(fitted_curve,
                                        fitted_curve + Nframe);
    std :: vector<double> tac_vec(tac_flatten, tac_flatten + Nframe);
    std :: vector<double> predicted_tac(fitted_curve, fitted_curve + Nframe);
    params.AIC  = this->computeAIC(tac_vec, predicted_tac, params.dof, wgt);
    params.BIC  = this->computeBIC(tac_vec, predicted_tac, params.dof, wgt);
    params.MASE = this->MASE(tac_vec, predicted_tac, wgt);
    params.chi2 = this->computeChi2(tac_vec, predicted_tac, wgt) / (Nframe-params.dof);
    params.loglik = this->computeLogLik(tac_vec, predicted_tac, wgt);
    params.r.resize(Nframe);
    for (size_t i = 0; i < Nframe; ++i) {
      params.r[i] = tac_vec[i] - predicted_tac[i];
    }
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
    params.vb = sens[0] ? fitted_params[0] : std::numeric_limits<double>::quiet_NaN();
    params.K1 = sens[1] ? fitted_params[1] : std::numeric_limits<double>::quiet_NaN();
    params.k2 = sens[2] ? fitted_params[2] : std::numeric_limits<double>::quiet_NaN();
    params.k3 = sens[3] ? fitted_params[3] : std::numeric_limits<double>::quiet_NaN();
    params.k4 = sens[4] ? fitted_params[4] : std::numeric_limits<double>::quiet_NaN();
    params.td = sens[5] ? fitted_params[5] : std::numeric_limits<double>::quiet_NaN();
    params.Ki = params.K1 * params.k3 / (params.k2 + params.k3);
    params.DV = params.K1/(params.k2 + 1e-16) * (1 + params.k3/(params.k4 + 1e-16));
    params.weights.assign(wt, wt + Nframe);
    params.dof = std::count(sens, sens + 6, true);
    std::vector<double> fittedTCMvalues(fitted_curve,
                                        fitted_curve + Nframe);
    std :: vector<double> tac_vec(tac_flatten, tac_flatten + Nframe);
    std :: vector<double> predicted_tac(fitted_curve, fitted_curve + Nframe);
    params.AIC  = this->computeAIC(tac_vec, predicted_tac, params.dof, wgt);
    params.BIC  = this->computeBIC(tac_vec, predicted_tac, params.dof, wgt);
    params.MASE = this->MASE(tac_vec, predicted_tac, wgt);
    params.chi2 = this->computeChi2(tac_vec, predicted_tac, wgt)/ (Nframe-params.dof);
    params.loglik = this->computeLogLik(tac_vec, predicted_tac, wgt);
    params.r.resize(Nframe);
    for (size_t i = 0; i < Nframe; ++i) {
      params.r[i] = tac_vec[i] - predicted_tac[i];
    }
  } else {
    std :: cerr << "Forbidden number of tissue compartment: " << n_tc << std :: endl;
  }

  // --- Cleanup ---
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

    double mse = (numfrm > 0 ? ss / static_cast<double>(numfrm) : 0.0);
    if (mse <= 0.0 || !std::isfinite(mse))
    {
        mse = std::numeric_limits<double>::min();
    }
    double AIC = numfrm * std::log(mse) + 2.0 * numpar_new;

    // Apply AICc correction if needed
    if (aicc && numfrm > (numpar_new + 1) && (static_cast<double>(numfrm) / numpar < 40.0))
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

double vtkSlicerKMAPLogic::computeBIC(const std::vector<double>& obs,
                                      const std::vector<double>& est,
                                      int numpar,
                                      const std::vector<double>* wgt)
{
    size_t n = obs.size();
    if (est.size() != n)
        throw std::invalid_argument("obs and est must have same length");

    // Default weights = 1
    std::vector<double> weights;
    if (wgt == nullptr)
        weights.assign(n, 1.0);
    else {
        if (wgt->size() != n)
            throw std::invalid_argument("weights size must match obs");
        weights = *wgt;
    }

    // Weighted SSE
    double ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = obs[i] - est[i];
        ss += weights[i] * diff * diff;
    }

    int p = numpar + 1; // intercept term
    double BIC = n * std::log(ss / static_cast<double>(n))
                 + p * std::log(static_cast<double>(n));

    return BIC;
}

double vtkSlicerKMAPLogic::computeR2(const std::vector<double>& obs,
                                     const std::vector<double>& est,
                                     const std::vector<double>* wgt)
{
    size_t n = obs.size();
    if (est.size() != n)
        throw std::invalid_argument("obs and est must have same length");

    std::vector<double> weights;
    if (wgt == nullptr)
        weights.assign(n, 1.0);
    else {
        if (wgt->size() != n)
            throw std::invalid_argument("weights size must match obs");
        weights = *wgt;
    }

    // Weighted mean of obs
    double wsum = 0.0, mean = 0.0;
    for (size_t i = 0; i < n; ++i) {
        wsum += weights[i];
        mean += weights[i] * obs[i];
    }
    mean /= wsum;

    // SST and SSE
    double sst = 0.0, sse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sst += weights[i] * (obs[i] - mean) * (obs[i] - mean);
        double diff = obs[i] - est[i];
        sse += weights[i] * diff * diff;
    }

    if (sst == 0.0)
        return 1.0; // all obs identical, perfect fit by definition
    return 1.0 - (sse / sst);
}

double vtkSlicerKMAPLogic::computeChi2(const std::vector<double>& y,
                                       const std::vector<double>& fitted,
                                       const std::vector<double>* wgt)
{
    size_t n = y.size();

    std::vector<double> weights;
    if (wgt == nullptr)
        weights.assign(n, 1.0);
    else {
        if (wgt->size() != n)
            throw std::invalid_argument("weights size must match obs");
        weights = *wgt;
    }

    if (y.size() != fitted.size() || y.size() != weights.size()) {
        throw std::invalid_argument("Vectors must have the same length");
    }

    double chi2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double r = y[i] - fitted[i];
        chi2 += weights[i] * r * r;
    }
    return chi2;
}

double vtkSlicerKMAPLogic::computeVuongP(const std::vector<double>& r1,
                                         const std::vector<double>& r2,
                                         const std::vector<double>* wgt,
                                         int k1, int k2,
                                         VuongCorrection corr,
                                         Tail tail
                                       )
{
    size_t N = r1.size();
    if (r1.size() != N || r2.size() != N)
        throw std::invalid_argument("Models do not have the same number of observations");

    // Weights: default to 1, then normalize to sum=1 (w~)
    std::vector<double> w(N, 1.0);
    if (wgt)
    {
      if (wgt->size() != N) {
        std::ostringstream oss;
        oss << "weights size (" << wgt->size() << ") must match residuals size (" << N << ")";
        throw std::invalid_argument(oss.str());
      }
      w = *wgt;
    }

    double sumw = 0.0;
    for (double wi : w) {
      if (wi <= 0.0) throw std::invalid_argument("All weights must be positive.");
      sumw += wi;
    }
    for (double& wi : w) wi /= sumw; // now sum(w) = 1

    // Pointwise log-likelihood differences under Gaussian with unit variance:
    // m_i = -0.5*r1^2 - (-0.5*r2^2) = 0.5*(r2^2 - r1^2)
    std::vector<double> m(N);
    for (size_t i = 0; i < N; ++i)
      m[i] = 0.5 * (r2[i]*r2[i] - r1[i]*r1[i]);

    // Weighted mean of m
    double mean_m = 0.0;
    for (size_t i = 0; i < N; ++i) mean_m += w[i] * m[i];

    // Unbiased weighted variance of m (with  reliability weights)
    // s2 = sum_i w_i (m_i - mean)^2 / (1 - sum_i w_i^2)
    double sumw2 = 0.0;
    for (double wi : w) sumw2 += wi*wi;

    double num = 0.0;
    for (size_t i = 0; i < N; ++i) {
      const double d = (m[i] - mean_m);
      num += w[i] * d * d;
    }
    const double denom = std::max(1e-16, 1.0 - sumw2); // guard
    const double s2 = num / denom;
    const double s  = std::sqrt(std::max(s2, 1e-16));

    // Effective sample size (Kish)
    const double n_eff = 1.0 / sumw2; // equals N if all weights equal

    // Penalty c_N
    const int dk = (k1 - k2);
    double cN = 0.0;
    if (corr == VuongCorrection::AIC) {
      cN = static_cast<double>(dk);
    } else if (corr == VuongCorrection::BIC) {
      cN = 0.5 * static_cast<double>(dk) * std::log(static_cast<double>(N));
      // If you're heavily weighting, you *may* prefer log(n_eff) here; keep N for standard practice.
    }

    // Corrected mean: subtract penalty per observation
    const double mean_corr = mean_m - (cN / static_cast<double>(N));

    // Vuong Z
    const double Z = std::sqrt(n_eff) * mean_corr / s;

    // p-value
    double p = 1.0;

    if (tail == Tail::TwoSided) {
      const double phi = norm_cdf(Z);
      p = 2.0 * std::min(phi, 1.0 - phi);
    } else if (tail == Tail::Model1Greater) {
      // H1: model1 better (Z large positive)
      p = 1.0 - norm_cdf(Z);
    } else { // Tail::Model2Greater
      // H1: model2 better (Z large negative)
      p = norm_cdf(Z);
    }

    return std::max(0.0, std::min(1.0, p));
}



double vtkSlicerKMAPLogic::computeLRTP(double logLik1,
                                       double logLik2,
                                       int df2, int df1
                                      )
{
  // Determine complexity ordering
  int dfDiff = df2 - df1;
  if (dfDiff == 0) {
    return std::numeric_limits<double>::quiet_NaN(); // cannot compute
  }

  // Ensure model2 is more complex
  if (dfDiff < 0) {
      std::swap(logLik1, logLik2);
      std::swap(df1, df2);
      dfDiff = -dfDiff;
  }

  double Lambda = -2.0 * (logLik1 - logLik2);
  if (Lambda < 0) Lambda = 0.0;

  double p = 1.0 - chi2_cdf(Lambda, dfDiff);
  return p;
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
  params.keep.clear();
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset)
      {
          double x = intCp[i] / (Cp[i] + 1e-16);
          double y = tac[i]   / (Cp[i] + 1e-16);
          outX.push_back(x);
          outY.push_back(y);
          outframe.push_back(i+1);
          if (wgt) {
            wgt_adj.push_back((*wgt)[i]);
            params.keep.push_back((*wgt)[i]==0. ? false : true);
          } else {
            params.keep.push_back(true);
          }
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
    std::vector<size_t> keepIndices;
    if (!wgt_adj.empty())
    {
        for (size_t iz = 0; iz < wgt_adj.size(); ++iz)
        {
            if (wgt_adj[iz] != 0.0)  // or use std::abs(wgt_adj[i]) < eps for floating-point
            {
                keepIndices.push_back(iz);

            } else {
              weights(iz) = 0.;
            }
        }
    }
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
      for (auto i: keepIndices) {
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
  params.r.resize(n);
  for (size_t i = 0; i < n; ++i) {
    fittedValues[i] = intercept + slope * outX[i];
    params.r[i] = outY[i] - fittedValues[i];
  }

  // AIC and MASE
  params.AIC = computeAIC(outY, fittedValues, 2, wgt ? &wgt_adj : nullptr);
  params.MASE = MASE(outY, fittedValues, wgt ? &wgt_adj : nullptr);
  params.R2 = computeR2(outY, fittedValues, wgt ? &wgt_adj : nullptr);
  params.chi2 = computeChi2(outY, fittedValues, wgt ? &wgt_adj : nullptr)/(n-2);

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

  params.dof = 2;
  if (wgt == nullptr)
  {
      params.weights.assign(n, 1.0);  // fills with ones
  }
  else
  {
      params.weights = wgt_adj;  // copy from your pre-filled vector
  }
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
  params.keep.clear();
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset)
      {
          double x = intCp[i] / (tac[i] + 1e-16);
          double y = intCt[i] / (tac[i] + 1e-16);
          outX.push_back(x);
          outY.push_back(y);
          outframe.push_back(i+1);
          if (wgt) {
            wgt_adj.push_back((*wgt)[i]);
            params.keep.push_back((*wgt)[i]==0. ? false : true);
          } else {
            params.keep.push_back(true);
          }
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
      std::vector<size_t> keepIndices;
      if (!wgt_adj.empty())
      {
          for (size_t iz = 0; iz < wgt_adj.size(); ++iz)
          {
              if (wgt_adj[iz] != 0.0)  // or use std::abs(wgt_adj[i]) < eps for floating-point
              {
                  keepIndices.push_back(iz);
              } else {
                weights(iz) = 0.;
              }
          }
      }
      Eigen::MatrixXd W = Eigen::MatrixXd::Identity(n, n);
      Eigen::VectorXd residuals(n), prev_coeff(2);
      for (int iter = 0; iter < max_iter; ++iter)
      {
          W.diagonal() = weights;
          coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
          if (iter > 0 && (coeff - prev_coeff).norm() < tol) break;
          prev_coeff = coeff;
          residuals = Yv - A * coeff;
          for (auto i: keepIndices) {
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
  params.r.resize(n);
  for (size_t i = 0; i < n; ++i) {
    fittedValues[i] = intercept + slope * outX[i];
    params.r[i] = outY[i] - fittedValues[i];
  }

  // AIC and MASE
  params.AIC = computeAIC(outY, fittedValues, 2, wgt ? &wgt_adj : nullptr);
  params.MASE = MASE(outY, fittedValues, wgt ? &wgt_adj : nullptr);
  params.R2 = computeR2(outY, fittedValues, wgt ? &wgt_adj : nullptr);
  params.chi2 = computeChi2(outY, fittedValues, wgt ? &wgt_adj : nullptr)/(n-2);

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

  params.dof = 2;
  if (wgt == nullptr)
  {
      params.weights.assign(n, 1.0);  // fills with ones
  }
  else
  {
      params.weights = wgt_adj;  // copy from your pre-filled vector
  }
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
  params.keep.clear();
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset)
      {
          double x = intCp[i] / (Cp[i] + 1e-16);
          double y = intCt[i] / (Cp[i] + 1e-16);
          outX.push_back(x);
          outY.push_back(y);
          outframe.push_back(i+1);
          if (wgt) {
            wgt_adj.push_back((*wgt)[i]);
            params.keep.push_back((*wgt)[i]==0. ? false : true);
          } else {
            params.keep.push_back(true);
          }
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
      std::vector<size_t> keepIndices;
      if (!wgt_adj.empty())
      {
          for (size_t iz = 0; iz < wgt_adj.size(); ++iz)
          {
              if (wgt_adj[iz] != 0.0)  // or use std::abs(wgt_adj[i]) < eps for floating-point
              {
                  keepIndices.push_back(iz);
              } else {
                weights(iz) = 0.;
              }
          }
      }
      Eigen::MatrixXd W = Eigen::MatrixXd::Identity(n, n);
      Eigen::VectorXd residuals(n), prev_coeff(2);
      for (int iter = 0; iter < max_iter; ++iter)
      {
          W.diagonal() = weights;
          coeff = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv);
          if (iter > 0 && (coeff - prev_coeff).norm() < tol) break;
          prev_coeff = coeff;
          residuals = Yv - A * coeff;
          for (auto i: keepIndices) {
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
  params.r.resize(n);
  for (size_t i = 0; i < n; ++i){
    fittedValues[i] = intercept + slope * outX[i];
    params.r[i] = outY[i] - fittedValues[i];
  }

  // AIC and MASE
  params.AIC = computeAIC(outY, fittedValues, 2, wgt ? &wgt_adj : nullptr);
  params.MASE = MASE(outY, fittedValues, wgt ? &wgt_adj : nullptr);
  params.R2 = computeR2(outY, fittedValues, wgt ? &wgt_adj : nullptr);
  params.chi2 = computeChi2(outY, fittedValues, wgt ? &wgt_adj : nullptr)/(n-2);

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

  params.dof = 2;
  if (wgt == nullptr)
  {
      params.weights.assign(n, 1.0);  // fills with ones
  }
  else
  {
      params.weights = wgt_adj;  // copy from your pre-filled vector
  }
}

void vtkSlicerKMAPLogic::Image2Flatten(
    vtkIdType petID,
    std::vector<std::vector<double>>& flatten_voxels_values,
    int (&dims)[3],
    int& numberOfTimepoints,
    QProgressBar* ProgressBar,
    QPushButton* stopButton,
    std::atomic<bool>& stopRequested
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
  // vtkMRMLScalarVolumeNode* ctNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(ctID));
  // if (!ctNode) {
  //   return;
  // }

  // Fetch PET
  vtkMRMLScalarVolumeNode* petNode = vtkMRMLScalarVolumeNode::SafeDownCast(shNode->GetItemDataNode(petID));
  if (!petNode) {
    return;
  }

  // Find corresponding sequence
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

  if (!sequencePETNode)
  {
    std::cerr << "Invalid input nodes!" << std::endl;
    return;
  }

  numberOfTimepoints = sequencePETNode->GetNumberOfDataNodes();

  // flatten_values.clear();
  // flatten_values.reserve(numberOfTimepoints * petNode->GetImageData()->GetNumberOfPoints());
  //
  // for (int i = 0; i < numberOfTimepoints; ++i)
  // {
  //   std::string indexValue = sequencePETNode->GetNthIndexValue(i);
  //   auto* PETVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
  //     sequencePETNode->GetDataNodeAtValue(indexValue));
  //   if (!PETVolume)
  //   {
  //     std::cerr << "Missing data for timepoint " << i << std::endl;
  //     continue;
  //   }
  //   vtkImageData* petImage = PETVolume->GetImageData();
  //   petImage->GetDimensions(dims);
  //   vtkDataArray* petArray = petImage->GetPointData()->GetScalars();
  //   if (!petArray)
  //     continue;
  //
  //   vtkIdType nVoxels = petArray->GetNumberOfTuples();
  //   for (vtkIdType v = 0; v < nVoxels; ++v)
  //   {
  //     flatten_values.push_back(petArray->GetComponent(v, 0));
  //   }
  // }

  auto* firstVol = vtkMRMLScalarVolumeNode::SafeDownCast(
    sequencePETNode->GetNthDataNode(0));
  if (!firstVol) return;

  vtkImageData* firstImage = firstVol->GetImageData();
  firstImage->GetDimensions(dims);

  vtkDataArray* firstArray = firstImage->GetPointData()->GetScalars();
  if (!firstArray) return;

  vtkIdType nVoxels = firstArray->GetNumberOfTuples();

  flatten_voxels_values.clear();
  flatten_voxels_values.resize(nVoxels);

  #ifdef HAVE_OPENMP
  int max_hw_threads = omp_get_num_procs();
  omp_set_num_threads(max_hw_threads);
  #pragma omp parallel for
  #endif
  for (int v = 0; v < nVoxels; ++v)
  {
    flatten_voxels_values[v].resize(numberOfTimepoints);
  }

  if (ProgressBar)
  {
      ProgressBar->setFormat("Flattening dPET (%p%)");
      ProgressBar->setVisible(true);
      ProgressBar->setMinimum(0);
      ProgressBar->setMaximum(100);
      ProgressBar->setValue(0);
      stopButton->setVisible(true);
      stopButton->show();
      qApp->processEvents();
  }
  for (int i = 0; i < numberOfTimepoints; ++i)
  {
    std::string indexValue = sequencePETNode->GetNthIndexValue(i);
    auto* PETVolume = vtkMRMLScalarVolumeNode::SafeDownCast(
      sequencePETNode->GetDataNodeAtValue(indexValue));
    if (!PETVolume)
    {
      std::cerr << "Missing data for timepoint " << i << std::endl;
      continue;
    }

    vtkImageData* petImage = PETVolume->GetImageData();
    vtkDataArray* petArray = petImage->GetPointData()->GetScalars();
    if (!petArray) continue;

    vtkIdType nVoxels = petArray->GetNumberOfTuples();


    #ifdef HAVE_OPENMP
    #pragma omp parallel for
    #endif
    for (int v = 0; v < nVoxels; ++v)
    {
      if (stopRequested) continue;
      flatten_voxels_values[v][i] = petArray->GetComponent(v, 0);
    }
    if (stopRequested) break;
    if (ProgressBar){
      ProgressBar->setValue(static_cast<double>(i + 1) / numberOfTimepoints*100.);
      // qApp->processEvents();
    }
  }

  if (ProgressBar)
  {
    ProgressBar->setVisible(false);
    stopButton->setVisible(false);
    // stopButton->deleteLater();
    if (stopRequested) {
      flatten_voxels_values.clear();
    }
    qApp->processEvents();
  }
}

vtkMRMLScalarVolumeNode* vtkSlicerKMAPLogic::Flatten2Image(
    const std::vector<double>& flatten_values,
    const int dims[3],
    const std::string& name
)
{
  vtkIdType nVoxels = static_cast<vtkIdType>(dims[0]) * dims[1] * dims[2];
  if (flatten_values.size() != static_cast<size_t>(nVoxels))
  {
    std::cerr << "Flatten2Image: size mismatch! Expected " << nVoxels
              << " values, got " << flatten_values.size() << std::endl;
    return nullptr;
  }

  vtkMRMLScene* scene = this->GetMRMLScene();

  if (!scene)
  {
      std::cerr
          << "Flatten2Image: MRML scene is null"
          << std::endl;
      return nullptr;
  }

  // Allocate vtkImageData
  vtkNew<vtkImageData> image;
  image->SetDimensions(dims);
  image->AllocateScalars(VTK_DOUBLE, 1);

  vtkDoubleArray* scalars = vtkDoubleArray::SafeDownCast(
      image->GetPointData()->GetScalars());
  if (!scalars)
  {
    std::cerr << "Flatten2Image: could not allocate double scalars" << std::endl;
    return nullptr;
  }

  double* scalarPtr = scalars->GetPointer(0);
  if (!scalarPtr)
  {
      std::cerr
          << "Flatten2Image: scalar pointer is null."
          << std::endl;
      return nullptr;
  }

  // Safe copy into the VTK array
  #ifdef HAVE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (vtkIdType i = 0; i < nVoxels; ++i)
  {
      scalarPtr[i] =
          flatten_values[static_cast<size_t>(i)];
  }

  // Create MRML volume node
  vtkMRMLScalarVolumeNode* volumeNode =
      vtkMRMLScalarVolumeNode::SafeDownCast(
          scene->AddNewNodeByClass(
              "vtkMRMLScalarVolumeNode",
              name.c_str()));
  if (!volumeNode)
  {
      std::cerr
          << "Flatten2Image: could not create MRML volume node"
          << std::endl;
      return nullptr;
  }

  volumeNode->SetAndObserveImageData(image);
  volumeNode->CreateDefaultDisplayNodes();

  return volumeNode;
}


void vtkSlicerKMAPLogic::Patlak4Img(
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
    QProgressBar* progressBar,
    int numThreads,
    QPushButton* stopButton
)
{
    #ifdef HAVE_OPENMP
    int max_hw_threads = omp_get_num_procs();
    if (numThreads > 0) {
        int n = std::min(numThreads, max_hw_threads);
        omp_set_num_threads(n);
    }
    #endif

    const double EPS = 1e-12;
    size_t N = Cp.size();
    if (framing.size() != N) return;
    size_t nVoxels = voxels.size();
    if (nVoxels == 0) return;

    // ---------- 1) Precompute frameScaled, timeAlong, intCp, invCp ----------
    std::vector<double> frameScaled(N);
    for (size_t i = 0; i < N; ++i) frameScaled[i] = framing[i] / framingNorm;

    std::vector<double> timeAlong(N);
    std::partial_sum(frameScaled.begin(), frameScaled.end(), timeAlong.begin());

    std::vector<double> intCp(N, 0.0);
    double acc = 0.0;
    for (size_t i = 1; i < N; ++i)
    {
        acc += 0.5 * (Cp[i] + Cp[i-1]) * frameScaled[i];
        intCp[i] = acc;
    }

    std::vector<double> invCp(N);
    for (size_t i = 0; i < N; ++i) invCp[i] = 1.0 / (Cp[i] + EPS);

    // ---------- 2) Filter frames >= timeOffset ----------
    std::vector<double> X_all;
    std::vector<int> keepIndex;
    for (size_t i = 0; i < N; ++i)
    {
        if (timeAlong[i] >= timeOffset)
        {
            X_all.push_back(intCp[i] * invCp[i]);
            keepIndex.push_back(static_cast<int>(i));
        }
    }
    size_t n = X_all.size();
    if (n < 2) return;

    // ---------- 3) Build design matrix A ----------
    Eigen::MatrixXd A(n, 2);
    A.col(0) = Eigen::VectorXd::Ones(n);
    for (size_t i = 0; i < n; ++i) A(i,1) = X_all[i];

    double meanX = 0.0, stdX = 1.0;
    if (standardize)
    {
        meanX = std::accumulate(X_all.begin(), X_all.end(), 0.0) / double(n);
        double sq = 0.0;
        for (size_t i = 0; i < n; ++i) sq += X_all[i] * X_all[i];
        stdX = std::sqrt(sq / double(n) - meanX * meanX);
        if (stdX < EPS) stdX = 1.0;
        for (size_t i = 0; i < n; ++i) A(i,1) = (X_all[i] - meanX) / stdX;
    }

    // Precompute pseudoinverse for unweighted OLS
    Eigen::Matrix2d AtA = A.transpose() * A;
    Eigen::Matrix2d AtA_inv = AtA.ldlt().solve(Eigen::Matrix2d::Identity());
    Eigen::MatrixXd pinv = AtA_inv * A.transpose(); // 2 x n

    // ---------- 4) Initialize output and progress ----------
    outputParams.clear();
    outputParams.resize(nVoxels);


    if (progressBar)
    {
        progressBar->setFormat("Fitting Patlak (%p%)");
        progressBar->setVisible(true);
        progressBar->setMinimum(0);
        progressBar->setMaximum(100);
        progressBar->setValue(0);
        stopButton->setVisible(true);
        stopButton->show();
        qApp->processEvents();
    }

    std::atomic<size_t> voxProcessed(0);

    // ---------- 5) Parallel voxel-wise processing ----------
    #ifdef HAVE_OPENMP
    #pragma omp parallel
    #endif
    {
        std::vector<double> Yvec(n);
        std::vector<double> fitted(n);
        std::vector<double> wgt_adj(n);
        Eigen::Vector2d coeff;
        Eigen::VectorXd coeff_vec(2);

        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int v = 0; v < nVoxels; ++v)
        {
            if (stopRequested) continue;
            MTGAParameters params;
            const std::vector<double>& tac = voxels[v];

            // Build Y and weights
            for (size_t ii = 0; ii < n; ++ii)
            {
                int idx = keepIndex[ii];
                Yvec[ii] = tac[idx] * invCp[idx];
                if (wgt_global) wgt_adj[ii] = (*wgt_global)[idx];
            }

            // Standardize Y if requested
            double meanY = 0.0, stdY = 1.0;
            if (standardize)
            {
                meanY = std::accumulate(Yvec.begin(), Yvec.end(), 0.0) / double(n);
                double sq = 0.0;
                for (size_t ii = 0; ii < n; ++ii) sq += Yvec[ii] * Yvec[ii];
                stdY = std::sqrt(sq / double(n) - meanY*meanY);
                if (stdY < EPS) stdY = 1.0;
                for (size_t ii = 0; ii < n; ++ii) Yvec[ii] = (Yvec[ii] - meanY) / stdY;
            }

            Eigen::Map<const Eigen::VectorXd> Yv_map(Yvec.data(), n);

            // ---------- Weighted / unweighted OLS ----------
            if (!robust && wgt_global)
            {
                Eigen::Map<const Eigen::VectorXd> Wv(wgt_adj.data(), n);
                Eigen::MatrixXd W = Wv.asDiagonal();
                coeff_vec = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv_map);
            }
            else
            {
                coeff_vec = pinv * Yv_map;
            }
            coeff(0) = coeff_vec(0); coeff(1) = coeff_vec(1);

            // ---------- Robust IRLS ----------
            if (robust)
            {
                Eigen::VectorXd wts = Eigen::VectorXd::Ones(n);
                if (wgt_global)
                {
                    for (size_t ii = 0; ii < n; ++ii)
                        wts(ii) = wgt_adj[ii] == 0.0 ? 0.0 : wgt_adj[ii];
                }
                Eigen::VectorXd prev_coeff = coeff;
                for (int iter = 0; iter < max_iter; ++iter)
                {
                    Eigen::MatrixXd W = wts.asDiagonal();
                    Eigen::Vector2d c = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv_map);
                    if ((c - prev_coeff).norm() < tol) { coeff = c; break; }
                    prev_coeff = c;
                    Eigen::VectorXd residuals = Yv_map - A * c;
                    for (size_t ii = 0; ii < n; ++ii)
                        wts(ii) = (std::abs(residuals(ii)) <= huber_tune) ? 1.0 : huber_tune / std::max(std::abs(residuals(ii)), EPS);
                    coeff = c;
                }
                for (size_t ii = 0; ii < n; ++ii) wgt_adj[ii] = wts(ii);
            }

            // ---------- Compute fitted values ----------
            for (size_t ii = 0; ii < n; ++ii)
            {
                double xval = (standardize ? ((X_all[ii] - meanX)/stdX) : X_all[ii]);
                fitted[ii] = coeff(0) + coeff(1) * xval;
            }

            // ---------- Compute statistics ----------
            params.AIC = computeAIC(Yvec, fitted, 2, wgt_global ? &wgt_adj : nullptr);
            params.MASE = MASE(Yvec, fitted, wgt_global ? &wgt_adj : nullptr);
            params.R2 = computeR2(Yvec, fitted, wgt_global ? &wgt_adj : nullptr);
            params.chi2 = computeChi2(Yvec, fitted, wgt_global ? &wgt_adj : nullptr)/(n-2);

            // ---------- De-standardize ----------
            double slope = coeff(1), intercept = coeff(0);
            if (standardize)
            {
                double devyoverdevx = stdY / stdX;
                slope = slope * devyoverdevx;
                intercept = meanY - slope * meanX + intercept;
                for (size_t ii = 0; ii < n; ++ii) fitted[ii] = fitted[ii] * stdY + meanY;
            }

            // ---------- Fill MTGAParameters ----------
            params.Ki = slope;
            params.Intercept = intercept;
            params.x = X_all;
            params.y.resize(n);
            for (size_t ii = 0; ii < n; ++ii) params.y[ii] = voxels[v][keepIndex[ii]] * invCp[keepIndex[ii]];
            params.fitted = fitted;
            params.frame.resize(n); for (size_t ii = 0; ii < n; ++ii) params.frame[ii] = keepIndex[ii] + 1;
            params.dof = 2;
            params.weights = (wgt_global ? std::vector<double>(wgt_adj.begin(), wgt_adj.begin()+n) : std::vector<double>(n, 1.0));
            params.r.resize(n); for (size_t ii = 0; ii < n; ++ii) params.r[ii] = params.y[ii] - params.fitted[ii];

            outputParams[v] = std::move(params);

            // ---------- Update progress bar safely ----------
            if (progressBar)
            {
                size_t done = ++voxProcessed;

                size_t updateInterval = std::max(static_cast<size_t>(1), std::min(nVoxels / 1000, static_cast<size_t>(10000)));

                if (done % updateInterval == 0 || done == nVoxels)
                {
                    int progress = static_cast<int>(100 * done / nVoxels);
                    if (QThread::currentThread() == progressBar->thread())
                    {
                        progressBar->setValue(progress);
                        // qApp->processEvents();
                    } else
                    {
                        QMetaObject::invokeMethod(progressBar, "setValue", Qt::QueuedConnection, Q_ARG(int, progress));
                    }
                }
            }

        }
    }
  if (progressBar)
  {
    progressBar->setVisible(false);
    stopButton->setVisible(false);
    // stopButton->deleteLater();
    if (stopRequested) {
      outputParams.clear();
    }
    qApp->processEvents();
  }
}

void vtkSlicerKMAPLogic::Logan4Img(
    const std::vector<std::vector<double>>& voxels, // TACs, [nVoxels][N]
    const std::vector<double>& Cp,
    const std::vector<double>& framing,
    const std::vector<double>* wgt_global, // nullptr if not used
    double timeOffset,
    double framingNorm,
    bool robust,
    bool standardize,
    double huber_tune,
    double tol,
    int max_iter,
    std::vector<MTGAParameters>& outputParams,
    std::atomic<bool>& stopRequested,
    QProgressBar* progressBar,
    int numThreads,
    QPushButton* stopButton
)
{
    #ifdef HAVE_OPENMP
    int max_hw_threads = omp_get_num_procs();
    if (numThreads > 0) {
        int n = std::min(numThreads, max_hw_threads);
        omp_set_num_threads(n);
    }
    #endif

    const double EPS = 1e-12;
    size_t N = Cp.size();
    if (framing.size() != N) return;
    size_t nVoxels = voxels.size();
    if (nVoxels == 0) return;

    // ---------- 1) Precompute frameScaled, timeAlong, intCp ----------
    std::vector<double> frameScaled(N);
    for (size_t i = 0; i < N; ++i) frameScaled[i] = framing[i] / framingNorm;

    std::vector<double> timeAlong(N);
    std::partial_sum(frameScaled.begin(), frameScaled.end(), timeAlong.begin());

    std::vector<double> intCp(N, 0.0);
    double acc = 0.0;
    for (size_t i = 1; i < N; ++i)
    {
        acc += 0.5 * (Cp[i] + Cp[i-1]) * frameScaled[i];
        intCp[i] = acc;
    }

    // ---------- 2) Determine frames >= timeOffset ----------
    std::vector<int> keepIndex;
    for (size_t i = 0; i < N; ++i)
    {
        if (timeAlong[i] >= timeOffset)
            keepIndex.push_back(static_cast<int>(i));
    }
    size_t n = keepIndex.size();
    if (n < 2) return;

    // ---------- 3) Prepare output ----------
    outputParams.clear();
    outputParams.resize(nVoxels);

    if (progressBar)
    {
        progressBar->setFormat("Fitting Logan (%p%)");
        progressBar->setVisible(true);
        progressBar->setMinimum(0);
        progressBar->setMaximum(100);
        progressBar->setValue(0);
        stopButton->setVisible(true);
        stopButton->show();
        qApp->processEvents();
    }

    std::atomic<size_t> voxProcessed(0);

    // ---------- 4) Parallel voxel-wise processing ----------
    #ifdef HAVE_OPENMP
    #pragma omp parallel
    #endif
    {
        std::vector<double> Xvec(n), Yvec(n), fitted(n), wgt_adj(n);
        Eigen::Vector2d coeff;
        Eigen::VectorXd coeff_vec(2);

        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int v = 0; v < nVoxels; ++v)
        {
            if (stopRequested) continue;
            MTGAParameters params;
            const std::vector<double>& tac = voxels[v];

            // ---------- Build X, Y ----------
            std::vector<double> intCt(N, 0.0);
            double accCt = 0.0;
            for (size_t i = 1; i < N; ++i)
            {
                accCt += 0.5 * (tac[i] + tac[i-1]) * frameScaled[i];
                intCt[i] = accCt;
            }

            for (size_t ii = 0; ii < n; ++ii)
            {
                int idx = keepIndex[ii];
                double denom = tac[idx] + EPS;
                Xvec[ii] = intCp[idx] / denom;
                Yvec[ii] = intCt[idx] / denom;
                if (wgt_global) wgt_adj[ii] = (*wgt_global)[idx];
            }

            // ---------- Standardization ----------
            double meanX=0.0, meanY=0.0, stdX=1.0, stdY=1.0;
            if (standardize)
            {
                meanX = std::accumulate(Xvec.begin(), Xvec.end(), 0.0) / double(n);
                meanY = std::accumulate(Yvec.begin(), Yvec.end(), 0.0) / double(n);
                double sqX=0.0, sqY=0.0;
                for (size_t ii=0; ii<n; ++ii){sqX += Xvec[ii]*Xvec[ii]; sqY += Yvec[ii]*Yvec[ii];}
                stdX = std::sqrt(sqX/double(n) - meanX*meanX);
                stdY = std::sqrt(sqY/double(n) - meanY*meanY);
                if (stdX < EPS) stdX = 1.0;
                if (stdY < EPS) stdY = 1.0;
                for (size_t ii=0; ii<n; ++ii){Xvec[ii] = (Xvec[ii]-meanX)/stdX; Yvec[ii] = (Yvec[ii]-meanY)/stdY;}
            }


            // ---------- Weighted / unweighted OLS ----------
            Eigen::MatrixXd A(n, 2);
            A.col(0) = Eigen::VectorXd::Ones(n);
            for (size_t i = 0; i < n; ++i) A(i,1) = Xvec[i];

            Eigen::Matrix2d AtA = A.transpose() * A;
            Eigen::Matrix2d AtA_inv = AtA.ldlt().solve(Eigen::Matrix2d::Identity());
            Eigen::MatrixXd pinv = AtA_inv * A.transpose(); // 2 x n

            Eigen::Map<const Eigen::VectorXd> Yv_map(Yvec.data(), n);
            if (!robust && wgt_global)
            {
                Eigen::Map<const Eigen::VectorXd> Wv(wgt_adj.data(), n);
                Eigen::MatrixXd W = Wv.asDiagonal();
                coeff_vec = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv_map);
            }
            else
            {
                coeff_vec = pinv * Yv_map;
            }
            coeff(0) = coeff_vec(0); coeff(1) = coeff_vec(1);


            if (robust)
            {
                Eigen::VectorXd wts = Eigen::VectorXd::Ones(n);
                if (wgt_global)
                {
                    for (size_t ii=0; ii<n; ++ii) wts(ii) = (wgt_adj[ii]==0.0 ? 0.0 : wgt_adj[ii]);
                }
                Eigen::VectorXd prev_coeff = coeff;
                for (int iter = 0; iter < max_iter; ++iter)
                {
                    Eigen::MatrixXd W = wts.asDiagonal();
                    Eigen::Vector2d c = (A.transpose() * W * A).ldlt().solve(A.transpose() * W * Yv_map);
                    if ((c - prev_coeff).norm() < tol) { coeff = c; break; }
                    prev_coeff = c;
                    Eigen::VectorXd residuals = Yv_map - A * c;
                    for (size_t ii = 0; ii < n; ++ii)
                        wts(ii) = (std::abs(residuals(ii)) <= huber_tune) ? 1.0 : huber_tune / std::max(std::abs(residuals(ii)), EPS);
                    coeff = c;
                }
                for (size_t ii = 0; ii < n; ++ii) wgt_adj[ii] = wts(ii);
            }

            // ---------- Compute fitted values ----------
            for (size_t ii=0; ii<n; ++ii) fitted[ii] = coeff(0) + coeff(1)*Xvec[ii];

            // ---------- Statistics ----------
            params.AIC = computeAIC(Yvec, fitted, 2, wgt_global ? &wgt_adj : nullptr);
            params.MASE = MASE(Yvec, fitted, wgt_global ? &wgt_adj : nullptr);
            params.R2 = computeR2(Yvec, fitted, wgt_global ? &wgt_adj : nullptr);
            params.chi2 = computeChi2(Yvec, fitted, wgt_global ? &wgt_adj : nullptr)/(n-2);

            // ---------- De-standardize ----------
            double slope = coeff(1), intercept = coeff(0);
            if (standardize)
            {
                double devyoverdevx = stdY / stdX;
                slope = slope * devyoverdevx;
                intercept = meanY - slope * meanX + intercept;
                for (size_t ii=0; ii<n; ++ii)
                {
                    fitted[ii] = fitted[ii]*stdY + meanY;
                    Xvec[ii] = Xvec[ii]*stdX + meanX;
                    Yvec[ii] = Yvec[ii]*stdY + meanY;
                }
            }

            // ---------- Fill MTGAParameters ----------
            params.DV = slope;
            params.Intercept = intercept;
            params.x.resize(n); for (size_t ii=0; ii<n; ++ii) params.x[ii] = Xvec[ii];
            params.y.resize(n); for (size_t ii=0; ii<n; ++ii) params.y[ii] = Yvec[ii];
            params.fitted = fitted;
            params.frame.resize(n); for (size_t ii=0; ii<n; ++ii) params.frame[ii] = keepIndex[ii]+1;
            params.dof = 2;
            params.weights = (wgt_global ? std::vector<double>(wgt_adj.begin(), wgt_adj.begin()+n) : std::vector<double>(n,1.0));
            params.r.resize(n); for (size_t ii=0; ii<n; ++ii) params.r[ii] = params.y[ii]-params.fitted[ii];

            outputParams[v] = std::move(params);

            // ---------- Progress bar ----------
            if (progressBar)
            {
                size_t done = ++voxProcessed;

                size_t updateInterval = std::max(static_cast<size_t>(1), std::min(nVoxels / 1000, static_cast<size_t>(10000)));

                if (done % updateInterval == 0 || done == nVoxels)
                {
                    int progress = static_cast<int>(100 * done / nVoxels);
                    if (QThread::currentThread() == progressBar->thread())
                    {
                        progressBar->setValue(progress);
                        // qApp->processEvents();
                    } else
                    {
                        QMetaObject::invokeMethod(progressBar, "setValue", Qt::QueuedConnection, Q_ARG(int, progress));
                    }
                }
            }
        }
    }

    if (progressBar)
    {
        progressBar->setVisible(false);
        stopButton->setVisible(false);
        if (stopRequested) {
          outputParams.clear();
        }
        qApp->processEvents();
    }
}

void vtkSlicerKMAPLogic::RE4Img(
    const std::vector<std::vector<double>>& voxels, // TACs, [nVoxels][N]
    const std::vector<double>& Cp,
    const std::vector<double>& framing,
    const std::vector<double>* wgt_global, // nullptr if not used
    double timeOffset,
    double framingNorm,
    bool robust,
    bool standardize,
    double huber_tune,
    double tol,
    int max_iter,
    std::vector<MTGAParameters>& outputParams,
    std::atomic<bool>& stopRequested,
    QProgressBar* progressBar,
    int numThreads,
    QPushButton* stopButton
)
{
    #ifdef HAVE_OPENMP
    int max_hw_threads = omp_get_num_procs();
    if (numThreads > 0) {
        int n = std::min(numThreads, max_hw_threads);
        omp_set_num_threads(n);
    }
    #endif

    const double EPS = 1e-12;
    size_t N = Cp.size();
    if (framing.size() != N) return;
    size_t nVoxels = voxels.size();
    if (nVoxels == 0) return;

    // ---------- 1) Precompute frameScaled, timeAlong, intCp ----------
    std::vector<double> frameScaled(N);
    for (size_t i = 0; i < N; ++i) frameScaled[i] = framing[i] / framingNorm;

    std::vector<double> timeAlong(N);
    std::partial_sum(frameScaled.begin(), frameScaled.end(), timeAlong.begin());

    std::vector<double> intCp(N, 0.0);
    double acc = 0.0;
    for (size_t i = 1; i < N; ++i)
    {
        acc += 0.5 * (Cp[i] + Cp[i-1]) * frameScaled[i];
        intCp[i] = acc;
    }

    std::vector<double> invCp(N);
    for (size_t i = 0; i < N; ++i) invCp[i] = 1.0 / (Cp[i] + EPS);

    // ---------- 2) Determine frames >= timeOffset ----------
    std::vector<double> X_all;
    std::vector<int> keepIndex;
    for (size_t i = 0; i < N; ++i)
    {
      if (timeAlong[i] >= timeOffset)
      {
          X_all.push_back(intCp[i] * invCp[i]);
          keepIndex.push_back(static_cast<int>(i));
      }
    }
    size_t n = X_all.size();
    if (n < 2) return;

    // ---------- 3) Build design matrix A ----------
    Eigen::MatrixXd A(n, 2);
    A.col(0) = Eigen::VectorXd::Ones(n);
    for (size_t i = 0; i < n; ++i) A(i,1) = X_all[i];

    double meanX = 0.0, stdX = 1.0;
    if (standardize)
    {
        meanX = std::accumulate(X_all.begin(), X_all.end(), 0.0) / double(n);
        double sq = 0.0;
        for (size_t i = 0; i < n; ++i) sq += X_all[i]*X_all[i];
        stdX = std::sqrt(sq/double(n) - meanX*meanX);
        if (stdX < EPS) stdX = 1.0;
        for (size_t i = 0; i < n; ++i) A(i,1) = (X_all[i]-meanX)/stdX;
    }

    // Precompute pseudoinverse for unweighted OLS
    Eigen::Matrix2d AtA = A.transpose() * A;
    Eigen::Matrix2d AtA_inv = AtA.ldlt().solve(Eigen::Matrix2d::Identity());
    Eigen::MatrixXd pinv = AtA_inv * A.transpose(); // 2 x n

    // ---------- 4) Prepare output ----------
    outputParams.clear();
    outputParams.resize(nVoxels);

    if (progressBar)
    {
        progressBar->setFormat("Fitting RE (%p%)");
        progressBar->setVisible(true);
        progressBar->setMinimum(0);
        progressBar->setMaximum(100);
        progressBar->setValue(0);
        stopButton->setVisible(true);
        stopButton->show();
        qApp->processEvents();
    }

    std::atomic<size_t> voxProcessed(0);

    // ---------- 5) Parallel voxel-wise processing ----------
    #ifdef HAVE_OPENMP
    #pragma omp parallel
    #endif
    {
        std::vector<double> Yvec(n);
        std::vector<double> fitted(n);
        std::vector<double> wgt_adj(n);
        Eigen::Vector2d coeff;
        Eigen::VectorXd coeff_vec(2);

        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int v = 0; v < nVoxels; ++v)
        {
            if (stopRequested) continue;
            MTGAParameters params;
            const std::vector<double>& tac = voxels[v];

            // ---------- Compute intCt (voxel-dependent) ----------
            std::vector<double> intCt(N, 0.0);
            double accCt = 0.0;
            for (size_t i = 1; i < N; ++i)
            {
                accCt += 0.5*(tac[i]+tac[i-1])*frameScaled[i];
                intCt[i] = accCt;
            }


            // ---------- Build Yvec and weights ----------
            for (size_t ii = 0; ii < n; ++ii)
            {
                int idx = keepIndex[ii];
                Yvec[ii] = intCt[idx] * invCp[idx];
                if (wgt_global) wgt_adj[ii] = (*wgt_global)[idx];
            }

            // ---------- Standardization ----------
            double meanY=0.0, stdY=1.0;
            if (standardize)
            {
                meanY = std::accumulate(Yvec.begin(), Yvec.end(), 0.0)/double(n);
                double sq = 0.0;
                for (size_t ii=0; ii<n; ++ii) sq += Yvec[ii]*Yvec[ii];
                stdY = std::sqrt(sq/double(n) - meanY*meanY);
                if (stdY < EPS) stdY = 1.0;
                for (size_t ii=0; ii<n; ++ii) Yvec[ii] = (Yvec[ii]-meanY)/stdY;
            }

            Eigen::Map<const Eigen::VectorXd> Yv_map(Yvec.data(), n);

            // ---------- Weighted / unweighted OLS ----------
            if (!robust && wgt_global)
            {
                Eigen::Map<const Eigen::VectorXd> Wv(wgt_adj.data(), n);
                Eigen::MatrixXd W = Wv.asDiagonal();
                coeff_vec = (A.transpose()*W*A).ldlt().solve(A.transpose()*W*Yv_map);
            }
            else
            {
                coeff_vec = pinv * Yv_map;
            }
            coeff(0) = coeff_vec(0); coeff(1) = coeff_vec(1);

            // ---------- Robust regression ----------
            if (robust)
            {
                Eigen::VectorXd wts = Eigen::VectorXd::Ones(n);
                if (wgt_global)
                {
                    for (size_t ii=0; ii<n; ++ii) wts(ii) = (wgt_adj[ii]==0.0 ? 0.0 : wgt_adj[ii]);
                }
                Eigen::VectorXd prev_coeff = coeff;
                for (int iter=0; iter<max_iter; ++iter)
                {
                    Eigen::MatrixXd W = wts.asDiagonal();
                    Eigen::Vector2d c = (A.transpose()*W*A).ldlt().solve(A.transpose()*W*Yv_map);
                    if ((c - prev_coeff).norm() < tol) { coeff = c; break; }
                    prev_coeff = c;
                    Eigen::VectorXd residuals = Yv_map - A*c;
                    for (size_t ii=0; ii<n; ++ii)
                        wts(ii) = (std::abs(residuals(ii)) <= huber_tune) ? 1.0 : huber_tune / std::max(std::abs(residuals(ii)), EPS);
                    coeff = c;
                }
                for (size_t ii=0; ii<n; ++ii) wgt_adj[ii] = wts(ii);
            }

            // ---------- Compute fitted values ----------
            for (size_t ii=0; ii<n; ++ii)
            {
                double xval = standardize ? (X_all[ii]-meanX)/stdX : X_all[ii];
                fitted[ii] = coeff(0) + coeff(1)*xval;
            }

            // ---------- Statistics ----------
            params.AIC = computeAIC(Yvec,fitted,2,wgt_global?&wgt_adj:nullptr);
            params.MASE = MASE(Yvec,fitted,wgt_global?&wgt_adj:nullptr);
            params.R2 = computeR2(Yvec,fitted,wgt_global?&wgt_adj:nullptr);
            params.chi2 = computeChi2(Yvec,fitted,wgt_global?&wgt_adj:nullptr)/(n-2);

            // ---------- De-standardize ----------
            double slope = coeff(1), intercept = coeff(0);
            if (standardize)
            {
                double devyoverdevx = stdY/stdX;
                slope *= devyoverdevx;
                intercept = meanY - slope*meanX + intercept;
                for (size_t ii=0; ii<n; ++ii)
                    fitted[ii] = fitted[ii]*stdY + meanY;
            }

            // ---------- Fill MTGAParameters ----------
            params.DV = slope;
            params.Intercept = intercept;
            params.x = X_all;
            params.y.resize(n);
            for (size_t ii=0; ii<n; ++ii) params.y[ii] = intCt[keepIndex[ii]]*invCp[keepIndex[ii]];
            params.fitted = fitted;
            params.frame.resize(n); for (size_t ii=0; ii<n; ++ii) params.frame[ii] = keepIndex[ii]+1;
            params.dof = 2;
            params.weights = (wgt_global ? std::vector<double>(wgt_adj.begin(),wgt_adj.begin()+n) : std::vector<double>(n,1.0));
            params.r.resize(n); for (size_t ii=0; ii<n; ++ii) params.r[ii] = params.y[ii]-params.fitted[ii];

            outputParams[v] = std::move(params);

            // ---------- Progress bar ----------
            if (progressBar)
            {
                size_t done = ++voxProcessed;

                size_t updateInterval = std::max(static_cast<size_t>(1), std::min(nVoxels / 1000, static_cast<size_t>(10000)));

                if (done % updateInterval == 0 || done == nVoxels)
                {
                    int progress = static_cast<int>(100 * done / nVoxels);
                    if (QThread::currentThread() == progressBar->thread())
                    {
                        progressBar->setValue(progress);
                        // qApp->processEvents();
                    } else {
                        QMetaObject::invokeMethod(progressBar, "setValue", Qt::QueuedConnection, Q_ARG(int, progress));
                    }
                }
            }
        }
    }

    if (progressBar)
    {
        progressBar->setVisible(false);
        stopButton->setVisible(false);
        if (stopRequested) {
          outputParams.clear();
        }
        qApp->processEvents();
    }
}


std::vector<double> vtkSlicerKMAPLogic::ExtractParameter(
    const std::vector<MTGAParameters>& outputParams,
    const std::string& field)
{
  std::vector<double> values;
  values.reserve(outputParams.size());

  for (const auto& param : outputParams)
  {
    if (field == "Ki") values.push_back(param.Ki);
    else if (field == "DV") values.push_back(param.DV);
    else if (field == "Intercept") values.push_back(param.Intercept);
    else if (field == "AIC") values.push_back(param.AIC);
    else if (field == "MASE") values.push_back(param.MASE);
    else if (field == "R2") values.push_back(param.R2);
    else if (field == "chi2") values.push_back(param.chi2);
    else values.push_back(0.0); // fallback
  }
  return values;
}

void vtkSlicerKMAPLogic::CreateMTGAParametricImages(
    const std::vector<MTGAParameters>& outputParams,
    const int dims[3],
    const std::vector<std::string>& fields,
    const std::string& modelID,
    vtkMRMLScalarVolumeNode* refNode,
    vtkMRMLSubjectHierarchyNode* refSH,
    vtkIdType refID
  )
{

  for (const auto& field : fields)
  {
    std::vector<double> flatten = this->ExtractParameter(outputParams, field);
    vtkMRMLScalarVolumeNode* node = this->Flatten2Image(flatten, dims, modelID + " - " + field);
    node->CopyOrientation(refNode);
    node->SetSpacing(refNode->GetSpacing());
    node->SetOrigin(refNode->GetOrigin());
    vtkIdType parentItemID = refSH->GetItemParent(refID);
    vtkIdType newItemID = refSH->GetItemByDataNode(node);
    refSH->SetItemParent(newItemID, parentItemID);
    if (!node)
    {
      std :: cerr << "Could not create image for " << modelID + " - " + field << std :: endl;
      return;
    }
  }
}


static volatile double sink = 0.0;

void bench_exp_strict()
{
  const int N = 20'000'000;
  double x = -0.001;

  auto t0 = std::chrono::high_resolution_clock::now();

  double acc = 0.0;
  for (int i = 0; i < N; ++i) {
    // vary input to prevent constant folding and SIMD tricks
    x += 1e-9 * (i & 1023);
    double y = std::exp(x);
    acc += y;
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  sink = acc; // force side-effect

  double sec = std::chrono::duration<double>(t1 - t0).count();
}

void vtkSlicerKMAPLogic::callTCMImg(
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
    std::atomic<bool>& stopRequested /*= false*/,
    const std::vector<double>* wgt_global /*= nullptr*/,
    int numThreads /*= 0 */,
    std::function<void(int)> progressCallback, /*= nullptr*/
    std::function<bool()> stopCallback /*= nullptr*/
)
{
    bench_exp_strict();
    const double EPS = 1e-12;
    const int Nframe = static_cast<int>(Cp.size());
    const int Nvox = static_cast<int>(voxels.size());
    if (Nvox == 0 || framing.size() != Cp.size()) return;

    // ---------- choose threads ----------
    #ifdef HAVE_OPENMP
    int max_hw_threads = omp_get_num_procs();
    if (numThreads > 0) {
        int n = std::min(numThreads, max_hw_threads);
        omp_set_num_threads(n);
    }
    #endif


    int N = 10000000;
    std::vector<double> a(N, 1.234567);
    auto t0 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for schedule(static)
    for(int i=0;i<N;++i){
        double x = a[i];
        // heavy-ish FP work to force vectorization
        for (int k=0;k<100;++k) {
            x = x*1.0000001 + 0.1234567;
            x = x / 1.0000003 + x*0.000001;
        }
        a[i] = x;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double sec = std::chrono::duration<double>(t1-t0).count();
    std::cout << "threads="<<numThreads<<" time="<<sec<<" s\n";



    // ---------- 1) Weights ----------
    std::vector<double> wt(Nframe, 1.0);
    if (wgt_global) wt = *wgt_global;
    std::vector<bool> keep(Nframe);
    for (int i = 0; i < Nframe; ++i) keep[i] = (wt[i] != 0.0);

    // ---------- 2) Cumulative frame times ----------
    std::vector<double> cumsum(Nframe, 0.0);
    double cum = 0.0;
    for (int i = 0; i < Nframe; ++i) {
        cum += framing[i];
        cumsum[i] = cum;
    }

    // ---------- 3) Scant and whole blood ----------
    std::vector<double*> scant(Nframe);
    for (int i = 0; i < Nframe; ++i) {
      scant[i] = new double[2];
    }
    std::vector<double> cwb(Nframe);
    for (int i = 0; i < Nframe; ++i) {
        scant[i][0] = (i==0)?0.0:cumsum[i-1];
        scant[i][1] = cumsum[i];
        double t = 0.5*(scant[i][0]+scant[i][1]);
        double pbr = pbrp[0]*exp(-pbrp[1]*t/60.0)+pbrp[2];
        // std :: cout << "Cp[i] = " << Cp[i] << std :: endl;
        cwb[i] = Cp[i]/pbr;
    }

    // ---------- 4) Fine-sample ----------
    long int N_cp = 0;
    double* Cp_new  = finesample2(scant, Cp,  N_cp, timestep, "linear");
    // for (int i = 0; i < N_cp; ++i) {
    //   std :: cout << "Cp_new[i] = " << Cp_new[i] << std :: endl;
    // }
    double* cwb_new = finesample2(scant, cwb, N_cp, timestep, "linear");

    // ---------- 5) Preallocate output ----------
    const int num_par = (n_tc == 1) ? 4 : 6;
    const int dof_fixed = static_cast<int>(std::count(sens, sens + num_par, true));
    outputParams.clear();
    outputParams.resize(Nvox);

    // ---------- 5b) Remove exact-zero TACs from fitting ----------
    auto markInvalidVoxel = [&](TCMParameters& params)
    {
        params.vb = -1.0;
        params.K1 = -1.0;
        params.k2 = -1.0;

        if (n_tc == 2)
        {
            params.k3 = -1.0;
            params.k4 = -1.0;
        }

        params.td = -1.0;
        params.Ki = -1.0;
        params.DV = -1.0;

        params.r.clear();
        params.weights = wt;
        params.keep = keep;
        params.dof = dof_fixed;

        params.AIC    = std::numeric_limits<double>::quiet_NaN();
        params.BIC    = std::numeric_limits<double>::quiet_NaN();
        params.MASE   = std::numeric_limits<double>::quiet_NaN();
        params.chi2   = std::numeric_limits<double>::quiet_NaN();
        params.loglik = std::numeric_limits<double>::quiet_NaN();
    };

    std::vector<int> fitVoxelIndices;
    fitVoxelIndices.reserve(Nvox);

    int zeroVoxelCount = 0;

    for (int v = 0; v < Nvox; ++v)
    {
        const auto& tac = voxels[v];

        const bool allZero = std::all_of(
            tac.begin(),
            tac.end(),
            [](double value)
            {
                return value == 0.0;
            });

        if (allZero)
        {
            markInvalidVoxel(outputParams[v]);
            ++zeroVoxelCount;
        }
        else
        {
            fitVoxelIndices.push_back(v);
        }
    }

    const int Nfit = static_cast<int>(fitVoxelIndices.size());

    std::cout
        << "TCM fitting: "
        << Nfit << " / " << Nvox
        << " voxels will be fitted; "
        << zeroVoxelCount << " zero TACs skipped."
        << std::endl;


    // ---------- 6) KMODEL_T ----------
    double * scant_flatten = new double[Nframe*2];
    for (int i=0; i<Nframe; ++i) {
      for (int j=0; j<2; ++j) {
        scant_flatten[i + j * Nframe] = scant[i][j];
      }
    }

    KMODEL_T km_template;
    km_template.dk = dk;
    km_template.td = timestep;
    km_template.cp = Cp_new;
    km_template.wb = cwb_new;
    km_template.num_frm = Nframe;
    km_template.num_vox = 1;  // voxel-wise
    km_template.scant = scant_flatten;
    km_template.tacfunc = (n_tc==1)? kconv_1tcm_tac : kconv_2tcm_tac;
    km_template.jacfunc = (n_tc==1)? kconv_1tcm_jac : kconv_2tcm_jac;

    // ---------- 7) Progress tracking ----------
    std::atomic<int> voxProcessed(0);
    const int progressUpdateInterval = std::max(1, Nfit / 200);

    // ---------- 8) Parallel voxel loop with per-thread scratch buffers ----------
    t0 = std::chrono::high_resolution_clock::now();
    #ifdef HAVE_OPENMP
    #pragma omp parallel
    #endif
    {
        // thread-local scratch buffers (allocated once per thread)
        std::vector<double> cfit_local(Nframe);
        std::vector<double> pinit_local(num_par);
        std::vector<int>    psens_local(num_par);
        std::vector<double> fitted_local(Nframe);
        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int fitIndex = 0; fitIndex < Nfit; ++fitIndex) {
            // if ((stopCallback && stopCallback()) || stopRequested) {
            //   // v = Nvox;
            //   continue;
            // }
            const int v = fitVoxelIndices[fitIndex];

            TCMParameters params;

            for (int i = 0; i < Nframe; ++i) cfit_local[i] = voxels[v][i];
            if (check_if_constant(cfit_local)) {
              // std :: cout << "Constant" << std :: endl;
              markInvalidVoxel(outputParams[v]);
            } else {

              for (int i = 0; i < num_par; ++i) pinit_local[i] = kinit[i];

              for (int i = 0; i < num_par; ++i) psens_local[i] = static_cast<int>(sens[i]);

              LevmarStats stx{};

              // ---------- Fit voxel ----------
              KMODEL_T km = km_template;
              kmap_levmar_stats(
                  cfit_local.data(), wt.data(), Nframe,
                  pinit_local.data(), num_par,
                  &km, tac_eval, jac_eval,
                  lb, ub, psens_local.data(), maxiter,
                  fitted_local.data(),
                  &stx
              );

              if (v == 0 || v == 99 || v == 999 || v == 9999) {
                std::cout
                  << "v=" << v
                  << " loops=" << stx.n_loops
                  << " acc=" << stx.it_accept
                  << " rej=" << stx.it_reject
                  << " guard=" << stx.hit_guard
                  << " nan_ct=" << stx.nan_ct
                  << " nan_st=" << stx.nan_st
                  << " last_rho=" << stx.last_rho
                  << " last_mu=" << stx.last_mu
                  << "\n";
              }
              // ---------- Fill TCMParameters ----------
              if (n_tc==1) {
                  params.vb = pinit_local[0];
                  params.K1 = pinit_local[1];
                  params.k2 = pinit_local[2];
                  params.td = pinit_local[3];
                  params.Ki = params.K1;
                  params.DV = params.K1/(params.k2+EPS);
              } else if (n_tc==2) {
                  params.vb = pinit_local[0];
                  params.K1 = pinit_local[1];
                  params.k2 = pinit_local[2];
                  params.k3 = pinit_local[3];
                  params.k4 = pinit_local[4];
                  params.td = pinit_local[5];
                  params.Ki = params.K1*params.k3/(params.k2+params.k3+EPS);
                  params.DV = params.K1/(params.k2+EPS)*(1.0+params.k3/(params.k4+EPS));
              }

              // Residuals
              params.r.resize(Nframe);
              for (int i=0;i<Nframe;++i) params.r[i] = cfit_local[i]-fitted_local[i];

              params.weights = wt;
              params.keep = keep;
              params.dof = dof_fixed;
              //
              // // // Statistics
              params.AIC  = this->computeAIC(cfit_local, fitted_local, params.dof, &wt);
              params.BIC  = this->computeBIC(cfit_local, fitted_local, params.dof, &wt);
              params.MASE = this->MASE(cfit_local, fitted_local, &wt);
              params.chi2 = this->computeChi2(cfit_local, fitted_local, &wt)/(Nframe-params.dof);
              params.loglik = this->computeLogLik(cfit_local, fitted_local, &wt);

              outputParams[v] = std::move(params);
            }
            // ---------- Progress update ----------
            if (progressCallback)
            {
                int done = ++voxProcessed;

                if (done % progressUpdateInterval == 0 || done == Nfit)
                {
                    int progress =
                        static_cast<int>(100LL * done / Nfit);

                    progressCallback(progress);
                }
            }

        }
    }
    t1 = std::chrono::high_resolution_clock::now();

    sec = std::chrono::duration<double>(t1-t0).count();
    std::cout << "elapsed time per call "<<sec/10000<<" s\n";

    delete[] Cp_new;
    delete[] cwb_new;
}

std::vector<double> vtkSlicerKMAPLogic::ExtractParameter(
    const std::vector<TCMParameters>& outputParams,
    const std::string& field)
{
  std::vector<double> values;
  values.reserve(outputParams.size());

  for (const auto& param : outputParams)
  {
    if (field == "K1") values.push_back(param.K1);
    else if (field == "k2") values.push_back(param.k2);
    else if (field == "k3") values.push_back(param.k3);
    else if (field == "k4") values.push_back(param.k4);
    else if (field == "vb") values.push_back(param.vb);
    else if (field == "td") values.push_back(param.td);
    else if (field == "Ki") values.push_back(param.Ki);
    else if (field == "DV") values.push_back(param.DV);
    else if (field == "AIC") values.push_back(param.AIC);
    else if (field == "MASE") values.push_back(param.MASE);
    else if (field == "BIC") values.push_back(param.BIC);
    else if (field == "chi2") values.push_back(param.chi2);
    else if (field == "loglik") values.push_back(param.loglik);
    else values.push_back(0.0); // fallback
  }
  return values;
}

void vtkSlicerKMAPLogic::CreateTCMParametricImages(
    const std::vector<TCMParameters>& outputParams,
    const int dims[3],
    const std::vector<std::string>& fields,
    const std::string& modelID,
    vtkMRMLScalarVolumeNode* refNode,
    vtkMRMLSubjectHierarchyNode* refSH,
    vtkIdType refID)
{

    const vtkIdType expectedVoxels =
        static_cast<vtkIdType>(dims[0]) *
        static_cast<vtkIdType>(dims[1]) *
        static_cast<vtkIdType>(dims[2]);

    if (outputParams.size() != static_cast<size_t>(expectedVoxels))
    {
        std::cerr
            << "CreateTCMParametricImages: output size mismatch. "
            << "Expected " << expectedVoxels
            << ", got " << outputParams.size()
            << std::endl;
        return;
    }

    if (!refNode)
    {
        std::cerr << "CreateTCMParametricImages: refNode is null" << std::endl;
        return;
    }

    if (!refSH)
    {
        std::cerr << "CreateTCMParametricImages: refSH is null" << std::endl;
        return;
    }

    vtkMRMLScene* scene = this->GetMRMLScene();
    if (!scene)
    {
        std::cerr << "CreateTCMParametricImages: MRML scene is null" << std::endl;
        return;
    }

    vtkIdType refItemID =
        refSH->GetItemByDataNode(refNode);

    if (refItemID ==
        vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
    {
        std::cerr
            << "Reference PET node has no subject hierarchy item."
            << std::endl;
        return;
    }

    vtkIdType parentItemID =
        refSH->GetItemParent(refItemID);

    for (const auto& field : fields)
    {
        std::vector<double> flatten =
            this->ExtractParameter(outputParams, field);

        vtkMRMLScalarVolumeNode* node =
            this->Flatten2Image(
                flatten,
                dims,
                modelID + " - " + field);

        if (!node)
        {
            std::cerr
                << "Could not create image for "
                << modelID << " - " << field
                << std::endl;
            return;
        }

        node->CopyOrientation(refNode);

        node->SetSpacing(refNode->GetSpacing());

        node->SetOrigin(refNode->GetOrigin());

        const vtkIdType newItemID =
            refSH->GetItemByDataNode(node);

        if (newItemID ==
            vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
        {
            std::cerr
                << "Could not find subject hierarchy item for "
                << field << std::endl;
        }
        else
        {
            refSH->SetItemParent(newItemID, parentItemID);
        }

    }

}
