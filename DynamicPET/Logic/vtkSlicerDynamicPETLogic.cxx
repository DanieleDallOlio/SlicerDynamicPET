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

// DynamicPET Logic includes
#include "vtkSlicerDynamicPETLogic.h"
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

#include <vtkImageThreshold.h>
#include <vtkImageConnectivityFilter.h>
#include <vtkImageDilateErode3D.h>
#include <vtkMatrix4x4.h>
#include <vtkMRMLSegmentationDisplayNode.h>
#include <vtkMRMLTransformNode.h>
#include <vtkAlgorithmOutput.h>

#include <vtkImageOpenClose3D.h>
#include <vtkIdTypeArray.h>

#include <array>
#include <cmath>
#include <numeric>
#include <limits>

// define M_PI in case of Win
#ifdef _WIN32
    #ifndef M_PI
        #define M_PI 3.14159265358979323846
    #endif
#endif


namespace
{

std::vector<double> GaussianSmooth1D(
    const std::vector<double>& input,
    double sigma)
{
    if (input.empty() || sigma <= 0.0)
    {
        return input;
    }

    const int radius =
        std::max(
            1,
            static_cast<int>(
                std::ceil(3.0 * sigma)));

    std::vector<double> kernel(
        2 * radius + 1);

    double kernelSum = 0.0;

    for (int k = -radius;
         k <= radius;
         ++k)
    {
        const double x =
            static_cast<double>(k) /
            sigma;

        const double value =
            std::exp(-0.5 * x * x);

        kernel[k + radius] =
            value;

        kernelSum += value;
    }

    for (double& value : kernel)
    {
        value /= kernelSum;
    }

    std::vector<double> output(
        input.size(),
        0.0);

    const int n =
        static_cast<int>(
            input.size());

    for (int i = 0; i < n; ++i)
    {
        double value = 0.0;

        for (int k = -radius;
             k <= radius;
             ++k)
        {
            const int index =
                std::max(
                    0,
                    std::min(
                        n - 1,
                        i + k));

            value +=
                kernel[k + radius] *
                input[index];
        }

        output[i] = value;
    }

    return output;
}


std::vector<int> FindSignificantPeaks(
    const std::vector<double>& values,
    double minimumRelativeHeight)
{
    std::vector<int> peaks;

    if (values.size() < 3)
    {
        return peaks;
    }

    const double maximum =
        *std::max_element(
            values.begin(),
            values.end());

    if (!(maximum > 0.0))
    {
        return peaks;
    }

    const double minimumHeight =
        minimumRelativeHeight *
        maximum;

    for (size_t i = 1;
         i + 1 < values.size();
         ++i)
    {
        if (values[i] >= minimumHeight &&
            values[i] > values[i - 1] &&
            values[i] >= values[i + 1])
        {
            peaks.push_back(
                static_cast<int>(i));
        }
    }

    return peaks;
}


std::vector<int> FindLocalMinima(
    const std::vector<double>& values)
{
    std::vector<int> minima;

    for (size_t i = 1;
         i + 1 < values.size();
         ++i)
    {
        if (values[i] <= values[i - 1] &&
            values[i] < values[i + 1])
        {
            minima.push_back(
                static_cast<int>(i));
        }
    }

    return minima;
}


int FindValleyBetween(
    const std::vector<double>& values,
    int leftPeak,
    int rightPeak)
{
    if (leftPeak < 0 ||
        rightPeak <= leftPeak + 1 ||
        rightPeak >=
            static_cast<int>(values.size()))
    {
        return -1;
    }

    int bestIndex =
        leftPeak + 1;

    double bestValue =
        values[bestIndex];

    for (int i = leftPeak + 2;
         i < rightPeak;
         ++i)
    {
        if (values[i] < bestValue)
        {
            bestValue =
                values[i];

            bestIndex =
                i;
        }
    }

    return bestIndex;
}


int ComputeOtsuHistogramThreshold(
    const std::vector<double>& histogram)
{
    if (histogram.empty())
    {
        return -1;
    }

    double total = 0.0;
    double totalWeightedIndex = 0.0;

    for (size_t i = 0;
         i < histogram.size();
         ++i)
    {
        total += histogram[i];

        totalWeightedIndex +=
            static_cast<double>(i) *
            histogram[i];
    }

    if (!(total > 0.0))
    {
        return -1;
    }

    double backgroundWeight = 0.0;
    double backgroundWeightedIndex = 0.0;

    double bestVariance = -1.0;
    int bestThreshold = -1;

    for (size_t i = 0;
         i + 1 < histogram.size();
         ++i)
    {
        backgroundWeight +=
            histogram[i];

        if (backgroundWeight <= 0.0)
        {
            continue;
        }

        const double foregroundWeight =
            total -
            backgroundWeight;

        if (foregroundWeight <= 0.0)
        {
            break;
        }

        backgroundWeightedIndex +=
            static_cast<double>(i) *
            histogram[i];

        const double backgroundMean =
            backgroundWeightedIndex /
            backgroundWeight;

        const double foregroundMean =
            (totalWeightedIndex -
             backgroundWeightedIndex) /
            foregroundWeight;

        const double delta =
            backgroundMean -
            foregroundMean;

        const double betweenClassVariance =
            backgroundWeight *
            foregroundWeight *
            delta *
            delta;

        if (betweenClassVariance >
            bestVariance)
        {
            bestVariance =
                betweenClassVariance;

            bestThreshold =
                static_cast<int>(i);
        }
    }

    return bestThreshold;
}


bool ComputeMultiscaleLogPETThreshold(
    const std::vector<double>& composite,
    double& threshold,
    bool& usedOtsuFallback)
{
    constexpr int numberOfBins = 512;

    // This small peak-height criterion prevents tiny histogram
    // fluctuations from becoming anatomical modes.
    constexpr double minimumPeakFraction = 0.005;

    const std::array<double, 5> sigmaLevels =
    {
        16.0,
        8.0,
        4.0,
        2.0,
        1.0
    };

    double minimumLog =
        std::numeric_limits<double>::infinity();

    double maximumLog =
        -std::numeric_limits<double>::infinity();

    size_t positiveFiniteCount = 0;

    // First pass: determine log-intensity range.
    for (double value : composite)
    {
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            continue;
        }

        const double logValue =
            std::log(value);

        minimumLog =
            std::min(
                minimumLog,
                logValue);

        maximumLog =
            std::max(
                maximumLog,
                logValue);

        ++positiveFiniteCount;
    }

    if (positiveFiniteCount < 2 ||
        !std::isfinite(minimumLog) ||
        !std::isfinite(maximumLog) ||
        maximumLog <= minimumLog)
    {
        return false;
    }

    const double binWidth =
        (maximumLog - minimumLog) /
        static_cast<double>(
            numberOfBins);

    if (!(binWidth > 0.0))
    {
        return false;
    }

    std::vector<double> histogram(
        numberOfBins,
        0.0);

    // Second pass: histogram of log-positive composite values.
    for (double value : composite)
    {
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            continue;
        }

        const double logValue =
            std::log(value);

        int bin =
            static_cast<int>(
                (logValue - minimumLog) /
                binWidth);

        bin =
            std::max(
                0,
                std::min(
                    numberOfBins - 1,
                    bin));

        histogram[bin] += 1.0;
    }

    // Always prepare a valid fallback.
    const int otsuBin =
        ComputeOtsuHistogramThreshold(
            histogram);

    if (otsuBin < 0)
    {
        return false;
    }

    const double otsuLogThreshold =
        minimumLog +
        (static_cast<double>(otsuBin) +
         0.5) *
        binWidth;

    int trackedValley = -1;

    std::vector<int>
        trackedValleys;

    for (double sigma :
         sigmaLevels)
    {
        const std::vector<double> smoothed =
            GaussianSmooth1D(
                histogram,
                sigma);

        const std::vector<int> peaks =
            FindSignificantPeaks(
                smoothed,
                minimumPeakFraction);

        if (trackedValley < 0)
        {
            if (peaks.size() < 2)
            {
                continue;
            }

            // The leftmost persistent mode is treated as
            // reconstructed background, the next as patient signal.
            const int backgroundPeak =
                peaks[0];

            const int tissuePeak =
                peaks[1];

            const int candidate =
                FindValleyBetween(
                    smoothed,
                    backgroundPeak,
                    tissuePeak);

            if (candidate < 0)
            {
                continue;
            }

            // Reject an almost-flat "valley".
            const double surroundingPeakHeight =
                std::min(
                    smoothed[backgroundPeak],
                    smoothed[tissuePeak]);

            if (!(smoothed[candidate] <
                  0.95 *
                  surroundingPeakHeight))
            {
                continue;
            }

            trackedValley =
                candidate;

            trackedValleys.push_back(
                candidate);

            continue;
        }

        // Follow the valley trajectory as the scale becomes finer.
        const std::vector<int> minima =
            FindLocalMinima(
                smoothed);

        int bestCandidate = -1;
        int bestDistance =
            std::numeric_limits<int>::max();

        const int maximumShift =
            std::max(
                4,
                static_cast<int>(
                    std::ceil(
                        2.5 * sigma)));

        for (int candidate : minima)
        {
            const int distance =
                std::abs(
                    candidate -
                    trackedValley);

            if (distance <=
                    maximumShift &&
                distance <
                    bestDistance)
            {
                bestDistance =
                    distance;

                bestCandidate =
                    candidate;
            }
        }

        if (bestCandidate >= 0)
        {
            trackedValley =
                bestCandidate;

            trackedValleys.push_back(
                bestCandidate);
        }
    }

    // Require persistence over several scales.
    if (trackedValleys.size() >= 3)
    {
        std::sort(
            trackedValleys.begin(),
            trackedValleys.end());

        const int medianBin =
            trackedValleys[
                trackedValleys.size() / 2];

        const double logThreshold =
            minimumLog +
            (static_cast<double>(medianBin) +
             0.5) *
            binWidth;

        threshold =
            std::exp(logThreshold);

        usedOtsuFallback =
            false;

        return std::isfinite(threshold) &&
               threshold > 0.0;
    }

    // No sufficiently stable scale-space trajectory.
    threshold =
        std::exp(
            otsuLogThreshold);

    usedOtsuFallback =
        true;

    return std::isfinite(threshold) &&
           threshold > 0.0;
}

} // end anonymous namespace

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
vtkStandardNewMacro(vtkSlicerDynamicPETLogic);

//----------------------------------------------------------------------------
vtkSlicerDynamicPETLogic::vtkSlicerDynamicPETLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerDynamicPETLogic::~vtkSlicerDynamicPETLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerDynamicPETLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerDynamicPETLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerDynamicPETLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerDynamicPETLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerDynamicPETLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerDynamicPETLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}


void vtkSlicerDynamicPETLogic::computeTAC(vtkIdType ctID,
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



void vtkSlicerDynamicPETLogic::setupSeg(vtkMRMLSegmentationNode* segNode)
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


VoxelStatistics vtkSlicerDynamicPETLogic::ComputeVoxelStatistics(vtkMRMLScalarVolumeNode* PETVolume, vtkImageData* labelmap, int labelValue)
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

double vtkSlicerDynamicPETLogic::boundaryLRTPvalue(double LR, int r_b, int r_i)
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

ModelComparisonResult vtkSlicerDynamicPETLogic::compareModels(
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

void vtkSlicerDynamicPETLogic::TAC(vtkMRMLSequenceNode* sequencePETNode,
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

double vtkSlicerDynamicPETLogic::computeLogLik(const std::vector<double>& y,
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


void vtkSlicerDynamicPETLogic::callTCM(std :: vector< std :: vector<double> > tac,
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



// void vtkSlicerDynamicPETLogic::getFittedTCM(double *& fitted_curve,
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

double vtkSlicerDynamicPETLogic::computeAIC(const std::vector<double>& obs,
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
double vtkSlicerDynamicPETLogic::MASE(const std::vector<double>& Actual,
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

double vtkSlicerDynamicPETLogic::computeBIC(const std::vector<double>& obs,
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

double vtkSlicerDynamicPETLogic::computeR2(const std::vector<double>& obs,
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

double vtkSlicerDynamicPETLogic::computeChi2(const std::vector<double>& y,
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

double vtkSlicerDynamicPETLogic::computeVuongP(const std::vector<double>& r1,
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



double vtkSlicerDynamicPETLogic::computeLRTP(double logLik1,
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

void vtkSlicerDynamicPETLogic::Patlak(const std::vector<double>& tac,
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

void vtkSlicerDynamicPETLogic::Logan(const std::vector<double>& tac,
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

void vtkSlicerDynamicPETLogic::RE(const std::vector<double>& tac,
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

vtkSmartPointer<vtkOrientedImageData>
vtkSlicerDynamicPETLogic::
CreateFullPETSupportMask(
    vtkMRMLScalarVolumeNode* referencePETNode)
{
    if (!referencePETNode)
    {
        return nullptr;
    }

    vtkMRMLTransformNode* referenceParentTransform =
        referencePETNode->
            GetParentTransformNode();

    vtkSmartPointer<vtkOrientedImageData> petOriented =
        vtkSmartPointer<vtkOrientedImageData>::Take(
            vtkSlicerSegmentationsModuleLogic::
                CreateOrientedImageDataFromVolumeNode(
                    referencePETNode,
                    referenceParentTransform));

    if (!petOriented)
    {
        return nullptr;
    }

    vtkSmartPointer<vtkOrientedImageData> mask =
        vtkSmartPointer<vtkOrientedImageData>::New();

    mask->SetExtent(
        petOriented->GetExtent());

    mask->AllocateScalars(
        VTK_UNSIGNED_CHAR,
        1);

    vtkNew<vtkMatrix4x4> imageToWorld;

    petOriented->GetImageToWorldMatrix(
        imageToWorld);

    mask->SetImageToWorldMatrix(
        imageToWorld);

    const vtkIdType numberOfVoxels =
        mask->GetNumberOfPoints();

    unsigned char* ptr =
        static_cast<unsigned char*>(
            mask->GetScalarPointer());

    if (!ptr)
    {
        return nullptr;
    }

    std::fill(
        ptr,
        ptr + numberOfVoxels,
        static_cast<unsigned char>(1));

    mask->Modified();

    return mask;
}

vtkSmartPointer<vtkOrientedImageData>
vtkSlicerDynamicPETLogic::CreateCTBodySupportMask(
    vtkMRMLScalarVolumeNode* ctNode,
    vtkMRMLScalarVolumeNode* referencePETNode,
    double ctThresholdHU,
    double bodyMarginMm,
    bool fillHoles)
{
    if (!ctNode || !referencePETNode)
    {
        std::cerr
            << "CreateCTBodySupportMask: missing CT or PET node."
            << std::endl;
        return nullptr;
    }

    // ---------------------------------------------------------------------
    // 1. Convert both MRML volumes to oriented image data in the same
    //    parent-transform coordinate system.
    // ---------------------------------------------------------------------

    vtkMRMLTransformNode* referenceParentTransform =
        referencePETNode->GetParentTransformNode();

    vtkSmartPointer<vtkOrientedImageData> ctOriented =
        vtkSmartPointer<vtkOrientedImageData>::Take(
            vtkSlicerSegmentationsModuleLogic::
                CreateOrientedImageDataFromVolumeNode(
                    ctNode,
                    referenceParentTransform));

    vtkSmartPointer<vtkOrientedImageData> petOriented =
        vtkSmartPointer<vtkOrientedImageData>::Take(
            vtkSlicerSegmentationsModuleLogic::
                CreateOrientedImageDataFromVolumeNode(
                    referencePETNode,
                    referenceParentTransform));

    if (!ctOriented || !petOriented)
    {
        std::cerr
            << "CreateCTBodySupportMask: "
               "could not create oriented image data."
            << std::endl;
        return nullptr;
    }

    // ---------------------------------------------------------------------
    // 2. Resample CT onto the PET voxel grid.
    //
    // CT is continuous data, therefore use linear interpolation.
    //
    // IMPORTANT: background is -1000 HU, not 0 HU.
    // Using 0 would accidentally classify areas outside the CT FOV as body
    // when applying a -500 HU threshold.
    // ---------------------------------------------------------------------

    vtkNew<vtkOrientedImageData> ctOnPETGrid;

    const bool resampleOK =
        vtkOrientedImageDataResample::
            ResampleOrientedImageToReferenceOrientedImage(
                ctOriented,
                petOriented,
                ctOnPETGrid,
                true,       // linear interpolation
                false,      // no padding beyond PET geometry
                nullptr,
                -1000.0);   // outside-CT value in HU

    if (!resampleOK)
    {
        std::cerr
            << "CreateCTBodySupportMask: CT-to-PET resampling failed."
            << std::endl;
        return nullptr;
    }

    // ---------------------------------------------------------------------
    // 3. CT body candidate:
    //
    //       HU >= threshold -> 1
    //       otherwise       -> 0
    // ---------------------------------------------------------------------

    vtkNew<vtkImageThreshold> threshold;

    threshold->SetInputData(ctOnPETGrid);
    threshold->ThresholdByUpper(ctThresholdHU);

    threshold->SetInValue(1);
    threshold->SetOutValue(0);

    threshold->SetOutputScalarType(VTK_UNSIGNED_CHAR);

    threshold->ReplaceInOn();
    threshold->ReplaceOutOn();

    threshold->Update();

    // ---------------------------------------------------------------------
    // 4. Keep the dominant connected body component.
    //
    // This removes most disconnected table/hardware/background structures.
    //
    // Note: if a table is physically connected to the body mask then a pure
    // connected-component approach cannot separate it. For our conservative
    // fitting-support purpose, keeping a few such false-positive voxels is
    // preferable to removing real anatomy.
    // ---------------------------------------------------------------------

    vtkNew<vtkImageConnectivityFilter> bodyConnectivity;

    bodyConnectivity->SetInputConnection(
        threshold->GetOutputPort());

    bodyConnectivity->SetScalarRange(
        1.0,
        1.0);

    bodyConnectivity->
        SetExtractionModeToLargestRegion();

    bodyConnectivity->
        SetLabelModeToConstantValue();

    bodyConnectivity->
        SetLabelConstantValue(1);

    bodyConnectivity->Update();

    vtkAlgorithmOutput* currentMaskPort =
        bodyConnectivity->GetOutputPort();

    // ---------------------------------------------------------------------
    // 5. Fill enclosed holes.
    //
    // Do this using connectivity rather than hand-written flood filling:
    //
    // body         = 1
    // background   = 0
    //
    // invert -> find largest connected background (outside air)
    //        -> invert that result
    //
    // Any background cavity that is NOT connected to external air becomes
    // foreground, thereby recovering lung/internal cavities.
    // ---------------------------------------------------------------------

    vtkNew<vtkImageThreshold> invertBody;

    vtkNew<vtkImageConnectivityFilter> exteriorBackground;

    vtkNew<vtkImageThreshold> filledBody;

    if (fillHoles)
    {
        invertBody->SetInputConnection(
            bodyConnectivity->GetOutputPort());

        // Original zero-valued background becomes one.
        invertBody->ThresholdByLower(0.0);
        invertBody->SetInValue(1);
        invertBody->SetOutValue(0);
        invertBody->SetOutputScalarType(
            VTK_UNSIGNED_CHAR);
        invertBody->ReplaceInOn();
        invertBody->ReplaceOutOn();

        exteriorBackground->SetInputConnection(
            invertBody->GetOutputPort());

        exteriorBackground->SetScalarRange(
            1.0,
            1.0);

        exteriorBackground->
            SetExtractionModeToLargestRegion();

        exteriorBackground->
            SetLabelModeToConstantValue();

        exteriorBackground->
            SetLabelConstantValue(1);

        // Invert exterior background:
        //
        // external air 1 -> 0
        // body          0 -> 1
        // internal hole 0 -> 1
        filledBody->SetInputConnection(
            exteriorBackground->GetOutputPort());

        filledBody->ThresholdByLower(0.0);
        filledBody->SetInValue(1);
        filledBody->SetOutValue(0);
        filledBody->SetOutputScalarType(
            VTK_UNSIGNED_CHAR);
        filledBody->ReplaceInOn();
        filledBody->ReplaceOutOn();

        filledBody->Update();

        currentMaskPort =
            filledBody->GetOutputPort();
    }

    // ---------------------------------------------------------------------
    // 6. Dilate by physical margin.
    //
    // vtkImageDilateErode3D uses voxel kernel dimensions, therefore convert
    // millimetres separately for X/Y/Z.
    // ---------------------------------------------------------------------

    vtkNew<vtkImageDilateErode3D> dilate;

    if (bodyMarginMm > 0.0)
    {
        const double* spacing =
            petOriented->GetSpacing();

        int radius[3] = {0, 0, 0};

        for (int axis = 0; axis < 3; ++axis)
        {
            const double safeSpacing =
                std::max(spacing[axis], 1e-6);

            radius[axis] =
                static_cast<int>(
                    std::ceil(
                        bodyMarginMm /
                        safeSpacing));
        }

        dilate->SetInputConnection(
            currentMaskPort);

        dilate->SetDilateValue(1);
        dilate->SetErodeValue(0);

        dilate->SetKernelSize(
            2 * radius[0] + 1,
            2 * radius[1] + 1,
            2 * radius[2] + 1);

        dilate->Update();

        currentMaskPort =
            dilate->GetOutputPort();
    }

    // ---------------------------------------------------------------------
    // 7. Convert back to vtkOrientedImageData and restore PET geometry.
    // ---------------------------------------------------------------------

    vtkAlgorithm* finalAlgorithm =
        currentMaskPort->GetProducer();

    finalAlgorithm->Update();

    vtkImageData* finalImage =
        vtkImageData::SafeDownCast(
            finalAlgorithm->GetOutputDataObject(0));

    if (!finalImage)
    {
        std::cerr
            << "CreateCTBodySupportMask: "
               "final binary image is invalid."
            << std::endl;
        return nullptr;
    }

    vtkSmartPointer<vtkOrientedImageData> result =
        vtkSmartPointer<vtkOrientedImageData>::New();

    result->DeepCopy(finalImage);

    vtkNew<vtkMatrix4x4> petImageToWorld;

    petOriented->GetImageToWorldMatrix(
        petImageToWorld);

    result->SetImageToWorldMatrix(
        petImageToWorld);

    return result;
}

vtkSmartPointer<vtkOrientedImageData>
vtkSlicerDynamicPETLogic::
CreatePETBodySupportMask(
    const std::vector<std::vector<double>>& voxels,
    const int dims[3],
    vtkMRMLScalarVolumeNode* referencePETNode,
    const std::vector<double>& frameDurations,
    PETCompositeMode compositeMode,
    double bodyMarginMm,
    bool fillHoles,
    double minComponentFractionOfLargest,
    double* thresholdOut,
    bool* usedOtsuFallbackOut)
{
    if (!referencePETNode ||
        voxels.empty())
    {
        return nullptr;
    }

    const size_t numberOfVoxels =
        voxels.size();

    const size_t numberOfFrames =
        voxels.front().size();

    const size_t expectedVoxels =
        static_cast<size_t>(dims[0]) *
        static_cast<size_t>(dims[1]) *
        static_cast<size_t>(dims[2]);

    if (numberOfVoxels != expectedVoxels ||
        numberOfFrames == 0)
    {
        std::cerr
            << "CreatePETBodySupportMask: "
               "invalid PET dimensions."
            << std::endl;

        return nullptr;
    }

    const bool durationWeighted =
        compositeMode ==
        PETCompositeMode::
            DurationWeightedSum;

    if (durationWeighted &&
        frameDurations.size() !=
            numberOfFrames)
    {
        std::cerr
            << "CreatePETBodySupportMask: "
               "duration vector does not match "
               "the number of PET frames."
            << std::endl;

        return nullptr;
    }

    // ------------------------------------------------------------------
    // 1. Composite PET.
    //
    // Maroy/Zbib default:
    //     sum_t PET(v,t)
    //
    // Optional SlicerDynamicPET variant:
    //     sum_t PET(v,t) * duration(t)
    // ------------------------------------------------------------------

    std::vector<double> composite(
        numberOfVoxels,
        0.0);

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (vtkIdType v = 0;
         v <
         static_cast<vtkIdType>(
             numberOfVoxels);
         ++v)
    {
        const auto& tac =
            voxels[
                static_cast<size_t>(v)];

        double sum = 0.0;

        const size_t usableFrames =
            std::min(
                numberOfFrames,
                tac.size());

        for (size_t frame = 0;
             frame < usableFrames;
             ++frame)
        {
            const double value =
                tac[frame];

            if (!std::isfinite(value))
            {
                continue;
            }

            if (durationWeighted)
            {
                sum +=
                    value *
                    frameDurations[frame];
            }
            else
            {
                sum += value;
            }
        }

        composite[
            static_cast<size_t>(v)] =
            sum;
    }

    // ------------------------------------------------------------------
    // 2. Maroy/Zbib-style multiscale threshold.
    //    Fall back to log-domain Otsu if scale-space tracking fails.
    // ------------------------------------------------------------------

    double threshold =
        std::numeric_limits<double>::
            quiet_NaN();

    bool usedOtsuFallback =
        false;

    if (!ComputeMultiscaleLogPETThreshold(
            composite,
            threshold,
            usedOtsuFallback))
    {
        std::cerr
            << "CreatePETBodySupportMask: "
               "could not determine a PET body threshold."
            << std::endl;

        return nullptr;
    }

    if (thresholdOut)
    {
        *thresholdOut =
            threshold;
    }

    if (usedOtsuFallbackOut)
    {
        *usedOtsuFallbackOut =
            usedOtsuFallback;
    }

    std::cout
        << "PET body-support threshold = "
        << threshold
        << " ("
        << (usedOtsuFallback
                ? "log-Otsu fallback"
                : "multiscale log-histogram")
        << ", "
        << (durationWeighted
                ? "duration-weighted sum"
                : "unweighted sum")
        << ")"
        << std::endl;

    // ------------------------------------------------------------------
    // 3. Obtain PET physical geometry.
    // ------------------------------------------------------------------

    vtkMRMLTransformNode*
        referenceParentTransform =
            referencePETNode->
                GetParentTransformNode();

    vtkSmartPointer<vtkOrientedImageData>
        petOriented =
            vtkSmartPointer<
                vtkOrientedImageData>::Take(
                vtkSlicerSegmentationsModuleLogic::
                    CreateOrientedImageDataFromVolumeNode(
                        referencePETNode,
                        referenceParentTransform));

    if (!petOriented)
    {
        return nullptr;
    }

    vtkNew<vtkOrientedImageData>
        candidateMask;

    candidateMask->SetExtent(
        petOriented->GetExtent());

    candidateMask->AllocateScalars(
        VTK_UNSIGNED_CHAR,
        1);

    vtkNew<vtkMatrix4x4>
        petImageToWorld;

    petOriented->
        GetImageToWorldMatrix(
            petImageToWorld);

    candidateMask->
        SetImageToWorldMatrix(
            petImageToWorld);

    unsigned char* maskPointer =
        static_cast<unsigned char*>(
            candidateMask->
                GetScalarPointer());

    if (!maskPointer)
    {
        return nullptr;
    }

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (vtkIdType v = 0;
         v <
         static_cast<vtkIdType>(
             numberOfVoxels);
         ++v)
    {
        const double value =
            composite[
                static_cast<size_t>(v)];

        maskPointer[v] =
            (std::isfinite(value) &&
             value >= threshold)
                ? static_cast<unsigned char>(1)
                : static_cast<unsigned char>(0);
    }

    candidateMask->Modified();

    // ------------------------------------------------------------------
    // 4. Small 3-D closing.
    //
    // vtkImageOpenClose3D performs proper binary closing using the VTK
    // implementation instead of hand-coded morphology.
    // ------------------------------------------------------------------

    vtkNew<vtkImageOpenClose3D>
        closing;

    closing->SetInputData(
        candidateMask);

    closing->SetKernelSize(
        3,
        3,
        3);

    closing->SetCloseValue(1);
    closing->SetOpenValue(0);

    closing->Update();

    // ------------------------------------------------------------------
    // 5. Measure connected-component sizes.
    // ------------------------------------------------------------------

    vtkNew<vtkImageConnectivityFilter>
        inspectComponents;

    inspectComponents->SetInputConnection(
        closing->GetOutputPort());

    inspectComponents->SetScalarRange(
        1.0,
        1.0);

    inspectComponents->
        SetExtractionModeToAllRegions();

    inspectComponents->
        SetLabelModeToSizeRank();

    inspectComponents->Update();

    vtkIdTypeArray* regionSizes =
        inspectComponents->
            GetExtractedRegionSizes();

    if (!regionSizes ||
        regionSizes->GetNumberOfValues() == 0)
    {
        return nullptr;
    }

    vtkIdType largestRegionSize = 0;

    for (vtkIdType i = 0;
         i <
         regionSizes->GetNumberOfValues();
         ++i)
    {
        largestRegionSize =
            std::max(
                largestRegionSize,
                regionSizes->GetValue(i));
    }

    const double safeFraction =
        std::max(
            0.0,
            std::min(
                1.0,
                minComponentFractionOfLargest));

    const vtkIdType minimumRegionSize =
        std::max<vtkIdType>(
            1,
            static_cast<vtkIdType>(
                std::ceil(
                    safeFraction *
                    static_cast<double>(
                        largestRegionSize))));

    // ------------------------------------------------------------------
    // 6. Retain all sufficiently large components.
    //
    // This is intentionally NOT "largest component only": disconnected
    // arms/hands can therefore survive.
    // ------------------------------------------------------------------

    vtkNew<vtkImageConnectivityFilter>
        bodyComponents;

    bodyComponents->SetInputConnection(
        closing->GetOutputPort());

    bodyComponents->SetScalarRange(
        1.0,
        1.0);

    bodyComponents->
        SetExtractionModeToAllRegions();

    bodyComponents->SetSizeRange(
        minimumRegionSize,
        std::numeric_limits<
            vtkIdType>::max());

    bodyComponents->
        SetLabelModeToConstantValue();

    bodyComponents->
        SetLabelConstantValue(1);

    bodyComponents->Update();

    vtkAlgorithmOutput* currentMaskPort =
        bodyComponents->
            GetOutputPort();

    // ------------------------------------------------------------------
    // 7. Fill enclosed cavities.
    // ------------------------------------------------------------------

    vtkNew<vtkImageThreshold>
        invertBody;

    vtkNew<vtkImageConnectivityFilter>
        exteriorBackground;

    vtkNew<vtkImageThreshold>
        filledBody;

    if (fillHoles)
    {
        invertBody->SetInputConnection(
            currentMaskPort);

        invertBody->ThresholdByLower(
            0.0);

        invertBody->SetInValue(1);
        invertBody->SetOutValue(0);

        invertBody->SetOutputScalarType(
            VTK_UNSIGNED_CHAR);

        invertBody->ReplaceInOn();
        invertBody->ReplaceOutOn();

        exteriorBackground->
            SetInputConnection(
                invertBody->
                    GetOutputPort());

        exteriorBackground->
            SetScalarRange(
                1.0,
                1.0);

        exteriorBackground->
            SetExtractionModeToLargestRegion();

        exteriorBackground->
            SetLabelModeToConstantValue();

        exteriorBackground->
            SetLabelConstantValue(1);

        filledBody->SetInputConnection(
            exteriorBackground->
                GetOutputPort());

        filledBody->ThresholdByLower(
            0.0);

        filledBody->SetInValue(1);
        filledBody->SetOutValue(0);

        filledBody->SetOutputScalarType(
            VTK_UNSIGNED_CHAR);

        filledBody->ReplaceInOn();
        filledBody->ReplaceOutOn();

        filledBody->Update();

        currentMaskPort =
            filledBody->GetOutputPort();
    }

    // ------------------------------------------------------------------
    // 8. Optional physical margin.
    // ------------------------------------------------------------------

    vtkNew<vtkImageDilateErode3D>
        dilate;

    if (bodyMarginMm > 0.0)
    {
        const double* spacing =
            petOriented->GetSpacing();

        int radius[3] =
        {
            0,
            0,
            0
        };

        for (int axis = 0;
             axis < 3;
             ++axis)
        {
            const double safeSpacing =
                std::max(
                    spacing[axis],
                    1e-6);

            radius[axis] =
                static_cast<int>(
                    std::ceil(
                        bodyMarginMm /
                        safeSpacing));
        }

        dilate->SetInputConnection(
            currentMaskPort);

        dilate->SetDilateValue(1);
        dilate->SetErodeValue(0);

        dilate->SetKernelSize(
            2 * radius[0] + 1,
            2 * radius[1] + 1,
            2 * radius[2] + 1);

        dilate->Update();

        currentMaskPort =
            dilate->GetOutputPort();
    }

    vtkAlgorithm* finalAlgorithm =
        currentMaskPort->
            GetProducer();

    finalAlgorithm->Update();

    vtkImageData* finalImage =
        vtkImageData::SafeDownCast(
            finalAlgorithm->
                GetOutputDataObject(0));

    if (!finalImage)
    {
        return nullptr;
    }

    vtkSmartPointer<
        vtkOrientedImageData> result =
            vtkSmartPointer<
                vtkOrientedImageData>::New();

    result->DeepCopy(
        finalImage);

    result->SetImageToWorldMatrix(
        petImageToWorld);

    return result;
}

vtkMRMLSegmentationNode*
vtkSlicerDynamicPETLogic::CreateOrUpdateBodySupportPreview(
    vtkOrientedImageData* mask,
    vtkMRMLScalarVolumeNode* referencePETNode,
    const std::string& nodeName)
{
    if (!mask || !referencePETNode)
    {
        return nullptr;
    }

    vtkMRMLScene* scene =
        this->GetMRMLScene();

    if (!scene)
    {
        return nullptr;
    }

    // Remove previous debug/preview node.
    vtkMRMLSegmentationNode* previous =
        vtkMRMLSegmentationNode::SafeDownCast(
            scene->GetFirstNodeByName(
                nodeName.c_str()));

    if (previous)
    {
        scene->RemoveNode(previous);
    }

    vtkMRMLSegmentationNode* segmentationNode =
        vtkMRMLSegmentationNode::SafeDownCast(
            scene->AddNewNodeByClass(
                "vtkMRMLSegmentationNode",
                nodeName.c_str()));

    if (!segmentationNode)
    {
        return nullptr;
    }

    segmentationNode->SetAttribute(
        "SlicerDynamicPET.InternalNode",
        "1");

    segmentationNode->
        SetReferenceImageGeometryParameterFromVolumeNode(
            referencePETNode);

    segmentationNode->CreateDefaultDisplayNodes();

    double color[3] =
    {
        0.2,
        0.9,
        0.3
    };

    const std::string segmentID =
        segmentationNode->
            AddSegmentFromBinaryLabelmapRepresentation(
                mask,
                "Fitting support",
                color);

    if (segmentID.empty())
    {
        scene->RemoveNode(
            segmentationNode);

        return nullptr;
    }

    vtkMRMLSegmentationDisplayNode* displayNode =
        vtkMRMLSegmentationDisplayNode::SafeDownCast(
            segmentationNode->GetDisplayNode());

    if (displayNode)
    {
        // Transparent body overlay + strong contour.
        displayNode->SetSegmentOpacity2DFill(
            segmentID,
            0.20);

        displayNode->SetSegmentOpacity2DOutline(
            segmentID,
            1.0);

        // Do not automatically generate/use a 3-D body surface.
        displayNode->SetSegmentVisibility3D(
            segmentID,
            false);
    }

    // Put the preview beside the source PET in Subject Hierarchy.
    vtkMRMLSubjectHierarchyNode* shNode =
        vtkMRMLSubjectHierarchyNode::
            GetSubjectHierarchyNode(scene);

    if (shNode)
    {
        const vtkIdType petItem =
            shNode->GetItemByDataNode(
                referencePETNode);

        const vtkIdType supportItem =
            shNode->GetItemByDataNode(
                segmentationNode);

        if (petItem !=
                vtkMRMLSubjectHierarchyNode::
                    INVALID_ITEM_ID &&
            supportItem !=
                vtkMRMLSubjectHierarchyNode::
                    INVALID_ITEM_ID)
        {
            shNode->SetItemParent(
                supportItem,
                shNode->GetItemParent(
                    petItem));
        }
    }

    return segmentationNode;
}

vtkSmartPointer<vtkOrientedImageData>
vtkSlicerDynamicPETLogic::
CombineBodySupportMasks(
    vtkOrientedImageData* firstMask,
    vtkOrientedImageData* secondMask,
    BodySupportCombination combination)
{
    if (!firstMask || !secondMask)
    {
        return nullptr;
    }

    if (!vtkOrientedImageDataResample::
            DoGeometriesMatch(
                firstMask,
                secondMask))
    {
        std::cerr
            << "CombineBodySupportMasks: "
               "mask geometries do not match."
            << std::endl;

        return nullptr;
    }

    vtkNew<vtkImageLogic> imageLogic;

    imageLogic->SetInput1Data(
        firstMask);

    imageLogic->SetInput2Data(
        secondMask);

    if (combination ==
        BodySupportCombination::Union)
    {
        imageLogic->SetOperationToOr();
    }
    else
    {
        imageLogic->SetOperationToAnd();
    }

    imageLogic->SetOutputTrueValue(1.0);
    imageLogic->Update();

    vtkSmartPointer<vtkOrientedImageData>
        result =
            vtkSmartPointer<
                vtkOrientedImageData>::New();

    result->DeepCopy(
        imageLogic->GetOutput());

    vtkNew<vtkMatrix4x4>
        imageToWorld;

    firstMask->GetImageToWorldMatrix(
        imageToWorld);

    result->SetImageToWorldMatrix(
        imageToWorld);

    return result;
}

void vtkSlicerDynamicPETLogic::Image2Flatten(
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

vtkMRMLScalarVolumeNode* vtkSlicerDynamicPETLogic::Flatten2Image(
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


void vtkSlicerDynamicPETLogic::Patlak4Img(
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
    std::function<void(int)> progressCallback)
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

        const int nFit =
            static_cast<int>(
                fitVoxelIndices.size());

        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int fitIndex = 0;
             fitIndex < nFit;
             ++fitIndex)
        {
            const int v =
                fitVoxelIndices[fitIndex];

            if (v < 0 ||
                static_cast<size_t>(v) >= nVoxels)
            {
                continue;
            }
            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
            }
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
                    if (stopRequested.load(
                            std::memory_order_relaxed))
                    {
                      break;
                    }
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

            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
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
            const size_t done =
                ++voxProcessed;

            const size_t updateInterval =
                std::max(
                    static_cast<size_t>(1),
                    static_cast<size_t>(
                        std::max(1, nFit / 200)));

            if (progressCallback &&
                (done % updateInterval == 0 ||
                 done == static_cast<size_t>(nFit)))
            {
                const int progress =
                    static_cast<int>(
                        100LL *
                        done /
                        static_cast<size_t>(nFit));

                progressCallback(progress);
            }

        }
    }
}

void vtkSlicerDynamicPETLogic::Logan4Img(
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
    const std::vector<int>& fitVoxelIndices,
    std::atomic<bool>& stopRequested,
    int numThreads,
    std::function<void(int)> progressCallback)
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

    std::atomic<size_t> voxProcessed(0);

    // ---------- 4) Parallel voxel-wise processing ----------
    #ifdef HAVE_OPENMP
    #pragma omp parallel
    #endif
    {
        std::vector<double> Xvec(n), Yvec(n), fitted(n), wgt_adj(n);
        Eigen::Vector2d coeff;
        Eigen::VectorXd coeff_vec(2);

        const int nFit =
            static_cast<int>(
                fitVoxelIndices.size());

        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int fitIndex = 0;
             fitIndex < nFit;
             ++fitIndex)
        {
            const int v =
                fitVoxelIndices[fitIndex];

            if (v < 0 ||
                static_cast<size_t>(v) >= nVoxels)
            {
                continue;
            }
            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
            }
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
                    if (stopRequested.load(
                            std::memory_order_relaxed))
                    {
                      break;
                    }
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

            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
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
            const size_t done =
                ++voxProcessed;

            const size_t updateInterval =
                std::max(
                    static_cast<size_t>(1),
                    static_cast<size_t>(
                        std::max(1, nFit / 200)));

            if (progressCallback &&
                (done % updateInterval == 0 ||
                 done == static_cast<size_t>(nFit)))
            {
                const int progress =
                    static_cast<int>(
                        100LL *
                        done /
                        static_cast<size_t>(nFit));

                progressCallback(progress);
            }
        }
    }

}

void vtkSlicerDynamicPETLogic::RE4Img(
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
    const std::vector<int>& fitVoxelIndices,
    std::atomic<bool>& stopRequested,
    int numThreads,
    std::function<void(int)> progressCallback)
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

        const int nFit =
            static_cast<int>(
                fitVoxelIndices.size());

        #ifdef HAVE_OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for (int fitIndex = 0;
             fitIndex < nFit;
             ++fitIndex)
        {
            const int v =
                fitVoxelIndices[fitIndex];

            if (v < 0 ||
                static_cast<size_t>(v) >= nVoxels)
            {
                continue;
            }
            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
            }
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
                    if (stopRequested.load(
                            std::memory_order_relaxed))
                    {
                      break;
                    }
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

            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
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
            const size_t done =
                ++voxProcessed;

            const size_t updateInterval =
                std::max(
                    static_cast<size_t>(1),
                    static_cast<size_t>(
                        std::max(1, nFit / 200)));

            if (progressCallback &&
                (done % updateInterval == 0 ||
                 done == static_cast<size_t>(nFit)))
            {
                const int progress =
                    static_cast<int>(
                        100LL *
                        done /
                        static_cast<size_t>(nFit));

                progressCallback(progress);
            }
        }
    }

}


std::vector<double> vtkSlicerDynamicPETLogic::ExtractParameter(
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

void vtkSlicerDynamicPETLogic::CreateMTGAParametricImages(
    const std::vector<MTGAParameters>& outputParams,
    const int dims[3],
    const std::vector<std::string>& fields,
    const std::string& modelID,
    vtkMRMLScalarVolumeNode* refNode,
    vtkMRMLSubjectHierarchyNode* refSH,
    vtkIdType refID
  )
{

  vtkMRMLScene* scene =
      this->GetMRMLScene();

  if (!scene)
  {
    std::cerr
        << "CreateMTGAParametricImages: "
           "MRML scene is null."
        << std::endl;
    return;
  }

  for (const auto& field : fields)
  {
    std::vector<double> flatten = this->ExtractParameter(outputParams, field);
    vtkMRMLScalarVolumeNode* existingNode = nullptr;

    const int numberOfVolumes =
        scene->GetNumberOfNodesByClass(
            "vtkMRMLScalarVolumeNode");

    for (int i = 0;
         i < numberOfVolumes;
         ++i)
    {
      vtkMRMLScalarVolumeNode* candidate =
          vtkMRMLScalarVolumeNode::SafeDownCast(
              scene->GetNthNodeByClass(
                  i,
                  "vtkMRMLScalarVolumeNode"));

      if (!candidate)
      {
        continue;
      }

      const char* resultType =
          candidate->GetAttribute(
              "SlicerDynamicPET.ResultType");

      const char* method =
          candidate->GetAttribute(
              "SlicerDynamicPET.Method");

      const char* model =
          candidate->GetAttribute(
              "SlicerDynamicPET.Model");

      const char* parameter =
          candidate->GetAttribute(
              "SlicerDynamicPET.Parameter");

      const char* sourceNodeID =
          candidate->GetAttribute(
              "SlicerDynamicPET.SourceNodeID");

      if (!resultType ||
          !method ||
          !model ||
          !parameter ||
          !sourceNodeID ||
          !refNode->GetID())
      {
        continue;
      }

      if (std::string(resultType) ==
              "ParametricMap" &&
          std::string(method) ==
              "MTGA" &&                // TCM below: use "TCM"
          std::string(model) ==
              modelID &&
          std::string(parameter) ==
              field &&
          std::string(sourceNodeID) ==
              refNode->GetID())
      {
        existingNode = candidate;
        break;
      }
    }

    if (existingNode)
    {
      scene->RemoveNode(existingNode);
    }

    vtkMRMLScalarVolumeNode* node =
        this->Flatten2Image(
            flatten,
            dims,
            modelID + " - " + field);

    if (!node)
    {
      std::cerr
          << "Could not create image for "
          << modelID + " - " + field
          << std::endl;
      return;
    }

    node->SetAttribute(
        "SlicerDynamicPET.ResultType",
        "ParametricMap");

    node->SetAttribute(
        "SlicerDynamicPET.Method",
        "MTGA");

    node->SetAttribute(
        "SlicerDynamicPET.Model",
        modelID.c_str());

    node->SetAttribute(
        "SlicerDynamicPET.Parameter",
        field.c_str());

    node->SetAttribute(
        "SlicerDynamicPET.SourceNodeID",
        refNode->GetID());

    node->CopyOrientation(refNode);
    node->SetSpacing(refNode->GetSpacing());
    node->SetOrigin(refNode->GetOrigin());

    vtkIdType parentItemID =
        refSH->GetItemParent(refID);

    vtkIdType newItemID =
        refSH->GetItemByDataNode(node);

    refSH->SetItemParent(
        newItemID,
        parentItemID);
  }
}


void vtkSlicerDynamicPETLogic::callTCMImg(
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
    const std::vector<int>& fitVoxelIndices,
    std::atomic<bool>& stopRequested,
    const std::vector<double>* wgt_global /*= nullptr*/,
    int numThreads /*= 0 */,
    std::function<void(int)> progressCallback, /*= nullptr*/
    std::function<bool()> stopCallback /*= nullptr*/
)
{
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


    // ---------- 1) Weights ----------
    std::vector<double> wt(Nframe, 1.0);
    if (wgt_global)
    {
        if (wgt_global->size() !=
            static_cast<size_t>(Nframe))
        {
            std::cerr
                << "callTCMImg: weight vector size mismatch."
                << std::endl;
            return;
        }

        wt = *wgt_global;
    }
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
        if (std::abs(pbr) < EPS)
        {
            pbr =
                (pbr < 0.0) ? -EPS : EPS;
        }
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
        params.weights.clear();
        params.keep.clear();
        params.dof = dof_fixed;

        params.AIC    = std::numeric_limits<double>::quiet_NaN();
        params.BIC    = std::numeric_limits<double>::quiet_NaN();
        params.MASE   = std::numeric_limits<double>::quiet_NaN();
        params.chi2   = std::numeric_limits<double>::quiet_NaN();
        params.loglik = std::numeric_limits<double>::quiet_NaN();
    };
    // All voxels start excluded.
    // Only voxels present in fitVoxelIndices will receive fitted values.
    for (TCMParameters& params : outputParams)
    {
        markInvalidVoxel(params);
    }

    const int Nfit =
        static_cast<int>(fitVoxelIndices.size());

    std::cout
        << "TCM fitting: "
        << Nfit << " / " << Nvox
        << " voxels eligible for fitting; "
        << (Nvox - Nfit)
        << " excluded by the common voxel mask."
        << std::endl;

    if (Nfit == 0)
    {
        if (progressCallback)
        {
            progressCallback(100);
        }

        // Clean allocations performed above.
        for (double* frame : scant)
        {
            delete[] frame;
        }

        delete[] Cp_new;
        delete[] cwb_new;

        return;
    }

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

    auto reportVoxelProcessed =
        [&]()
        {
            const int done = ++voxProcessed;

            if (progressCallback &&
                (done % progressUpdateInterval == 0 ||
                 done == Nfit))
            {
                const int progress =
                    static_cast<int>(
                        100LL * done / Nfit);

                progressCallback(progress);
            }
        };

    std::atomic<int> guardHitCount(0);

    // ---------- 8) Parallel voxel loop with per-thread scratch buffers ----------
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
            if (stopRequested.load(std::memory_order_relaxed) ||
                (stopCallback && stopCallback()))
            {
                continue;
            }

            const int v =
                fitVoxelIndices[fitIndex];

            if (v < 0 || v >= Nvox)
            {
                continue;
            }

            TCMParameters params;

            for (int i = 0; i < Nframe; ++i) cfit_local[i] = voxels[v][i];
            if (check_if_constant(cfit_local)) {
              reportVoxelProcessed();
              continue;
            }

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

            if (stopRequested.load(std::memory_order_relaxed) ||
                (stopCallback && stopCallback()))
            {
              continue;
            }

            if (stx.hit_guard > 0)
            {
                ++guardHitCount;
                reportVoxelProcessed();
                continue;
            }

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
            reportVoxelProcessed();
        }
    }

    if (guardHitCount.load() > 0)
    {
        std::cout
            << "TCM fitting: "
            << guardHitCount.load()
            << " voxels failed to converge "
               "(LM iteration guard reached)."
            << std::endl;
    }


    for (double* frame : scant)
    {
        delete[] frame;
    }

    delete[] scant_flatten;
    delete[] Cp_new;
    delete[] cwb_new;
}

std::vector<double> vtkSlicerDynamicPETLogic::ExtractParameter(
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

void vtkSlicerDynamicPETLogic::CreateTCMParametricImages(
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

        vtkMRMLScalarVolumeNode* existingNode = nullptr;

        const int numberOfVolumes =
            scene->GetNumberOfNodesByClass(
                "vtkMRMLScalarVolumeNode");

        for (int i = 0;
             i < numberOfVolumes;
             ++i)
        {
          vtkMRMLScalarVolumeNode* candidate =
              vtkMRMLScalarVolumeNode::SafeDownCast(
                  scene->GetNthNodeByClass(
                      i,
                      "vtkMRMLScalarVolumeNode"));

          if (!candidate)
          {
            continue;
          }

          const char* resultType =
              candidate->GetAttribute(
                  "SlicerDynamicPET.ResultType");

          const char* method =
              candidate->GetAttribute(
                  "SlicerDynamicPET.Method");

          const char* model =
              candidate->GetAttribute(
                  "SlicerDynamicPET.Model");

          const char* parameter =
              candidate->GetAttribute(
                  "SlicerDynamicPET.Parameter");

          const char* sourceNodeID =
              candidate->GetAttribute(
                  "SlicerDynamicPET.SourceNodeID");

          if (!resultType ||
              !method ||
              !model ||
              !parameter ||
              !sourceNodeID ||
              !refNode->GetID())
          {
            continue;
          }

          if (std::string(resultType) ==
                  "ParametricMap" &&
              std::string(method) ==
                  "TCM" &&
              std::string(model) ==
                  modelID &&
              std::string(parameter) ==
                  field &&
              std::string(sourceNodeID) ==
                  refNode->GetID())
          {
            existingNode = candidate;
            break;
          }
        }

        if (existingNode)
        {
          scene->RemoveNode(existingNode);
        }

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

        node->SetAttribute(
            "SlicerDynamicPET.ResultType",
            "ParametricMap");

        node->SetAttribute(
            "SlicerDynamicPET.Method",
            "TCM");

        node->SetAttribute(
            "SlicerDynamicPET.Model",
            modelID.c_str());

        node->SetAttribute(
            "SlicerDynamicPET.Parameter",
            field.c_str());

        node->SetAttribute(
            "SlicerDynamicPET.SourceNodeID",
            refNode->GetID());

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
