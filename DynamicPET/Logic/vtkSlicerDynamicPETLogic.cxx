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
#include <kmaplib.h>
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
#include <algorithm>
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

  bool ComputePchipSlopes(
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
      std::vector<double> h(n - 1), delta(n - 1);
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
              return 0.0;
          if (d0 * d1 < 0.0 && std::abs(m) > 3.0 * std::abs(d0))
              return 3.0 * d0;
          return m;
      };
      slopes.front() = endpointSlope(h[0], h[1], delta[0], delta[1]);
      slopes.back() = endpointSlope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
      return true;
  }

  double EvaluatePchip(
      const std::vector<double>& x,
      const std::vector<double>& y,
      const std::vector<double>& slopes,
      double t)
  {
      if (x.size() < 2 || x.size() != y.size() || slopes.size() != x.size() ||
          t < x.front() || t > x.back())
      {
          return std::numeric_limits<double>::quiet_NaN();
      }
      if (t == x.back())
      {
          return y.back();
      }
      const auto upper = std::upper_bound(x.begin(), x.end(), t);
      const size_t right = static_cast<size_t>(std::distance(x.begin(), upper));
      const size_t left = right - 1;
      const double h = x[right] - x[left];
      const double u = (t - x[left]) / h;
      const double u2 = u * u;
      const double u3 = u2 * u;
      const double h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
      const double h10 = u3 - 2.0 * u2 + u;
      const double h01 = -2.0 * u3 + 3.0 * u2;
      const double h11 = u3 - u2;
      return std::max(0.0,
          h00 * y[left] + h10 * h * slopes[left] +
          h01 * y[right] + h11 * h * slopes[right]);
  }

  void BuildFrameRepresentativeCurve(
      double** scant,
      const std::vector<double>& frameValues,
      std::vector<double>& times,
      std::vector<double>& values)
  {
      const size_t n = frameValues.size();
      times.clear();
      values.clear();
      if (n == 0)
      {
          return;
      }
      times.reserve(n + 2);
      values.reserve(n + 2);
      times.push_back(scant[0][0]);
      values.push_back(frameValues.front());
      for (size_t i = 0; i < n; ++i)
      {
          times.push_back(0.5 * (scant[i][0] + scant[i][1]));
          values.push_back(frameValues[i]);
      }
      times.push_back(scant[n - 1][1]);
      values.push_back(frameValues.back());
  }

  double*
  FineSampleExplicitInputFunction(
      const std::vector<double>& timesSec,
      const std::vector<double>& values,
      double endTimeSec,
      double timeStepSec,
      const std::string& interpolationType,
      long int& numberOfSamples)
  {
      if (timesSec.size() < 2 ||
          timesSec.size() != values.size() ||
          timeStepSec <= 0.0 ||
          timesSec.front() > 0.0 ||
          timesSec.back() < endTimeSec)
      {
          throw std::invalid_argument(
              "Invalid explicit input-function sampling.");
      }

      numberOfSamples =
          static_cast<long int>(
              endTimeSec /
              timeStepSec);

      if (numberOfSamples <= 0)
      {
          throw std::invalid_argument(
              "Invalid TCM integration grid.");
      }

      double* output =
          new double[numberOfSamples];

      std::vector<double> pchipSlope;
      if (interpolationType == "pchip" &&
          !ComputePchipSlopes(timesSec, values, pchipSlope))
      {
          delete[] output;
          throw std::invalid_argument("Invalid PCHIP input-function sampling.");
      }

      for (long int i = 0;
           i < numberOfSamples;
           ++i)
      {
          const double t =
              (static_cast<double>(i) + 0.5) *
              timeStepSec;

          if (t < timesSec.front() ||
              t > timesSec.back())
          {
              delete[] output;

              throw std::runtime_error(
                  "Input-function extrapolation "
                  "would be required.");
          }

          const auto upper =
              std::upper_bound(
                  timesSec.begin(),
                  timesSec.end(),
                  t);

          if (upper == timesSec.end())
          {
              output[i] =
                  values.back();

              continue;
          }

          if (upper == timesSec.begin())
          {
              output[i] =
                  values.front();

              continue;
          }

          const size_t right =
              static_cast<size_t>(
                  std::distance(
                      timesSec.begin(),
                      upper));

          const size_t left =
              right - 1;

          if (interpolationType == "const")
          {
              output[i] = values[left];
          }
          else if (interpolationType == "pchip")
          {
              output[i] = EvaluatePchip(timesSec, values, pchipSlope, t);
          }
          else
          {
              const double t1 =
                  timesSec[left];

              const double t2 =
                  timesSec[right];

              const double y1 =
                  values[left];

              const double y2 =
                  values[right];

              output[i] =
                  std::max(
                      0.0,
                      y1 +
                      (y2 - y1) *
                      (t - t1) /
                      (t2 - t1));
          }
      }

      return output;
  }

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


struct FengFitContext
{
    std::vector<double> timesMin;
    std::vector<double> frameStartMin;
    std::vector<double> frameEndMin;
    bool frameAverages{false};
};

static double fengSafeExp(double x)
{
    return std::exp(std::max(-40.0, std::min(40.0, x)));
}

static void fengInternalToPhysical(
    const double* u,
    FengParameters& p)
{
    p.tau = u[0];
    p.A1 = fengSafeExp(u[1]);
    p.A2 = fengSafeExp(u[2]);
    p.A3 = fengSafeExp(u[3]);

    p.lambda3 = fengSafeExp(u[6]);
    p.lambda2 = p.lambda3 + fengSafeExp(u[5]);
    p.lambda1 = p.lambda2 + fengSafeExp(u[4]);
}

static double fengI0(double x, double lambda)
{
    if (x <= 0.0)
    {
        return 0.0;
    }

    const double z = lambda * x;
    if (std::abs(z) < 1e-6)
    {
        return x *
            (1.0 - 0.5 * z + z * z / 6.0);
    }

    return -std::expm1(-z) / lambda;
}

static double fengI1(double x, double lambda)
{
    if (x <= 0.0)
    {
        return 0.0;
    }

    const double z = lambda * x;
    if (std::abs(z) < 1e-5)
    {
        return x * x *
            (0.5 - z / 3.0 + z * z / 8.0);
    }

    return
        (1.0 - std::exp(-z) * (1.0 + z)) /
        (lambda * lambda);
}

static double fengValueMin(
    double timeMin,
    const FengParameters& p)
{
    const double x = timeMin - p.tau;
    if (x < 0.0)
    {
        return 0.0;
    }

    const double value =
        (p.A1 * x - p.A2 - p.A3) *
            std::exp(-p.lambda1 * x) +
        p.A2 * std::exp(-p.lambda2 * x) +
        p.A3 * std::exp(-p.lambda3 * x);

    return value;
}

static double fengCumulativeMin(
    double timeMin,
    const FengParameters& p)
{
    const double x = std::max(0.0, timeMin - p.tau);

    return
        p.A1 * fengI1(x, p.lambda1) -
        (p.A2 + p.A3) * fengI0(x, p.lambda1) +
        p.A2 * fengI0(x, p.lambda2) +
        p.A3 * fengI0(x, p.lambda3);
}

static void fengEvaluateInternal(
    double* u,
    void* param,
    double* out)
{
    auto* ctx = static_cast<FengFitContext*>(param);
    FengParameters physical;
    fengInternalToPhysical(u, physical);

    const size_t n = ctx->frameAverages
        ? ctx->frameStartMin.size()
        : ctx->timesMin.size();

    for (size_t i = 0; i < n; ++i)
    {
        if (ctx->frameAverages)
        {
            const double start = ctx->frameStartMin[i];
            const double end = ctx->frameEndMin[i];
            const double duration = end - start;
            out[i] = duration > 0.0
                ? (fengCumulativeMin(end, physical) -
                   fengCumulativeMin(start, physical)) / duration
                : 0.0;
        }
        else
        {
            out[i] = fengValueMin(ctx->timesMin[i], physical);
        }
    }
}

static void fengJacobianInternal(
    double* u,
    void* param,
    double* out,
    int* psens,
    double* jac)
{
    auto* ctx = static_cast<FengFitContext*>(param);
    const size_t n = ctx->frameAverages
        ? ctx->frameStartMin.size()
        : ctx->timesMin.size();

    fengEvaluateInternal(u, param, out);

    std::array<double, 7> plus{};
    std::array<double, 7> minus{};
    for (int j = 0; j < 7; ++j)
    {
        plus[static_cast<size_t>(j)] = u[j];
        minus[static_cast<size_t>(j)] = u[j];
    }

    std::vector<double> yPlus(n, 0.0);
    std::vector<double> yMinus(n, 0.0);

    int column = 0;
    for (int j = 0; j < 7; ++j)
    {
        if (!psens[j])
        {
            continue;
        }

        for (int k = 0; k < 7; ++k)
        {
            plus[static_cast<size_t>(k)] = u[k];
            minus[static_cast<size_t>(k)] = u[k];
        }

        const double step =
            1e-5 * (1.0 + std::abs(u[j]));
        plus[static_cast<size_t>(j)] += step;
        minus[static_cast<size_t>(j)] -= step;

        fengEvaluateInternal(plus.data(), param, yPlus.data());
        fengEvaluateInternal(minus.data(), param, yMinus.data());

        for (size_t i = 0; i < n; ++i)
        {
            jac[i + static_cast<size_t>(column) * n] =
                (yPlus[i] - yMinus[i]) / (2.0 * step);
        }
        ++column;
    }
}


struct ParentFractionFitContext
{
    std::vector<double> timesMin;
    ParentFractionModel model{ParentFractionModel::Hill};
    double maximumTimeMin{1.0};
};

static double parentFractionSafeExp(double x)
{
    return std::exp(std::max(-30.0, std::min(30.0, x)));
}

static double parentFractionSigmoid(double x)
{
    if (x >= 0.0)
    {
        const double z = std::exp(-std::min(40.0, x));
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(std::max(-40.0, x));
    return z / (1.0 + z);
}

static double parentFractionLogit(double p)
{
    const double clipped = std::max(1e-6, std::min(1.0 - 1e-6, p));
    return std::log(clipped / (1.0 - clipped));
}

static int parentFractionParameterCount(ParentFractionModel model)
{
    switch (model)
    {
      case ParentFractionModel::Hill:
        return 3;
      case ParentFractionModel::ExtendedHill:
        return 5;
      case ParentFractionModel::ExponentialPlateau:
        return 2;
      default:
        return 0;
    }
}

static void parentFractionInternalToPhysical(
    const double* u,
    const ParentFractionFitContext& ctx,
    ParentFractionFitParameters& p)
{
    p = ParentFractionFitParameters{};

    if (ctx.model == ParentFractionModel::Hill)
    {
        p.A = parentFractionSigmoid(u[0]);
        p.B = 1.0 + parentFractionSafeExp(u[1]);
        p.C = parentFractionSafeExp(u[2]);
    }
    else if (ctx.model == ParentFractionModel::ExtendedHill)
    {
        p.D = parentFractionSigmoid(u[0]);
        const double lowerRatio = parentFractionSigmoid(u[1]);
        p.A = p.D * lowerRatio;
        p.B = 1.0 + parentFractionSafeExp(u[2]);
        p.C = parentFractionSafeExp(u[3]);
        p.E = std::max(0.0, std::min(ctx.maximumTimeMin, u[4]));
    }
    else if (ctx.model == ParentFractionModel::ExponentialPlateau)
    {
        p.plateau = parentFractionSigmoid(u[0]);
        p.rate = parentFractionSafeExp(u[1]);
    }
}

static double parentFractionValueMin(
    double timeMin,
    ParentFractionModel model,
    const ParentFractionFitParameters& p)
{
    const double t = std::max(0.0, timeMin);
    double value = std::numeric_limits<double>::quiet_NaN();

    if (model == ParentFractionModel::Hill)
    {
        if (!(p.C > 0.0) || !(p.B >= 1.0))
            return value;
        if (t <= 0.0)
            return 1.0;
        const double tb = std::pow(t, p.B);
        value = 1.0 - (1.0 - p.A) * tb / (p.C + tb);
    }
    else if (model == ParentFractionModel::ExtendedHill)
    {
        if (!(p.C > 0.0) || !(p.B >= 1.0) ||
            !(p.D >= 0.0) || !(p.D <= 1.0) ||
            !(p.A >= 0.0) || !(p.A <= p.D))
            return value;
        if (t <= p.E)
        {
            value = p.D;
        }
        else
        {
            const double x = t - p.E;
            const double xb = std::pow(x, p.B);
            value = p.D - (p.D - p.A) * xb / (p.C + xb);
        }
    }
    else if (model == ParentFractionModel::ExponentialPlateau)
    {
        if (!(p.plateau >= 0.0) || !(p.plateau <= 1.0) || !(p.rate > 0.0))
            return value;
        value = p.plateau + (1.0 - p.plateau) * std::exp(-p.rate * t);
    }

    if (!std::isfinite(value))
        return value;
    return std::max(0.0, std::min(1.0, value));
}

static void parentFractionEvaluateInternal(
    double* u,
    void* param,
    double* out)
{
    auto* ctx = static_cast<ParentFractionFitContext*>(param);
    ParentFractionFitParameters p;
    parentFractionInternalToPhysical(u, *ctx, p);
    for (size_t i = 0; i < ctx->timesMin.size(); ++i)
    {
        out[i] = parentFractionValueMin(ctx->timesMin[i], ctx->model, p);
    }
}

static void parentFractionJacobianInternal(
    double* u,
    void* param,
    double* out,
    int* psens,
    double* jac)
{
    auto* ctx = static_cast<ParentFractionFitContext*>(param);
    const size_t n = ctx->timesMin.size();
    const int pCount = parentFractionParameterCount(ctx->model);

    parentFractionEvaluateInternal(u, param, out);

    std::vector<double> plus(static_cast<size_t>(pCount), 0.0);
    std::vector<double> minus(static_cast<size_t>(pCount), 0.0);
    std::vector<double> yPlus(n, 0.0);
    std::vector<double> yMinus(n, 0.0);

    int column = 0;
    for (int j = 0; j < pCount; ++j)
    {
        if (!psens[j])
            continue;

        for (int k = 0; k < pCount; ++k)
        {
            plus[static_cast<size_t>(k)] = u[k];
            minus[static_cast<size_t>(k)] = u[k];
        }

        const double step = 1e-5 * (1.0 + std::abs(u[j]));
        plus[static_cast<size_t>(j)] += step;
        minus[static_cast<size_t>(j)] -= step;

        parentFractionEvaluateInternal(plus.data(), param, yPlus.data());
        parentFractionEvaluateInternal(minus.data(), param, yMinus.data());

        for (size_t i = 0; i < n; ++i)
        {
            jac[i + static_cast<size_t>(column) * n] =
                (yPlus[i] - yMinus[i]) / (2.0 * step);
        }
        ++column;
    }
}

} // end anonymous namespace

static double DynamicPETBinaryLabelmapVolumeMm3(
    vtkMRMLScalarVolumeNode* referenceVolume,
    vtkImageData* labelmap,
    int labelValue = 1)
{
    if (!referenceVolume || !labelmap ||
        !labelmap->GetPointData() ||
        !labelmap->GetPointData()->GetScalars())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double spacing[3] = {1.0, 1.0, 1.0};
    referenceVolume->GetSpacing(spacing);
    const double voxelVolume =
        std::abs(spacing[0] * spacing[1] * spacing[2]);

    vtkDataArray* labels =
        labelmap->GetPointData()->GetScalars();
    const vtkIdType n = labelmap->GetNumberOfPoints();
    vtkIdType count = 0;
    for (vtkIdType i = 0; i < n; ++i)
    {
        if (static_cast<int>(labels->GetComponent(i, 0)) == labelValue)
        {
            ++count;
        }
    }
    return voxelVolume * static_cast<double>(count);
}

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

static double robustMADScale(
    const Eigen::VectorXd& residuals,
    const std::vector<double>& baseWeights)
{
    std::vector<double> activeResiduals;
    activeResiduals.reserve(
        static_cast<size_t>(residuals.size()));

    for (Eigen::Index i = 0; i < residuals.size(); ++i)
    {
        const double baseWeight =
            (static_cast<size_t>(i) < baseWeights.size())
            ? baseWeights[static_cast<size_t>(i)]
            : 1.0;

        if (baseWeight > 0.0 && std::isfinite(residuals(i)))
        {
            activeResiduals.push_back(residuals(i));
        }
    }

    if (activeResiduals.empty())
    {
        return 1.0;
    }

    const double residualMedian =
        median(activeResiduals);

    std::vector<double> absoluteDeviations;
    absoluteDeviations.reserve(activeResiduals.size());

    for (const double residual : activeResiduals)
    {
        absoluteDeviations.push_back(
            std::abs(residual - residualMedian));
    }

    // For Gaussian residuals, MAD / 0.67448975 estimates sigma.
    double scale =
        1.482602218505602 *
        median(absoluteDeviations);

    if (!std::isfinite(scale) || scale <= 1e-12)
    {
        double weightedSquareSum = 0.0;
        double weightSum = 0.0;

        for (Eigen::Index i = 0; i < residuals.size(); ++i)
        {
            const double baseWeight =
                (static_cast<size_t>(i) < baseWeights.size())
                ? baseWeights[static_cast<size_t>(i)]
                : 1.0;

            if (baseWeight <= 0.0 || !std::isfinite(residuals(i)))
            {
                continue;
            }

            weightedSquareSum +=
                baseWeight * residuals(i) * residuals(i);
            weightSum += baseWeight;
        }

        if (weightSum > 0.0)
        {
            scale = std::sqrt(weightedSquareSum / weightSum);
        }
    }

    if (!std::isfinite(scale) || scale <= 1e-12)
    {
        scale = 1.0;
    }

    return scale;
}

static bool solveHuberIRLS(
    const Eigen::MatrixXd& design,
    const Eigen::VectorXd& observations,
    const std::vector<double>& baseWeights,
    double huberTune,
    double tolerance,
    int maxIterations,
    Eigen::VectorXd& coefficients,
    std::vector<double>& finalWeights)
{
    const Eigen::Index n = observations.size();

    if (n <= 0 || design.rows() != n || design.cols() <= 0)
    {
        return false;
    }

    const double tune =
        std::max(huberTune, 1e-6);

    std::vector<double> base(
        static_cast<size_t>(n),
        1.0);

    if (!baseWeights.empty())
    {
        if (baseWeights.size() != static_cast<size_t>(n))
        {
            return false;
        }
        base = baseWeights;
    }

    auto solveWithWeights =
        [&](const std::vector<double>& combined,
            Eigen::VectorXd& solution) -> bool
        {
            Eigen::VectorXd diagonal(n);
            for (Eigen::Index i = 0; i < n; ++i)
            {
                const double weight =
                    combined[static_cast<size_t>(i)];
                diagonal(i) =
                    (std::isfinite(weight) && weight > 0.0)
                    ? weight
                    : 0.0;
            }

            if (diagonal.sum() <= 0.0)
            {
                return false;
            }

            const Eigen::MatrixXd normalMatrix =
                design.transpose() *
                diagonal.asDiagonal() *
                design;

            const Eigen::VectorXd rightHandSide =
                design.transpose() *
                diagonal.asDiagonal() *
                observations;

            solution =
                normalMatrix.ldlt().solve(rightHandSide);

            return solution.allFinite();
        };

    std::vector<double> combinedWeights(
        static_cast<size_t>(n),
        0.0);

    for (Eigen::Index i = 0; i < n; ++i)
    {
        combinedWeights[static_cast<size_t>(i)] =
            std::max(0.0, base[static_cast<size_t>(i)]);
    }

    if (!solveWithWeights(combinedWeights, coefficients))
    {
        return false;
    }

    for (int iteration = 0;
         iteration < std::max(1, maxIterations);
         ++iteration)
    {
        const Eigen::VectorXd residuals =
            observations - design * coefficients;

        const double scale =
            robustMADScale(residuals, base);

        for (Eigen::Index i = 0; i < n; ++i)
        {
            const double baseWeight =
                std::max(0.0, base[static_cast<size_t>(i)]);

            if (baseWeight <= 0.0)
            {
                combinedWeights[static_cast<size_t>(i)] = 0.0;
                continue;
            }

            const double standardizedResidual =
                std::abs(residuals(i)) / scale;

            const double robustWeight =
                (standardizedResidual <= tune)
                ? 1.0
                : tune / std::max(standardizedResidual, 1e-12);

            combinedWeights[static_cast<size_t>(i)] =
                baseWeight * robustWeight;
        }

        Eigen::VectorXd nextCoefficients;
        if (!solveWithWeights(
                combinedWeights,
                nextCoefficients))
        {
            return false;
        }

        const double difference =
            (nextCoefficients - coefficients).norm();

        const double reference =
            1.0 + coefficients.norm();

        coefficients = nextCoefficients;

        if (difference <= tolerance * reference)
        {
            break;
        }
    }

    finalWeights = combinedWeights;
    return true;
}

static void markBoundIfNeeded(
    unsigned int& flags,
    double value,
    double lower,
    double upper,
    bool sensitive,
    unsigned int lowerFlag,
    unsigned int upperFlag)
{
    if (!sensitive ||
        !std::isfinite(value) ||
        !std::isfinite(lower) ||
        !std::isfinite(upper) ||
        upper < lower)
    {
        return;
    }

    const double span =
        std::max(std::abs(upper - lower), 1e-8);

    const double tolerance =
        std::max(1e-8, 1e-3 * span);

    if (std::abs(value - lower) <= tolerance)
    {
        flags |= lowerFlag;
    }

    if (std::abs(value - upper) <= tolerance)
    {
        flags |= upperFlag;
    }
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
bool vtkSlicerDynamicPETLogic::FitFengInputFunction(
    const std::vector<double>& timesSec,
    const std::vector<double>& values,
    const std::vector<double>* frameStartSec,
    const std::vector<double>* frameEndSec,
    bool observationsAreFrameAverages,
    FengParameters& params,
    std::vector<double>& fittedObservationValues,
    std::string* errorMessage)
{
    params = FengParameters{};
    fittedObservationValues.clear();

    const size_t n = values.size();
    if (n < 8)
    {
        if (errorMessage)
        {
            *errorMessage =
                "Feng fitting requires at least 8 retained observations.";
        }
        return false;
    }

    FengFitContext context;
    context.frameAverages = observationsAreFrameAverages;

    if (observationsAreFrameAverages)
    {
        if (!frameStartSec || !frameEndSec ||
            frameStartSec->size() != n ||
            frameEndSec->size() != n)
        {
            if (errorMessage)
            {
                *errorMessage =
                    "Feng frame-average fitting requires matching frame start/end arrays.";
            }
            return false;
        }

        context.frameStartMin.resize(n);
        context.frameEndMin.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            if (!std::isfinite((*frameStartSec)[i]) ||
                !std::isfinite((*frameEndSec)[i]) ||
                (*frameEndSec)[i] <= (*frameStartSec)[i])
            {
                if (errorMessage)
                {
                    *errorMessage =
                        "Invalid PET frame interval supplied to Feng fitting.";
                }
                return false;
            }
            context.frameStartMin[i] = (*frameStartSec)[i] / 60.0;
            context.frameEndMin[i] = (*frameEndSec)[i] / 60.0;
        }
    }
    else
    {
        if (timesSec.size() != n)
        {
            if (errorMessage)
            {
                *errorMessage =
                    "Feng point-sampled fitting requires one time per observation.";
            }
            return false;
        }

        context.timesMin.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            if (!std::isfinite(timesSec[i]))
            {
                if (errorMessage)
                {
                    *errorMessage =
                        "Non-finite sample time supplied to Feng fitting.";
                }
                return false;
            }
            if (i > 0 && timesSec[i] <= timesSec[i - 1])
            {
                if (errorMessage)
                {
                    *errorMessage =
                        "Feng sample times must be strictly increasing.";
                }
                return false;
            }
            context.timesMin[i] = timesSec[i] / 60.0;
        }
    }

    double peak = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
        {
            if (errorMessage)
            {
                *errorMessage =
                    "Non-finite input-function value supplied to Feng fitting.";
            }
            return false;
        }
        peak = std::max(peak, value);
    }

    if (!(peak > 0.0))
    {
        if (errorMessage)
        {
            *errorMessage =
                "Feng fitting requires a positive input-function peak.";
        }
        return false;
    }

    const double maximumTimeMin = observationsAreFrameAverages
        ? context.frameEndMin.back()
        : context.timesMin.back();

    const double tauLower = -1.0;
    const double tauUpper = std::max(
        0.25,
        std::min(5.0, maximumTimeMin));

    const std::array<std::array<double, 3>, 3> lambdaStarts{{
        {{4.0, 0.5, 0.05}},
        {{10.0, 1.0, 0.10}},
        {{2.0, 0.25, 0.02}}
    }};
    const std::array<double, 3> tauStarts{{0.0, -0.05, 0.10}};

    double bestSSE = std::numeric_limits<double>::infinity();
    std::array<double, 7> bestInternal{};
    std::vector<double> bestFitted(n, 0.0);

    for (double tau0 : tauStarts)
    {
        for (const auto& lambdas : lambdaStarts)
        {
            std::array<double, 7> u{{
                std::max(tauLower, std::min(tauUpper, tau0)),
                std::log(std::max(1e-12, 2.0 * peak)),
                std::log(std::max(1e-12, 0.65 * peak)),
                std::log(std::max(1e-12, 0.35 * peak)),
                std::log(std::max(1e-8, lambdas[0] - lambdas[1])),
                std::log(std::max(1e-8, lambdas[1] - lambdas[2])),
                std::log(std::max(1e-8, lambdas[2]))
            }};

            std::array<double, 7> lower{{
                tauLower, -30.0, -30.0, -30.0,
                -20.0, -20.0, -20.0
            }};
            std::array<double, 7> upper{{
                tauUpper, 30.0, 30.0, 30.0,
                10.0, 10.0, 10.0
            }};
            std::array<int, 7> sensitive{{1, 1, 1, 1, 1, 1, 1}};
            std::vector<double> weights(n, 1.0);
            std::vector<double> predicted(n, 0.0);

            std::vector<double> observations = values;

            kmap_levmar(
                observations.data(),
                weights.data(),
                static_cast<int>(n),
                u.data(),
                7,
                &context,
                fengEvaluateInternal,
                fengJacobianInternal,
                lower.data(),
                upper.data(),
                sensitive.data(),
                300,
                predicted.data());

            fengEvaluateInternal(
                u.data(),
                &context,
                predicted.data());

            double sse = 0.0;
            bool finite = true;
            for (size_t i = 0; i < n; ++i)
            {
                if (!std::isfinite(predicted[i]))
                {
                    finite = false;
                    break;
                }
                const double residual = values[i] - predicted[i];
                sse += residual * residual;
            }

            if (finite && sse < bestSSE)
            {
                bestSSE = sse;
                bestInternal = u;
                bestFitted = predicted;
            }
        }
    }

    if (!std::isfinite(bestSSE))
    {
        if (errorMessage)
        {
            *errorMessage =
                "Feng optimization did not produce a finite solution.";
        }
        return false;
    }

    fengInternalToPhysical(bestInternal.data(), params);
    params.SSE = bestSSE;
    fittedObservationValues = std::move(bestFitted);
    return true;
}

//----------------------------------------------------------------------------
double vtkSlicerDynamicPETLogic::EvaluateFengInputFunction(
    double timeSec,
    const FengParameters& params) const
{
    return fengValueMin(timeSec / 60.0, params);
}

//----------------------------------------------------------------------------
double vtkSlicerDynamicPETLogic::AverageFengInputFunction(
    double frameStartSec,
    double frameEndSec,
    const FengParameters& params) const
{
    if (!(frameEndSec > frameStartSec))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double startMin = frameStartSec / 60.0;
    const double endMin = frameEndSec / 60.0;
    return
        (fengCumulativeMin(endMin, params) -
         fengCumulativeMin(startMin, params)) /
        (endMin - startMin);
}

//----------------------------------------------------------------------------
bool vtkSlicerDynamicPETLogic::FitParentFraction(
    const std::vector<double>& timesSec,
    const std::vector<double>& values,
    ParentFractionModel model,
    ParentFractionFitParameters& params,
    std::vector<double>& fittedObservationValues,
    std::string* errorMessage)
{
    params = ParentFractionFitParameters{};
    fittedObservationValues.clear();

    if (model == ParentFractionModel::Linear)
    {
        if (errorMessage)
            *errorMessage = "Linear parent-fraction interpolation does not require parametric fitting.";
        return false;
    }

    if (timesSec.size() != values.size())
    {
        if (errorMessage)
            *errorMessage = "Parent-fraction fitting requires one time for every fraction value.";
        return false;
    }

    const int pCount = parentFractionParameterCount(model);
    const size_t minimumObservations =
        model == ParentFractionModel::ExtendedHill ? 6u :
        model == ParentFractionModel::Hill ? 4u : 3u;

    if (values.size() < minimumObservations)
    {
        if (errorMessage)
        {
            *errorMessage =
                "Not enough retained parent-fraction observations for the selected model.";
        }
        return false;
    }

    ParentFractionFitContext context;
    context.model = model;
    context.timesMin.resize(timesSec.size());

    for (size_t i = 0; i < timesSec.size(); ++i)
    {
        if (!std::isfinite(timesSec[i]) || !std::isfinite(values[i]) ||
            values[i] < 0.0 || values[i] > 1.0)
        {
            if (errorMessage)
                *errorMessage = "Parent-fraction samples must be finite and between 0 and 1.";
            return false;
        }
        if (i > 0 && timesSec[i] <= timesSec[i - 1])
        {
            if (errorMessage)
                *errorMessage = "Parent-fraction sample times must be strictly increasing.";
            return false;
        }
        context.timesMin[i] = timesSec[i] / 60.0;
    }

    context.maximumTimeMin = std::max(1e-3, context.timesMin.back());

    std::vector<std::vector<double>> starts;
    std::vector<double> lower(static_cast<size_t>(pCount), -30.0);
    std::vector<double> upper(static_cast<size_t>(pCount), 30.0);

    const double first = std::max(1e-4, std::min(0.9999, values.front()));
    const double last = std::max(1e-4, std::min(0.9999, values.back()));
    const double tScale = std::max(
        1e-3,
        context.timesMin[context.timesMin.size() / 2]);

    if (model == ParentFractionModel::Hill)
    {
        lower = {-12.0, -6.0, -20.0};
        upper = { 12.0,  3.0,  20.0};
        for (double b0 : {1.5, 2.0, 3.0})
        {
            const double c0 = std::max(1e-6, std::pow(tScale, b0));
            starts.push_back({
                parentFractionLogit(last),
                std::log(std::max(1e-6, b0 - 1.0)),
                std::log(c0)});
        }
    }
    else if (model == ParentFractionModel::ExtendedHill)
    {
        lower = {-12.0, -12.0, -6.0, -20.0, 0.0};
        upper = { 12.0,  12.0,  3.0,  20.0, context.maximumTimeMin};
        const double d0 = std::max(first, last + 1e-4);
        const double ratio0 = std::max(1e-4, std::min(0.9999, last / d0));
        for (double b0 : {1.5, 2.0, 3.0})
        {
            const double c0 = std::max(1e-6, std::pow(tScale, b0));
            for (double e0 : {0.0, 0.05 * context.maximumTimeMin})
            {
                starts.push_back({
                    parentFractionLogit(d0),
                    parentFractionLogit(ratio0),
                    std::log(std::max(1e-6, b0 - 1.0)),
                    std::log(c0),
                    std::min(context.maximumTimeMin, e0)});
            }
        }
    }
    else if (model == ParentFractionModel::ExponentialPlateau)
    {
        lower = {-12.0, -20.0};
        upper = { 12.0,  10.0};
        const double baseRate = 1.0 / context.maximumTimeMin;
        for (double multiplier : {0.25, 1.0, 4.0, 10.0})
        {
            starts.push_back({
                parentFractionLogit(last),
                std::log(std::max(1e-8, multiplier * baseRate))});
        }
    }

    double bestSSE = std::numeric_limits<double>::infinity();
    std::vector<double> bestInternal;
    std::vector<double> bestFitted(values.size(), 0.0);

    for (std::vector<double> u : starts)
    {
        std::vector<int> sensitive(static_cast<size_t>(pCount), 1);
        std::vector<double> weights(values.size(), 1.0);
        std::vector<double> predicted(values.size(), 0.0);
        std::vector<double> observations = values;

        kmap_levmar(
            observations.data(),
            weights.data(),
            static_cast<int>(observations.size()),
            u.data(),
            pCount,
            &context,
            parentFractionEvaluateInternal,
            parentFractionJacobianInternal,
            lower.data(),
            upper.data(),
            sensitive.data(),
            300,
            predicted.data());

        parentFractionEvaluateInternal(u.data(), &context, predicted.data());

        double sse = 0.0;
        bool finite = true;
        for (size_t i = 0; i < values.size(); ++i)
        {
            if (!std::isfinite(predicted[i]))
            {
                finite = false;
                break;
            }
            const double residual = values[i] - predicted[i];
            sse += residual * residual;
        }

        if (finite && sse < bestSSE)
        {
            bestSSE = sse;
            bestInternal = u;
            bestFitted = predicted;
        }
    }

    if (!std::isfinite(bestSSE) || bestInternal.empty())
    {
        if (errorMessage)
            *errorMessage = "Parent-fraction optimization did not produce a finite solution.";
        return false;
    }

    parentFractionInternalToPhysical(bestInternal.data(), context, params);
    params.SSE = bestSSE;
    params.numberOfObservations = static_cast<int>(values.size());
    fittedObservationValues = std::move(bestFitted);
    return true;
}

//----------------------------------------------------------------------------
double vtkSlicerDynamicPETLogic::EvaluateParentFraction(
    double timeSec,
    ParentFractionModel model,
    const ParentFractionFitParameters& params) const
{
    if (!std::isfinite(timeSec) || timeSec < 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    return parentFractionValueMin(timeSec / 60.0, model, params);
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
  // Fetch CT when available. TAC extraction itself remains PET-centered;
  // CT is optional and is used only to quantify the same ROI on the CT grid.
  vtkMRMLScalarVolumeNode* ctNode = nullptr;
  if (ctID != vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID)
  {
    ctNode = vtkMRMLScalarVolumeNode::SafeDownCast(
        shNode->GetItemDataNode(ctID));
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

          if (stats.count == 1)
          {
            stats.min = val;
            stats.max = val;
            max_ijk[0] = x;
            max_ijk[1] = y;
            max_ijk[2] = z;
          }
          else
          {
            stats.min = std::min(stats.min, val);
            if (val > stats.max)
            {
              stats.max = val;
              max_ijk[0] = x;
              max_ijk[1] = y;
              max_ijk[2] = z;
            }
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

  if (stats.count == 0)
  {
    stats.keep = false;
    stats.empty = true;
    stats.mean = std::numeric_limits<double>::quiet_NaN();
    stats.median = std::numeric_limits<double>::quiet_NaN();
    stats.min = std::numeric_limits<double>::quiet_NaN();
    stats.max = std::numeric_limits<double>::quiet_NaN();
    stats.stddev = std::numeric_limits<double>::quiet_NaN();
    stats.q1 = std::numeric_limits<double>::quiet_NaN();
    stats.q3 = std::numeric_limits<double>::quiet_NaN();
    stats.iqr = std::numeric_limits<double>::quiet_NaN();
    stats.peak = std::numeric_limits<double>::quiet_NaN();
    stats.peakStddev = std::numeric_limits<double>::quiet_NaN();
    stats.peakCount = 0;
    return stats;
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
  std::vector<double> peakValues;
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
    peakValues.push_back(val);
  }

  if (countPeak > 0)
  {
    stats.peak = sumPeak / countPeak;
    stats.peakCount = countPeak;

    double peakVariance = 0.0;
    for (const double value : peakValues)
    {
      const double delta = value - stats.peak;
      peakVariance += delta * delta;
    }

    stats.peakStddev =
        std::sqrt(peakVariance / countPeak);
  }
  else
  {
    stats.peak = std::numeric_limits<double>::quiet_NaN();
    stats.peakStddev = std::numeric_limits<double>::quiet_NaN();
    stats.peakCount = 0;
  }

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

  // Liver DBIF changes the input-function construction itself (arterial plus
  // portal contribution). Even though it shares several kinetic parameter
  // names with 2TCM-family models, it is not a nested restriction of them.
  // Comparisons involving Liver DBIF are therefore explicitly non-nested
  // and use Vuong on the common tissue-response residual scale.
  if (modelA == "Liver DBIF" || modelB == "Liver DBIF")
  {
    if (m1.r.size() != m2.r.size() ||
        m1.weights.size() != m2.weights.size() ||
        m1.weights.size() != m1.r.size())
    {
      throw std::invalid_argument(
          "Liver DBIF Vuong comparison requires matching residual/weight vectors.");
    }

    std::vector<double> wgt(m1.weights.size(), 1.0);
    for (size_t i = 0; i < wgt.size(); ++i)
    {
      wgt[i] = 0.5 * (m1.weights[i] + m2.weights[i]);
    }

    res.type = "Vuong";
    res.p_value = this->computeVuongP(
        m1.r,
        m2.r,
        &wgt,
        m1.dof,
        m2.dof,
        VuongCorrection::BIC,
        Tail::TwoSided);
    return res;
  }

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
    const std::string segmentID = id.toStdString();

    // Always allocate the full TAC first. A newly added segment may be absent
    // from early segmentation-sequence frames; writing an invalid frame into an
    // unallocated vector was the cause of the move-between-segmentations crash.
    segmentTACs[segmentID].assign(
        numberOfTimepoints,
        VoxelStatistics{});

    bool foundSegment = false;

    for (int i = 0; i < numberOfTimepoints; ++i)
    {
      const std::string indexValue =
          sequencePETNode->GetNthIndexValue(i);

      vtkMRMLSegmentationNode* segmentationNode =
          vtkMRMLSegmentationNode::SafeDownCast(
              segSequenceNode->GetDataNodeAtValue(indexValue));

      if (!segmentationNode ||
          !segmentationNode->GetSegmentation())
      {
        continue;
      }

      vtkSegment* segment =
          segmentationNode->GetSegmentation()->GetSegment(segmentID);

      if (!segment)
      {
        continue;
      }

      segmentTACsnames[segmentID] =
          segment->GetName() ? segment->GetName() : segmentID;
      foundSegment = true;
      break;
    }

    if (!foundSegment)
    {
      std::cerr
          << "Segment not found in any segmentation timepoint: "
          << segmentID << std::endl;
      segmentTACsnames.erase(segmentID);
    }
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

    // MRML/segmentation export is intentionally serialized.
    // Heavy numerical work may be parallelized, but scene-dependent MRML
    // segmentation operations must not run concurrently here.
    for (int s = 0; s < segmentsID.size(); ++s)
    {
      if (stopRequested) continue;
      const std::string& segmentID = segmentsID[s].toStdString();
      VoxelStatistics stats;

      vtkSegmentation* segmentation = segmentationNode->GetSegmentation();
      if (!segmentation || !segmentation->GetSegment(segmentID))
      {
        std::cerr
            << "Segment not present at timepoint " << i
            << ": " << segmentID << std::endl;
        stats.keep = false;
        stats.empty = true;
        segmentTACs[segmentID][i] = stats;
        continue;
      }

      vtkNew<vtkStringArray> segmentArray;
      segmentArray->InsertNextValue(segmentID);

      vtkSmartPointer<vtkOrientedImageData> labelmap =
          vtkSmartPointer<vtkOrientedImageData>::New();

      vtkSlicerSegmentationsModuleLogic::
          GenerateMergedLabelmapInReferenceGeometry(
              segmentationNode,
              PETVolume,
              segmentArray,
              vtkSegmentation::EXTENT_UNION_OF_EFFECTIVE_SEGMENTS,
              labelmap);

      if (!labelmap ||
          !labelmap->GetPointData() ||
          !labelmap->GetPointData()->GetScalars())
      {
        std::cerr
            << "Failed to generate labelmap for segment: "
            << segmentID
            << " at timepoint " << i
            << std::endl;
        stats.keep = false;
        stats.empty = true;
      }
      else
      {
        stats = ComputeVoxelStatistics(
            PETVolume,
            labelmap,
            1);

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


void vtkSlicerDynamicPETLogic::callTCM(
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
    const std::string& interpolationType,
    const std::vector<double>* nativePlasmaTimesSec,
    const std::vector<double>* nativePlasmaValues,
    const std::vector<double>* nativeWholeBloodTimesSec,
    const std::vector<double>* nativeWholeBloodValues,
    const std::vector<double>* parentFractionTimesSec,
    const std::vector<double>* parentFractionValues,
    bool plasmaIsParent,
    double acquisitionStartSec)
{
  const int nth = 1;

  // Basic validation
  if (containsNaN(tac)) return error_nan("TAC");
  if (containsNaN(Cp)) return error_nan("Cp");
  if (containsNaN(Cwb)) return error_nan("Cwb");
  if (containsNaN(framing)) return error_nan("framing");

  if (tac.size() != framing.size()) return error_size("TAC", "framing", tac.size(), framing.size());
  if (tac.size() != Cp.size()) return error_size("TAC", "Cp", tac.size(), Cp.size());
  if (tac.size() != Cwb.size()) return error_size("TAC", "Cwb", tac.size(), Cwb.size());

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

  double** scant =
      new double*[Nframe];

  for (long int i = 0;
       i < Nframe;
       ++i)
  {
      scant[i] = new double[2];
      scant[i][0] =
          acquisitionStartSec +
          ((i == 0)
           ? 0.0
           : cumsum[i - 1]);
      scant[i][1] =
          acquisitionStartSec + cumsum[i];
  }

  const double scanEndSec =
      acquisitionStartSec + cumsum[Nframe - 1];

  long int N_cp = 0;
  long int N_wb = 0;

  double* Cp_new = nullptr;
  double* cwb_new = nullptr;

  const bool useNativePlasma =
      nativePlasmaTimesSec &&
      nativePlasmaValues &&
      nativePlasmaTimesSec->size() >= 2 &&
      nativePlasmaTimesSec->size() ==
          nativePlasmaValues->size();

  if (useNativePlasma)
  {
      Cp_new =
          FineSampleExplicitInputFunction(
              *nativePlasmaTimesSec,
              *nativePlasmaValues,
              scanEndSec,
              timestep,
              interpolationType,
              N_cp);
  }
  else
  {
      if (interpolationType == "pchip")
      {
          std::vector<double> frameValues(Nframe);
          for (long int i = 0; i < Nframe; ++i) frameValues[i] = Cp[i][0];
          std::vector<double> times, values;
          BuildFrameRepresentativeCurve(scant, frameValues, times, values);
          Cp_new = FineSampleExplicitInputFunction(times, values, scanEndSec, timestep, "pchip", N_cp);
      }
      else
      {
          Cp_new = finesample(scant, Cp, Nframe, N_cp, timestep, interpolationType);
      }
  }

  const bool useNativeWholeBlood =
      nativeWholeBloodTimesSec &&
      nativeWholeBloodValues &&
      nativeWholeBloodTimesSec->size() >= 2 &&
      nativeWholeBloodTimesSec->size() ==
          nativeWholeBloodValues->size();

  if (useNativeWholeBlood)
  {
      cwb_new =
          FineSampleExplicitInputFunction(
              *nativeWholeBloodTimesSec,
              *nativeWholeBloodValues,
              scanEndSec,
              timestep,
              interpolationType,
              N_wb);
  }
  else
  {
      if (interpolationType == "pchip")
      {
          std::vector<double> frameValues(Nframe);
          for (long int i = 0; i < Nframe; ++i) frameValues[i] = Cwb[i][0];
          std::vector<double> times, values;
          BuildFrameRepresentativeCurve(scant, frameValues, times, values);
          cwb_new = FineSampleExplicitInputFunction(times, values, scanEndSec, timestep, "pchip", N_wb);
      }
      else
      {
          cwb_new = finesample(scant, Cwb, Nframe, N_wb, timestep, interpolationType);
      }
  }

  if (N_cp != N_wb)
  {
      delete[] Cp_new;
      delete[] cwb_new;

      throw std::runtime_error(
          "Plasma and whole-blood fine-sampling grids differ.");
  }

  const bool applyParentFraction =
      !plasmaIsParent &&
      parentFractionTimesSec &&
      parentFractionValues &&
      parentFractionTimesSec->size() >= 2 &&
      parentFractionTimesSec->size() ==
          parentFractionValues->size();

  if (applyParentFraction)
  {
      long int N_parent = 0;

      double* parentFractionNew =
          FineSampleExplicitInputFunction(
              *parentFractionTimesSec,
              *parentFractionValues,
              scanEndSec,
              timestep,
              "linear",
              N_parent);

      if (N_parent != N_cp)
      {
          delete[] parentFractionNew;
          delete[] Cp_new;
          delete[] cwb_new;

          throw std::runtime_error(
              "Parent-fraction and plasma fine-sampling grids differ.");
      }

      for (long int i = 0;
           i < N_cp;
           ++i)
      {
          Cp_new[i] *= parentFractionNew[i];
      }

      delete[] parentFractionNew;
  }

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

    params.boundFlags = TCM_BOUND_NONE;
    markBoundIfNeeded(params.boundFlags, params.vb, lb[0], ub[0], sens[0], TCM_BOUND_VB_LOWER, TCM_BOUND_VB_UPPER);
    markBoundIfNeeded(params.boundFlags, params.K1, lb[1], ub[1], sens[1], TCM_BOUND_K1_LOWER, TCM_BOUND_K1_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k2, lb[2], ub[2], sens[2], TCM_BOUND_K2_LOWER, TCM_BOUND_K2_UPPER);
    markBoundIfNeeded(params.boundFlags, params.td, lb[3], ub[3], sens[3], TCM_BOUND_TD_LOWER, TCM_BOUND_TD_UPPER);

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

    params.boundFlags = TCM_BOUND_NONE;
    markBoundIfNeeded(params.boundFlags, params.vb, lb[0], ub[0], sens[0], TCM_BOUND_VB_LOWER, TCM_BOUND_VB_UPPER);
    markBoundIfNeeded(params.boundFlags, params.K1, lb[1], ub[1], sens[1], TCM_BOUND_K1_LOWER, TCM_BOUND_K1_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k2, lb[2], ub[2], sens[2], TCM_BOUND_K2_LOWER, TCM_BOUND_K2_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k3, lb[3], ub[3], sens[3], TCM_BOUND_K3_LOWER, TCM_BOUND_K3_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k4, lb[4], ub[4], sens[4], TCM_BOUND_K4_LOWER, TCM_BOUND_K4_UPPER);
    markBoundIfNeeded(params.boundFlags, params.td, lb[5], ub[5], sens[5], TCM_BOUND_TD_LOWER, TCM_BOUND_TD_UPPER);

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


void vtkSlicerDynamicPETLogic::callLiverTCM(
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
    const std::string& interpolationType,
    const std::vector<double>* nativeWholeBloodTimesSec,
    const std::vector<double>* nativeWholeBloodValues,
    double acquisitionStartSec)
{
    constexpr double EPS = 1e-16;
    constexpr int numberOfParameters = 8;

    if (Nframe <= 0 || timestep <= 0.0)
    {
        throw std::invalid_argument(
            "Invalid liver DBIF fitting dimensions or integration timestep.");
    }

    if (tac.size() != static_cast<size_t>(Nframe) ||
        Cwb.size() != static_cast<size_t>(Nframe) ||
        framing.size() != static_cast<size_t>(Nframe))
    {
        throw std::invalid_argument(
            "Liver DBIF TAC, whole-blood input, and framing sizes must match.");
    }

    if (containsNaN(tac)) return error_nan("TAC");
    if (containsNaN(Cwb)) return error_nan("Cwb");
    if (containsNaN(framing)) return error_nan("framing");

    std::vector<double> weights(
        static_cast<size_t>(Nframe),
        1.0);

    if (wgt)
    {
        if (wgt->size() != static_cast<size_t>(Nframe))
        {
            throw std::runtime_error(
                "Weight vector length does not match Nframe.");
        }

        weights = *wgt;
    }

    params.keep.resize(
        static_cast<size_t>(Nframe));

    for (long int i = 0; i < Nframe; ++i)
    {
        params.keep[static_cast<size_t>(i)] =
            weights[static_cast<size_t>(i)] != 0.0;
    }

    std::vector<double> cumulative(
        static_cast<size_t>(Nframe),
        0.0);

    double accumulatedTime = 0.0;

    double** scant =
        new double*[Nframe];

    for (long int i = 0; i < Nframe; ++i)
    {
        if (framing[static_cast<size_t>(i)].empty() ||
            framing[static_cast<size_t>(i)][0] <= 0.0)
        {
            for (long int j = 0; j < i; ++j)
            {
                delete[] scant[j];
            }

            delete[] scant;

            throw std::invalid_argument(
                "Liver DBIF frame durations must be positive.");
        }

        accumulatedTime +=
            framing[static_cast<size_t>(i)][0];

        cumulative[static_cast<size_t>(i)] =
            accumulatedTime;

        scant[i] = new double[2];

        scant[i][0] =
            acquisitionStartSec +
            ((i == 0)
             ? 0.0
             : cumulative[static_cast<size_t>(i - 1)]);

        scant[i][1] =
            acquisitionStartSec +
            cumulative[static_cast<size_t>(i)];
    }

    const double scanEndSec =
        acquisitionStartSec + cumulative.back();

    long int numberOfFineSamples = 0;
    double* arterialInputFine = nullptr;

    const bool useNativeWholeBlood =
        nativeWholeBloodTimesSec &&
        nativeWholeBloodValues &&
        nativeWholeBloodTimesSec->size() >= 2 &&
        nativeWholeBloodTimesSec->size() ==
            nativeWholeBloodValues->size();

    if (useNativeWholeBlood)
    {
        arterialInputFine =
            FineSampleExplicitInputFunction(
                *nativeWholeBloodTimesSec,
                *nativeWholeBloodValues,
                scanEndSec,
                timestep,
                interpolationType,
                numberOfFineSamples);
    }
    else
    {
        // Segment-derived IDIFs intentionally retain the established
        // frame-aware fine-sampling pathway.
        std::vector<std::vector<double>>
            wholeBloodFrameValues = Cwb;

        arterialInputFine =
            finesample(
                scant,
                wholeBloodFrameValues,
                Nframe,
                numberOfFineSamples,
                timestep,
                interpolationType);
    }

    if (!arterialInputFine ||
        numberOfFineSamples <= 0)
    {
        for (long int i = 0; i < Nframe; ++i)
        {
            delete[] scant[i];
        }

        delete[] scant;

        throw std::runtime_error(
            "Could not fine-sample the arterial whole-blood input "
            "for the liver DBIF model.");
    }

    std::vector<double> scantFlattened(
        static_cast<size_t>(Nframe) * 2);

    for (long int i = 0; i < Nframe; ++i)
    {
        scantFlattened[static_cast<size_t>(i)] =
            scant[i][0];

        scantFlattened[static_cast<size_t>(i + Nframe)] =
            scant[i][1];
    }

    std::vector<double> tacValues(
        static_cast<size_t>(Nframe));

    for (long int i = 0; i < Nframe; ++i)
    {
        if (tac[static_cast<size_t>(i)].empty())
        {
            delete[] arterialInputFine;

            for (long int j = 0; j < Nframe; ++j)
            {
                delete[] scant[j];
            }

            delete[] scant;

            throw std::invalid_argument(
                "Liver DBIF TAC contains an empty frame.");
        }

        tacValues[static_cast<size_t>(i)] =
            tac[static_cast<size_t>(i)][0];
    }

    std::vector<double> fittedParameters(
        numberOfParameters);

    std::vector<int> parameterSensitivity(
        numberOfParameters);

    for (int i = 0; i < numberOfParameters; ++i)
    {
        fittedParameters[static_cast<size_t>(i)] =
            kinit[i];

        parameterSensitivity[static_cast<size_t>(i)] =
            sens[i] ? 1 : 0;
    }

    KMODEL_T model;
    model.dk = dk;
    model.td = timestep;

    // kconv_liver_* names this argument "ca". It is the
    // arterial/aortic whole-blood input from which the portal component
    // and effective dual input are generated.
    model.cp = arterialInputFine;

    // The generic callback signature also contains wb. Supply a valid
    // curve without changing the established KMAP API.
    model.wb = arterialInputFine;

    model.num_frm =
        static_cast<int>(Nframe);

    model.num_vox = 1;

    model.scant =
        scantFlattened.data();

    model.tacfunc =
        kconv_liver_tac;

    model.jacfunc =
        kconv_liver_jac;

    fitted_curve =
        new double[Nframe];

    kmap_levmar(
        tacValues.data(),
        weights.data(),
        static_cast<int>(Nframe),
        fittedParameters.data(),
        numberOfParameters,
        &model,
        tac_eval,
        jac_eval,
        lb,
        ub,
        parameterSensitivity.data(),
        maxiter,
        fitted_curve);

    params.vb = sens[0]
        ? fittedParameters[0]
        : std::numeric_limits<double>::quiet_NaN();

    params.K1 = sens[1]
        ? fittedParameters[1]
        : std::numeric_limits<double>::quiet_NaN();

    params.k2 = sens[2]
        ? fittedParameters[2]
        : std::numeric_limits<double>::quiet_NaN();

    params.k3 = sens[3]
        ? fittedParameters[3]
        : std::numeric_limits<double>::quiet_NaN();

    params.k4 = sens[4]
        ? fittedParameters[4]
        : std::numeric_limits<double>::quiet_NaN();

    params.ka = sens[5]
        ? fittedParameters[5]
        : std::numeric_limits<double>::quiet_NaN();

    params.fa = sens[6]
        ? fittedParameters[6]
        : std::numeric_limits<double>::quiet_NaN();

    params.td = sens[7]
        ? fittedParameters[7]
        : std::numeric_limits<double>::quiet_NaN();

    params.boundFlags = TCM_BOUND_NONE;
    markBoundIfNeeded(params.boundFlags, params.vb, lb[0], ub[0], sens[0], TCM_BOUND_VB_LOWER, TCM_BOUND_VB_UPPER);
    markBoundIfNeeded(params.boundFlags, params.K1, lb[1], ub[1], sens[1], TCM_BOUND_K1_LOWER, TCM_BOUND_K1_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k2, lb[2], ub[2], sens[2], TCM_BOUND_K2_LOWER, TCM_BOUND_K2_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k3, lb[3], ub[3], sens[3], TCM_BOUND_K3_LOWER, TCM_BOUND_K3_UPPER);
    markBoundIfNeeded(params.boundFlags, params.k4, lb[4], ub[4], sens[4], TCM_BOUND_K4_LOWER, TCM_BOUND_K4_UPPER);
    markBoundIfNeeded(params.boundFlags, params.ka, lb[5], ub[5], sens[5], TCM_BOUND_KA_LOWER, TCM_BOUND_KA_UPPER);
    markBoundIfNeeded(params.boundFlags, params.fa, lb[6], ub[6], sens[6], TCM_BOUND_FA_LOWER, TCM_BOUND_FA_UPPER);
    markBoundIfNeeded(params.boundFlags, params.td, lb[7], ub[7], sens[7], TCM_BOUND_TD_LOWER, TCM_BOUND_TD_UPPER);

    params.Ki =
        params.K1 *
        params.k3 /
        (params.k2 +
         params.k3 +
         EPS);

    params.DV =
        params.K1 /
        (params.k2 + EPS) *
        (1.0 +
         params.k3 /
         (params.k4 + EPS));

    params.weights =
        weights;

    params.dof =
        std::count(
            sens,
            sens + numberOfParameters,
            true);

    std::vector<double> predicted(
        fitted_curve,
        fitted_curve + Nframe);

    params.AIC =
        this->computeAIC(
            tacValues,
            predicted,
            params.dof,
            wgt);

    params.BIC =
        this->computeBIC(
            tacValues,
            predicted,
            params.dof,
            wgt);

    params.MASE =
        this->MASE(
            tacValues,
            predicted,
            wgt);

    if (Nframe > params.dof)
    {
        params.chi2 =
            this->computeChi2(
                tacValues,
                predicted,
                wgt) /
            static_cast<double>(
                Nframe -
                params.dof);
    }
    else
    {
        params.chi2 =
            std::numeric_limits<double>::quiet_NaN();
    }

    params.loglik =
        this->computeLogLik(
            tacValues,
            predicted,
            wgt);

    params.r.resize(
        static_cast<size_t>(Nframe));

    for (long int i = 0; i < Nframe; ++i)
    {
        params.r[static_cast<size_t>(i)] =
            tacValues[static_cast<size_t>(i)] -
            predicted[static_cast<size_t>(i)];
    }

    delete[] arterialInputFine;

    for (long int i = 0; i < Nframe; ++i)
    {
        delete[] scant[i];
    }

    delete[] scant;
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

    int numpar_new = numpar + 1; // fitted kinetic parameters + residual-variance scale

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

    int p = numpar + 1; // fitted kinetic parameters + residual-variance scale
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
                                int max_iter,
                                double initialPlasmaIntegral
                              )
{
  size_t N = tac.size();
  std::vector<double> outX, outY, fittedValues;
  std::vector<int> outframe;

  // normalize framing to minutes
  std::vector<double> frameScaled(N);
  for (size_t i = 0; i < N; ++i)
      frameScaled[i] = framing[i] / framingNorm;

  // PET values are frame averages. Associate each observation with its
  // frame midpoint and integrate frame-average values directly. Full-frame
  // areas are known exactly; only the current half-frame is approximated.
  std::vector<double> timeAlong(N, 0.0);
  std::vector<double> intCp(N, 0.0);
  double elapsed = 0.0;
  double accumulatedCp = initialPlasmaIntegral;
  for (size_t i = 0; i < N; ++i)
  {
      timeAlong[i] = elapsed + 0.5 * frameScaled[i];
      intCp[i] = accumulatedCp + 0.5 * Cp[i] * frameScaled[i];
      accumulatedCp += Cp[i] * frameScaled[i];
      elapsed += frameScaled[i];
  }

  // build X, Y with time filter
  std::vector<double> wgt_adj;
  params.keep.clear();
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset &&
          (!wgt || (*wgt)[i] > 0.0))
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
      stdX = std::sqrt(std::max(0.0, std::inner_product(outX.begin(), outX.end(), outX.begin(), 0.0) / n - meanX*meanX));
      stdY = std::sqrt(std::max(0.0, std::inner_product(outY.begin(), outY.end(), outY.begin(), 0.0) / n - meanY*meanY));
      if (stdX < 1e-12) stdX = 1.0;
      if (stdY < 1e-12) stdY = 1.0;

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
  std::vector<double> baseFitWeights(n, 1.0);
  if (!wgt_adj.empty())
  {
      baseFitWeights = wgt_adj;
  }

  std::vector<double> finalFitWeights = baseFitWeights;

  if (robust)
  {
      if (!solveHuberIRLS(
              A,
              Yv,
              baseFitWeights,
              huber_tune,
              tol,
              max_iter,
              coeff,
              finalFitWeights))
      {
          throw std::runtime_error(
              "Huber robust regression failed.");
      }
  }
  else if (!wgt_adj.empty())
  {
      Eigen::Map<const Eigen::VectorXd> Wv(
          wgt_adj.data(),
          n);
      const Eigen::MatrixXd W = Wv.asDiagonal();
      coeff =
          (A.transpose() * W * A)
              .ldlt()
              .solve(A.transpose() * W * Yv);
  }
  else
  {
      coeff =
          (A.transpose() * A)
              .ldlt()
              .solve(A.transpose() * Yv);
  }

  double intercept = coeff(0);
  double slope = coeff(1);

  fittedValues.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
      fittedValues[i] =
          intercept + slope * outX[i];
  }

  // Keep diagnostics on the standardized graphical-analysis scale when
  // standardization is enabled. The reported kinetic slope/intercept and
  // plotted X/Y values are converted back to their original coordinates.
  const std::vector<double> diagnosticY = outY;
  const std::vector<double> diagnosticFitted = fittedValues;

  // Convert the regression parameters and displayed vectors back to
  // the original graphical-analysis scale.
  if (std)
  {
      const double slopeRaw =
          slope * stdY / stdX;
      const double interceptRaw =
          meanY + stdY * intercept - slopeRaw * meanX;

      slope = slopeRaw;
      intercept = interceptRaw;

      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = outX[i] * stdX + meanX;
          outY[i] = outY[i] * stdY + meanY;
          fittedValues[i] =
              fittedValues[i] * stdY + meanY;
      }
  }

  // Huber weights are adaptive fitting weights, not a measurement-variance
  // model.  For AIC/R2/MASE/chi-square diagnostics evaluate the final robust
  // solution using the original frame weights/mask.
  const std::vector<double>* diagnosticWeights =
      baseFitWeights.empty()
      ? nullptr
      : &baseFitWeights;

  params.AIC =
      computeAIC(
          diagnosticY,
          diagnosticFitted,
          2,
          diagnosticWeights);
  params.MASE =
      MASE(
          diagnosticY,
          diagnosticFitted,
          diagnosticWeights);
  params.R2 =
      computeR2(
          diagnosticY,
          diagnosticFitted,
          diagnosticWeights);
  params.chi2 =
      (n > 2)
      ? computeChi2(
            diagnosticY,
            diagnosticFitted,
            diagnosticWeights) /
            static_cast<double>(n - 2)
      : std::numeric_limits<double>::quiet_NaN();

  params.r.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
      params.r[i] =
          diagnosticY[i] - diagnosticFitted[i];
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
  params.weights =
      robust
      ? finalFitWeights
      : baseFitWeights;
}

void vtkSlicerDynamicPETLogic::RelativePatlak(
    const std::vector<double>& tac,
    const std::vector<double>& Cp,
    const std::vector<double>& framing,
    MTGAParameters& params,
    const std::vector<double>* wgt,
    const double timeOffset,
    const double framingNorm,
    bool robust,
    bool std,
    double huber_tune,
    double tol,
    int max_iter,
    size_t dataStartIndex)
{
    const size_t N = tac.size();
    if (Cp.size() != N || framing.size() != N || dataStartIndex >= N)
    {
        throw std::invalid_argument("Invalid Relative Patlak input dimensions.");
    }

    // Relative Patlak (Zuo, Qi & Wang, PMB 2018) restarts the plasma
    // integral at t*.  Here t* is the selected MTGA start frame, constrained
    // to be no earlier than the first jointly available tissue/input frame.
    size_t relativeStartIndex = N;
    double frameStart = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        const double frameMid =
            frameStart + 0.5 * framing[i] / framingNorm;
        if (i >= dataStartIndex && frameMid + 1e-12 >= timeOffset)
        {
            relativeStartIndex = i;
            break;
        }
        frameStart += framing[i] / framingNorm;
    }

    if (relativeStartIndex >= N || N - relativeStartIndex < 2)
    {
        throw std::invalid_argument(
            "Relative Patlak requires at least two frames at or after the selected start time.");
    }

    std::vector<double> tacSub(
        tac.begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), tac.end());
    std::vector<double> cpSub(
        Cp.begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), Cp.end());
    std::vector<double> framingSub(
        framing.begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), framing.end());

    std::vector<double> weightSub;
    const std::vector<double>* weightPtr = nullptr;
    if (wgt)
    {
        if (wgt->size() != N)
        {
            throw std::invalid_argument("Relative Patlak weight vector size mismatch.");
        }
        weightSub.assign(
            wgt->begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), wgt->end());
        weightPtr = &weightSub;
    }

    // After subsetting, t*=0 on the relative clock and the cumulative plasma
    // integral in Patlak() therefore implements integral_{t*}^{t} Cp(tau)dtau.
    this->Patlak(tacSub, cpSub, framingSub, params, weightPtr,
                 0.0, framingNorm, robust, std,
                 huber_tune, tol, max_iter, 0.0);

    for (int& frame : params.frame)
    {
        frame += static_cast<int>(relativeStartIndex);
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

  // PET values are frame averages. Associate each observation with its
  // frame midpoint and integrate frame-average values directly. Full-frame
  // areas are known exactly; only the current half-frame is approximated.
  std::vector<double> timeAlong(N, 0.0);
  std::vector<double> intCp(N, 0.0);
  double elapsed = 0.0;
  double accumulatedCp = 0.0;
  for (size_t i = 0; i < N; ++i)
  {
      timeAlong[i] = elapsed + 0.5 * frameScaled[i];
      intCp[i] = accumulatedCp + 0.5 * Cp[i] * frameScaled[i];
      accumulatedCp += Cp[i] * frameScaled[i];
      elapsed += frameScaled[i];
  }

  // Cumulative tissue area at each frame midpoint, using the same
  // frame-average convention as for the input function.
  std::vector<double> intCt(N, 0.0);
  double accumulatedCt = 0.0;
  for (size_t i = 0; i < N; ++i)
  {
      intCt[i] = accumulatedCt + 0.5 * tac[i] * frameScaled[i];
      accumulatedCt += tac[i] * frameScaled[i];
  }

  // Build X, Y with time filter
  std::vector<double> wgt_adj;
  params.keep.clear();
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset &&
          (!wgt || (*wgt)[i] > 0.0))
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
      stdX = std::sqrt(std::max(0.0, std::inner_product(outX.begin(), outX.end(), outX.begin(), 0.0) / n - meanX * meanX));
      stdY = std::sqrt(std::max(0.0, std::inner_product(outY.begin(), outY.end(), outY.begin(), 0.0) / n - meanY * meanY));
      if (stdX < 1e-12) stdX = 1.0;
      if (stdY < 1e-12) stdY = 1.0;

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
  std::vector<double> baseFitWeights(n, 1.0);
  if (!wgt_adj.empty())
  {
      baseFitWeights = wgt_adj;
  }

  std::vector<double> finalFitWeights = baseFitWeights;

  if (robust)
  {
      if (!solveHuberIRLS(
              A,
              Yv,
              baseFitWeights,
              huber_tune,
              tol,
              max_iter,
              coeff,
              finalFitWeights))
      {
          throw std::runtime_error(
              "Huber robust regression failed.");
      }
  }
  else if (!wgt_adj.empty())
  {
      Eigen::Map<const Eigen::VectorXd> Wv(
          wgt_adj.data(),
          n);
      const Eigen::MatrixXd W = Wv.asDiagonal();
      coeff =
          (A.transpose() * W * A)
              .ldlt()
              .solve(A.transpose() * W * Yv);
  }
  else
  {
      coeff =
          (A.transpose() * A)
              .ldlt()
              .solve(A.transpose() * Yv);
  }

  double intercept = coeff(0);
  double slope = coeff(1);

  fittedValues.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
      fittedValues[i] =
          intercept + slope * outX[i];
  }

  // Keep diagnostics on the standardized graphical-analysis scale when
  // standardization is enabled. The reported kinetic slope/intercept and
  // plotted X/Y values are converted back to their original coordinates.
  const std::vector<double> diagnosticY = outY;
  const std::vector<double> diagnosticFitted = fittedValues;

  // Convert the regression parameters and displayed vectors back to
  // the original graphical-analysis scale.
  if (std)
  {
      const double slopeRaw =
          slope * stdY / stdX;
      const double interceptRaw =
          meanY + stdY * intercept - slopeRaw * meanX;

      slope = slopeRaw;
      intercept = interceptRaw;

      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = outX[i] * stdX + meanX;
          outY[i] = outY[i] * stdY + meanY;
          fittedValues[i] =
              fittedValues[i] * stdY + meanY;
      }
  }

  // Huber weights are adaptive fitting weights, not a measurement-variance
  // model.  For AIC/R2/MASE/chi-square diagnostics evaluate the final robust
  // solution using the original frame weights/mask.
  const std::vector<double>* diagnosticWeights =
      baseFitWeights.empty()
      ? nullptr
      : &baseFitWeights;

  params.AIC =
      computeAIC(
          diagnosticY,
          diagnosticFitted,
          2,
          diagnosticWeights);
  params.MASE =
      MASE(
          diagnosticY,
          diagnosticFitted,
          diagnosticWeights);
  params.R2 =
      computeR2(
          diagnosticY,
          diagnosticFitted,
          diagnosticWeights);
  params.chi2 =
      (n > 2)
      ? computeChi2(
            diagnosticY,
            diagnosticFitted,
            diagnosticWeights) /
            static_cast<double>(n - 2)
      : std::numeric_limits<double>::quiet_NaN();

  params.r.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
      params.r[i] =
          diagnosticY[i] - diagnosticFitted[i];
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
  params.weights =
      robust
      ? finalFitWeights
      : baseFitWeights;
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

  // PET values are frame averages. Associate each observation with its
  // frame midpoint and integrate frame-average values directly. Full-frame
  // areas are known exactly; only the current half-frame is approximated.
  std::vector<double> timeAlong(N, 0.0);
  std::vector<double> intCp(N, 0.0);
  double elapsed = 0.0;
  double accumulatedCp = 0.0;
  for (size_t i = 0; i < N; ++i)
  {
      timeAlong[i] = elapsed + 0.5 * frameScaled[i];
      intCp[i] = accumulatedCp + 0.5 * Cp[i] * frameScaled[i];
      accumulatedCp += Cp[i] * frameScaled[i];
      elapsed += frameScaled[i];
  }

  // Cumulative tissue area at each frame midpoint, using the same
  // frame-average convention as for the input function.
  std::vector<double> intCt(N, 0.0);
  double accumulatedCt = 0.0;
  for (size_t i = 0; i < N; ++i)
  {
      intCt[i] = accumulatedCt + 0.5 * tac[i] * frameScaled[i];
      accumulatedCt += tac[i] * frameScaled[i];
  }

  // Build X, Y with time filter
  std::vector<double> wgt_adj;
  params.keep.clear();
  for (size_t i = 0; i < N; ++i)
  {
      if (timeAlong[i] >= timeOffset &&
          (!wgt || (*wgt)[i] > 0.0))
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
      stdX = std::sqrt(std::max(0.0, std::inner_product(outX.begin(), outX.end(), outX.begin(), 0.0) / n - meanX * meanX));
      stdY = std::sqrt(std::max(0.0, std::inner_product(outY.begin(), outY.end(), outY.begin(), 0.0) / n - meanY * meanY));
      if (stdX < 1e-12) stdX = 1.0;
      if (stdY < 1e-12) stdY = 1.0;

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
  std::vector<double> baseFitWeights(n, 1.0);
  if (!wgt_adj.empty())
  {
      baseFitWeights = wgt_adj;
  }

  std::vector<double> finalFitWeights = baseFitWeights;

  if (robust)
  {
      if (!solveHuberIRLS(
              A,
              Yv,
              baseFitWeights,
              huber_tune,
              tol,
              max_iter,
              coeff,
              finalFitWeights))
      {
          throw std::runtime_error(
              "Huber robust regression failed.");
      }
  }
  else if (!wgt_adj.empty())
  {
      Eigen::Map<const Eigen::VectorXd> Wv(
          wgt_adj.data(),
          n);
      const Eigen::MatrixXd W = Wv.asDiagonal();
      coeff =
          (A.transpose() * W * A)
              .ldlt()
              .solve(A.transpose() * W * Yv);
  }
  else
  {
      coeff =
          (A.transpose() * A)
              .ldlt()
              .solve(A.transpose() * Yv);
  }

  double intercept = coeff(0);
  double slope = coeff(1);

  fittedValues.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
      fittedValues[i] =
          intercept + slope * outX[i];
  }

  // Keep diagnostics on the standardized graphical-analysis scale when
  // standardization is enabled. The reported kinetic slope/intercept and
  // plotted X/Y values are converted back to their original coordinates.
  const std::vector<double> diagnosticY = outY;
  const std::vector<double> diagnosticFitted = fittedValues;

  // Convert the regression parameters and displayed vectors back to
  // the original graphical-analysis scale.
  if (std)
  {
      const double slopeRaw =
          slope * stdY / stdX;
      const double interceptRaw =
          meanY + stdY * intercept - slopeRaw * meanX;

      slope = slopeRaw;
      intercept = interceptRaw;

      for (size_t i = 0; i < n; ++i)
      {
          outX[i] = outX[i] * stdX + meanX;
          outY[i] = outY[i] * stdY + meanY;
          fittedValues[i] =
              fittedValues[i] * stdY + meanY;
      }
  }

  // Huber weights are adaptive fitting weights, not a measurement-variance
  // model.  For AIC/R2/MASE/chi-square diagnostics evaluate the final robust
  // solution using the original frame weights/mask.
  const std::vector<double>* diagnosticWeights =
      baseFitWeights.empty()
      ? nullptr
      : &baseFitWeights;

  params.AIC =
      computeAIC(
          diagnosticY,
          diagnosticFitted,
          2,
          diagnosticWeights);
  params.MASE =
      MASE(
          diagnosticY,
          diagnosticFitted,
          diagnosticWeights);
  params.R2 =
      computeR2(
          diagnosticY,
          diagnosticFitted,
          diagnosticWeights);
  params.chi2 =
      (n > 2)
      ? computeChi2(
            diagnosticY,
            diagnosticFitted,
            diagnosticWeights) /
            static_cast<double>(n - 2)
      : std::numeric_limits<double>::quiet_NaN();

  params.r.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
      params.r[i] =
          diagnosticY[i] - diagnosticFitted[i];
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
  params.weights =
      robust
      ? finalFitWeights
      : baseFitWeights;
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

void vtkSlicerDynamicPETLogic::RelativeRE(
    const std::vector<double>& tac,
    const std::vector<double>& Cp,
    const std::vector<double>& framing,
    MTGAParameters& params,
    const std::vector<double>* wgt,
    const double timeOffset,
    const double framingNorm,
    bool robust,
    bool std,
    double huber_tune,
    double tol,
    int max_iter,
    size_t dataStartIndex)
{
    const size_t N = tac.size();
    if (Cp.size() != N || framing.size() != N || dataStartIndex >= N)
    {
        throw std::invalid_argument("Invalid Relative RE input dimensions.");
    }

    size_t relativeStartIndex = N;
    double frameStart = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        const double frameMid = frameStart + 0.5 * framing[i] / framingNorm;
        if (i >= dataStartIndex && frameMid + 1e-12 >= timeOffset)
        {
            relativeStartIndex = i;
            break;
        }
        frameStart += framing[i] / framingNorm;
    }
    if (relativeStartIndex >= N || N - relativeStartIndex < 2)
    {
        throw std::invalid_argument("Relative RE requires at least two frames at or after the selected equilibrium start.");
    }

    std::vector<double> tacSub(tac.begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), tac.end());
    std::vector<double> cpSub(Cp.begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), Cp.end());
    std::vector<double> framingSub(framing.begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), framing.end());
    std::vector<double> weightSub;
    const std::vector<double>* weightPtr = nullptr;
    if (wgt)
    {
        if (wgt->size() != N)
        {
            throw std::invalid_argument("Relative RE weight vector size mismatch.");
        }
        weightSub.assign(wgt->begin() + static_cast<std::ptrdiff_t>(relativeStartIndex), wgt->end());
        weightPtr = &weightSub;
    }

    this->RE(tacSub, cpSub, framingSub, params, weightPtr, 0.0, framingNorm,
             robust, std, huber_tune, tol, max_iter);
    for (int& frame : params.frame)
    {
        frame += static_cast<int>(relativeStartIndex);
    }
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
    std::function<void(int)> progressCallback,
    double initialPlasmaIntegral,
    size_t dataStartIndex)
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
    if (framing.size() != N || dataStartIndex >= N) return;
    size_t nVoxels = voxels.size();
    if (nVoxels == 0) return;

    // ---------- 1) Precompute frameScaled, timeAlong, intCp, invCp ----------
    std::vector<double> frameScaled(N);
    for (size_t i = 0; i < N; ++i) frameScaled[i] = framing[i] / framingNorm;

    // PET observations are frame averages and are associated with frame
    // midpoints. Preserve prior full-frame areas, then add half of the
    // current frame to evaluate the cumulative integral at its midpoint.
    if (dataStartIndex >= N) return;
    std::vector<double> timeAlong(N, 0.0);
    std::vector<double> intCp(N, std::numeric_limits<double>::quiet_NaN());
    double elapsed = 0.0;
    for (size_t i = 0; i < dataStartIndex; ++i)
    {
        elapsed += frameScaled[i];
    }
    double accumulatedCp = initialPlasmaIntegral;
    for (size_t i = dataStartIndex; i < N; ++i)
    {
        timeAlong[i] = elapsed + 0.5 * frameScaled[i];
        intCp[i] = accumulatedCp + 0.5 * Cp[i] * frameScaled[i];
        accumulatedCp += Cp[i] * frameScaled[i];
        elapsed += frameScaled[i];
    }

    std::vector<double> invCp(N);
    for (size_t i = 0; i < N; ++i) invCp[i] = 1.0 / (Cp[i] + EPS);

    // ---------- 2) Filter frames >= timeOffset ----------
    std::vector<double> X_all;
    std::vector<int> keepIndex;
    for (size_t i = dataStartIndex; i < N; ++i)
    {
        if (timeAlong[i] >= timeOffset &&
            (!wgt_global || (*wgt_global)[i] > 0.0))
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
                stdY = std::sqrt(std::max(0.0, sq / double(n) - meanY*meanY));
                if (stdY < EPS) stdY = 1.0;
                for (size_t ii = 0; ii < n; ++ii) Yvec[ii] = (Yvec[ii] - meanY) / stdY;
            }

            Eigen::Map<const Eigen::VectorXd> Yv_map(Yvec.data(), n);

            std::vector<double> baseFitWeights(n, 1.0);
            if (wgt_global)
            {
                baseFitWeights = wgt_adj;
            }
            std::vector<double> finalFitWeights = baseFitWeights;

            if (robust)
            {
                if (!solveHuberIRLS(
                        A,
                        Yv_map,
                        baseFitWeights,
                        huber_tune,
                        tol,
                        max_iter,
                        coeff_vec,
                        finalFitWeights))
                {
                    continue;
                }
            }
            else if (wgt_global)
            {
                Eigen::Map<const Eigen::VectorXd> Wv(
                    wgt_adj.data(),
                    n);
                const Eigen::MatrixXd W = Wv.asDiagonal();
                coeff_vec =
                    (A.transpose() * W * A)
                        .ldlt()
                        .solve(A.transpose() * W * Yv_map);
            }
            else
            {
                coeff_vec = pinv * Yv_map;
            }
            coeff(0) = coeff_vec(0);
            coeff(1) = coeff_vec(1);

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

            // Keep residual/model-selection diagnostics in standardized
            // graphical-analysis space when requested. Reported kinetic
            // parameters and plotted vectors are converted back below.
            const std::vector<double> diagnosticY = Yvec;
            const std::vector<double> diagnosticFitted = fitted;

            // ---------- De-standardize ----------
            double slope = coeff(1), intercept = coeff(0);
            if (standardize)
            {
                const double slopeRaw = slope * stdY / stdX;
                const double interceptRaw =
                    meanY + stdY * intercept - slopeRaw * meanX;
                slope = slopeRaw;
                intercept = interceptRaw;
                for (size_t ii = 0; ii < n; ++ii)
                {
                    fitted[ii] = fitted[ii] * stdY + meanY;
                }
            }

            // ---------- Fill MTGAParameters ----------
            params.Ki = slope;
            params.Intercept = intercept;
            params.x = X_all;
            params.y.resize(n);
            for (size_t ii = 0; ii < n; ++ii)
            {
                params.y[ii] =
                    voxels[v][keepIndex[ii]] *
                    invCp[keepIndex[ii]];
            }
            params.fitted = fitted;
            params.frame.resize(n);
            for (size_t ii = 0; ii < n; ++ii)
            {
                params.frame[ii] = keepIndex[ii] + 1;
            }
            params.dof = 2;
            params.weights =
                robust
                ? finalFitWeights
                : baseFitWeights;
            params.r.resize(n);
            for (size_t ii = 0; ii < n; ++ii)
            {
                params.r[ii] =
                    diagnosticY[ii] - diagnosticFitted[ii];
            }

            const std::vector<double>* diagnosticWeights =
                &baseFitWeights;
            params.AIC = computeAIC(diagnosticY, diagnosticFitted, 2, diagnosticWeights);
            params.MASE = MASE(diagnosticY, diagnosticFitted, diagnosticWeights);
            params.R2 = computeR2(diagnosticY, diagnosticFitted, diagnosticWeights);
            params.chi2 =
                (n > 2)
                ? computeChi2(diagnosticY, diagnosticFitted, diagnosticWeights) /
                    static_cast<double>(n - 2)
                : std::numeric_limits<double>::quiet_NaN();

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

void vtkSlicerDynamicPETLogic::RelativePatlak4Img(
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
    std::function<void(int)> progressCallback,
    size_t dataStartIndex)
{
    if (Cp.empty() || framing.size() != Cp.size() || dataStartIndex >= Cp.size())
    {
        return;
    }

    // Relative Patlak restarts only the plasma integral at the first
    // jointly available late-time frame. The regression itself may begin
    // later according to the user-selected MTGA start.
    this->Patlak4Img(
        voxels, Cp, framing, wgt_global, timeOffset, framingNorm,
        robust, standardize, huber_tune, tol, max_iter, outputParams,
        fitVoxelIndices, stopRequested, numThreads, progressCallback,
        0.0, dataStartIndex);
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

    // PET observations are frame averages and are associated with frame
    // midpoints. Preserve prior full-frame areas, then add half of the
    // current frame to evaluate the cumulative integral at its midpoint.
    std::vector<double> timeAlong(N, 0.0);
    std::vector<double> intCp(N, 0.0);
    double elapsed = 0.0;
    double accumulatedCp = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        timeAlong[i] = elapsed + 0.5 * frameScaled[i];
        intCp[i] = accumulatedCp + 0.5 * Cp[i] * frameScaled[i];
        accumulatedCp += Cp[i] * frameScaled[i];
        elapsed += frameScaled[i];
    }

    // ---------- 2) Determine frames >= timeOffset ----------
    std::vector<int> keepIndex;
    for (size_t i = 0; i < N; ++i)
    {
        if (timeAlong[i] >= timeOffset &&
            (!wgt_global || (*wgt_global)[i] > 0.0))
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
            double accumulatedCt = 0.0;
            for (size_t i = 0; i < N; ++i)
            {
                intCt[i] = accumulatedCt + 0.5 * tac[i] * frameScaled[i];
                accumulatedCt += tac[i] * frameScaled[i];
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
                stdX = std::sqrt(std::max(0.0, sqX/double(n) - meanX*meanX));
                stdY = std::sqrt(std::max(0.0, sqY/double(n) - meanY*meanY));
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

            std::vector<double> baseFitWeights(n, 1.0);
            if (wgt_global)
            {
                baseFitWeights = wgt_adj;
            }
            std::vector<double> finalFitWeights = baseFitWeights;

            if (robust)
            {
                if (!solveHuberIRLS(
                        A,
                        Yv_map,
                        baseFitWeights,
                        huber_tune,
                        tol,
                        max_iter,
                        coeff_vec,
                        finalFitWeights))
                {
                    continue;
                }
            }
            else if (wgt_global)
            {
                Eigen::Map<const Eigen::VectorXd> Wv(
                    wgt_adj.data(),
                    n);
                const Eigen::MatrixXd W = Wv.asDiagonal();
                coeff_vec =
                    (A.transpose() * W * A)
                        .ldlt()
                        .solve(A.transpose() * W * Yv_map);
            }
            else
            {
                coeff_vec = pinv * Yv_map;
            }
            coeff(0) = coeff_vec(0);
            coeff(1) = coeff_vec(1);

            if (stopRequested.load(
                    std::memory_order_relaxed))
            {
              continue;
            }

            // ---------- Compute fitted values ----------
            for (size_t ii=0; ii<n; ++ii) fitted[ii] = coeff(0) + coeff(1)*Xvec[ii];

            // Keep residual/model-selection diagnostics in standardized
            // graphical-analysis space when requested. Reported kinetic
            // parameters and plotted vectors are converted back below.
            const std::vector<double> diagnosticY = Yvec;
            const std::vector<double> diagnosticFitted = fitted;

            // ---------- De-standardize ----------
            double slope = coeff(1), intercept = coeff(0);
            if (standardize)
            {
                const double slopeRaw = slope * stdY / stdX;
                const double interceptRaw =
                    meanY + stdY * intercept - slopeRaw * meanX;
                slope = slopeRaw;
                intercept = interceptRaw;
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
            params.weights = robust ? finalFitWeights : baseFitWeights;
            params.r.resize(n); for (size_t ii=0; ii<n; ++ii) params.r[ii] = diagnosticY[ii] - diagnosticFitted[ii];

            const std::vector<double>* diagnosticWeights = &baseFitWeights;
            params.AIC = computeAIC(diagnosticY, diagnosticFitted, 2, diagnosticWeights);
            params.MASE = MASE(diagnosticY, diagnosticFitted, diagnosticWeights);
            params.R2 = computeR2(diagnosticY, diagnosticFitted, diagnosticWeights);
            params.chi2 = (n > 2)
                ? computeChi2(diagnosticY, diagnosticFitted, diagnosticWeights) / static_cast<double>(n - 2)
                : std::numeric_limits<double>::quiet_NaN();

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
    std::function<void(int)> progressCallback,
    size_t integralStartIndex)
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
    if (framing.size() != N || integralStartIndex >= N) return;
    size_t nVoxels = voxels.size();
    if (nVoxels == 0) return;

    // ---------- 1) Precompute frameScaled, timeAlong, intCp ----------
    std::vector<double> frameScaled(N);
    for (size_t i = 0; i < N; ++i) frameScaled[i] = framing[i] / framingNorm;

    // PET observations are frame averages and are associated with frame
    // midpoints. Preserve prior full-frame areas, then add half of the
    // current frame to evaluate the cumulative integral at its midpoint.
    std::vector<double> timeAlong(N, 0.0);
    std::vector<double> intCp(N, std::numeric_limits<double>::quiet_NaN());
    double elapsed = 0.0;
    for (size_t i = 0; i < integralStartIndex; ++i)
    {
        timeAlong[i] = elapsed + 0.5 * frameScaled[i];
        elapsed += frameScaled[i];
    }
    double accumulatedCp = 0.0;
    for (size_t i = integralStartIndex; i < N; ++i)
    {
        timeAlong[i] = elapsed + 0.5 * frameScaled[i];
        intCp[i] = accumulatedCp + 0.5 * Cp[i] * frameScaled[i];
        accumulatedCp += Cp[i] * frameScaled[i];
        elapsed += frameScaled[i];
    }

    std::vector<double> invCp(N);
    for (size_t i = 0; i < N; ++i) invCp[i] = 1.0 / (Cp[i] + EPS);

    // ---------- 2) Determine frames >= timeOffset ----------
    std::vector<double> X_all;
    std::vector<int> keepIndex;
    for (size_t i = integralStartIndex; i < N; ++i)
    {
      if (timeAlong[i] >= timeOffset &&
          (!wgt_global || (*wgt_global)[i] > 0.0))
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
            std::vector<double> intCt(N, std::numeric_limits<double>::quiet_NaN());
            double accumulatedCt = 0.0;
            for (size_t i = integralStartIndex; i < N; ++i)
            {
                intCt[i] = accumulatedCt + 0.5 * tac[i] * frameScaled[i];
                accumulatedCt += tac[i] * frameScaled[i];
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
                stdY = std::sqrt(std::max(0.0, sq/double(n) - meanY*meanY));
                if (stdY < EPS) stdY = 1.0;
                for (size_t ii=0; ii<n; ++ii) Yvec[ii] = (Yvec[ii]-meanY)/stdY;
            }

            Eigen::Map<const Eigen::VectorXd> Yv_map(Yvec.data(), n);

            std::vector<double> baseFitWeights(n, 1.0);
            if (wgt_global)
            {
                baseFitWeights = wgt_adj;
            }
            std::vector<double> finalFitWeights = baseFitWeights;

            if (robust)
            {
                if (!solveHuberIRLS(
                        A,
                        Yv_map,
                        baseFitWeights,
                        huber_tune,
                        tol,
                        max_iter,
                        coeff_vec,
                        finalFitWeights))
                {
                    continue;
                }
            }
            else if (wgt_global)
            {
                Eigen::Map<const Eigen::VectorXd> Wv(
                    wgt_adj.data(),
                    n);
                const Eigen::MatrixXd W = Wv.asDiagonal();
                coeff_vec =
                    (A.transpose() * W * A)
                        .ldlt()
                        .solve(A.transpose() * W * Yv_map);
            }
            else
            {
                coeff_vec = pinv * Yv_map;
            }
            coeff(0) = coeff_vec(0);
            coeff(1) = coeff_vec(1);

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

            // Keep residual/model-selection diagnostics in standardized
            // graphical-analysis space when requested. Reported kinetic
            // parameters and plotted vectors are converted back below.
            const std::vector<double> diagnosticY = Yvec;
            const std::vector<double> diagnosticFitted = fitted;

            // ---------- De-standardize ----------
            double slope = coeff(1), intercept = coeff(0);
            if (standardize)
            {
                const double slopeRaw = slope * stdY / stdX;
                const double interceptRaw =
                    meanY + stdY * intercept - slopeRaw * meanX;
                slope = slopeRaw;
                intercept = interceptRaw;
                for (size_t ii=0; ii<n; ++ii)
                {
                    fitted[ii] = fitted[ii]*stdY + meanY;
                }
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
            params.weights = robust ? finalFitWeights : baseFitWeights;
            params.r.resize(n); for (size_t ii=0; ii<n; ++ii) params.r[ii] = diagnosticY[ii] - diagnosticFitted[ii];

            const std::vector<double>* diagnosticWeights = &baseFitWeights;
            params.AIC = computeAIC(diagnosticY, diagnosticFitted, 2, diagnosticWeights);
            params.MASE = MASE(diagnosticY, diagnosticFitted, diagnosticWeights);
            params.R2 = computeR2(diagnosticY, diagnosticFitted, diagnosticWeights);
            params.chi2 = (n > 2)
                ? computeChi2(diagnosticY, diagnosticFitted, diagnosticWeights) / static_cast<double>(n - 2)
                : std::numeric_limits<double>::quiet_NaN();

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


void vtkSlicerDynamicPETLogic::RelativeRE4Img(
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
    std::function<void(int)> progressCallback,
    size_t dataStartIndex)
{
    const size_t N = Cp.size();
    if (N == 0 || framing.size() != N || dataStartIndex >= N) return;

    size_t relativeStartIndex = N;
    double frameStart = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        const double frameMid = frameStart + 0.5 * framing[i] / framingNorm;
        if (i >= dataStartIndex && frameMid + 1e-12 >= timeOffset)
        {
            relativeStartIndex = i;
            break;
        }
        frameStart += framing[i] / framingNorm;
    }
    if (relativeStartIndex >= N || N - relativeStartIndex < 2) return;

    // Relative RE is the published RE transform with both cumulative
    // integrals restarted at the selected reversible-equilibrium frame.
    // Reuse the same voxelwise kernel without copying the voxel matrix.
    this->RE4Img(
        voxels, Cp, framing, wgt_global, timeOffset, framingNorm,
        robust, standardize, huber_tune, tol, max_iter, outputParams,
        fitVoxelIndices, stopRequested, numThreads, progressCallback,
        relativeStartIndex);
}

std::vector<double> vtkSlicerDynamicPETLogic::ExtractParameter(
    const std::vector<MTGAParameters>& outputParams,
    const std::string& field)
{
  std::vector<double> values;
  values.reserve(outputParams.size());

  for (const auto& param : outputParams)
  {
    if (field == "Ki" || field == "KiPrime") values.push_back(param.Ki);
    else if (field == "DV" || field == "DVPrime") values.push_back(param.DV);
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
    const std::vector<double>* wgt_global /*= nullptr*/,
    int numThreads /*= 0 */,
    std::function<void(int)> progressCallback,
    std::function<bool()> stopCallback,
    const std::string& interpolationType,
    const std::vector<double>* nativePlasmaTimesSec,
    const std::vector<double>* nativePlasmaValues,
    const std::vector<double>* nativeWholeBloodTimesSec,
    const std::vector<double>* nativeWholeBloodValues,
    const std::vector<double>* parentFractionTimesSec,
    const std::vector<double>* parentFractionValues,
    bool plasmaIsParent,
    double acquisitionStartSec
    )
{
    constexpr double EPS = 1e-16;
    const int Nframe = static_cast<int>(Cp.size());
    const int Nvox = static_cast<int>(voxels.size());
    if (Nvox == 0 ||
        framing.size() != Cp.size() ||
        Cwb.size() != Cp.size())
    {
        return;
    }

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

    // ---------- 3) Scan intervals ----------
    std::vector<double*> scant(Nframe);
    for (int i = 0; i < Nframe; ++i)
    {
        scant[i] = new double[2];
        scant[i][0] = acquisitionStartSec + ((i == 0) ? 0.0 : cumsum[i - 1]);
        scant[i][1] = acquisitionStartSec + cumsum[i];
    }

    // ---------- 4) Fine-sample plasma and whole blood independently ----------
    long int N_cp = 0;
    long int N_wb = 0;

    double* Cp_new = nullptr;
    double* cwb_new = nullptr;

    const bool useNativePlasma =
        nativePlasmaTimesSec &&
        nativePlasmaValues &&
        nativePlasmaTimesSec->size() >= 2 &&
        nativePlasmaTimesSec->size() ==
            nativePlasmaValues->size();

    if (useNativePlasma)
    {
        Cp_new =
            FineSampleExplicitInputFunction(
                *nativePlasmaTimesSec,
                *nativePlasmaValues,
                acquisitionStartSec + cumsum.back(),
                timestep,
                interpolationType,
                N_cp);
    }
    else
    {
        if (interpolationType == "pchip")
        {
            std::vector<double> times, values;
            BuildFrameRepresentativeCurve(scant.data(), Cp, times, values);
            Cp_new = FineSampleExplicitInputFunction(times, values, acquisitionStartSec + cumsum.back(), timestep, "pchip", N_cp);
        }
        else
        {
            Cp_new = finesample2(scant, Cp, N_cp, timestep, interpolationType);
        }
    }

    const bool useNativeWholeBlood =
        nativeWholeBloodTimesSec &&
        nativeWholeBloodValues &&
        nativeWholeBloodTimesSec->size() >= 2 &&
        nativeWholeBloodTimesSec->size() ==
            nativeWholeBloodValues->size();

    if (useNativeWholeBlood)
    {
        cwb_new =
            FineSampleExplicitInputFunction(
                *nativeWholeBloodTimesSec,
                *nativeWholeBloodValues,
                acquisitionStartSec + cumsum.back(),
                timestep,
                interpolationType,
                N_wb);
    }
    else
    {
        if (interpolationType == "pchip")
        {
            std::vector<double> times, values;
            BuildFrameRepresentativeCurve(scant.data(), Cwb, times, values);
            cwb_new = FineSampleExplicitInputFunction(times, values, acquisitionStartSec + cumsum.back(), timestep, "pchip", N_wb);
        }
        else
        {
            cwb_new = finesample2(scant, Cwb, N_wb, timestep, interpolationType);
        }
    }

    if (N_cp != N_wb)
    {
        delete[] Cp_new;
        delete[] cwb_new;

        throw std::runtime_error(
            "Plasma and whole-blood fine-sampling grids differ.");
    }

    const bool applyParentFraction =
        !plasmaIsParent &&
        parentFractionTimesSec &&
        parentFractionValues &&
        parentFractionTimesSec->size() >= 2 &&
        parentFractionTimesSec->size() ==
            parentFractionValues->size();

    if (applyParentFraction)
    {
        long int N_parent = 0;

        double* parentFractionNew =
            FineSampleExplicitInputFunction(
                *parentFractionTimesSec,
                *parentFractionValues,
                acquisitionStartSec + cumsum.back(),
                timestep,
                "linear",
                N_parent);

        if (N_parent != N_cp)
        {
            delete[] parentFractionNew;
            delete[] Cp_new;
            delete[] cwb_new;

            throw std::runtime_error(
                "Parent-fraction and plasma fine-sampling grids differ.");
        }

        for (long int i = 0;
             i < N_cp;
             ++i)
        {
            Cp_new[i] *= parentFractionNew[i];
        }

        delete[] parentFractionNew;
    }

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
        params.boundFlags = TCM_BOUND_NONE;
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
    std::atomic<int> boundHitVoxelCount(0);

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

            params.boundFlags = TCM_BOUND_NONE;
            markBoundIfNeeded(params.boundFlags, params.vb, lb[0], ub[0], sens[0], TCM_BOUND_VB_LOWER, TCM_BOUND_VB_UPPER);
            markBoundIfNeeded(params.boundFlags, params.K1, lb[1], ub[1], sens[1], TCM_BOUND_K1_LOWER, TCM_BOUND_K1_UPPER);
            markBoundIfNeeded(params.boundFlags, params.k2, lb[2], ub[2], sens[2], TCM_BOUND_K2_LOWER, TCM_BOUND_K2_UPPER);
            if (n_tc == 2)
            {
                markBoundIfNeeded(params.boundFlags, params.k3, lb[3], ub[3], sens[3], TCM_BOUND_K3_LOWER, TCM_BOUND_K3_UPPER);
                markBoundIfNeeded(params.boundFlags, params.k4, lb[4], ub[4], sens[4], TCM_BOUND_K4_LOWER, TCM_BOUND_K4_UPPER);
                markBoundIfNeeded(params.boundFlags, params.td, lb[5], ub[5], sens[5], TCM_BOUND_TD_LOWER, TCM_BOUND_TD_UPPER);
            }
            else
            {
                markBoundIfNeeded(params.boundFlags, params.td, lb[3], ub[3], sens[3], TCM_BOUND_TD_LOWER, TCM_BOUND_TD_UPPER);
            }

            if (params.boundFlags != TCM_BOUND_NONE)
            {
                ++boundHitVoxelCount;
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

    if (boundHitVoxelCount.load() > 0)
    {
        std::cout
            << "TCM fitting: "
            << boundHitVoxelCount.load()
            << " fitted voxels reached or approached at least one active parameter bound."
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
