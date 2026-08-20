#include <algorithm>
#include <utils.h>

bool file_exists (const std :: string & filename)
{
//#if have_filesystem == 1
//
//  return std :: filesystem :: exists(std :: filesystem :: path(filename));
//
//#else

  if (FILE *file = fopen(filename.c_str(), "r"))
  {
    fclose(file);
    return true;
  }
  return false;

//#endif
}

std :: vector < std :: string > split (const std :: string & txt, const std :: string & del)
{
  std :: vector < std :: string > token;

  std :: size_t pos = txt.find_first_of(del);
  std :: size_t start = 0;
  std :: size_t end = txt.size();

  while (pos != std :: string :: npos)
  {
    if (pos) token.push_back(txt.substr(start, pos));
    start += pos + 1;
    pos = txt.substr(start, end).find_first_of(del);
  }

  if (start != end) token.push_back(txt.substr(start, pos));

  return token;
}

bool containsNaN(const std::vector<std::vector<double>>& matrix) {
    return std::any_of(matrix.begin(), matrix.end(), [](const std::vector<double>& row) {
        return std::any_of(row.begin(), row.end(), [](double value) {
            return std::isnan(value);
        });
    });
}

double* finesample(double** & scant,
                   std :: vector < std :: vector <double> > & blood,
                   long int N,
                   long int & Nint,
                   long int res,
                   const std::string & itype)
{
    // Create the time vectors
    double* ts = new double[N];
    double* te = new double[N];
    double* tm = new double[N];
    for (long int i = 0L; i < N; ++i) {
        ts[i] = scant[i][0L] / res;
        te[i] = scant[i][1L] / res;
        tm[i] = (ts[i] + te[i]) * 0.5;
    }

    // Create a new time vector with step size 1 (resolution = 1 second)
    std::vector<double> tt;
    for (double t = 0.5; t <= te[N-1]; t += 1.0) {
        tt.push_back(t);
    }
    Nint = tt.size();

    if (ts[0L]<1.) {
      ts[0L] = 1.;
    }

    // Create the Ca vector (output)
    double * Ca = new double[Nint];
    double t, t1, t2, y1, y2;
    int idx1, idx2;
    if (itype == "linear") {
        // Linear interpolation
        idx1 = 0;
        idx2 = 1;
        for (int i = 0; i < Nint; ++i) {
            // Interpolate using linear interpolation
            t = tt[i];

            while (t > tm[idx2] && idx2 < N-1) {
              idx1 += 1;
              idx2 += 1;
            }

            // for (int j = 0; j < N; ++j) {
            //     if (tm[j] <= t )
            //
            //     if (tm[j] <= t && (tm[j] > tm[idx1] || j == 0)) {
            //         idx1 = j;
            //     }
            //     if (tm[j] >= t && (tm[j] < tm[idx2] || j == 0)) {
            //         idx2 = j;
            //     }
            // }
            // Perform linear interpolation between idx1 and idx2
            t1 = tm[idx1];
            t2 = tm[idx2];
            y1 = blood[idx1][0L];
            y2 = blood[idx2][0L];
            Ca[i] = std::max(0.0, y1 + (y2 - y1) * (t - t1) / (t2 - t1));
        }
    }
    else if (itype == "const") {
        // Constant values between each scan
        for (int i = 0; i < N; ++i) {
            for (t = ts[i]; t <= te[i]; ++t) {
                int idx = std::round(t) - 1;
                if (idx < Nint) {
                    Ca[idx] = blood[i][0L];
                }
            }
        }
    } else {
        // Handle other interpolation types (similar to "linear")
        error_interpolation(itype);
    }

    return Ca;
}

double* finesample2(
    const std::vector<double*>& scant,  // [Nframe][2] start/end times
    const std::vector<double>& blood,               // [Nframe]
    long int& Nint,
    long int res,
    const std::string& itype)
{
    const long int N = static_cast<long int>(scant.size());

    // Frame start, end, mid
    std::vector<double> ts(N), te(N), tm(N);
    for (long int i = 0; i < N; ++i) {
        ts[i] = scant[i][0] / res;
        te[i] = scant[i][1] / res;
        tm[i] = 0.5 * (ts[i] + te[i]);
    }

    // Fine time vector (step = 1 second)
    std::vector<double> tt;
    for (double t = 0.5; t <= te[N-1]; t += 1.0) {
        tt.push_back(t);
    }
    Nint = static_cast<long int>(tt.size());

    if (ts[0] < 1.0) {
        ts[0] = 1.0;
    }

    // Output vector
    double* Ca = new double[Nint];

    if (itype == "linear") {
        // Linear interpolation
        long int idx1 = 0, idx2 = 1;
        for (long int i = 0; i < Nint; ++i) {
            double t = tt[i];
            while (t > tm[idx2] && idx2 < N-1) {
                idx1++;
                idx2++;
            }
            double t1 = tm[idx1];
            double t2 = tm[idx2];
            double y1 = blood[idx1];
            double y2 = blood[idx2];
            Ca[i] = std::max(0.0, y1 + (y2 - y1) * (t - t1) / (t2 - t1));
        }
    }
    else if (itype == "const") {
        // Step-wise constant values
        std::fill(Ca, Ca+Nint, 0.0);
        for (long int i = 0; i < N; ++i) {
            for (double t = ts[i]; t <= te[i]; ++t) {
                int idx = static_cast<int>(std::round(t)) - 1;
                if (idx >= 0 && idx < Nint) {
                    Ca[idx] = blood[i];
                }
            }
        }
    }
    else {
        error_interpolation(itype);
    }

    return Ca;
}
