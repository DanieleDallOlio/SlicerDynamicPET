#include <utils.hpp>
#include <cmd.h>
#include <kmodels.h>

int main (int argc, char *argv[]) {

  std :: string tac_file,
                blood_file,
                framing_file,
                delimiter;
  int nth;

  parse_demo_tcm_realdata(argc,
                          argv,
                          tac_file,
                          blood_file,
                          framing_file,
                          nth,
                          delimiter
                         );


  std :: vector< std :: vector<double> > tac, Cp, framing;
  long int Nframe;
  long int Nvox;
  read_matrix(tac, Nframe, Nvox, tac_file, delimiter);
  // print_matrix(tac, Nframe, Nvox);

  read_matrix(Cp, Nframe, Nvox, blood_file, delimiter);
  // print_matrix(Cp, Nframe, Nvox);

  read_matrix(framing, Nframe, Nvox, framing_file, delimiter);
  // print_matrix(dt, Nframe, Nvox);

  if (containsNaN(tac)) {
    error_nan("TAC");
  }
  if (containsNaN(Cp)) {
    error_nan("Cp");
  }
  if (containsNaN(framing)) {
    error_nan("framing");
  }

  if (tac.size() != framing.size()) {
    error_size("TAC", "framing", tac.size(), framing.size());
  }

  if (tac.size() != Cp.size()) {
    error_size("TAC", "Cp", tac.size(), Cp.size());
  }


  double dk = 0.;
  double timestep = 1.;
  double pbrp[] = {1., 0., 0.};
  int maxiter = 100;
  double lb[] = {0., 0., 0., 0., 0., 0.};
  double ub[] = {1., 2., 2., 2., 0.5, 20.};
  bool sens[] = {true, true, true, true, true, true};
  double kinit[] = {0.1, 0.1, 0.1, 0.1, 0.01, 0.};
  double *wt = new double[Nframe];
  for (long int i = 0L; i < Nframe; ++i) wt[i] = 1.;

  double *cumsum = new double[Nframe];
  double cumsum_iter = 0.;
  for (long int i = 0L; i < Nframe; ++i) {
      cumsum_iter += framing[i][0L];
      cumsum[i] = cumsum_iter;
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

  double *fitted_params = new double[6*Nvox];
  double *fitted_curve   = new double[Nframe*Nvox];

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
  std :: cout << "fitted_params" << std ::endl;
  print_vec(fitted_params, 6*Nvox);
  std :: cout << "fitted_curve" << std ::endl;
  print_vec(fitted_curve, Nframe*Nvox);

  return 0;
}
