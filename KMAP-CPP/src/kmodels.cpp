#include <kmodels.h>

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This file implements the fitting of a one-tissue kinetic model (1TCM) using
// the Levenberg-Marquardt algorithm with OpenMP for parallel processing
// within the MATLAB environment.
//
// Usage:
// kfit_1tcm_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
//
// Input parameters:
// - tac: Time activity curve (TAC) data.
// - w: Weights for the TAC data.
// - scant: Scan time data.
// - blood: Blood data.
// - wblood: Whole blood data.
// - dk: Decay constant.
// - pinit: Initial parameters for the model.
// - lb: Lower bounds for the parameters.
// - ub: Upper bounds for the parameters.
// - psens: Sensitivity matrix for the parameters.
// - maxit: Maximum number of iterations for the fitting algorithm.
// - td: Time duration for the scan.
// - nth: number of threads (in case of parallelization)
//
// Output:
// - p: Estimated parameters.
// - c: Fitted curve.
//
//
// This will produce a MEX file named 'kfit_1tcm_mex_omp', which you can call from MATLAB
// as kfit_1tcm_mex_omp(...) with the same arguments as described above.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void kfit_1tcm_mex_omp(double * & tac,
                       long int & num_frm,
                       long int & num_vox,
                       long int & nw,
                       double * & w1,
                       double * & scant,
                       double * & cp,
                       double * & wb,
                       const double dk,
                       const double * pinit,
                       const int & num_par,
                       const int & np,
                       double * plb,
                       double * pub,
                       const bool * temp,
                       const int maxit,
                       const double & td,
                       const int &  nth,
                       double * & p,
                       double * & c
                      )
{

    double *w;
    int i, j;
    int psens[4];
    KMODEL_T km;

    // Set up the kinetic model parameters
    km.dk = dk;
    km.td = td;
    km.cp = cp;
    km.wb = wb;
    km.num_frm = num_frm;
    km.num_vox = 1; // Process one voxel at a time
    km.scant = scant;
    km.tacfunc = kconv_1tcm_tac; // TAC function for 1TCM model
    km.jacfunc = kconv_1tcm_jac; // Jacobian function for 1TCM model

    for (i = 0; i < num_par; i++) {
        psens[i] = (int)temp[i];
    }

    // Initialize parameter p
    if (np == 1) {
        for (i = 0; i < num_par; i++) {
            for (j = 0; j < num_vox; j++) {
                p[i + j * num_par] = pinit[i];
            }
        }
    } else if (np == num_vox) {
        for (i = 0; i < num_par; i++) {
            for (j = 0; j < num_vox; j++) {
                p[i + j * num_par] = pinit[i + j * num_par];
            }
        }
    }

    // Prepare weights
    if (nw == 1) {
        w = (double*) malloc(sizeof(double) * num_frm * num_vox);
        for (i = 0; i < num_frm; i++) {
            for (j = 0; j < num_vox; j++) {
                w[i + j * num_frm] = w1[i];
            }
        }
    } else {
        w = w1;
    }

    i = nth;


    // #ifdef _OPENMP
    // #pragma omp parallel
    {
    // #endif
        int m;
        double *cj, *wj, *pj, *cfit;

        // Allocate memory for each thread
        cj = (double*) malloc(sizeof(double) * num_frm);
        wj = (double*) malloc(sizeof(double) * num_frm);
        pj = (double*) malloc(sizeof(double) * num_par);
        cfit = (double*) malloc(sizeof(double) * num_frm);

        // #ifdef _OPENMP
        // #pragma omp for nowait
        // #endif
        // Voxel-wise fitting
        for (j = 0; j < num_vox; j++) {

            // Copy TAC and weights for the current voxel
            for (m = 0; m < num_frm; m++) {
                cj[m] = tac[m + j * num_frm];
                wj[m] = w[m + j * num_frm];
            }

            // Copy initial parameters for the current voxel
            for (m = 0; m < num_par; m++) {
                pj[m] = p[m + j * num_par];
            }

            // Perform Levenberg-Marquardt fitting
            kmap_levmar(cj, wj, num_frm, pj, num_par, &km, tac_eval, jac_eval, plb, pub,
                        psens, maxit, cfit);

            // Copy fitted curve and parameters back to the output
            for (m = 0; m < num_frm; m++) {
                c[m + j * num_frm] = cfit[m];
            }
            for (m = 0; m < num_par; m++) {
                p[m + j * num_par] = pj[m];
            }
        }

        // Free memory allocated for each thread
        if (cj) free(cj);
        if (wj) free(wj);
        if (pj) free(pj);
        if (cfit) free(cfit);
    // #ifdef _OPENMP
    }
    // #endif

    // Free the weight array if it was dynamically allocated
    if ((w) && (nw == 1)) {
        free(w);
    }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This file implements the fitting of a two-tissue kinetic model (2TCM) using
// the Levenberg-Marquardt algorithm with OpenMP for parallel processing
// within the MATLAB environment.
//
// Usage:
// kfit_2tcm_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
//
//
// This will produce a MEX file named 'kfit_2tcm_mex_omp', which you can call from MATLAB
// as kfit_2tcm_mex_omp(...) with the same arguments as described above.

// Input parameters:
// - tac: Time activity curve (TAC) data.
// - w: Weights for the TAC data.
// - scant: Scan time data.
// - blood: Blood data.
// - wblood: Whole blood data.
// - dk: Decay constant.
// - pinit: Initial parameters for the model.
// - lb: Lower bounds for the parameters.
// - ub: Upper bounds for the parameters.
// - psens: Sensitivity matrix for the parameters.
// - maxit: Maximum number of iterations for the fitting algorithm.
// - td: Time duration for the scan.
//
// Output:
// - p: Estimated parameters.
// - c: Fitted curve.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void kfit_2tcm_mex_omp(double * & tac,
                       long int & num_frm,
                       long int & num_vox,
                       long int & nw,
                       double * & w1,
                       double * & scant,
                       double * & cp,
                       double * & wb,
                       const double dk,
                       const double * pinit,
                       const int & num_par,
                       const int & np,
                       double * plb,
                       double * pub,
                       const bool * temp,
                       const int maxit,
                       const double & td,
                       const int &  nth,
                       double * & p,
                       double * & c
                      )
{

    int i, j;
    double *w;
    int psens[6];
    KMODEL_T km;

    // Set up the kinetic model parameters
    km.dk = dk;
    km.td = td;
    km.cp = cp;
    km.wb = wb;
    km.num_frm = num_frm;
    km.num_vox = 1; // Process one voxel at a time
    km.scant = scant;
    km.tacfunc = kconv_2tcm_tac; // TAC function for 2TCM model
    km.jacfunc = kconv_2tcm_jac; // Jacobian function for 2TCM model

    for (i = 0; i < num_par; i++) {
        psens[i] = (int)temp[i];
    }

    // Initialize parameter p
    if (np == 1) {
        for (i = 0; i < num_par; i++) {
            for (j = 0; j < num_vox; j++) {
                p[i + j * num_par] = pinit[i];
            }
        }
    } else if (np == num_vox) {
        for (i = 0; i < num_par; i++) {
            for (j = 0; j < num_vox; j++) {
                p[i + j * num_par] = pinit[i + j * num_par];
            }
        }
    }

    // Prepare weights
    if (nw == 1) {
        w = (double*) malloc(sizeof(double) * num_frm * num_vox);
        for (i = 0; i < num_frm; i++) {
            for (j = 0; j < num_vox; j++) {
                w[i + j * num_frm] = w1[i];
            }
        }
    } else {
        w = w1;
    }

    i = nth;

    // #ifdef _OPENMP
    // #pragma omp parallel
    {
    // #endif
        int m;
        double *cj, *wj, *pj, *cfit;

        // Allocate memory for each thread
        cj = (double*) malloc(sizeof(double) * num_frm);
        wj = (double*) malloc(sizeof(double) * num_frm);
        pj = (double*) malloc(sizeof(double) * num_par);
        cfit = (double*) malloc(sizeof(double) * num_frm);

        // #ifdef _OPENMP
        // #pragma omp for nowait
        // #endif
        // Voxel-wise fitting
        for (j = 0; j < num_vox; j++) {

            // Copy TAC and weights for the current voxel
            for (m = 0; m < num_frm; m++) {
                cj[m] = tac[m + j * num_frm];
                wj[m] = w[m + j * num_frm];
            }

            // Copy initial parameters for the current voxel
            for (m = 0; m < num_par; m++) {
                pj[m] = p[m + j * num_par];
            }

            // Perform Levenberg-Marquardt fitting
            kmap_levmar(cj, wj, num_frm, pj, num_par, &km, tac_eval, jac_eval, plb, pub,
                        psens, maxit, cfit);

            // Copy fitted curve and parameters back to the output
            for (m = 0; m < num_frm; m++) {
                c[m + j * num_frm] = cfit[m];
            }
            for (m = 0; m < num_par; m++) {
                p[m + j * num_par] = pj[m];
            }
        }

        // Free memory allocated for each thread
        if (cj) free(cj);
        if (wj) free(wj);
        if (pj) free(pj);
        if (cfit) free(cfit);
    // #ifdef _OPENMP
    }
    // #endif

    // Free the weight array if it was dynamically allocated
    if ((w) && (nw == 1)) {
        free(w);
    }
}
