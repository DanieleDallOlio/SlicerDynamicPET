#ifndef __kmodels_h__
#define __kmodels_h__

#include "kmaplib.h"


/*
 * kmodels.h
 *
 * This file contains the function declarations for each single compartment model.
 */

 /*
  * kfit_1tcm_mex_omp
  *
  * Fitting of a one-tissue kinetic model (1TCM)
*/
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
                      );

/*
  * kfit_2tcm_mex_omp
  *
  * Fitting of a two-tissue kinetic model (1TCM)
*/
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
                      );

#endif // __kmodels_h__
