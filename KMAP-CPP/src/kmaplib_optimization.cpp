/*
 * kmaplib_optimization.cpp
 *
 * This file contains the optimization functions used for parameter estimation in
 * kinetic modeling. The Levenberg-Marquardt algorithm and associated helper
 * functions are implemented here.
 */

#include "kmaplib.h"
// #include "mex.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

//------------------------------------------------------------------------------
// kmap_levmar
//------------------------------------------------------------------------------
// The Levenberg-Marquardt algorithm to solve the nonlinear least squares problem
// subject to linear bounds.
// This function is used to estimate the parameters that best fit the model to the data.
void kmap_levmar(double *y, double *w, int num_frm, double *pinit, int num_p,
                 void *param, void (*func)(double *, void *, double *),
                 void (*jacf)(double *, void *, double *, int *, double *),
                 double *plb, double *pub, int *psens, int maxit, double *ct)
{
   int      i, j, it, n;
   double   mu, v, tau, rho, maxaii;
   int      num_par;
   double   *f, *fnew, *r;
   double   *h;
   double   *pnew, *pest, *plb_sens, *pub_sens;
   double   *st;
   double   etol = 1e-9;
   double   F, Fnew, Lnew;
   int      subit = 100;
   double   tmp;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (j = 0; j < num_p; j++)
      if (psens[j] == 1) ++num_par;

   // Memory allocation
   f = (double*) malloc(sizeof(double) * num_frm);
   fnew = (double*) malloc(sizeof(double) * num_frm);
   r = (double*) malloc(sizeof(double) * num_frm);
   h = (double*) malloc(sizeof(double) * num_par);
   pnew = (double*) malloc(sizeof(double) * num_par);
   pest = (double*) malloc(sizeof(double) * num_par);
   plb_sens = (double*) malloc(sizeof(double) * num_par);
   pub_sens = (double*) malloc(sizeof(double) * num_par);
   st = (double*) malloc(sizeof(double) * num_frm * num_par);

   // Initialize the parameter estimates and their bounds
   getkin(pinit, psens, num_par, pest);
   getkin(plb, psens, num_par, plb_sens);
   getkin(pub, psens, num_par, pub_sens);

   // Initial TAC and sensitivity matrix
   (*jacf)(pinit, param, ct, psens, st);

   // Compute the initial residual and the cost value
   for (i = 0; i < num_frm; i++) {
      f[i] = y[i] - ct[i];
   }
   F = vecnormw(w, f, num_frm) * 0.5;

   // Initialize parameters for the search process
   v = 2.0;
   tau = 1.0e-3;
   maxaii = 0;
   for (j = 0; j < num_par; j++) {
      tmp = vecnorm2(st + num_frm * j, num_frm);
      tmp *= tmp;
      if (maxaii < tmp)
         maxaii = tmp;
   }
   mu = tau * maxaii;

   // Iterative optimization loop
   it = 0; n = 0;
   while (it < maxit) {
      // Estimate new parameters
      for (j = 0; j < num_par; j++)
         pnew[j] = pest[j];

      // Coordinate least squares step
      boundpls_cd(f, w, st, mu, num_frm, num_par, plb_sens, pub_sens, subit, pnew, r);

      // Check for convergence or update mu
      for (i = 0; i < num_par; i++) {
         h[i] = pnew[i] - pest[i];
      }
      if (vecnorm2(h, num_par) <= etol * (vecnorm2(pest, num_par) + etol))
         break;

      // Calculate the ratio rho
      setkin(pnew, num_par, psens, pinit);
      (*func)(pinit, param, ct);
      for (i = 0; i < num_frm; i++)
         fnew[i] = y[i] - ct[i];
      Fnew = vecnormw(w, fnew, num_frm) * 0.5;
      Lnew = vecnormw(w, r, num_frm) * 0.5;
      rho  = (F - Fnew) / (F - Lnew);

      // Accept the step or find a new one
      if (rho > 0) {
         // Update estimates
         for (j = 0; j < num_par; j++) {
            pest[j] = pnew[j];
         }
         (*jacf)(pinit, param, ct, psens, st);
         for (i = 0; i < num_frm; i++) {
            f[i] = fnew[i];
         }
         F = Fnew;
         ++it;

         // Update mu
         tmp = 2 * rho - 1;
         tmp = 1 - tmp * tmp * tmp;
         if (tmp < 1.0 / 3.0)
            tmp = 1.0 / 3.0;
         mu *= tmp;
         v = 2.0;
      } else {
         mu *= v;
         v *= 2.0;
      }
      ++n;
      if (n > 1000) {
         printf("stopped: maximum iteration number exceeds 1000!\n");
         setkin(pest, num_par, psens, pinit);
         (*func)(pinit, param, ct);
         break;
      }
   }

   // Free allocated memory
   free(f);
   free(fnew);
   free(r);
   free(h);
   free(pnew);
   free(pest);
   free(plb_sens);
   free(pub_sens);
   free(st);
}

//------------------------------------------------------------------------------
// boundpls_cd
//------------------------------------------------------------------------------
// Coordinate descent algorithm to solve the penalized least squares problem
// subject to box bounds. Minimizes || y - A * (x - x_n) ||^2.
void boundpls_cd(double *y, double *w, double *a, double alpha, int num_y,
                 int num_x, double *xlb, double *xub, int maxit, double *x,
                 double *r)
{
   int      i, j, ij;
   double   err0, err, dj, xj;
   double   tmp1, tmp2;
   double   etol = 1.0e-9;
   int      it;

   // Initialize the residual
   for (i = 0; i < num_y; i++){
      r[i] = y[i];
   }
   err0 = vecnorm2(y, num_y);

   // std :: cout << "alpha = " << alpha << std :: endl;
   // std :: cout << "num_y = " << num_y << std :: endl;
   // std :: cout << "num_x = " << num_x << std :: endl;
   // std :: cout << "maxit = " << maxit << std :: endl;

   // for (j = 0; j < num_y; j++) {
   //   std :: cout << j << " - y[j] = " << y[j] << std :: endl;
   //   std :: cout << j << " - r[j] = " << r[j] << std :: endl;
   // }
   //
   // for (j = 0; j < num_x; j++) {
   //   std :: cout << j << " - w[j] = " << w[j] << std :: endl;
   //   std :: cout << j << " - a[j] = " << a[j] << std :: endl;
   //   std :: cout << j << " - xlb[j] = " << xlb[j] << std :: endl;
   //   std :: cout << j << " - xub[j] = " << xub[j] << std :: endl;
   // }
   //
   // for (j = 0; j < num_x; j++) {
   //   std :: cout << "pre - " << j << " - x[j] = " << x[j] << std :: endl;
   // }

   // Iterative optimization
   for (it = 0; it < maxit; it++) {
      // Check for convergence
      err = vecnorm2(r, num_y);
      // std :: cout << "it = " << it << ", err = " << err << std :: endl;
      if (err / err0 < etol) {
        // std :: cout << "boundpls_cd break" << std :: endl;
        break;
      }

      // Coordinate-wise optimization
      for (j = 0; j < num_x; j++) {
         tmp1 = tmp2 = 0;
         for (i = 0; i < num_y; i++) {
            ij = i + j * num_y;
            tmp1 += w[i] * a[ij] * r[i];
            tmp2 += w[i] * a[ij] * a[ij];
         }
         tmp2 += alpha * tmp2;
         if (tmp2 == 0)
            dj = 0;
         else
            dj = tmp1 / tmp2;
         xj = x[j] + dj;

         // Apply bounds to the estimate
         if (xj < xlb[j])
            xj = xlb[j];
         else if (xj > xub[j])
            xj = xub[j];

         // Update x and r
         dj = xj - x[j];
         x[j] = xj;
         for (i = 0; i < num_y; i++)
            r[i] -= a[i + j * num_y] * dj;
      }
   }
   // for (j = 0; j < num_x; j++) {
   //   std :: cout << "post - " << j << " - x[j] = " << x[j] << std :: endl;
   // }
}

void kmap_levmar_stats(
    double *y, double *w, int num_frm, double *pinit, int num_p,
    void *param, void (*func)(double *, void *, double *),
    void (*jacf)(double *, void *, double *, int *, double *),
    double *plb, double *pub, int *psens, int maxit, double *ct,
    LevmarStats *stx
)
{
   int      i, j, it, n;
   double   mu, v, tau, rho, maxaii;
   int      num_par;
   double   *f, *fnew, *r;
   double   *h;
   double   *pnew, *pest, *plb_sens, *pub_sens;
   double   *st;
   double   etol = 1e-9;
   double   F, Fnew, Lnew;
   int      subit = 100;
   double   tmp;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (j = 0; j < num_p; j++)
      if (psens[j] == 1) ++num_par;

   // Allocate using your existing malloc/free style (minimal changes)
   f        = (double*) malloc(sizeof(double) * num_frm);
   fnew     = (double*) malloc(sizeof(double) * num_frm);
   r        = (double*) malloc(sizeof(double) * num_frm);
   h        = (double*) malloc(sizeof(double) * num_par);
   pnew     = (double*) malloc(sizeof(double) * num_par);
   pest     = (double*) malloc(sizeof(double) * num_par);
   plb_sens = (double*) malloc(sizeof(double) * num_par);
   pub_sens = (double*) malloc(sizeof(double) * num_par);
   st       = (double*) malloc(sizeof(double) * num_frm * num_par);

   // Initialize
   getkin(pinit, psens, num_par, pest);
   getkin(plb, psens, num_par, plb_sens);
   getkin(pub, psens, num_par, pub_sens);

   // Initial jacobian and ct
   // for (i = 0; i < num_frm; i++) {
   //   std :: cout <<"Before initial jacobian " <<  i << " - ct[i] = " << ct[i] << std :: endl;
   // }
   (*jacf)(pinit, param, ct, psens, st);
   if (stx && kmap_isnan(st[0])) stx->nan_st++;

   // std :: cout << "num_frm = " << num_frm << std :: endl;
   for (i = 0; i < num_frm; i++) {
      f[i] = y[i] - ct[i];
      // std :: cout << i << " - f[i] = " << f[i] << std :: endl;
      // std :: cout << i << " - y[i] = " << y[i] << std :: endl;
      // std :: cout << i << " - ct[i] = " << ct[i] << std :: endl;
      // std :: cout << i << " - w[i] = " << w[i] << std :: endl;
   }
   F = vecnormw(w, f, num_frm) * 0.5;
   // std :: cout << "F = " << F << std :: endl;

   v = 2.0;
   tau = 1.0e-3;
   maxaii = 0;
   for (j = 0; j < num_par; j++) {
      tmp = vecnorm2(st + num_frm * j, num_frm);
      // std :: cout << j << " - tmp = " << tmp << std :: endl;
      tmp *= tmp;
      if (maxaii < tmp)
         maxaii = tmp;
   }
   // std :: cout << "maxaii = " << maxaii << std :: endl;
   mu = tau * maxaii;
   // std :: cout << "mu = " << mu << std :: endl;

   it = 0; n = 0;
   // std :: cout << "etol = " << etol << std :: endl;
   // std :: cout << "num_par = " << num_par << std :: endl;
   while (it < maxit) {
      if (stx) stx->n_loops++; // TODO: count loops; do I get the same numbers?

      // Copy current estimate
      for (j = 0; j < num_par; j++)
         pnew[j] = pest[j];

      // Coordinate descent step
      boundpls_cd(f, w, st, mu, num_frm, num_par, plb_sens, pub_sens, subit, pnew, r);

      // h = pnew - pest
      for (i = 0; i < num_par; i++) {
         h[i] = pnew[i] - pest[i];
      }
      // std :: cout << "vecnorm2(h, num_par) = " << vecnorm2(h, num_par) << std :: endl;
      // std :: cout << "vecnorm2(pest, num_par) = " << vecnorm2(pest, num_par) << std :: endl;
      if (vecnorm2(h, num_par) <= etol * (vecnorm2(pest, num_par) + etol)) {
        // std :: cout << "first break" << std :: endl;
        break;
      }

      // Evaluate new parameters
      setkin(pnew, num_par, psens, pinit);
      (*func)(pinit, param, ct);
      if (stx && kmap_isnan(ct[0])) stx->nan_ct++;

      for (i = 0; i < num_frm; i++)
         fnew[i] = y[i] - ct[i];

      Fnew = vecnormw(w, fnew, num_frm) * 0.5;
      Lnew = vecnormw(w, r, num_frm) * 0.5;
      rho  = (F - Fnew) / (F - Lnew);

      if (stx) { stx->last_rho = rho; stx->last_mu = mu; }

      if (rho > 0) {
         if (stx) stx->it_accept++;

         for (j = 0; j < num_par; j++) {
            pest[j] = pnew[j];
         }

         (*jacf)(pinit, param, ct, psens, st);
         if (stx && kmap_isnan(st[0])) stx->nan_st++;

         for (i = 0; i < num_frm; i++) {
            f[i] = fnew[i];
         }
         F = Fnew;
         ++it;

         tmp = 2 * rho - 1;
         tmp = 1 - tmp * tmp * tmp;
         if (tmp < 1.0 / 3.0)
            tmp = 1.0 / 3.0;
         mu *= tmp;
         v = 2.0;
      } else {
         if (stx) stx->it_reject++;
         mu *= v;
         v *= 2.0;
      }

      ++n;
      if (n > 1000) {
         if (stx) stx->hit_guard++;

         setkin(pest, num_par, psens, pinit);
         (*func)(pinit, param, ct);
         if (stx && kmap_isnan(ct[0])) stx->nan_ct++;
         std :: cout << "second break" << std :: endl;
         break;
      }
   }

   free(f);
   free(fnew);
   free(r);
   free(h);
   free(pnew);
   free(pest);
   free(plb_sens);
   free(pub_sens);
   free(st);
}
