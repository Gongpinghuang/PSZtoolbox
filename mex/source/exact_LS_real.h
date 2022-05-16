// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include <mex.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#ifndef FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
    #define FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
#endif
#include <blas.h>
#include <lapack.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#ifndef EXACT_LS_REAL
#define EXACT_LS_REAL

// -------------------------------------------------------------------------
void exact_LS_real(double *h,
                   double *d,
                   const int L,
                   const int M,
                   const int Ih,
                   const int Ig,
                   const int Id,
                   const double beta,
                   double *g);
// -------------------------------------------------------------------------
void matC_x_herm_times_y_step(fftw_complex *x,
                              fftw_complex *y,
                              fftw_complex *z,
                              const ptrdiff_t M,
                              const ptrdiff_t stepx,
                              const ptrdiff_t stepy);
// -------------------------------------------------------------------------
void matR_chol(double *X,const ptrdiff_t N);
// -------------------------------------------------------------------------
void matR_solve_sys(double *L,const ptrdiff_t N,double *b);
// -------------------------------------------------------------------------
void matR_scale(double *X,const ptrdiff_t N,const double alpha);
// -------------------------------------------------------------------------


#endif