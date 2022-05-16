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

#ifndef EXACT_LS_COMPLEX
#define EXACT_LS_COMPLEX

// -------------------------------------------------------------------------
void exact_LS_complex(fftw_complex *h,
                      fftw_complex *d,
                      const int L,
                      const int M,
                      const int Ih,
                      const int Ig,
                      const int Id,
                      const double beta,
                      fftw_complex *g);
// -------------------------------------------------------------------------
void matC_x_herm_times_y_step(fftw_complex *x,
                              fftw_complex *y,
                              fftw_complex *z,
                              const ptrdiff_t M,
                              const ptrdiff_t stepx,
                              const ptrdiff_t stepy);
//.............................................................................
void matC_chol(fftw_complex *X,const ptrdiff_t N);
//.............................................................................
void matC_solve_sys(fftw_complex *U,
                    const ptrdiff_t N,
                    fftw_complex *b);

#endif