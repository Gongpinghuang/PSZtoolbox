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

#ifndef MPI
#define MPI 3.14159265358979323846
#endif

#ifndef GDFT_FB_SC_DEC
#define GDFT_FB_SC_DEC

// -------------------------------------------------------------------------
void gdft_fb_sb_dec(fftw_complex *x_an,
                    int Ix_an,
                    int Nx,
                    double *hp,
                    const int K,
                    const int R,
                    const int Ip,
                    fftw_complex *x_k);
// -------------------------------------------------------------------------
void matR_chol(double *X,const ptrdiff_t N);
//
void matR_chol_banded(double *X,const ptrdiff_t N,const ptrdiff_t Nsubdiag);
// -------------------------------------------------------------------------
void matR_solve_2sys(double *L,const ptrdiff_t N,double *b);
// -------------------------------------------------------------------------
void matR_solve_2sys_banded(double *L,
                            const ptrdiff_t N,
                            const ptrdiff_t Nsubdiag,
                            double *b);
// -------------------------------------------------------------------------
void matC_scale(fftw_complex *X,const ptrdiff_t N,const double alpha);
#endif