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

#ifndef GDFT_FB_SYNTHESIS
#define GDFT_FB_SYNTHESIS

// -------------------------------------------------------------------------
void gdft_fb_synthesis(fftw_complex *x,
                       int Ix,
                       int Nx,
                       int K,
                       int R,
                       int Ip,
                       int poly_Nnon0,
                       double *poly_coeff,
                       int *poly_tap,
                       int *poly_branch,
                       double *y);

#endif