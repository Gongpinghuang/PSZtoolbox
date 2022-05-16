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

#ifndef SUPERFAST_LS_COMPLEX
#define SUPERFAST_LS_COMPLEX

// -------------------------------------------------------------------------
void superfast_LS_complex(fftw_complex *h,
                          fftw_complex *d,
                          const int L,
                          const int M,
                          const int Ih,
                          const int Ig,
                          const int Id,
                          const int Nfft,
                          const int P,
                          const double beta,
                          fftw_complex *gout);
// -------------------------------------------------------------------------
void matC_copy_step(fftw_complex *X,
                    fftw_complex *Y,
                    const ptrdiff_t N,
                    const ptrdiff_t stepX,
                    const ptrdiff_t stepY);
// -------------------------------------------------------------------------
void matC_X_times_X_herm(fftw_complex *X, 
                         const size_t* sizeX,
                         fftw_complex *Z);
// -------------------------------------------------------------------------
void matC_X_herm_times_Y(fftw_complex *X, 
                         const size_t* sizeX,
                         fftw_complex *Y, 
                         const size_t* sizeY,
                         fftw_complex *Z);
// -------------------------------------------------------------------------
void matC_X_herm_times_y_step(fftw_complex *X, 
                              const size_t *sizeX,
                              fftw_complex *y, 
                              const ptrdiff_t Ly,
                              fftw_complex *z,
                              const ptrdiff_t stepZ);
// -------------------------------------------------------------------------
void matC_X_times_y_step(fftw_complex *X, 
                        const size_t* sizeX,
                        fftw_complex *y, 
                        const ptrdiff_t Ly,
                        const ptrdiff_t stepY,
                        fftw_complex *z,
                        const ptrdiff_t stepZ);
// -------------------------------------------------------------------------
void matC_linsolve(fftw_complex *A,
                   const size_t* sizeA,
                   fftw_complex *B,
                   const size_t* sizeB);
// -------------------------------------------------------------------------
void matC_addEye(fftw_complex *X, const size_t* sizeX,const double alpha);
// -------------------------------------------------------------------------
void matC_scale(fftw_complex *X,const ptrdiff_t N,const double alpha);

#endif