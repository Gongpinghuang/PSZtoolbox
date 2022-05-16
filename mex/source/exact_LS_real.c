// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// C implementation to compute the exact solution for the real-valued Least
// Squares (LS) problem related with wPM using the Cholesky 
// decomposition. The correlation matrices and vectors are computed using 
// the FFT.
// ................................. INPUTS ................................
// - h:       RIR, -> [Ih x M x L].
// - d:       Target response, -> [Id x M].
// - L:       Number of loudspeakers.
// - M:       Number of control points.
// - Ih:      Length of RIR.
// - Ig:      Length of filters.
// - Id:      Length of target.
// - beta:    Regularization factor, beta>0.
// ................................. OUTPUTS ...............................
// - g:       Computed filters, -> [Ig x L].
// .........................................................................
#include "exact_LS_real.h"

// -------------------------------------------------------------------------
void exact_LS_real(double *h,
                   double *d,
                   const int L,
                   const int M,
                   const int Ih,
                   const int Ig,
                   const int Id,
                   const double beta,
                   double *g){
	
    // ----------------------- FFT-related parameters ----------------------
    // FFT size
    int Nfft              = 2*((int)ceil((double)Id/(double)2));
    // Number of positive frequency bins
    int Kfft              = (Nfft/2)+1;
    // Define FFT plans
    double *x             = calloc(Nfft,sizeof(double)); 
    fftw_complex *X       = calloc(Kfft,sizeof(fftw_complex)); 
    fftw_plan fftPlan     = fftw_plan_dft_r2c_1d(Nfft,
          	                                     x,
          	                                     X,
          	                                     FFTW_ESTIMATE); 
    fftw_plan fftPlan_inv = fftw_plan_dft_c2r_1d(Nfft,
          									     X,
          									     x,
          									     FFTW_ESTIMATE); 

    // ----------------------- Compute FFT of inputs -----------------------
    // Control point and loudspeaker index
    int m,l;
    // Allocate required variables  
    fftw_complex *hf = fftw_malloc(M*L*Kfft*sizeof(fftw_complex)); 
    fftw_complex *df = fftw_malloc(M*Kfft*sizeof(fftw_complex)); 
    // For each loudspeaker and control point in the RIR...
    for(int i=0;i<M*L;i++){
  	  // Compute FFT of RIR
    	// Control point 
  	  m = (int) floor(((double)i)/((double)L));
  	  // Control point loudspeaker index
  	  l = i%L;
      // Copy the non zero elements of the RIR
      memcpy(x,h+i*Ih,Ih*sizeof(double));
      // Compute FFT
      fftw_execute_dft_r2c(fftPlan,x,hf+i*Kfft);
    }
    // For each control point in the target...
    for(int m=0;m<M;m++){
      // Compute FFT of target
      // Copy the non zero elements of the target
      memcpy(x,d+m*Id,Id*sizeof(double));
      // Compute FFT
      fftw_execute_dft_r2c(fftPlan,x,df+m*Kfft);
  	}

    // -------------------- Compute correlation matrix ---------------------
    // It is important to highlight that only the upper-triangular part of 
    // the matrix is computed, since it is symmetric
    // Initialize 
    double *R          = calloc(L*Ig*L*Ig,sizeof(double)); 
    double *Rt_i       = calloc(2*Ig-1,sizeof(double)); 

    // For each unique combination of loudspeakers...
    for(int l1=0;l1<L;l1++){
      for(int l2=l1;l2<L;l2++){

          // ................. Compute freq. domain corr. ..................
          // For each positive freq. bin...
          for(int k=0;k<Kfft;k++)
          {
              // Mean correlation over all control points for the current 
              // combination of loudspeakers
              matC_x_herm_times_y_step(hf+l1*Kfft+k,
                                       hf+l2*Kfft+k,
                                       X+k,
                                       M,
                                       Kfft*L,
                                       Kfft*L);
          }

          // ................. Compute time domain corr. ...................
          // Compute IFFT
          fftw_execute(fftPlan_inv);
          // Keep only the required correlation values
          memcpy(Rt_i+Ig-1,x,Ig*sizeof(double));
          for(int n=0;n<Ig-1;n++){
            Rt_i[n] = x[Nfft-Ig+1+n];
          }
          // Add regularization
          if(l1==l2){
            Rt_i[Ig-1] = Rt_i[Ig-1]+beta*Nfft;
          }

          //................. Map to the toeplitz matrix ...................
          // Row and column idx of the first element in the sub-matrix for 
          // the current combination
          int row = l1*Ig;
          int col = l2*Ig;
          // For all the columns in the sub-matrix...
          for(int n=0;n<Ig;n++){
            // Copy values to the current column
            memcpy(R+col*L*Ig+row,Rt_i+Ig-1-n,Ig*sizeof(double));
            // Increase column index
            col++;
          }

        }
    }
    free(Rt_i);

    // -------------------- Compute correlation vector ---------------------
    // For each loudspeaker...
    for(int l=0;l<L;l++){

      // ................... Compute freq. domain corr. ....................
      // For each positive freq. bin...
      for(int k=0;k<Kfft;k++)
      {
          // Mean cross-correlation over all control points for the current 
          // loudspeaker
          matC_x_herm_times_y_step(hf+l*Kfft+k,
                                   df+k,
                                   X+k,
                                   M,
                                   Kfft*L,
                                   Kfft);
      }

      // ................... Compute time domain corr. .....................
      // Compute IFFT
      fftw_execute(fftPlan_inv);
      // Map required indices
      memcpy(g+l*Ig,x,Ig*sizeof(double));
    }

    // -------------------------- Compute Filters --------------------------
    // Compute Cholesky decomposition 
    matR_chol(R,L*Ig);
    // Solve the system of linear equations
    matR_solve_sys(R,L*Ig,g);

    // ---------------------------- Free-memory ----------------------------
	  free(x);
	  free(X);
    fftw_destroy_plan(fftPlan);
    fftw_destroy_plan(fftPlan_inv);
    free(R);
}

// -------------------------------------------------------------------------
// Perform the product of M elements of two colum vectors x and y as 
// z = x(0:stepx:M*stepx-1)^H*y(0:stepy:M*stepy-1)
void matC_x_herm_times_y_step(fftw_complex *x,
                              fftw_complex *y,
                              fftw_complex *z,
                              const ptrdiff_t M,
                              const ptrdiff_t stepx,
                              const ptrdiff_t stepy){

  zdotc((doublecomplex*)z,&M,(double*)x,&stepx,(double*)y,&stepy);
}

//.............................................................................
// Compute the Cholesky dec. of a  symmetric positive-definite REAL matrix X.
void matR_chol(double *X,const ptrdiff_t N){

    // Initialize info variable for the factorization
    ptrdiff_t info_Chol;
    const char upTri = 'U';
    // Compute Cholesky factorization
    dpotrf(&upTri,&N,X,&N,&info_Chol);
}

//.............................................................................
// Solve the linear system Ax=b given the upper triangular factor of the 
// Cholesky decomposition of matrix A.
void matR_solve_sys(double *U,
                    const ptrdiff_t N,
                    double *b){

    // Initialize info variable for the inverse computation
    ptrdiff_t info_sys;
    const char upTri       = 'U';
    const ptrdiff_t long_1 = 1;
    // Solve system
    dpotrs(&upTri,&N,&long_1,U,&N,b,&N,&info_sys);

}

// -------------------------------------------------------------------------
// Scale N values of X, i.e., X(0:N-1) = alpha*X(0:N-1)
void matR_scale(double *X,const ptrdiff_t N,const double alpha){
  const ptrdiff_t long_1   = 1;
  dscal(&N,&alpha,(double*)X,&long_1);
}





