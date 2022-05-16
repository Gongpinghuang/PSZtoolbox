// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// C implementation of the superfast solver for complex-valued LS problems
// proposed in:
//
// M. Poletti et al, "A Superfast Toeplitz Matrix Inversion Method for 
// Single- and Multi-Channel Inverse Filters and its Application to Room
// Equalization," in  IEEE/ACM Transactions on Audio Speech and Language 
// Processing, 29, 3144–3157, 2021. 
// https://doi.org/10.1109/TASLP.2021.3120650
// ................................. INPUTS ................................
// - h:       RIR, -> [Ih x M x L].
// - d:       Target response, -> [Id x M].
// - L:       Number of loudspeakers.
// - M:       Number of control points.
// - Ih:      Length of RIR.
// - Ig:      Length of filters.
// - Id:      Length of target.
// - Nfft:    FFT size.
// - P:       Approximation order.
// - beta:    Regularization factor, beta>0.
// ................................. OUTPUTS ...............................
// - g:       Computed filters, -> [Ig x L].
// .........................................................................
#include "superfast_LS_complex.h"

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
                          fftw_complex *g){
	
    // ----------------------- FFT-related parameters ----------------------
    // Check
    if(Nfft<Id){
    	mexErrMsgTxt("superfast_LS_complex.c: Incorrect FFT size");
    }
    // FFT scaling
    double fft_sc         = 1.0/Nfft;
    // Define FFT plans
    fftw_complex *x       = calloc(Nfft,sizeof(fftw_complex)); 
    fftw_complex *X       = calloc(Nfft,sizeof(fftw_complex)); 
    fftw_plan fftPlan     = fftw_plan_dft_1d(Nfft,
    	                                       x,
    	                                       X,
    	                                       FFTW_FORWARD,
    	                                       FFTW_ESTIMATE); 
    fftw_plan fftPlan_inv = fftw_plan_dft_1d(Nfft,
                      									     X,
                      									     x,
                      									     FFTW_BACKWARD,
                      									     FFTW_ESTIMATE); 

    // ----------------------- Compute FFT of inputs -----------------------
    // Control point and loudspeaker index
    int m,l;
    // Allocate required variables  
    fftw_complex *hf = fftw_malloc(M*L*Nfft*sizeof(fftw_complex)); 
    fftw_complex *df = fftw_malloc(M*Nfft*sizeof(fftw_complex)); 
    // For each loudspeaker and control point in the RIR...
     for(int i=0;i<M*L;i++){
  	  // Compute FFT of RIR
  	  // Control point 
	  m = (int) floor(((double)i)/((double)L));
	  // Control point loudspeaker index
	  l = i%L;
      // Copy the non zero elements of the RIR
      memcpy(x,h+i*Ih,Ih*sizeof(fftw_complex));
      // Compute FFT
      fftw_execute(fftPlan);
      matC_copy_step(X,hf+M*l+m,Nfft,1,M*L);
    }
    // For each control point in the target...
    for(int m=0;m<M;m++){
      // Compute FFT of target
      // Copy the non zero elements of the target
      memcpy(x,d+m*Id,Id*sizeof(fftw_complex));
      // Compute FFT
      fftw_execute(fftPlan);
      matC_copy_step(X,df+m,Nfft,1,M);
	}

    // ----------- Compute filters using freq. domain formulation ----------
    // Allocate required variables
    fftw_complex *gf         = fftw_malloc(Nfft*L*sizeof(fftw_complex)); 
    fftw_complex *Rf         = fftw_malloc(M*M*sizeof(fftw_complex)); 
    fftw_complex *Hf_ps      = fftw_malloc(L*M*sizeof(fftw_complex));  
    fftw_complex *Hf_ps_x_Hf = fftw_malloc(L*L*Nfft*sizeof(fftw_complex));  
    // Size vectors
    size_t size_M_x_L[2]     = {M,L}; 
    size_t size_L_x_M[2]     = {L,M}; 
    size_t size_L_x_L[2]     = {L,L}; 
    size_t size_M_x_M[2]     = {M,M}; 

    // For each freq. bin...
    for(int f=0;f<Nfft;f++){
      // Compute Rf = Hf*Hf' for the f-th bin
      matC_X_times_X_herm(hf+M*L*f,size_M_x_L,Rf);
      // Add regularization to Rf for the f-th bin
      matC_addEye(Rf,size_M_x_M,beta);
      // Compute Hf^pseudo = Hf'*(Hf*Hf'+beta*I)^{-1} 
      memcpy(Hf_ps,hf+M*L*f,M*L*sizeof(fftw_complex));
      matC_linsolve(Rf,size_M_x_M,Hf_ps,size_M_x_L);
      // Compute Hf^pseudo*Hf
      matC_X_herm_times_Y(Hf_ps,
                          size_M_x_L,
                          hf+M*L*f,
                          size_M_x_L,
                          Hf_ps_x_Hf+L*L*f);
      // Compute filters coefficients for the f-th bin
      matC_X_herm_times_y_step(Hf_ps, 
                               size_M_x_L,
                               df+M*f, 
                               M,
                               gf+f,
                               Nfft);
    }
    // Allocate time-domain filters
    fftw_complex *g_fd  = calloc(Nfft*L,sizeof(fftw_complex)); 
    // For each loudspeaker...
    for(int l=0;l<L;l++){
      // Compute filters time domain 
      fftw_execute_dft(fftPlan_inv,gf+Nfft*l,g_fd+Nfft*l);
      // Correct fft scaling
      matC_scale(g_fd+l*Nfft,Nfft,fft_sc);
    }
    // Free memory
    fftw_free(Rf);
    fftw_free(Hf_ps);
    fftw_free(hf);
    fftw_free(df);
    fftw_free(gf);

    // -------------------------- Correct filters --------------------------
    // Allocate residual components for correction
    fftw_complex *rt     = calloc(Nfft*L,sizeof(fftw_complex)); 
    fftw_complex *rf_aux = fftw_malloc(Nfft*L*sizeof(fftw_complex));
    fftw_complex *rf     = fftw_malloc(Nfft*L*sizeof(fftw_complex));

    // Obtain initial time-domain residuals for this iteration
    // For each loudspeaker...
    for(int l=0;l<L;l++){
      memcpy(rt+l*Nfft,g_fd+l*Nfft,Nfft*sizeof(fftw_complex));
    }

    // Iterate until the selected approximation order is reached...
    for(int p=0;p<P;p++){


        // ................. Compute freq. domain residuals ................
        // For each loudspeaker...
        for(int l=0;l<L;l++){
            // Set to 0 Ig firts samples
            memset(rt+Nfft*l,0.0,Ig*sizeof(fftw_complex));
            // Compute FFT
            fftw_execute_dft(fftPlan,rt+Nfft*l,rf_aux+Nfft*l);
        }

        // ................. Correct freq. domain residuals ................ 
        // For each loudspeaker...
        for(int f=0;f<Nfft;f++){
            matC_X_times_y_step(Hf_ps_x_Hf+L*L*f, 
                                size_L_x_L,
                                rf_aux+f, 
                                L,
                                Nfft,
                                rf+f,
                                Nfft);
        }

        // .................. Correct time-domain filters .................. 
        // For each loudspeaker...
        for(int l=0;l<L;l++){
            // Obtain initial time-domain residuals 
            fftw_execute_dft(fftPlan_inv,rf+Nfft*l,rt+Nfft*l);
            // Correct fft scaling
      		  matC_scale(rt+l*Nfft,Nfft,fft_sc);
            // Correct coefficients
            for(int n=0;n<Ig;n++){
            	g_fd[Nfft*l+n]+=rt[Nfft*l+n];
            }
        }
    }

    // -------------------------- Truncate filters -------------------------
    // Truncate filters and store in output
    // For each loudspeaker...
    for(int l=0;l<L;l++){
     	// Keep first Ig samples
      memcpy(g+l*Ig,g_fd+l*Nfft,Ig*sizeof(fftw_complex));
    }

    // ---------------------------- Free-memory ----------------------------
    free(g_fd);
    fftw_free(Hf_ps_x_Hf);
    fftw_free(rf);
	  fftw_free(rf_aux);
	  free(x);
	  free(X);
    fftw_destroy_plan(fftPlan);
    fftw_destroy_plan(fftPlan_inv);

}

// -------------------------------------------------------------------------
// Copy N elements of X to Y with step, i.e., 
// Y(0:stepY:stepY*(N-1)) = X(0:stepX:stepX*(N-1)
void matC_copy_step(fftw_complex *X,
                    fftw_complex *Y,
                    const ptrdiff_t N,
                    const ptrdiff_t stepX,
                    const ptrdiff_t stepY){
  zcopy(&N,(double*)X,&stepX,(double*)Y,&stepY);
}

// -------------------------------------------------------------------------
// Product of two matrices: Z = X*X', where the result Z is hermitian and it
// is stored only in the upper triangular part
void matC_X_times_X_herm(fftw_complex *X, 
                         const size_t* sizeX,
                         fftw_complex *Z){

    const ptrdiff_t n   = sizeX[0];
    const ptrdiff_t f   = sizeX[1];
    const double alpha  = 1.0;
    const double omega  = 0.0;
    const char UpTriang = 'U';
    const char NonT     = 'N';

    zherk(&UpTriang,
          &NonT,
          &n,
          &f,
          &alpha,
          (double*)X,
          &n,
          &omega,
          (double*)Z,
          &n);
}

// -------------------------------------------------------------------------
// Product of two matrices: Z = X'*Y
void matC_X_herm_times_Y(fftw_complex *X, 
                         const size_t* sizeX,
                         fftw_complex *Y, 
                         const size_t* sizeY,
                         fftw_complex *Z){

   const ptrdiff_t m        = sizeX[1];
   const ptrdiff_t p        = sizeX[0];
   const ptrdiff_t n        = sizeY[1];
   const char Herm          = 'C';
   const char NonT          = 'N';
   const fftw_complex alpha = 1.0+0.0*I;
   const fftw_complex omega = 0.0+0.0*I;

   zgemm(&Herm,
         &NonT,
         &m,
         &n,
         &p,
         (double*)&alpha,
         (double*)X,
         &p,
         (double*)Y,
         &p,
         (double*)&omega,
         (double*)Z,
         &m);

}

// -------------------------------------------------------------------------
// Product of a matrix and a vector: z(0:stepZ:end) = X'*y
void matC_X_herm_times_y_step(fftw_complex *X, 
                              const size_t *sizeX,
                              fftw_complex *y, 
                              const ptrdiff_t Ly,
                              fftw_complex *z,
                              const ptrdiff_t stepZ){

  const ptrdiff_t m        = sizeX[0];
  const ptrdiff_t n        = sizeX[1];
  const char Herm          = 'C';
  const fftw_complex alpha = 1.0+0.0*I;
  const fftw_complex omega = 0.0+0.0*I;
  const ptrdiff_t long_1   = 1;

  zgemv(&Herm,
        &m,
        &n,
        (double*)&alpha,
        (double*)X,
        &m,
        (double*)y,
        &long_1,
        (double*)&omega,
        (double*)z,
        &stepZ);
}

// -------------------------------------------------------------------------
// Product of a matrix and a vector: z(0:stepZ:end) = X*y(0:stepY:end)
void matC_X_times_y_step(fftw_complex *X, 
                         const size_t* sizeX,
                         fftw_complex *y, 
                         const ptrdiff_t Ly,
                         const ptrdiff_t stepY,
                         fftw_complex *z,
                         const ptrdiff_t stepZ){

  const ptrdiff_t m        = sizeX[0];
  const ptrdiff_t n        = sizeX[1];
  const fftw_complex alpha = 1.0+0.0*I;
  const fftw_complex omega = 0.0+0.0*I;
  const char NonT          = 'N';

  zgemv(&NonT,
      &m,
      &n,
      (double*)&alpha,
      (double*)X,
      &m,
      (double*)y,
      &stepY,
      (double*)&omega,
      (double*)z,
      &stepZ);
}

// -------------------------------------------------------------------------
// Add an an indentity matrix to a square matrix, i.e., X = X+alpha*I
void matC_addEye(fftw_complex *X, const size_t* sizeX,const double alpha){
  int idx        = 0;
  int aux        = (1+sizeX[0]);
  fftw_complex a = alpha+I*0.0;
  for(int i=0;i<sizeX[0];i++){
    X[idx] = X[idx]+a;
    idx    += aux;
  }
}

// -------------------------------------------------------------------------
// Compute the solution of the linear system A*X=B, where A is a hermitian
// positive-definite matrix. 
void matC_linsolve(fftw_complex *A,
                   const size_t* sizeA,
                   fftw_complex *B,
                   const size_t* sizeB){

    const ptrdiff_t m   = sizeA[0];
    const ptrdiff_t p   = sizeB[0];
    const ptrdiff_t n   = sizeB[1];
    const char UpTriang = 'U';
    ptrdiff_t info_Chol;
    ptrdiff_t info_sys;
    // Compute Cholesky factorization
    zpotrf(&UpTriang,&m,(double*)A,&m,&info_Chol);
    // Solve system
    zpotrs(&UpTriang,&m,&n,(double*)A,&m,(double*)B,&p,&info_sys);
}

// -------------------------------------------------------------------------
// Scale N values of X, i.e., X(0:N-1) = alpha*X(0:N-1)
void matC_scale(fftw_complex *X,const ptrdiff_t N,const double alpha){
  const ptrdiff_t long_1   = 1;
  zdscal(&N,&alpha,(double*)X,&long_1);
}
