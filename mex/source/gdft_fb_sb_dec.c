// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// C implementation to compute the subband components of a set of Nx FIR 
// signals x, given its analysis components. The decomposition is performed
// using the method proposed in:
//
// V. Moles-Cases et al, "Personal Sound Zones by Subband Filtering and Time
// Domain Optimization," in IEEE/ACM Transactions on Audio, Speech, and 
// Language Processing,  2020, doi: 10.1109/TASLP.2020.3023628.
// ................................. INPUTS ................................
// - x_an:    Analysis components of the signal to decompose, 
//            -> [Ix_an x Nx x K/2].
// - Ix_an:   Length of the subband components.
// - Nx:      Number of signals to decompose.
// - hp:      Prototype filter.
// - K:       Number of subbands of the GDFT filter bank.
// - R:       Resampling factor of the GDFT filter bank.
// - Ip:      Length of the prototype filter.
// ................................. OUTPUTS ...............................
// - x_k:     Subband components, -> [Ix_k x Nx x K/2].
// .........................................................................
#include "gdft_fb_sb_dec.h"

// -------------------------------------------------------------------------
void gdft_fb_sb_dec(fftw_complex *x_an,
                    int Ix_an,
                    int Nx,
                    double *hp,
                    const int K,
                    const int R,
                    const int Ip,
                    fftw_complex *x_k){
    
    // ------------------------ General parameters -------------------------
    // Length of the downsampled prototype filter
    int Ip_ds = ceil(((double)Ip)/((double) R));
    // Length of the subband components
    int Ix_k  = Ix_an-Ip_ds+1;
    // FFT size
    int Nfft             = ceil(Ix_an/2.0)*2;;
    // FFT scaling
    double fft_sc        = 1.0/sqrt((double) Nfft);
    // Determine if the corr. matrix is stored in full or banded format
    bool isbanded        = true;
    if(Ip_ds>=Ix_k){
        mexPrintf("aaaaaa\n");
        isbanded         = false;
    }

    // --------------------- Compute correlation matrix --------------------

    // ..................... Downsample prototype filt. ....................
    // Allocate array for downsampled prototype filter
    double *hp_ds = calloc(Nfft,sizeof(double)); 
    // Downsample
    for(int n=0;n<Ip_ds;n++){
        hp_ds[n] = hp[n*R];
    }
    
    // ................. FFT of Downsampled prototype filt. ................
    // Number of positive freq. bins
    int Kfft             = (Nfft/2)+1;
    // Allocate array for the fft of the downsampled prototype filter
    fftw_complex *Hp_ds  = calloc(Kfft,sizeof(fftw_complex)); 
    // Define FFT plans
    fftw_plan fftPlan_r  = fftw_plan_dft_r2c_1d(Nfft,hp_ds,Hp_ds,FFTW_ESTIMATE);
    // Compute FFT of dowsampled prototype filter
    fftw_execute(fftPlan_r);
    // Free
    fftw_destroy_plan(fftPlan_r);
    free(hp_ds);
    // Scale
    matC_scale(Hp_ds,Kfft,fft_sc);

    // .......................... Compute matrix ...........................
    // Allocate array for frequency correlation
    fftw_complex *Rpf    = calloc(Kfft,sizeof(fftw_complex)); 
    // For each positive freq. bin...
    for(int u=0;u<Kfft;u++)
    {
        Rpf[u] = Hp_ds[u]*conj(Hp_ds[u]);
    }
    // Allocate array for time correlation
    double *Rpt          = calloc(Nfft,sizeof(double)); 
    // Initialize IFFT plan  
    fftw_plan ifftPlan_r = fftw_plan_dft_c2r_1d(Nfft,Rpf,Rpt,FFTW_ESTIMATE);
    // Compute time correlation
    fftw_execute(ifftPlan_r);
    // Map correlation matrix
    double *Rp;
    if(isbanded){
       // Allocate array correlation matrix 
       Rp = calloc(Ip_ds*Ix_k,sizeof(double));  
       // Map to toeplitz matrix in banded format (the lower triangular part)
       for(int n=0;n<Ix_k;n++){
          // Copy values to the current column
          memcpy(Rp+n*Ip_ds,Rpt,Ip_ds*sizeof(double));
       }
    }else{
       // Allocate array correlation matrix
       Rp = calloc(Ix_k*Ix_k,sizeof(double)); 
       // Map to the toeplitz matrix (only the lower triangular part)
       for(int n=0;n<Ix_k;n++){
          // Copy values to the current column
          memcpy(Rp+n*Ix_k+n,Rpt,(Ix_k-n)*sizeof(double));
       }
    }
               
    
    // Free
    fftw_destroy_plan(ifftPlan_r);
    free(Rpf);
    free(Rpt);

    // ------------------------ Factorize corr. matrix -----------------------
    // Cholesky factorization
    if(isbanded){
        matR_chol_banded(Rp,Ix_k,Ip_ds-1);
    }else{
        matR_chol(Rp,Ix_k);
    }

    // ---------------------- Compute subband components ---------------------
    // Allocate array for analysis component 
    fftw_complex *X_an    = calloc(Nfft,sizeof(fftw_complex)); 
    // Allocate array for cross-correlation
    fftw_complex *r       = calloc(Nfft,sizeof(fftw_complex)); 
    // Allocate arrays for real and imaginary components of cross-correlation
    // Allocate array for cross-correlation
    double *r_real_imag   = calloc(Ix_k*2,sizeof(double)); 
    // Initialize FFT plans required to compute cross-correlation
    fftw_plan fftPlan_c   = fftw_plan_dft_1d(Nfft,X_an,X_an,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan ifftPlan_c  = fftw_plan_dft_1d(Nfft,r,r, FFTW_BACKWARD, FFTW_ESTIMATE);
    // Scale again
    matC_scale(Hp_ds,Kfft,fft_sc);


    // For each subband...
    for(int k=0;k<K/2;k++){
        // Auxiliary factor
        double sc = 2.0*MPI*((double)R)/((double)K)*(((double)k)+0.5);

        // For each signal...
        for(int s=0;s<Nx;s++){

            // Indices
            int idx_an = k*Ix_an*Nx+s*Ix_an;
            int idx_k  = k*Ix_k*Nx+s*Ix_k;

            // ................. Shift analysis components ..................
            // For each sample in the analysis components...
            for(int n=0;n<Ix_an;n++){
                x_an[idx_an+n] *= cexp(-I*sc*n);
            } 

            // ................. FFT of analysis components .................
            // Copy the analysis component for the current signal and subband
            memset(X_an+Ix_an,0,(Nfft-Ix_an)*sizeof(fftw_complex));
            memcpy(X_an,x_an+idx_an,Ix_an*sizeof(fftw_complex));
            // Compute its fft
            fftw_execute(fftPlan_c);

            // ..................... Cross-correlation .....................
            // Compute freq. domain cross-correlation for positive spectrum
            for(int u=0;u<Kfft;u++)
            {
                r[u] = conj(Hp_ds[u])*X_an[u];
            }
            // Compute freq. domain cross-correlation for negative spectrum
            for(int u=Kfft;u<Nfft;u++)
            {
                r[u] = Hp_ds[Nfft-u]*X_an[u];
            }
            // Compute ifft
            fftw_execute(ifftPlan_c);

            // Split real and imaginary components
            for(int n=0;n<Ix_k;n++){
                r_real_imag[n]       = creal(r[n]);
                r_real_imag[n+Ix_k]  = cimag(r[n]);
            }

            // ........................ Solve system ........................
            // Solver system for real and imaginary components
            if(isbanded){
                matR_solve_2sys_banded(Rp,Ix_k,Ip_ds-1,r_real_imag);
            }else{
                matR_solve_2sys(Rp,Ix_k,r_real_imag);
            }

            // ...................... Merge and shift .......................
            for(int n=0;n<Ix_k;n++){
                x_k[idx_k+n] = (r_real_imag[n]+I*r_real_imag[n+Ix_k])*cexp(I*sc*n);
            }
        }
    }

    // Free
    fftw_destroy_plan(fftPlan_c);
    fftw_destroy_plan(ifftPlan_c);
    free(Hp_ds);
    free(Rp);
    free(X_an);
    free(r);
    free(r_real_imag);

}


// -------------------------------------------------------------------------
// Compute the Cholesky dec. of a  symmetric positive-definite real matrix 
// X (only the lower triangular components are required)
void matR_chol(double *X,const ptrdiff_t N){

    // Initialize info variable for the factorization
    ptrdiff_t info_Chol;
    const char lowTri = 'L';
    // Compute Cholesky factorization
    dpotrf(&lowTri,&N,X,&N,&info_Chol);
}

// -------------------------------------------------------------------------
// Compute the Cholesky dec. of a  symmetric positive-definite real matrix 
// X that is stored banded format (only the lower triangular components are 
// required)
void matR_chol_banded(double *X,const ptrdiff_t N,const ptrdiff_t Nsubdiag){

    // Initialize info variable for the factorization
    ptrdiff_t info_Chol;
    const char lowTri = 'L';
    ptrdiff_t ldab = Nsubdiag+1;
    // Compute Cholesky factorization
    dpbtrf(&lowTri,&N,&Nsubdiag,X,&ldab,&info_Chol);
}

// -------------------------------------------------------------------------
// Solve the linear system Ax=[b1 b2] given the lower triangular factor of  
// the Cholesky decomposition of matrix A.
void matR_solve_2sys(double *L,const ptrdiff_t N,double *b){

    // Initialize info variable for the inverse computation
    ptrdiff_t info_sys;
    const char lowTri      = 'L';
    const ptrdiff_t long_2 = 2;
    // Solve system
    dpotrs(&lowTri,&N,&long_2,L,&N,b,&N,&info_sys);
}

// -------------------------------------------------------------------------
// Solve the linear system Ax=[b1 b2] given the lower triangular factor of  
// the Cholesky decomposition of matrix A stored in banded format
void matR_solve_2sys_banded(double *L,
                            const ptrdiff_t N,
                            const ptrdiff_t Nsubdiag,
                            double *b){

    // Initialize info variable for the inverse computation
    ptrdiff_t info_sys;
    const char lowTri      = 'L';
    const ptrdiff_t long_2 = 2;
    ptrdiff_t ldab         = Nsubdiag+1;
    // Solve system
    dpbtrs(&lowTri,&N,&Nsubdiag,&long_2,L,&ldab,b,&N,&info_sys);
}

// -------------------------------------------------------------------------
// Scale N values of X, i.e., X(0:N-1) = alpha*X(0:N-1)
void matC_scale(fftw_complex *X,const ptrdiff_t N,const double alpha){
  const ptrdiff_t long_1   = 1;
  zdscal(&N,&alpha,(double*)X,&long_1);
}





