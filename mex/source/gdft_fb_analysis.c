// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// C implementation of the polyphase network for the analysis stage  of a 
// GDFT filter bank proposed in:
// 
// S.Weiss,"On adaptive filtering in oversampled subbands," Ph.D 
// dissertation, 1998.
// ................................. INPUTS ................................
// - x:           Input signals, -> [Ix x Nx].
// - Ix:          Length of the input signals.
// - Nx:          Number of input signals.
// - K:           Number of subbands of the GDFT filter bank.
// - R:           Resampling factor of the GDFT filter bank.
// - Ip:          Length of the prototype filter.
// - poly_Nnon0:  Number of non-zero polyphase components. 
// - poly_coeff:  Value of the non-zero polyphase components,
//                -> [poly_Nnon0 x 1].
// - poly_tap:    Tap index of the non-zero polyphase components,
//                -> [poly_Nnon0 x 1].
// - poly_branch: Branch index of the non-zero polyphase components,
//                -> [poly_Nnon0 x 1].
// ................................. OUTPUTS ...............................
// - y:           Filtered signals, -> [K x Iy x Nx].
// .........................................................................
#include "gdft_fb_analysis.h"

// -------------------------------------------------------------------------
void gdft_fb_analysis(double *x,
					  int Ix,
					  int Nx,
					  int K,
					  int R,
					  int Ip,
					  int poly_Nnon0,
					  double *poly_coeff,
					  int *poly_tap,
					  int *poly_branch,
					  fftw_complex *y){

	// ---------------------------- Initialize -----------------------------
	// Required Frequency shift
    fftw_complex *freqShift =  fftw_malloc(K*sizeof(fftw_complex));  
    for(int k=0;k<K;k++){
    	freqShift[k] =  cexp(I*2*MPI*0.5*k/K);
    }
	// Lenght of the output signal
	int Iy                 = (int) ceil((double)(Ix+Ip-1)/(double)R); 
	// Allocate auxiliary variable
	fftw_complex *y_aux    =  calloc(K*Iy*Nx,sizeof(fftw_complex)); 
    // FFT related variables 
    fftw_complex *fft_in   = fftw_malloc(K*sizeof(fftw_complex));
    fftw_complex *fft_out  = fftw_malloc(K*sizeof(fftw_complex));
    fftw_plan fftPlan      = fftw_plan_dft_1d(K,
    	                                      fft_in,
    	                                      fft_out,
    	                                      FFTW_BACKWARD,
    	                                      FFTW_ESTIMATE);

    // ------------------------------ Filter -------------------------------
    // Initialize
    int idx_x,idx_y;
    // For each signal...
    for(int s=0;s<Nx;s++){
    	// For each output sample...
    	for(int n=0;n<Iy;n++){
	         // For all the non-zero taps in the polyphase network...
	        for(int i=0;i<poly_Nnon0;i++){
	            // Index of x that is multiplied by the i-th non-zero 
	            // polyphase component
	            idx_x = s*Ix+R*(n-poly_tap[i])-(poly_branch[i]%R);
	            // Filter only if valid indices are obtained
	            if(idx_x>=s*Ix && idx_x<((s+1)*Ix)){
	            	// Index of the output signal in which the product is 
	            	// added
	            	idx_y = s*Iy*K+(n*K)+(poly_branch[i]%K);
	            	// Multiply and add
	            	y_aux[idx_y] += x[idx_x]*poly_coeff[i];
	            }
	        }  
	        // Frequency shift
	        int idx_y = s*Iy*K+(n*K);
	        for(int k=0;k<K;k++){
	        	y_aux[idx_y+k] = y_aux[idx_y+k]*freqShift[k];
	        }
	        // Apply IFFT
	        fftw_execute_dft(fftPlan,y_aux+idx_y,y+idx_y); 
    	}
    }

    // --------------------------- Free Memory -----------------------------
    fftw_free(freqShift);
    fftw_free(fft_in);
    fftw_free(fft_out);
    free(y_aux);
    fftw_destroy_plan(fftPlan);
}























