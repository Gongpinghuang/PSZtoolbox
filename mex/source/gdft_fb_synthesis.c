// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// C implementation of the polyphase network for the synthesis stage  of a 
// GDFT filter bank proposed in:
// 
// S.Weiss,"On adaptive filtering in oversampled subbands," Ph.D 
// dissertation, 1998.
// ................................. INPUTS ................................
// - x:           Input signals, -> [K x Ix x Nx]. 
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
// - y:           Filtered signals, -> [Iy x Nx].
// .........................................................................
#include "gdft_fb_synthesis.h"

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
					   double *y){

	// ---------------------------- Initialize -----------------------------
	// Required Frequency shift
    fftw_complex *freqShift =  fftw_malloc(K*sizeof(fftw_complex));  
    for(int k=0;k<K;k++){
    	freqShift[k] =  cexp(-I*2*MPI*0.5*k/K);
    }
	// Lenght of the output signal
	int Iy                 = Ix*R+Ip-1; 

	// Allocate auxiliary variable
	fftw_complex *x_aux    =  calloc(K*Ix*Nx,sizeof(fftw_complex)); 
    // FFT related variables 
    fftw_complex *fft_in   = fftw_malloc(K*sizeof(fftw_complex));
    fftw_complex *fft_out  = fftw_malloc(K*sizeof(fftw_complex));
    fftw_plan fftPlan      = fftw_plan_dft_1d(K,
    	                                      fft_in,
    	                                      fft_out,
    	                                      FFTW_FORWARD,
    	                                      FFTW_ESTIMATE);
    // Find the maximum tap index
    int max_tap = 0;
    for(int t=0;t<poly_Nnon0;t++){
    	if(poly_tap[t]>max_tap){
    		max_tap = poly_tap[t];
    	}
    }

    // ------------------------------ Filter -------------------------------
    // Initialize
    int idx_x,idx_y,sample_idx;
    int offset   = (max_tap+1)*R-Ip;
    // For each signal...
    for(int s=0;s<Nx;s++){
    	// For each input sample...
    	for(int n=0;n<Ix+max_tap;n++){
    		if(n<Ix){
	    		idx_x = s*Ix*K+n*K;
	    		// Apply FFT
		        fftw_execute_dft(fftPlan,x+idx_x,x_aux+idx_x); 
		        // Frequency shift
		        for(int k=0;k<K;k++){
		        	x_aux[idx_x+k] = x_aux[idx_x+k]*freqShift[k];
		        }
    		}
	         // For all the non-zero taps in the polyphase network...
	        for(int i=0;i<poly_Nnon0;i++){
	        	// Index of x that is multiplied by the i-th non-zero 
	            // polyphase component
	            idx_x = s*Ix*K+(n-poly_tap[i])*K+(poly_branch[i]%K);
	            // Index of the output signal in which the product is 
	            // added
	            idx_y = s*Iy+n*R+R-1-(poly_branch[i]%R)-offset;
	            // Filter only if valid indices are obtained
            	if(idx_y>=s*Iy && idx_y<(s+1)*Iy &&
            	               idx_x>=s*Ix*K && idx_x<(s+1)*Ix*K){
            	   	// Multiply and add
	                y[idx_y] += x_aux[idx_x]*poly_coeff[i];
	            }
	        }  
    	}
    }

    // --------------------------- Free Memory -----------------------------
    fftw_free(freqShift);
    fftw_free(fft_in);
    fftw_free(fft_out);
    free(x_aux);
    fftw_destroy_plan(fftPlan);
}























