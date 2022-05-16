// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// MEX implementation to compute the subband components of a set of Nx FIR 
// signals x, given its analysis components. The decomposition is performed
// using the method proposed in:
//
// V. Moles-Cases et al, "Personal Sound Zones by Subband Filtering and Time
// Domain Optimization," in IEEE/ACM Transactions on Audio, Speech, and 
// Language Processing,  2020, doi: 10.1109/TASLP.2020.3023628.
// -------------------------------------------------------------------------
#include "gdft_fb_sb_dec.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,const mxArray *prhs[])
{
    // ------------------------ CHECK INPUT/OUTPUTS ------------------------
    // Check for proper number of arguments
    if(nrhs!=4) {
        mexErrMsgTxt("gdft_fb_sb_dec_mex.c: 4 inputs are required.");
    } 
    
    if(nlhs!=1) {
        mexErrMsgTxt("gdft_fb_sb_dec_mex.c: 1 output required.");
    }

    // --------------------------- READ INPUTS -----------------------------
    // Variagles required to read inputs 
    int Ix_an,Nx,K,R,Ip,Ix_k;
    double *hp;
    fftw_complex *x_an;
    size_t nDim_hp,nDim_x_an; 
    const size_t *size_hp,*size_x_an; 
       
    // 1st input (Analsys components [Ix_an x Nx x K/2]) 
    x_an      = (fftw_complex*) mxGetComplexDoubles(prhs[0]);
    nDim_x_an = mxGetNumberOfDimensions(prhs[0]);
    size_x_an = mxGetDimensions(prhs[0]); 
    // 2nd input (prototype filter [Id x M]) 
    hp        = (double*) mxGetDoubles(prhs[1]);
    nDim_hp   = mxGetNumberOfDimensions(prhs[1]);
    size_hp   = mxGetDimensions(prhs[1]); 
    // 3rd input (number of subbands) 
    K         = (int) mxGetScalar(prhs[2]);
    // 4th input (resampling factor) 
    R         = (int) mxGetScalar(prhs[3]);
    // Related parameters
    Ix_an     = size_x_an[0];
    Nx        = size_x_an[1];
    Ip        = size_hp[0];
    Ix_k      = Ix_an-(int)ceil((double)Ip/(double)R)+1;

    // Check parameters
    if(nDim_x_an!=3){
        mexErrMsgTxt("gdft_fb_sb_dec_mex.c: Incorrect size for x_an"); 
    }
    if(size_x_an[2]!=K/2){
        mexErrMsgTxt("gdft_fb_sb_dec_mex.c: Incorrect size for x_an"); 
    }
    if(nDim_hp!=2){
        mexErrMsgTxt("gdft_fb_sb_dec_mex.c: Incorrect size for hp"); 
    }
    if(size_hp[1]!=1){
        mexErrMsgTxt("gdft_fb_sb_dec_mex.c: Incorrect size for hp"); 
    }
 
    // --------------------------- DEFINE OUTPUT ---------------------------
    // Output (estimated filter [Ig x L])
    fftw_complex *x_k;  
    // Set output dimensions
    size_t nDim_x_k = 3;
    const size_t size_x_k[3] = {Ix_k,Nx,K/2};
    // Create output array
    plhs[0] = mxCreateNumericArray(nDim_x_k,
    	                           (const mwSize*) size_x_k,
    	                           mxDOUBLE_CLASS,mxCOMPLEX);
    x_k     =  (fftw_complex*) mxGetComplexDoubles(plhs[0]);     

    // ------------------------------ COMPUTE ------------------------------ 
    // Compute subband components
    gdft_fb_sb_dec(x_an,Ix_an,Nx,hp,K,R,Ip,x_k);

}

 

                

