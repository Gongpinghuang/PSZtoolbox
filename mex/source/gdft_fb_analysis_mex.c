// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// MEX implementation of the polyphase network for the analysis stage  of a 
// GDFT filter bank proposed in:
// 
// S.Weiss,"On adaptive filtering in oversampled subbands," Ph.D 
// dissertation, 1998.
// -------------------------------------------------------------------------
#include "gdft_fb_analysis.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,const mxArray *prhs[])
{
    // ------------------------ CHECK INPUT/OUTPUTS ------------------------
    // Check for proper number of arguments
    if(nrhs!=8) {
        mexErrMsgTxt("gdft_fb_analysis_mex.c: 8 inputs are required.");
    } 
    
    if(nlhs!=1) {
        mexErrMsgTxt("gdft_fb_analysis_mex.c: 1 output required.");
    }

    // --------------------------- READ INPUTS -----------------------------
    // Variagles required to read inputs 
    int Ix,Nx,K,R,Ip,poly_Nnon0;
    double *x,*poly_coeff,*poly_tap_in,*poly_branch_in;
    size_t nDim_x,nDim_p1,nDim_p2,nDim_p3; 
    const size_t *size_x,*size_p1,*size_p2,*size_p3; 
       
    // 1st input (input signals [Ix x Nx]) 
    x              = (double*) mxGetDoubles(prhs[0]);
    nDim_x         = mxGetNumberOfDimensions(prhs[0]);
    size_x         = mxGetDimensions(prhs[0]); 
    // 2nd input (number of subbands of the filter bank) 
    K              = (int) mxGetScalar(prhs[1]);
    // 3nd input (resampling factor of the filter bank) 
    R              = (int) mxGetScalar(prhs[2]);
    // 4th input (length of the prototype filter) 
    Ip             = (int) mxGetScalar(prhs[3]);
    // 5th input (number of non-zero elements on the polyphase network) 
    poly_Nnon0     = (int) mxGetScalar(prhs[4]);
    // 6st input (non-zero polyphase coefficients [poly_Nnon0 x 1]) 
    poly_coeff     = (double*) mxGetDoubles(prhs[5]);
    nDim_p1        = mxGetNumberOfDimensions(prhs[5]);
    size_p1        = mxGetDimensions(prhs[5]); 
    // 7st input (tap idx of the non-zero polyphase coefficients [poly_Nnon0 x 1]) 
    poly_tap_in    = (double*) mxGetDoubles(prhs[6]);
    nDim_p2        = mxGetNumberOfDimensions(prhs[6]);
    size_p2        = mxGetDimensions(prhs[6]); 
    // 8st input (branch idx of the non-zero polyphase coefficients [poly_Nnon0 x 1]) 
    poly_branch_in = (double*) mxGetDoubles(prhs[7]);
    nDim_p3        = mxGetNumberOfDimensions(prhs[7]);
    size_p3        = mxGetDimensions(prhs[7]); 
    // Length of the input signal
    Ix             = (int) size_x[0];
    // Number of input signals
    Nx             = (int) size_x[1];

    // Check dimensions
    if(nDim_p1!=2 || nDim_p2!=2 || nDim_p3!=2){
       mexErrMsgTxt("gdft_fb_analysis_mex.c: Incorrect dimensions for poly_coeff, poly_tap, or poly_branch");
    }
    if(size_p1[0] != poly_Nnon0 || size_p2[0] != poly_Nnon0 || size_p3[0] != poly_Nnon0){
       mexErrMsgTxt("gdft_fb_analysis_mex.c: Incorrect dimensions for poly_coeff, poly_tap, or poly_branch");
    }
    if(size_p1[1] != 1 || size_p2[1] != 1 || size_p3[1] != 1){
       mexErrMsgTxt("gdft_fb_analysis_mex.c: Incorrect dimensions for poly_coeff, poly_tap, or poly_branch");
    }

    // Cast input vectors to int
    int poly_tap[poly_Nnon0], poly_branch[poly_Nnon0];  
    for(int i=0;i<poly_Nnon0;i++){
        poly_tap[i]    = (int) poly_tap_in[i];
        poly_branch[i] = (int) poly_branch_in[i];
    } 

    // --------------------------- DEFINE OUTPUT ---------------------------
    // Output (filtered signal [K x Iy x Nx])
    fftw_complex *y;  
    int Iy = (int) ceil((double)(Ix+Ip-1)/(double)R); 
    // Set output dimensions
    size_t nDim_y = 3;
    const size_t size_y[3] = {K,Iy,Nx};
    // Create output array
    plhs[0] = mxCreateNumericArray(nDim_y,
    	                           (const mwSize*) size_y,
    	                           mxDOUBLE_CLASS,mxCOMPLEX);
    y       = (fftw_complex*) mxGetComplexDoubles(plhs[0]);   

    // ------------------------------- FILTER ------------------------------ 
    gdft_fb_analysis(x,
                     Ix,
                     Nx,
                     K,
                     R,
                     Ip,
                     poly_Nnon0,
                     poly_coeff,
                     poly_tap,
                     poly_branch,
                     y);
}

 

                

