// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022                   
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// MEX implementation of the superfast solver for real-valued LS problems
// proposed in:
//
// M. Poletti et al, "A Superfast Toeplitz Matrix Inversion Method for 
// Single- and Multi-Channel Inverse Filters and its Application to Room
// Equalization," in  IEEE/ACM Transactions on Audio Speech and Language 
// Processing, 29, 3144–3157, 2021. 
// https://doi.org/10.1109/TASLP.2021.3120650
// -------------------------------------------------------------------------
#include "superfast_LS_real.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,const mxArray *prhs[])
{
    // ------------------------ CHECK INPUT/OUTPUTS ------------------------
    // Check for proper number of arguments
    if(nrhs!=8) {
        mexErrMsgTxt("superfast_LS_real_mex.c: 8 inputs are required.");
    } 
    
    if(nlhs!=1) {
        mexErrMsgTxt("superfast_LS_real_mex.c: 1 output required.");
    }

    // --------------------------- READ INPUTS -----------------------------
    // Variagles required to read inputs 
    int L,M,Ig,Ih,Id,Nfft,P;
    double beta;
    double *d,*h;
    size_t nDim_h,nDim_d; 
    const size_t *size_h,*size_d; 
       
    // 1st input (RIR [Ih x L x M]) 
    h         = mxGetDoubles(prhs[0]);
    nDim_h    = mxGetNumberOfDimensions(prhs[0]);
    size_h    = mxGetDimensions(prhs[0]); 
    // 2nd input (target response [Id x M]) 
    d         = mxGetDoubles(prhs[1]);
    nDim_d    = mxGetNumberOfDimensions(prhs[1]);
    size_d    = mxGetDimensions(prhs[1]); 
    // 3rd input (number of loudspeakers) 
    L         = (int) mxGetScalar(prhs[2]);
    // 4th input (number of control points) 
    M         = (int) mxGetScalar(prhs[3]);
    // 5th input (filter length) 
    Ig        = (int) mxGetScalar(prhs[4]);
    // 6th input (regularization factor) 
    beta      = mxGetScalar(prhs[5]);
    // 7th input (solver selected) 
    Nfft      = (int) mxGetScalar(prhs[6]);
    // 8th input (approximation order)
    P         = (int) mxGetScalar(prhs[7]);
    // Length of the RIR
    Ih        = (int) size_h[0];
    // Length of the target
    Id        = (int) size_d[0];

    // Check parameters
    if(nDim_h!=3 && nDim_h!=2){
        mexErrMsgTxt("superfast_LS_real_mex.c: Incorrect size for h or d"); 
    }
    if(size_h[1]!=L || size_h[2]!=M || size_d[1]!=M){
        mexErrMsgTxt("superfast_LS_real_mex.c: Incorrect size for h or d"); 
    }
    if(beta<=0.0){
        mexErrMsgTxt("superfast_LS_real_mex.c: Incorrect value for beta"); 
    }
    if(P<0){
        mexErrMsgTxt("superfast_LS_real_mex.c: Incorrect approx. order"); 
    }
 
    // --------------------------- DEFINE OUTPUT ---------------------------
    // Output (estimated filter [Ig x L])
    double *g;  
    // Set output dimensions
    size_t nDim_g = 2;
    const size_t size_g[2] = {Ig*L,1};
    // Create output array
    plhs[0] = mxCreateNumericArray(nDim_g,
    	                           (const mwSize*) size_g,
    	                           mxDOUBLE_CLASS,mxREAL);
    g       =  mxGetDoubles(plhs[0]);     

    // ------------------------------ COMPUTE ------------------------------ 
    // Solve LS problem
    superfast_LS_real(h,d,L,M,Ih,Ig,Id,Nfft,P,beta,g);
}

 

                

