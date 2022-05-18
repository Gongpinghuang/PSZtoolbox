# PSZ-toolbox
[MATLAB](https://www.mathworks.com/products/matlab.html) toolbox to evaluate the performance of **Personal Sound Zones (PSZ)** systems. The use of the toolbox is exemplified with different examples included in the `examples` folder. In particular, the toolbox allows to evaluate the performance of the filters computed with the different formulations of the weighted Pressure Matching (wPM) algorithm proposed in:

* Chang, J.-H., & Jacobsen, F. (2013). Experimental validation of sound field control with a circular double-layer array of loudspeakers. _The Journal of the Acoustical Society of America_, 133, 2046. https://doi.org/10.1121/1.4792486
* Simon Galvez, M. F., Elliott, S. J., & Cheer, J. (2015). Time Domain Optimization of Filters Used in a Loudspeaker Array for Personal Audio. _IEEE/ACM Transactions on Audio, Speech, and Language Processing_, 23(11), 1869–1878. https://doi.org/10.1109/TASLP.2015.2456428
* Molés-Cases, V., Piñero, G., De Diego, M., & Gonzalez, A. (2020). Personal Sound Zones by Subband Filtering and Time Domain Optimization. _IEEE/ACM Transactions on Audio Speech and Language Processing_, 28, 2684–2696. https://doi.org/10.1109/TASLP.2020.3023628

Also, the toolbox includes an implementation of the superfast solver for Least Squares (LS) problems proposed in:
* Poletti, M. A., & Teal, P. (2021). A Superfast Toeplitz Matrix Inversion Method for Single- and Multi-Channel Inverse Filters and its Application to Room Equalization. _IEEE/ACM Transactions on Audio Speech and Language Processing_, 29, 3144–3157. https://doi.org/10.1109/TASLP.2021.3120650

Moreover, the toolbox includes the polyphase implementation for Generalized Discrete Fourtier Transform (GDFT) filter banks proposed in:

* Harteneck, M., Weiss, S., & Stewart, R. W. (1999). Design of near perfect reconstruction oversampled filter banks for subband adaptive filters. _IEEE Transactions on Circuits and Systems II: Analog and Digital Signal Processing_, 46(8), 1081–1085. https://doi.org/10.1109/82.782056

This toolbox has been used to obtain the results included in:
* Molés-Cases, V., Piñero, G., & Gonzalez, A. (2022). On the Performance of Personal Sound Zones Systems with Subband Filtering. Submitted to  _IEEE/ACM Transactions on Audio Speech and Language Processing_. 

## Folders
The toolbox is organized in the following folders:
* `algorithms`: source files required to compute the filters of the PSZ system.
* `examples`: scripts to exemplify the use of the toolbox.
* `filterbank`: source files related with filter banks.
* `mex`: source files for the MEX functions.
* `rir`: folder where the Room Impulse Responses (RIR) of the system are stored.
* `utilities`: source files with different utilities for the toolbox.

## Room Impulse Responses
The Room Impulse Responses (RIR) of the system must be stored in a `.mat` file within the `rir` folder. The `.mat` file must contain an structure called `RIR` with the following fields:
* `RIR.h_ctrl`: Array of size Ih x L x M with the RIRs between the L loudspeakers and the M control points.
* `RIR.h_val`: Array of size Ih x L x M, with the RIRs between the L loudspeakers and the M validation points.
* `RIR.fs`: Sampling frequency of the RIR.

An example of such an structure is inclunded in `rir/2021_6_30_13_12_27/RIR.mat`.

## MEX functions
All the functionalities of the toolbox are implemented with MATLAB, however, some functionalities are also implemented with C language and ran from MATLAB using MEX functions to speed up the simulations. To use the MEX implementation set `mexFlag = true;`. To build the MEX functions you must have installed one of the [C compilers supported by MATLAB](https://es.mathworks.com/support/requirements/supported-compilers.html). The provided MEX files use the [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/) routines of the [Intel(R) Math Kernel Library](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html), which are provided by MATLAB, and the external [FFTW library](https://www.fftw.org/).

### **Setup FFTW for Microsoft Windows** 
The toolbox includes all the required files of the FFTW v.3.3.5 library for platforms using Microsoft Windows x64. These files have been obtained by downloading a pre-compiled version of the library from http://www.fftw.org/install/windows.html, and by following the steps in the previous link to create the import libraries with the [Library Manager of Microsoft Visual Studio](https://docs.microsoft.com/en-us/cpp/build/reference/lib-reference?redirectedfrom=MSDN&view=msvc-170).

### **Setup FFTW for Mac OS X**
Download and install MacPorts from: https://www.macports.org/install.php. Open the terminal and type `sudo port install fftw-3` to install FFTW. To build the MEX functions, modify the file `mex/buildmex.m` such that 
* Variable `FFTWpath` contains the path to the installed FFTW library (usually '/opt/local/lib/').
* Variable `FFTWheader` contains the path to the header files of the library (usually '/opt/local/include/').
* Variable `FFTWlibfile` contains the name of the `.a` file for the library.

