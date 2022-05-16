% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to perform Fractional-Octave Smoothing of transfer functions, 
% based on:
% 
% Hatziantoniou, P. D. et al, "Generalized Fractional-Octave Smoothing of
% Audio and Acoustic Responses," J. Audio Eng. Soc., 2000. 
% http://www.aes.org/e-lib/browse.cfm?elib=12070
% -------------------------------- Inputs ---------------------------------
% - X:    Transfer function to smooth, -> either [Nfft x P] or  [Kfft x P],
%         where P is the number of transfer functions to smooth.
% - frac: Fraction of octave.
% - Nfft: Total number of frequency bins
% - fs:   Sampling frequency.
% -------------------------------- Outputs --------------------------------
% - Y:    Averaged transfer function -> either [Nfft x P] or [Kfft x P]
% -------------------------------------------------------------------------
function Y = octaveBandAveraging(X,frac,Nfft,fs)
    
    % ----------------------------- Check ---------------------------------
    if size(X,1) == Nfft
        Kmax = Nfft;
    elseif size(X,1) == (Nfft/2)+1
        Kmax = (Nfft/2)+1;
        X    = [X;flip(conj(X(2:end-1,:)))];
    else
        error('octaveBandAveraging: Incorrect dimensions for X');
    end

    % ----------- Compute the window required for each freq. bin ----------
    % All frequency bins
    k       = 0:Nfft-1;
    % Replace negative frequency bins with complex conjugate property
    idx     = find(k>Nfft/2);
    k(idx)  = Nfft-k(idx);
    % Obtain the frequency value for each bin
    f       = k*fs/Nfft;
    % Determine the maximum and minimum frequency of a frac-octave window
    fU      = f*2^(0.5*frac);
    fL      = f*0.5^(0.5*frac);
    % Determine the number of samples for each window
    m       = floor(0.5*(fU-fL)*Nfft/fs)+1; 
    % Total number of different windows
    M       = max(m);
    % Generate all the required windows
    Non0idx = cell(M,1);
    for mi = sort(unique(m))
        % Index where the Window is non-zero
        Non0idx{mi} = [1:mi+1,Nfft-mi+1:Nfft];
    end
    
    % ---------------------------- Averaging ------------------------------
    % Initialize the smoothed magnitude and phase components
    Y    = zeros(Kmax,size(X,2));
    idx  = 1:Nfft;
    % For each bin...
    for k=0:Kmax-1
        % Index of the frequency bins that we average for the current bin
        idxC     = idx(Non0idx{m(k+1)});
        Y(k+1,:) = mean(X(idxC,:),1);
        % Rotate the frequency bins to match the window positions
        idx      = circshift(idx,-1);
    end
end



