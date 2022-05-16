% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the required metrics to evaluate the performance of a
% PSZ system.
% -------------------------------- Inputs ---------------------------------
% - h:          RIRs in the validation points, -> [Ih x L x M]. 
% - g:          Broadband/Subband filters for the system,
%                  -> [Ig x L]/[Ig_k x L x K/2]
% - s:          Input signal to the system, -> [Is x 1].
% - idx_b:      IDX of the val. points in the bright zone, -> [Mb x 1]
% - idx_d:      IDX of the val. pointsin  the dark zone, -> [Md x 1]
% - Nval:       Number of frequency bins.
% - fs:         Sampling frequency.
% - FB:         Handle of an object of the class gdftFB_Class.
% - mod_delay:  Modelling delay.
% - LSPref:     Index of the reference loudspeaker.
% - smoothFlag: Flag indicating if freq. averaging is applied
% - frac:       Fraction of the octave band used for the averaging. 
% -------------------------------- Outputs --------------------------------
% - metrics:    Structure containing the computed metrics:
%                 * AC_dB:   Acoustic contrast in dB, -> [Kfft x 1]
%                 * MSE_dB:  Mean Square Error in dB of x with  
%                            respect to d in the bright zone, -> [Kfft x 1]
%                 * AE_dB:   Array Effort in dB, -> [Kfft x 1]
% -------------------------------------------------------------------------
function met = computeMetrics(h,...
                              g,...
                              s,...
                              idx_b,...
                              idx_d,...
                              Nval,...
                              fs,...
                              FB,...
                              mod_delay,...
                              LSPref,...
                              smoothFlag,...
                              frac)
    
    % ---------------------- Compute target response ----------------------
    % Total delay for the target
    if ~isempty(FB)
       delay = mod_delay+FB.delay; 
    else
       delay = mod_delay;
    end
    % Initialize target IR
    d = zeros(size(h,1)+size(g,1)-1,size(h,3));
    % Set target for the bright zone as the RIR for the ref LSP
    d(delay+1:delay+size(h,1),idx_b) = squeeze(h(:,LSPref,idx_b));
    % Set target for the dark zone as a null response
    d(delay+1:delay+size(h,1),idx_d) = 0;
    % Filter input signal with the target using FFTs
    Nfft = size(d,1)+length(s)-1;
    d    = ifft(fft(d,Nfft,1).*fft(s(:),Nfft,1),Nfft,1);
          
    % ------------- Compute signal captured in the microphones ------------
    % Compute the signal fed to the loudspeakers
    if ~isempty(FB)
       y    = FB.filter(s(:),g);         
    else
       Nfft = size(g,1)+length(s)-1;
       y    = ifft(fft(g,Nfft,1).*fft(s(:),Nfft,1),Nfft,1);
    end
    % Compute signal captured in the microphones
    Nfft = size(h,1)+size(y,1)-1;
    x    = squeeze(ifft(sum(fft(y,Nfft,1).*fft(h,Nfft,1),2),Nfft,1));
    
    % -------------------- Compute frequency responses --------------------
    % Check if the number of validation frequency points is correct
    assert(Nval>size(x,1),'Incorrect number of validation freq. points');
     % FFT size
    Nfft = Nval;
    % Number of frequency bins in the positive spectrum
    Kfft = (Nfft/2)+1;
    % Compute FFTs
    H    = fft(h,Nfft,1);
    X    = fft(x,Nfft,1);
    D    = fft(d,Nfft,1);
    Y    = fft(y,Nfft,1);

    %-------------------------- Compute Metrics ---------------------------
    % Compute Acoustic Contrast
    met.AC_dB  = 10*log10(mean(abs(X(1:Kfft,idx_b)).^2,2)./...
                          mean(abs(X(1:Kfft,idx_d)).^2,2));
    % Compute MSE in the bright zone
    met.MSE_dB = 10*log10(mean(abs(X(1:Kfft,idx_b)-D(1:Kfft,idx_b)).^2,2));
    % Compute reference energy
    Er         = sum(abs(X(1:Kfft,idx_b)).^2,2)./...
                 sum(abs(H(1:Kfft,LSPref,idx_b)).^2,3);
    % Compute Array Effort
    met.AE_dB  = 10*log10(sum(abs(Y(1:Kfft,:)).^2,2)./Er);
 
    %--------------------------- Smooth Metrics ---------------------------
    if smoothFlag
       % List of all the metrics in the structure 
       name = fieldnames(met);
       % Smooth all metrics
       for i=1:numel(name)
           met.(name{i}) = octaveBandAveraging(met.(name{i}),...
                                               frac,...
                                               Nfft,...
                                               fs);
       end
    end
    
    % Save validation frequencies
    met.fval = (0:Kfft-1)*fs/Nfft;
end

