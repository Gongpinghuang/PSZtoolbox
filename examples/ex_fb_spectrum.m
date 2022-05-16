% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to plot the spectrum of the analysis filter a two different
% GDFT filter bank configurations
% -------------------------------------------------------------------------
% Change path
cd('../');
% Initialize
init();
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sampling frequency
fs_proc       = 44100/7; 
% Number of subbands
K             = [16,16];
% Resampling factor
R             = [10,16];
% Prototype filter length
Ip            = 51;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT size
Nfft = 10*Ip;
Kfft = (Nfft/2)+1;
f    = (0:Kfft-1)*fs_proc/Nfft;
% Initialize
Ha   = cell(length(R),1);

% For each fb configuration...
for i=1:length(R)
    % Load filter bank
    FB = gdftFB_Class(K(i),R(i),Ip,false); 
    % Compute spectrum
    Ha{i} = fft(FB.ha/sqrt(R(i)),Nfft,1);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each configuration...
for i=1:length(R)
    % Init figure
    figure('Color',[1 1 1]);
    % For each positive subband
    for k=1:K(i)/2
        plot(f,20*log10(abs(Ha{i}(1:Kfft,k))),'LineWidth',1.5);
        hold on;
    end
    grid;
    xlabel('Frequency [Hz]','Interpreter','Latex','FontSize',16);
    ylabel('Magnitude [dB]','Interpreter','Latex','FontSize',16);
    xlim([0,fs_proc/2]);
    ylim([-70,2]);
end

dispPSZ('Done.',5);
