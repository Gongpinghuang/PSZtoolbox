% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to compute the metrics of wPM-S when the superfast solver 
% with subband-dependent approximation orders is used to compute the 
% subband filters
% -------------------------------------------------------------------------
% Change path
cd('../');
% Initialize
init();

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% .......................... General parameters ...........................
% Flag that indicates if mex functions are used
mexFlag    = false;
% Sampling frequency for the simulation
fs_proc    = 44100/7;
% Frequency range to display
freqRange  = [125,fs_proc/2];   
% Number of frequencies for which the metrics are computed
Nval       = 2^14;
% Frequencies for which the metrics are computed
fval       = (0:(Nval/2))*fs_proc/Nval;
% Flag inficating if the freq. metrics are averaged
smoothFlag = true;
% Fraction of octave band used for the averaging
frac       = 1/5;
% Input signal to the PSZ system (white noise of length 1 seconds)
rng('default');
s          = randn(fs_proc,1);
    
% .......................... Scenario parameters ..........................
 % Path to the RIR
RIRpath    = 'RIR/2021_6_30_13_12_27/RIR.mat';
% Index of the control/validation points in the bright zone
idx_b      = 1:16;
% Index of the control/validation points in the dark zone
idx_d      = 17:32;

% .......................... Algorithm parameters .........................
% Algorithm 
algorithm  = 'wPM-S';
% Length of the equivalent broadband filters
Ig         = 2500;
% Modelling delay
mod_delay  = 100;
% Index of the reference loudspealer
LSPref     = 4;
% Regularization level relative to the mean energy of the RIR
beta_rel   = 0.03; 
% Length of the prototype filter of the filter bank
Ip         = 51;
% Number of subbands of the filter bank
K          = 16;
% Resampling factor of the filter bank
R          = 10;
% Solver used to compute the filters
solver     = {'exact','superfast'       ,'superfast'        ,'superfast'             };
% Approximation order for the superfast solver
P          = {[]     ,100               ,2000               ,[2000,100*ones(1,K/2-1)]};
% Legend for plots
lgnStr     = {'Exact','Superfast Pk=100','Superfast Pk=2000','Superfast Pk-SBD'};
% Number of different evaluations 
Neval      = numel(solver);
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          BUILD MEX FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build mex functions (if required)
if mexFlag
    buildmex();
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                LOAD RIR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load RIR        
load([RIRpath]);
% Store RIR in control points
h.ctrl   = RIR.h_ctrl;
% Store RIR in validation points
h.val    = RIR.h_val;
% Check that the control and validation RIRs have the same dimensions
assert(all(size(h.ctrl)==size(h.val)),'Incorrect RIR sizes');
% Check that the idx of the control/validation points in the bright and
% dark zone are correct
assert(all([idx_b(:);idx_d(:)]<=size(h.ctrl,3)),...
       'Incorrect indices for the bright and dark zone'); 
% Resample RIR
h        = resampleRIR(h,RIR.fs,fs_proc); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                SIMULATE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
dispPSZ('Simulating',1);
% Initialize metrics
met = cell(Neval,1);

% For each evaluation...
for i=1:Neval
    % Display
    dispPSZ(['Solver: ',lgnStr{i}],2);

    % .......................... Compute filters ..........................
    % Setup filter bank
    FB = gdftFB_Class(K,R,Ip,mexFlag); 
    % Compute filters            
    g = computeFilters(algorithm,...
                       Ig,...
                       mod_delay,...
                       LSPref,...
                       beta_rel,...
                       solver{i},...
                       P{i},...
                       FB,...
                       mexFlag,...
                       h.ctrl,...
                       idx_b,...
                       idx_d);
    
    % .......................... Compute metrics ..........................
    met{i} = computeMetrics(h.val,...
                            g,...
                            s,...
                            idx_b,...
                            idx_d,...
                            Nval,...
                            fs_proc,...
                            FB,...
                            mod_delay,...
                            LSPref,...
                            smoothFlag,...
                            frac);                        
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID of the metrics to plot
met2plot = {'AC_dB','MSE_dB', 'AE_dB'};
% Label for y axis            
Ylabel   = {'AC [dB]','MSE [dB]','AE [dB]'};            
% Line type
lineS    = {'-','-.',':','--'};

% ---------------------------- Plot each metric ---------------------------
% For each metric...
for p = 1:numel(met2plot)
    % Initialize figure
    figure('Color',[1 1 1]);
    % For each evaluation...
    for i=1:Neval
        % Plot
        semilogx(fval,met{i}.(met2plot{p}),lineS{i},'LineWidth',2);
        hold on;
    end
    xlabel('Frequency [Hz]','Interpreter','Latex','FontSize',16);
    ylabel(Ylabel{p},'Interpreter','Latex','FontSize',16);
    xlim(freqRange);
    lgn = legend(lgnStr);
    set(lgn,'Interpreter','Latex','Location','Best');
    grid;
end

dispPSZ('Done.',5);




