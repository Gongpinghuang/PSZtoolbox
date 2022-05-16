% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to compute the average metrics over 125-1500 Hz for a PSZ 
% system using the wPM-S algorithm and different filter bank configurations
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
freqRange  = [125,1500];
% Number of frequencies for which the metrics are computed
Nval       = 2^14;
% Frequencies for which the metrics are computed
fval       = (0:(Nval/2))*fs_proc/Nval;
% IDX of the frequencies in the selected range
f_idx      = find(fval>freqRange(1) & fval<freqRange(2));  
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
% The exact solver is used in all subbands
solver     = 'exact';
% Length of the prototype filter of the filter bank
Ip         = 11:10:211;
% Number of subbands of the filter bank
K          = {8 ,8 ,8 ,16,16,16,24,24,24};
% Resampling factor of the filter bank
R          = {4 ,6 ,8 ,8 ,10,16,12,16,24};

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
%   SIMULATE THE PERFORMANCE OF WPM-S FOR THE DIFFERENT FB CONFIGURATIONS                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
dispPSZ('Simulating for wPM-S',1);

% ----------------------- Processing for FB config ------------------------
% Initialize
met_avg.AC_dB   = zeros([length(Ip),numel(R)]);
met_avg.MSE_dB  = zeros([length(Ip),numel(R)]);
lgnStr  = cell(numel(R),1);

% For each fb configuration...
for i=1:numel(R)
    % For each filter length...
    for p=1:length(Ip)
        % Display
        dispPSZ(['Ip: ',int2str(Ip(p)),' K: ',int2str(K{i}),' R: ',int2str(R{i})],2);

        % ....................... Compute filters .........................
        % Setup filter bank
        FB = gdftFB_Class(K{i},R{i},Ip(p),mexFlag); 
        % Compute filters            
        g = computeFilters(algorithm,...
                           Ig,...
                           mod_delay,...
                           LSPref,...
                           beta_rel,...
                           solver,...
                           [],...
                           FB,...
                           mexFlag,...
                           h.ctrl,...
                           idx_b,...
                           idx_d);

        % ........................ Compute metrics ........................
        met = computeMetrics(h.val,...
                                g,...
                                s,...
                                idx_b,...
                                idx_d,...
                                Nval,...
                                fs_proc,...
                                FB,...
                                mod_delay,...
                                LSPref,...
                                false,...
                                []);
        % Compute average metrics 
        met_avg.AC_dB(p,i)   = mean(met.AC_dB(f_idx));
        met_avg.MSE_dB(p,i)  = mean(met.MSE_dB(f_idx));
    end
    lgnStr{i} = ['K=',int2str(K{i}),', N=',int2str(R{i})];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 SIMULATE THE REFERENCE PERFORMANCE OF WPM-T         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
dispPSZ('Simulating for wPM-T',1);

% ----------------------- Processing for FB config ------------------------

% ............................ Compute filters ............................
% Compute filters            
g = computeFilters('wPM-T',...
                   Ig,...
                   mod_delay,...
                   LSPref,...
                   beta_rel,...
                   'exact',...
                   [],...
                   [],...
                   mexFlag,...
                   h.ctrl,...
                   idx_b,...
                   idx_d);

% ............................ Compute metrics ............................
met = computeMetrics(h.val,...
                     g,...
                     s,...
                     idx_b,...
                     idx_d,...
                     Nval,...
                     fs_proc,...
                     [],...
                     mod_delay,...
                     LSPref,...
                     false,...
                     []);
% Compute average metrics 
met_avg_ref.AC_dB  = mean(met.AC_dB(f_idx));
met_avg_ref.MSE_dB = mean(met.MSE_dB(f_idx));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID of the metrics to plot
met2plot = {'AC_dB','MSE_dB'};
% Label for y axis            
Ylabel       = {'AC [dB]','MSE [dB]'};            
% Determine line type
cRGB         = get(gca,'ColorOrder');close all;
cRGB         = repelem(cRGB(1:3,:),3,1);
lineS        = repmat({'-','--',':'},1,3);
marker       = repmat({'o','<','x'},1,3);

% ---------------------------- Plot each metric ---------------------------
% For each metric...
for p = 1:numel(met2plot)
    % Create figure
    figure('Color',[1 1 1]);

    % Plot reference metrics
    plot(Ip,met_avg_ref.(met2plot{p})*ones(size(Ip)),'-',...
             'Color',[0,0,0],'LineWidth',5,'HandleVisibility','off');
    hold on;

    % For each configuration
    for i=1:numel(R)
        % Plot
        plot(Ip,met_avg.(met2plot{p})(:,i),'LineStyle',lineS{i},...
               'Marker',marker{i},'Color',cRGB(i,:),'LineWidth',2);
        hold on;
    end
    xlabel('Prototype Filter Length [samples]','Interpreter','Latex','FontSize',16);
    ylabel(Ylabel{p},'Interpreter','Latex','FontSize',16);
    xlim([Ip(1),Ip(end)]);
    lgn = legend(lgnStr);
    set(lgn,'Interpreter','Latex','NumColumns',3,'Location','Best');
    grid;
end

dispPSZ('Done.',5);


