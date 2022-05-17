% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to compute the average metrics of wPM-S in the different 
% subbands when the superfast solver is used to compute the subband filters
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
% Length of the prototype filter of the filter bank
Ip         = 51;
% Number of subbands of the filter bank
K          = 16;
% Resampling factor of the filter bank
R          = 10;
% Regularization level relative to the mean energy of the RIR
beta_rel   = 0.03;
% The exact solver is used in all subbands
solver     = {'exact'    ,'superfast'};
% Approximation order for the superfast solver
P          = round(logspace(0,4,20));
% Number of solvers to evaluate
Nsolver        = numel(solver);
% Number of evaluation per solver
Neval_x_solver = numel(P);

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
dispPSZ('Simulating for wPM-S',1);

% --------------------------- Compute for wPM-S ---------------------------
% Initialize metrics
met_av_x_sb = cell(Nsolver,1);
% Initialize NMSE
NMSEfilt_dB = zeros(Neval_x_solver,K/2);
% Initialize filters
g           = cell(Nsolver,1);
% Setup filter bank
FB          = gdftFB_Class(K,R,Ip,mexFlag); 
% Min and maximum frequency for each subband in the filter bank
fmin        = (0:(K/2)-1)*fs_proc/K;
fmax        = fmin+(fs_proc/K);
fmin(1)     = 125;
    
% For each solver...
for i=1:Nsolver
    
    met_av_x_sb{i}.AC_dB  = zeros(length(Neval_x_solver),K/2);
    met_av_x_sb{i}.MSE_dB = zeros(length(Neval_x_solver),K/2);
    
    % Display
    dispPSZ(['Solver: ',solver{i}],2);
    
    % For each evaluation...
    for p=1:Neval_x_solver
         if ~(p>1 && i==1)
             dispPSZ(['Order ',int2str(P(p))],3);
             
             % ..................... Compute filters ......................
             g{i} = computeFilters(algorithm,...
                                   Ig,...
                                   mod_delay,...
                                   LSPref,...
                                   beta_rel,...
                                   solver{i},...
                                   P(p),...
                                   FB,...
                                   mexFlag,...
                                   h.ctrl,...
                                   idx_b,...
                                   idx_d);
                            
              % ..................... Compute metrics .....................
              met_l = computeMetrics(h.val,...
                                     g{i},...
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
               % For each subband
               for k=1:K/2
                    % Compute average AC and MSE in the subband
                    f_idx  = find(fval>=fmin(k) & fval<=fmax(k));
                    met_av_x_sb{i}.AC_dB(p,k)  = mean(met_l.AC_dB(f_idx));
                    met_av_x_sb{i}.MSE_dB(p,k) = mean(met_l.MSE_dB(f_idx));
                    % Compute NMSE with respect to the optimal filters
                    NMSEfilt_dB(p,k) = ...
                    20*log10(norm(reshape(g{i}(:,:,k)-g{1}(:,:,k),[],1))/...
                             norm(reshape(g{1}(:,:,k),[],1)));
               end
         else
                met_av_x_sb{i}.AC_dB(p,:)  = met_av_x_sb{i}.AC_dB(1,:);
                met_av_x_sb{i}.MSE_dB(p,:) = met_av_x_sb{i}.MSE_dB(1,:);
         end
    end
end

% ---------------- Compute reference performance for wPM-T ----------------
% Display
dispPSZ('Simulating for wPM-T',1);

% ............................ Compute filters ............................
gt = computeFilters('wPM-T',...
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
met_l = computeMetrics(h.val,...
                      gt,...
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
% Init                  
met_av_x_sb_ref.AC_dB  = zeros(K/2,1); 
met_av_x_sb_ref.MSE_dB = zeros(K/2,1);   
% For each subband
for k=1:K/2
    % Compute average AC and MSE in the subband
    f_idx  = find(fval>=fmin(k) & fval<=fmax(k));
    met_av_x_sb_ref.AC_dB(k)  = mean(met_l.AC_dB(f_idx));
    met_av_x_sb_ref.MSE_dB(k) = mean(met_l.MSE_dB(f_idx));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PLOT ACCURACY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of subbands to plot
Nk = K/2;
% Legend
for k=1:Nk
   lgnStr{k} = ['$k=',int2str(k),'$']; 
end
% Line Style
lineS      = {'-','-','-.','-.','--','--',':',':'};
marker     = {'x','o','+','diamond','square','*','<','>'};
markerSize = [8,5,7,5,6,6,5,5];
cRGB       = copper(Nk);

figure('Color',[1 1 1]);
for k=1:Nk
    semilogx(P,NMSEfilt_dB(:,k),'LineWidth',2,'Color',cRGB(k,:),...
             'LineStyle',lineS{k},'Marker',marker{k},...
             'MarkerSize',markerSize(k));
    hold on;
end
xlabel('Approximation Order Pk','Interpreter','Latex','FontSize',16);
ylabel('NMSE [dB]','Interpreter','Latex','FontSize',16);
lgn = legend(lgnStr);
set(lgn,'Interpreter','Latex');
grid;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PLOT PSZ PERFORMANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID of the metrics to plot
met2plot = {'AC_dB','MSE_dB'};
% Label for y axis            
Ylabel   = {'AC [dB]','MSE [dB]'};  
% Number of subbands to plot
Nk = 4;
% Legend
lgnStr = cell(Nsolver*Nk,1);
for k=1:Nk
   for i=1:Nsolver
      lgnStr{(k-1)*Nsolver+i} = [solver{i},' - $k=',int2str(k),'$']; 
   end
end
% Plot style
cRGB    = lines(8);
lineS   = {'-',':'};
marker  = {'none','x';'none','o';'none','+'; 'none','<'};       
markerS = [8,8;6,6;8,8;6,6];
lineW   = {3,2};

% For each metric...
for p = 1:numel(met2plot)
    % Initialize figure
    figure('Color',[1 1 1]);
    % For each subband...
    for k=1:Nk
        % For each solver...
        for i=1:Nsolver
            % Plot
            semilogx(P,met_av_x_sb{i}.(met2plot{p})(:,k),lineS{i},...
                     'LineWidth',lineW{i},'Color',cRGB(k,:),...
                     'LineStyle',lineS{i},'Marker',marker{k,i},...
                     'MarkerSize',markerS(k,i));
            hold on;
        end
    end
    xlabel('Aproximation Order Pk','Interpreter','Latex','FontSize',16);
    ylabel(Ylabel{p},'Interpreter','Latex','FontSize',16);
    lgn = legend(lgnStr);
    set(lgn,'Interpreter','Latex','Location','Best');
    grid;
end

dispPSZ('Done.',5);