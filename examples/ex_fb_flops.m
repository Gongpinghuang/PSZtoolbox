% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to compute the number of Floating Point Operations (FLOPs)
% required to compute the filters for wPM-S using different filter bank
% configurations.
% -------------------------------------------------------------------------
% Change path
cd('../');
% Initialize
init();
 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% .......................... General parameters ...........................
% Sampling frequency for the simulation
fs_proc    = 44100/7;

% .......................... Scenario parameters ..........................
 % Path to the RIR
RIRpath    = 'RIR/2021_6_30_13_12_27/RIR.mat';

% .......................... Algorithm parameters .........................
% Algorithm 
algorithm  = 'wPM-S';
% Solver
solver     = 'exact';
% Length of the equivalent broadband filters
Ig         = 2500;
% Number of subbands
K          = [8 ,8 ,8 ,16,16,16,24,24,24];
% Resampling factor
R          = [4 ,6 ,8 ,8 ,10,16,12,16,24];
% Prototype filter length
Ip         = 11:10:211;

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
% Resample RIR
h        = resampleRIR(h,RIR.fs,fs_proc); 
% Store RIR length
Ih       = size(h.ctrl,1);
% Store number of loudspeakers
L        = size(h.ctrl,2);
% Store number of control points
M        = size(h.ctrl,3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     COMPUTE COMPUTATIONAL COMPLEXITY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
FLOPs  = zeros(length(Ip),length(R));
lgnStr = cell(length(R),1);

% For each filter bank configuration...
for n=1:length(R)
    % For each filter length...
    for p=1:length(Ip)
        
        % Display
        dispPSZ(['Ip: ',int2str(Ip(p)),' K: ',int2str(K(n)),...
                                                 ' N: ',int2str(R(n))],2);
        % Setup filter bank
        FB = gdftFB_Class(K(n),R(n),Ip(p),false);
        % Compute costs
        FLOPs(p,n)  = wPMs_Class.getSolverFLOPs(solver,...
                                                Ig,...
                                                Ih,...
                                                L,...
                                                M,...
                                                K(n),...
                                                R(n),...
                                                Ip(n),...
                                                []);
    end  
    
    lgnStr{n} = ['K=',int2str(K(n)),', R=',int2str(R(n))];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine plot colors
cRGB   = repelem(get(gca,'ColorOrder'),3,1);
close all;
% Determine line type
lineS  = repmat({'-','--',':'},1,3);
marker = repelem({'o','<','x'},3);
% Init figure
figure('Color',[1 1 1]);
% Plot for each configuration
for n=1:length(R)
    % Plot
    semilogy(Ip,FLOPs(:,n)/1e9,[lineS{n},marker{n}],'Color',cRGB(n,:),'LineWidth',2);
    hold on;
end
grid;
xlabel('Prototype Filter Length [samples]','Interpreter','Latex','FontSize',16);
ylabel('GFLOPs','Interpreter','Latex','FontSize',16);
xlim([Ip(1),Ip(end)]);
lgn = legend(lgnStr);
set(lgn,'Interpreter','Latex','NumColumns',4,'Location','Best');
xlim([30,200]);
dispPSZ('Done.',5);

