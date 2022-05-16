% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to compute the number of Floating Point Operations (FLOPs)
% required to compute the filters for wPM-S using different configurations
% for the superfast solver.
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
% Sampling frequency for the simulation
fs_proc    = 44100/7;

% .......................... Scenario parameters ..........................
 % Path to the RIR
RIRpath    = 'RIR/2021_6_30_13_12_27/RIR.mat';

% .......................... Algorithm parameters .........................
% Algorithm 
algorithm  = 'wPM-S';
% Length of the equivalent broadband filters
Ig         = 400:100:5000;
% Length of the prototype filter of the filter bank
Ip         = 51;
% Number of subbands of the filter bank
K          = 16;
% Resampling factor of the filter bank
R          = 10;
% The exact solver is used in all subbands
solver     = {'exact','superfast','superfast'};
% Approximation order
P          = {[]     ,2000               ,[2000,100*ones(1,K/2-1)]}; 
% Legend for plots
lgnStr     = {'Exact','Superfast Pk=2000','Superfast Pk-SBD'      };
% Number of different evaluations 
Neval      = numel(solver);

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
% Display
dispPSZ('Computing computational complexity',1);
% Initialize
FLOPs = zeros(length(Ig),Neval);
% For each evaluation...
for i=1:Neval
    % Display
    dispPSZ(['Solver: ',lgnStr{i}],2);
    % For each filter length...
    for p=1:length(Ig)
        % Display
        dispPSZ(['Filter length ',int2str(Ig(p))],3);
        % Compute FLOPs
        FLOPs(p,i) = wPMs_Class.getSolverFLOPs(solver{i},...
                                               Ig(p),...
                                               Ih,...
                                               L,...
                                               M,...
                                               K,...
                                               R,...
                                               Ip,...
                                               P{i});
    end
end
   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init figure
figFLOPs = figure('Color',[1 1 1]);
% Plot
semilogy(ceil(Ig/R),FLOPs/1e9,'LineWidth',2);
grid;
xlim(ceil([Ig(1),Ig(end)]/R));
xlabel('Subband Filter length [samples]','Interpreter','Latex','FontSize',14);
ylabel('GFLOPs','Interpreter','Latex','FontSize',14);
lgn = legend(lgnStr);
set(lgn,'Interpreter','Latex','FontSize',14,'Location','Best');

dispPSZ('Done.',5);
