% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used to compute the Aliast-To-Signal Ratio (ASR) and the 
% Reconstruction Error (RE) for different configurations of a GDFT filter
% bank. 
% -------------------------------------------------------------------------
% Change path
cd('../');
% Initialize
init();

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of subbands
K             = [8 ,8 ,8 ,16,16,16,24,24,24];
% Resampling factor
R             = [4 ,6 ,8 ,8 ,10,16,12,16,24];
% Prototype filter length
Ip            = 11:10:211;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
ASR_dB   = zeros(length(Ip),length(R));
RE_dB    = zeros(length(Ip),length(R));
lgnStr   = cell(length(R),1);

% For each filter bank configuration...
for n=1:length(R)
    % For each filter length...
    for p=1:length(Ip)
        
        % Display
        dispPSZ(['Ip: ',int2str(Ip(p)),...
                 ' K: ',int2str(K(n)),...
                 ' R: ',int2str(R(n))],2);
        % Setup filter bank
        FB = gdftFB_Class(K(n),R(n),Ip(p),false);
        % Store metrics
        ASR_dB(p,n) = FB.ASR_dB;
        RE_dB(p,n)  = FB.RE_dB;
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

% -------------------------------- Plot RE --------------------------------
figRE = figure('Color',[1 1 1]);
for n=1:length(R)
    % Plot
    plot(Ip,RE_dB(:,n),[lineS{n},marker{n}],'Color',cRGB(n,:),'LineWidth',2);
    hold on;
end
grid;
xlabel('Prototype Filter Length [samples]','Interpreter','Latex','FontSize',16);
ylabel('RE [dB]','Interpreter','Latex','FontSize',16);
xlim([Ip(1),Ip(end)]);
lgn = legend(lgnStr);
set(lgn,'Interpreter','Latex','NumColumns',4,'Location','Best');
xlim([30,200]);
ylim([-100,0]);

% ------------------------------- Plot Alias ------------------------------
figASR = figure('Color',[1 1 1]);
for n=1:length(R)
    % Plot
    plot(Ip,ASR_dB(:,n),[lineS{n},marker{n}],'Color',cRGB(n,:),'LineWidth',2);
    hold on;
end
grid;
xlabel('Prototype Filter Length [samples]','Interpreter','Latex','FontSize',16);
ylabel('ASR [dB]','Interpreter','Latex','FontSize',16);
xlim([Ip(1),Ip(end)]);
lgn = legend(lgnStr);
set(lgn,'Interpreter','Latex','NumColumns',4,'Location','Best');
xlim([30,200]);
ylim([-100,0]);

dispPSZ('Done.',5);
