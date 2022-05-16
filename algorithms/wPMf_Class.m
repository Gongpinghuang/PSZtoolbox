classdef wPMf_Class < handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class to compute the filters for a Personal Sound Zones (PSZ) system
    % the Weighted Pressure Matching with Frequency Domain formulation 
    % (wPM-F) algorithm presented in:
    % 
    % J. H. Chang et al, "Experimental validation of sound field control 
    % with a circular double-layer array of loudspeakers," in Journal of 
    % the Acousitca Society of America, 2013, 
    % doi: 10.1121/1.4792486.
    % ---------------------------------------------------------------------
    properties
        % ...................... General properties .......................
        
        L;               % Number of loudspeakers.
        M;               % Number of control points.
        Ig;              % Length of the filters.
        Ih;              % Length of the RIR.
        Id;              % Length of the target response.                 
        beta;            % Regularization factor.
        tag='wPM-F';     % Tag to identify the algorithm.

        % ....................... Impulse responses .......................
        
        h=[];            % RIR, -> [Ih x L x M].
        d=[];            % Target,-> [Id x M x Q].

        % ............................ Flags ..............................  

        issetup = false; % Flag indicating if the object is setup.
    end
    
    methods
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,g,ProcTime] = wPMf_Class(Ig,...
                                               mod_delay,...
                                               LSPref,...
                                               beta_rel,...
                                               h,...
                                               idx_b,...
                                               idx_d)
            % -------------------------------------------------------------
            % Constructor with the required parameters. 
            % --------------------------- Inputs --------------------------
            % - Ig:         Lenght of the filters.
            % - mod_delay:  Modelling delay.
            % - LSPref:     Index of the reference loudspeaker.
            % - beta_rel:   Relative regularization factor.
            % - h:          RIR between all LSP and control points, 
            %               -> [Ih x L x M]
            % - idx_b:      IDX of the control points in the bright zone,
            %               -> [Mb x 1]
            % - idx_b:      IDX of the control points in the dark zone,
            %               -> [Md x 1]            
            % -------------------------- Outputs --------------------------
            % - obj:        Handle to object.   
            % - g:          PSZ filters,-> [Ig x L].
            % - ProcTime:   Processing Time used to compute the filters. 
            % -------------------------------------------------------------

            
            % ----------------------- Set properties ----------------------
            obj.L        = size(h,2);
            obj.M        = size(h,3);
            obj.beta     = beta_rel*mean(mean(sum(h.^2,1)));
            obj.Ih       = size(h,1);
            obj.Ig       = Ig;
            obj.Id       = obj.Ig+obj.Ih-1;
            obj.h        = h;
     
            % -------------------- Generate Target IR ---------------------
            % Initialize target
            obj.d  = zeros(obj.Id,obj.M);
            % Set target IR for the bright zone as the RIR for the ref LSP
            obj.d(mod_delay+1:mod_delay+size(h,1),idx_b) = ...
                                                squeeze(h(:,LSPref,idx_b));
            % Set target IR for the dark zone as a null response
            obj.d(mod_delay+1:mod_delay+size(h,1),idx_d) = 0;
            
            % -------------------------------------------------------------
            % Object is setup
            obj.issetup   = true;
                           
            % ---------------------- Compute filters ----------------------
            tic
            % Compute filters
            g = obj.compFilt();

            ProcTime = toc;
        
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g = compFilt(obj)
            % -------------------------------------------------------------
            % Method to compute the filters using an exact solver, i.e.,
            % using the expressiong gf=(Hf'*Hf+beta*I)^{-1}*Hf'*df. 
            % -------------------------- Outputs --------------------------
            % - g:   Filters, -> [Ig x L].    
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMf_Class/compFilt_exact: ',...
                                 'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Required FFT size 
            Nfft  = ceil((obj.Id)/2)*2;
            % Number of positive freq. bins
            Kfft  = (Nfft/2)+1;
            % FFT of the RIR
            hf    = permute(fft(obj.h,Nfft,1),[3,2,1]);
            % FFT of the target
            df    = permute(fft(obj.d,Nfft,1),[2,1]);
            
            % ----------------------- Compute filters ---------------------     
            % Init
            Rf    = zeros(obj.L,obj.L);
            rf    = zeros(obj.L,1);
            gf    = zeros(Nfft,obj.L);
            
            % ............. Compute optimal freq. coefficients ............
            % For each positive freq. bin
            for f=1:Kfft
               % Compute regularizaed correlation matrix 
               Rf(:) = hf(:,:,f)'*hf(:,:,f)+obj.beta*eye(obj.L);
               % Compute correlation vector
               rf(:) = hf(:,:,f)'*df(:,f);
               % Obtain optimal filters by solving the linear system
               gf(f,:) = linsolve(Rf,rf,struct('SYM',true,'POSDEF',true));                                
            end
            
            % ......... Compute impulse responses of the filters ..........
            % Compute IFFT
            g = ifft(gf,Nfft,1,'symmetric');
            % Truncate using a rectangular window
            g = g(1:obj.Ig,:);
        end

    end
    
    methods(Static)

        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function FLOPs = getSolverFLOPs(Ig,Ih,L,M)
            % -------------------------------------------------------------
            % Method to compute the number of Floating-Point Operations 
            % (FLOPS) required to compute the filters. The FLOPs counts can
            % be found in:
            %
            % V. Moles-Cases, "Filter Optimization for Personal Sound Zones 
            % Systems," PhD thesis, 2022, Available in September 2022 in 
            % UPV repository.
            % --------------------------- Inputs --------------------------
            % - Ig:          Filter length.
            % - Ih:          RIR length.
            % - L:           Number of loudspeakers.
            % - M:           Number of control points.
            % -------------------------- Outputs --------------------------
            % - FLOPs:       Total number of FLOPs.              
            % -------------------------------------------------------------
            
            % FFT size 
            Nfft     = Ih+Ig-1;
            % FLOPs required to compute the RTF
            FLOPs    = (L*M+M+L)*5*Nfft*log2(Nfft);
            % FLOPs required to compute the correlation matrices
            FLOPs    = FLOPs+((2/3)*L^3+2*M*L^2+6*M*L+(11/2)*L^2+(5/6)*L)*Nfft;
        end
        
    end
end

