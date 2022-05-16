classdef wPMs_Class < handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class to compute the subband filters for a Personal Sound Zones (PSZ)
    % system using the Weighted Pressure Matching with Subband Domain
    % formulation (wPM-S) algorithm presented in:
    % 
    % V. Moles-Cases et al, "Personal Sound Zones by Subband Filtering and 
    % Time Domain Optimization," in IEEE/ACM Transactions on Audio, Speech, 
    % and Language Processing, 2020, doi: 10.1109/TASLP.2020.3023628.
    % ---------------------------------------------------------------------
    
    properties
        
        % ...................... General properties .......................

        L;               % Number of loudspeakers.
        M;               % Number of control points.
        Ig_k;            % Length of the subband filters.
        Ih_k;            % Length of the subband components of the RIR.
        Id_k;            % Length of the subband components of the target.
        solver;          % Solver used to compute the filters, either 
                         % 'exact' or 'superfast'.
        P=0;             % Approximation order for the superfast solver.
        beta;            % Regularization factor.
        FB;              % Object of the type gdftFB_Class with the FB 
                         % methods and attributes.        
        tag='wPM-S';     % Tag to identify the algorithm.
        
        % ................... Subband Impulse responses ...................
        
        h_k=[];          % Subband components of the RIR, 
                         % -> [Ih_k x L x M x K/2].
        d_k=[];          % Subband components of the target, 
                         % -> [Id_k x M x K/2].

        % ............................ Flags .............................. 
        
        mexFlag = false; % Flag indicating if mex functions are used.
        issetup = false; % Flag indicating if the object is setup.
    end
    
    methods
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,g_k,ProcTime] = wPMs_Class(Ig,...
                                                 mod_delay,...
                                                 LSPref,...
                                                 beta_rel,...
                                                 solver,...
                                                 P,...
                                                 FB,...
                                                 mexFlag,...
                                                 h,...
                                                 idx_b,...
                                                 idx_d)
            % -------------------------------------------------------------
            % Constructor with the required parameters. 
            % --------------------------- Inputs --------------------------
            % - Ig:         Lenght of the equivalent broadband filters.
            % - mod_delay:  Modelling delay.
            % - LSPref:     Index of the reference loudspeaker.
            % - beta_rel:   Relative regularization factor.
            % - solver:     Solver used for computing the filters.
            % - P:          Approximation order for the superfast solver.
            % - FB:         Handle of an object of the class gdftFB_Class.
            % - mexFlag:    Flag indicating if mex functions are used.
            % - h:          RIR between all LSP and control points, 
            %               -> [Ih x L x M]
            % - idx_b:      IDX of the control points in the bright zone,
            %               -> [Mb x 1]
            % - idx_b:      IDX of the control points in the dark zone,
            %               -> [Md x 1]            
            % -------------------------- Outputs --------------------------
            % - obj:        Handle to the object.   
            % - g_k:        Subband PSZ filters,-> [Ig_k x L x K/2].
            % - ProcTime:   Processing Time used to compute the filters. 
            % -------------------------------------------------------------

            % ----------------------- Set properties ----------------------
            obj.L        = size(h,2);
            obj.M        = size(h,3);
            obj.solver   = solver;  
            obj.mexFlag  = mexFlag;
            obj.FB       = FB; 
            obj.beta     = beta_rel*mean(mean(sum(h.^2,1)));
            obj.Ih_k     = obj.FB.getSBlength(size(h,1));
            obj.Ig_k     = ceil(Ig/obj.FB.R);
            obj.Id_k     = obj.Ig_k+obj.Ih_k-1;
            if strcmp(obj.solver,'superfast')
                assert(length(P)==1 | length(P)==obj.FB.K/2,...
                           'wPMs_Class: The size of P must be 1 or K/2');
                % Store an approximation order for each subband       
                obj.P    = repmat(P(:),(obj.FB.K/2)/length(P),1);
            end
     
            % ---------------- Generate Broadband Target IR ---------------
            % Length required for the broadband target such that a length
            % Id_k is obtained for the subband target
            Id = (obj.Id_k+ceil(obj.FB.Ip/obj.FB.R)-2)*obj.FB.R-obj.FB.Ip+2;
            % Initialize broadband target
            d  = zeros(Id,obj.M);
            % Set target IR for the bright zone as the RIR for the ref LSP
            d(mod_delay+1:mod_delay+size(h,1),idx_b) = ...
                                                squeeze(h(:,LSPref,idx_b));
            % Set target IR for the dark zone as a null response
            d(mod_delay+1:mod_delay+size(h,1),idx_d) = 0;

            tic
            % ------------------- Subband decomposition -------------------
            % Subband decomposition of the RIR
            obj.h_k = reshape(obj.FB.subband_dec(h(:,:)),...
                                                 obj.Ih_k,obj.L,obj.M,[]);
            % Subband decomposition of the target response
            obj.d_k = reshape(obj.FB.subband_dec(d(:,:)),...
                                                 obj.Id_k,obj.M,[]);
            
            % -------------------------------------------------------------
            % Object is setup
            obj.issetup   = true;
                           
            % ---------------------- Compute filters ----------------------
            % Compute filters according to solver
            if strcmp(obj.solver,'exact')
                if obj.mexFlag
                    g_k = obj.compFilt_exact_mex();
                else
                    g_k = obj.compFilt_exact_matlab();
                end
            elseif strcmp(obj.solver,'superfast')
                if obj.mexFlag
                    g_k = obj.compFilt_superfast_mex();
                else
                    g_k = obj.compFilt_superfast_matlab();
                end
            else
                error('wPMs_Class: Not supported solver');
            end
          
            ProcTime = toc;
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g_k = compFilt_exact_matlab(obj)
            % -------------------------------------------------------------
            % Matlab implementation to compute the filters using an exact 
            % solver, i.e., using g_k=(H_k'*H_k+beta_k*I)^{-1}*H_k'*d_k.
            % The exact method is based on using the Choleksy decomposition
            % and uses the FFT to compute H_k'*H_k, and H_k'*d_k. A
            % detalied description of this method can be found in:
            %
            % V. Moles-Cases, "Filter Optimization for Personal Sound Zones 
            % Systems," PhD thesis, 2022, Available in September 2022 in 
            % UPV repository.
            % -------------------------- Outputs --------------------------
            % - g_k:  Subband filters, -> [Ig_k x L x K/2].    
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMs_Class/compFilt_exact: ',...
                                 'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Initialize filters
            g_k   = zeros(obj.Ig_k,obj.L,obj.FB.K/2);
            % Required FFT size to compute H_k'*H_k and H_k'*d_k 
            Nfft  = obj.Ih_k+obj.Ig_k-1;
            % FFT of the RIR
            Hf_k  = fft(obj.h_k,Nfft,1);
            % FFT of the target
            Df_k  = fft(obj.d_k,Nfft,1);
            % Idx of the Ig_k positive lags of the correlation
            lag_p = 1:obj.Ig_k;
            % Idx of the Ig_k-1 negative lags of the correlation
            lag_n = 2:obj.Ig_k;
            
            % ----------------------- Compute filters ---------------------        
            % For each active subband...
            for k = 1:obj.FB.K/2

                % ........ Compute correlation matrix R_k=H_k'*H_k ........
                % Compute correlation of h_k using the FFT
                R_k = ifft(sum(conj(Hf_k(:,:,:,k)).*...
                               permute(Hf_k(:,:,:,k),[1,4,3,2]),3),Nfft,1);
                % Keep lag_p and lag_n 
                R_k = cat(3,permute(R_k(lag_p,:,:,:),[2,4,1,3]),...
                            permute(conj(R_k(lag_n,:,:,:)),[4,2,1,3]));
                % Transform to a cell with 2*Ig_k-1 elements of size LxL       
                R_k = num2cell(R_k,[1,2]);
                % Form the block-toeplitz correlation matrix of h_k
                R_k = cell2mat(R_k(toeplitz(lag_p,[1,lag_n+obj.Ig_k-1])));

                % ........ Compute correlation vector r_k=H_k'*d_k ........
                % Compute cross-correlation using the FFT
                r_k = ifft(sum(conj(Hf_k(:,:,:,k)).*...
                               reshape(Df_k(:,:,k),Nfft,1,[]),3),Nfft,1);
                % Keep lag_p           
                r_k = reshape(r_k(lag_p,:).',[],1);

                % ................ Compute optimal filters ................
                % Obtain optimal filters by solving the linear system
                g_aux = linsolve(R_k+obj.beta*eye(obj.L*obj.Ig_k),...
                                 r_k,...
                                 struct('SYM',true,'POSDEF',true));

                % ......................... Store .........................     
                g_k(:,:,k) = reshape(g_aux,[],obj.Ig_k).';
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g_k = compFilt_exact_mex(obj)
            % -------------------------------------------------------------
            % MEX implementation to compute the filters using an exact 
            % solver, i.e., using g_k=(H_k'*H_k+beta_k*I)^{-1}*H_k'*d_k.
            % The exact method is based on using the Choleksy decomposition
            % and uses the FFT to compute H_k'*H_k, and H_k'*d_k. A
            % detalied description of this method can be found in:
            %
            % V. Moles-Cases, "Filter Optimization for Personal Sound Zones 
            % Systems," PhD thesis, 2022, Available in September 2022 in 
            % UPV repository.
            % -------------------------- Outputs --------------------------
            % - g_k:  Subband filters, -> [Ig_k x L x K/2].    
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMs_Class/compFilt_exact: ',...
                                 'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Initialize filters
            g_k   = zeros(obj.Ig_k,obj.L,obj.FB.K/2);
            
            % ----------------------- Compute filters ---------------------        
            % For each active subband...
            for k = 1:obj.FB.K/2

               % Compute using mex function
                 gaux = exact_LS_complex_mex(obj.h_k(:,:,:,k),...
                                             obj.d_k(:,:,k),...
                                             obj.L,...
                                             obj.M,...
                                             obj.Ig_k,...
                                             obj.beta);             
                % Store
                g_k(:,:,k) = reshape(gaux,obj.Ig_k,[]);    
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g_k = compFilt_superfast_matlab(obj)
            % -------------------------------------------------------------
            % Matlab implementation to compute an approximation of the 
            % subband filters using the superfast solver proposed in:
            %
            % M. Poletti et al, "A Superfast Toeplitz Matrix Inversion
            % Method for Single- and Multi-Channel Inverse Filters and its
            % Application to Room Equalization," in  IEEE/ACM Transactions
            % on Audio Speech and Language Processing, 29, 3144�3157, 2021.
            % https://doi.org/10.1109/TASLP.2021.3120650
            % -------------------------- Outputs --------------------------
            % - g_k:   Subband filters, -> [Ig_k x L x K/2].     
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMs_Class/compFilt_superfast_matlab: ',...
                                'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Initialize filters
            g_k          = zeros(obj.Ig_k,obj.L,obj.FB.K/2);
            % Required FFT size for the superfast solver
            Nfft         = obj.findFFTsize(obj.Id_k,11);
            % FFT of the RIR
            hf_k         = permute(fft(obj.h_k,Nfft,1),[3,2,1,4]);
            % FFT of the target
            df_k         = permute(fft(obj.d_k,Nfft,1),[2,1,3]);
            % Init auxilairy matrices
            Hf_ps        = zeros(obj.L,obj.M);
            Hf_ps_x_Hf   = zeros(Nfft,obj.L,obj.L);
            gf           = zeros(Nfft,obj.L);
            
            % ----------------------- Compute filters ---------------------        
            % For each active subband...
            for k = 1:obj.FB.K/2

                % .... Compute filters using freq. domain formulation ....
                % For each freq. bin...
                for f=1:Nfft
                    % Compute Hf^pseudo = Hf^H*(Hf*Hf^H+beta*I)^{-1} 
                    Hf_ps(:) = ...
                        linsolve(hf_k(:,:,f,k)*hf_k(:,:,f,k)'+...
                                 obj.beta*eye(obj.M),...
                                 hf_k(:,:,f,k),...
                                 struct('SYM',true,'POSDEF',true))';
                    % Compute Hf^pseudo*Hf
                    Hf_ps_x_Hf(f,:,:) = ...
                              reshape(Hf_ps*hf_k(:,:,f,k),1,obj.L,obj.L);
                    % Compute optimal freq. filter coefficients    
                    gf(f,:)  = Hf_ps*df_k(:,f,k);
                end
                % Compute IFFT
                g_fd_full = ifft(gf,Nfft,1);
                % Truncate filters
                g_fd      = g_fd_full(1:obj.Ig_k,:);

                % ........... Compute and apply correction terms ..........
                % Init time and freq. domain residuals
                rt = g_fd_full;
                rf = zeros(Nfft,obj.L);
                % Iterate until reach the selected approximation...
                for p=0:obj.P(k)
                    % Truncate residual
                    rt(1:obj.Ig_k,:) = 0;
                    % Compute its FFT
                    rf(:)            = fft(rt,Nfft,1);
                    % Correct freq. domain residuals
                    rf(:)            = ...
                                sum(Hf_ps_x_Hf.*reshape(rf,[],1,obj.L),3);
                    % Compute time-domain corrected residuals
                    rt(:)            = ifft(rf,Nfft,1);
                    % Correct filters
                    g_fd(:)          = g_fd+rt(1:obj.Ig_k,:);
                end
                % Store
                g_k(:,:,k) = g_fd;
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g_k = compFilt_superfast_mex(obj)
            % -------------------------------------------------------------
            % MEX implementation to compute an approximation of the 
            % subband filters using the superfast solver proposed in:
            %
            % M. Poletti et al, "A Superfast Toeplitz Matrix Inversion
            % Method for Single- and Multi-Channel Inverse Filters and its
            % Application to Room Equalization," in  IEEE/ACM Transactions
            % on Audio Speech and Language Processing, 29, 3144�3157, 2021.
            % https://doi.org/10.1109/TASLP.2021.3120650
            % -------------------------- Outputs --------------------------
            % - g_k:   Subband filters, -> [Ig_k x L x K/2].  
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMs_Class/compFilt_superfast_mex: ',...
                                 'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Initialize filters
            g_k  = zeros(obj.Ig_k,obj.L,obj.FB.K/2);
            % Required FFT size to compute H'*H and H'*d using the FFT
            Nfft = obj.findFFTsize(obj.Id_k,11);

            % ----------------------- Compute filters ---------------------        
            % For each active subband...
            for k = 1:obj.FB.K/2
                 % Compute using mex function
                 gaux = superfast_LS_complex_mex(obj.h_k(:,:,:,k),...
                                                 obj.d_k(:,:,k),...
                                                 obj.L,...
                                                 obj.M,...
                                                 obj.Ig_k,...
                                                 obj.beta,...
                                                 Nfft,...
                                                 obj.P(k));             
                % Store
                g_k(:,:,k) = reshape(gaux,obj.Ig_k,[]);    
            end
        end

    end
    
    methods(Static)

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Nfft = findFFTsize(Nfft_in,prime_max)
            % -------------------------------------------------------------
            % Method to compute the FFT size that fulfills the following
            % properties: 1) Nfft>=Nfft_in; 2) The prime factors of Nfft 
            % are equal or smaller than fac_max.
            % --------------------------- Inputs --------------------------
            % - Nfft_in:   Minimum FFT size.
            % - prime_max: Maximum allowed prime factor for Nfft.
            % -------------------------- Outputs --------------------------
            % - Nfft:      Computed FFT size.
            % -------------------------------------------------------------
            
            % Make sure that 2 is a factor
            Nfft_in     = ceil(Nfft_in/2)*2;
            % Valid factors
            prime_valid = primes(prime_max);
            % Auxiliary vector for each valid factor
            Nfft_r      = zeros(size(prime_valid));
            
            % For each valid factor...
            for r=1:length(prime_valid)
                % Init
                Nfft_r(r)   = Nfft_in;
                prime_fact  = factor(Nfft_r(r));
                % Iter until the maximum factor is below the threshold
                while any(prime_fact>prime_max)
                    % Find the factors that dont fulfill the threshold
                    idx1       = find(prime_fact>prime_max);
                    % Find the factors that fulfill the threshold
                    idx2       = setdiff(1:length(prime_fact),idx1);
                    % Compute new factors
                    new_factor = ceil(prod(prime_fact(idx1))/...
                                      prime_valid(r))*prime_valid(r);
                    % Update Nfft_r with the new factors
                    Nfft_r(r)  = prod(prime_fact(idx2))*new_factor;
                    % Compute factors of new Nfft_r
                    prime_fact = factor(Nfft_r(r));
                end
            end
            % Find the smalles FFT size
            Nfft = min(Nfft_r); 
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function FLOPs = getSolverFLOPs(solver,Ig,Ih,L,M,K,R,Ip,P)
            % -------------------------------------------------------------
            % Method to compute the number of Floating-Point Operations 
            % (FLOPS) required by one of the available solvers to compute 
            % the subband filters. A detail count of the FLOPs for this 
            % method can be found in found in the document wPM-S_flops.pdf
            % provided with this toolbox.
            % --------------------------- Inputs --------------------------
            % - solver:      String idincating the selected solver.
            % - Ig:          Filter length.
            % - Ih:          RIR length.
            % - L:           Number of loudspeakers.
            % - M:           Number of control points.
            % - K:           Number of subbands.
            % - R:           Resampling factor.
            % - Ip:          Prototype filter length.
            % - P:           Approximation order.
            % -------------------------- Outputs --------------------------
            % - FLOPs:       Total number of FLOPs.              
            % -------------------------------------------------------------

            % Length of subband filters
            Ig_k = ceil(Ig/R);
            % Length of the subband components of the RIR
            Ih_k = gdftFB_Class.getSBlength_static(Ih,Ip,R);
            % Length of the subband components of the target
            Id_k = Ig_k+Ih_k-1;
            % Length required for the broadband target such that a length
            % Id_k is obtained for the subband target
            Id = (Id_k+ceil(Ip/R)-2)*R-Ip+2;
            % Approximation order
            if ~isempty(P)
                assert(length(P)==1 | length(P)==K/2,...
               'wPMs_Class/getSolverFLOPs: The size of P must be 1 or K/2');
                % Store an approximation order for each subband       
                P    = repmat(P(:),(K/2)/length(P),1);
            end
            % FLOPs required to compute subband components of the RIR
            FLOPs_hk = gdftFB_Class.getsbDecFLOPs(Ih,M*L,K,R,Ip);
            % FLOPs required to compute subband components of the TR
            FLOPs_dk = gdftFB_Class.getsbDecFLOPs(Id,M,K,R,Ip);
            % FLOPs required to compute the subband filters
            if strcmp(solver,'exact')
                % FFT size 
                Nfft     = Id_k;
                % 1st order components
                FLOPs1   = (4*M*L*(L+3))*Nfft+(5/3)*L*Ig_k;
                % 2nd order components
                FLOPs2   = 11*(L*Ig_k).^2;
                 % 3rd order components
                FLOPs3   = (4/3)*(L*Ig_k)^3;
                % Logarithmic components
                FLOPslog = (M*L+M+(1/2)*L^2+(3/2)*L)*5*Nfft*log2(Nfft);
                % Total FLOPS
                FLOP_gk  = (K/2)*(FLOPs1+FLOPs2+FLOPs3+FLOPslog);
            elseif strcmp(solver,'superfast')
                % FFT size
                Nfft      = Id_k;
                % Logarithmic components
                FLOPslogP = 10*L*P(:)*Nfft*log2(Nfft);
                FLOPslog  = (M*L+M+3*L)*5*Nfft*log2(Nfft);
                % 1st order components
                FLOPs1P   =  8*L^2*P(:)*Nfft+2*L*P(:)*Ig_k;
                FLOPs1    = ((4/3)*L^3+11*L^2+(5/3)*L+4*M*L*(5*L+3))*Nfft+...
                             2*L*Ig_k;
                % Total FLOPS as the sum of the FLOPs for each subband
                FLOP_gk   = sum(FLOPs1P+FLOPs1+FLOPslog+FLOPslogP);
            else
                error('wPMsd_Class\getSolverFLOPS: Not supported solver');   
            end
            % Total number of FLOPs
            FLOPs = FLOPs_hk+FLOPs_dk+FLOP_gk;
        end
    end
end

