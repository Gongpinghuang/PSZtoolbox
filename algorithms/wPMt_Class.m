classdef wPMt_Class < handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class to compute the filters for a Personal Sound Zones (PSZ) system
    % the Weighted Pressure Matching with Time Domain formulation (wPM-T)
    % algorithm presented in:
    % 
    % M. F. Simon Galvez et al, "Time Domain Optimization of Filters Used 
    % in a Loudspeaker Array for Personal Audio," in IEEE/ACM Transactions 
    % on Audio, Speech, and Language Processing, 2015, 
    % doi: 10.1109/TASLP.2015.2456428.
    % ---------------------------------------------------------------------
    properties
        % ...................... General properties .......................
        
        L;               % Number of loudspeakers.
        M;               % Number of control points.
        Ig;              % Length of the filters.
        Ih;              % Length of the RIR.
        Id;              % Length of the target response.
        solver;          % Solver used to compute the filters, either 
                         % 'exact' or 'superfast'.
        P=0;             % Approximation order for the superfast solver.                 
        beta;            % Regularization factor.
        tag='wPM-T';     % Tag to identify the algorithm.

        % ....................... Impulse responses .......................
        
        h=[];            % RIR, -> [Ih x L x M].
        d=[];            % Target,-> [Id x M].

        % ............................ Flags ..............................  
        
        mexFlag = false; % Flag indicating if mex functions are used.
        issetup = false; % Flag indicating if the object is setup.
    end
    
    methods
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,g,ProcTime] = wPMt_Class(Ig,...
                                               mod_delay,...
                                               LSPref,...
                                               beta_rel,...
                                               solver,...
                                               P,...
                                               mexFlag,...
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
            % - solver:     Solver used for computing the filters.
            % - P:          Approximation order for the superfast solver.
            % - mexFlag:    Flag indicating if mex functions are used
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
            obj.solver   = solver;  
            obj.mexFlag  = mexFlag;
            obj.beta     = beta_rel*mean(mean(sum(h.^2,1)));
            obj.Ih       = size(h,1);
            obj.Ig       = Ig;
            obj.Id       = obj.Ig+obj.Ih-1;
            obj.h        = h;
            if strcmp(obj.solver,'superfast') 
                obj.P    = P;
            end
     
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
            % Compute filters according to solver
            if strcmp(obj.solver,'exact')
                if obj.mexFlag
                    g = obj.compFilt_exact_mex();
                else
                    g = obj.compFilt_exact_matlab();
                end
            elseif strcmp(obj.solver,'superfast')
                if obj.mexFlag
                    g = obj.compFilt_superfast_mex();
                else
                    g = obj.compFilt_superfast_matlab();
                end
            else
                error('wPMt_Class: Not supported solver');
            end
          
            ProcTime = toc;
        
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g = compFilt_exact_matlab(obj)
            % -------------------------------------------------------------
            % Matlab implementation to compute the filters using an exact 
            % solver, i.e., using the expressiong g=(H'*H+beta*I)^-1*H'*d. 
            % The exact method is based on using the Choleksy decomposition 
            % and uses the FFT to compute H'*H, and H'*d. A detalied 
            % description of this method can be found in:
            %
            % V. Moles-Cases, "Filter Optimization for Personal Sound Zones 
            % Systems," PhD thesis, 2022, Available in September 2022 in 
            % UPV repository.
            % -------------------------- Outputs --------------------------
            % - g:   Filters, -> [Ig x L].    
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMt_Class/compFilt_exact: ',...
                                 'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Required FFT size to compute H'*H and H'*d using FFTs, we
            % force that the size is even
            Nfft  = ceil((obj.Ih+obj.Ig-1)/2)*2;
            % FFT of the RIR
            Hf    = fft(obj.h,Nfft,1);
            % FFT of the target
            Df    = fft(obj.d,Nfft,1);
            % Keep positive freq. bins
            Hf    = Hf(1:(Nfft/2)+1,:,:);
            Df    = Df(1:(Nfft/2)+1,:);
            % Idx of the Ig_k positive lags of the correlation
            lag_p = 1:obj.Ig;
            % Idx of the Ig_k-1 negative lags of the correlation
            lag_n = 2:obj.Ig;
            
            % ----------------------- Compute filters ---------------------      
            
            % ............. Compute correlation matrix R=H'*H .............
            % Compute correlation of h using the FFT
            R = ifft(sum(conj(Hf).*permute(Hf,[1,4,3,2]),3),...
                        Nfft,1,'symmetric');
            % Keep lag_p and lag_n 
            R = cat(3,permute(R(lag_p,:,:,:),[2,4,1,3]),...
                      permute(conj(R(lag_n,:,:,:)),[4,2,1,3]));
            % Transform to a cell with 2*Ig-1 elements of size LxL       
            R = num2cell(R,[1,2]);
            % Form the block-toeplitz correlation matrix of h
            R = cell2mat(R(toeplitz(lag_p,[1,lag_n+obj.Ig-1])));

            % ............. Compute correlation vector r=H'*d .............
            % Compute cross-correlation using the FFT
            r = ifft(sum(conj(Hf).*reshape(Df,(Nfft/2)+1,1,[]),3),...
                         Nfft,1,'symmetric');
            % Keep lag_p           
            r = reshape(r(lag_p,:).',[],1);

            % .................. Compute optimal filters ..................
            % Obtain optimal filters by solving the linear system
            g_aux = linsolve(R+obj.beta*eye(obj.L*obj.Ig),...
                             r,...
                             struct('SYM',true,'POSDEF',true));

            % ........................... Store ...........................     
            g     = reshape(g_aux,[],obj.Ig).';
        end
        
                %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g = compFilt_exact_mex(obj)
            % -------------------------------------------------------------
            % MEX implementation to compute the filters using an exact 
            % solver, i.e., using the expressiong g=(H'*H+beta*I)^-1*H'*d. 
            % The exact method is based on using the Choleksy decomposition 
            % and uses the FFT to compute H'*H, and H'*d. A detalied 
            % description of this method can be found in:
            %
            % V. Moles-Cases, "Filter Optimization for Personal Sound Zones 
            % Systems," PhD thesis, 2022, Available in September 2022 in 
            % UPV repository.
            % -------------------------- Outputs --------------------------
            % - g:   Filters, -> [Ig x L].    
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMt_Class/compFilt_exact_mex: ',...
                                 'The object is not setup']);
                             
            % ----------------------- Compute filters ---------------------        
            % Compute using mex function
            g = exact_LS_real_mex(obj.h,...
                                  obj.d,...
                                  obj.L,...
                                  obj.M,...
                                  obj.Ig,...
                                  obj.beta);   
            % Reshape                  
            g = reshape(g,obj.Ig,[]);                      
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g = compFilt_superfast_matlab(obj)
            % -------------------------------------------------------------
            % Matlab implementation to compute an approximation of the 
            % filters using the superfast solver proposed in:
            %
            % M. Poletti et al, "A Superfast Toeplitz Matrix Inversion
            % Method for Single- and Multi-Channel Inverse Filters and its
            % Application to Room Equalization," in  IEEE/ACM Transactions
            % on Audio Speech and Language Processing, 29, 3144�3157, 2021.
            % https://doi.org/10.1109/TASLP.2021.3120650
            % -------------------------- Outputs --------------------------
            % - g:   Filters, -> [Ig x L x K/2].     
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMt_Class/compFilt_superfast_matlab: ',...
                                'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Required FFT size for the superfast solver
            Nfft         = obj.findFFTsize(obj.Id,11);
            % Number of positive freq. bins
            Kfft         = (Nfft/2)+1;
            % FFT of the RIR
            hf           = permute(fft(obj.h,Nfft,1),[3,2,1,4]);
            % FFT of the target
            df           = permute(fft(obj.d,Nfft,1),[2,1,3]);
            % Init auxilairy matrices
            Hf_ps        = zeros(obj.L,obj.M);
            Hf_ps_x_Hf   = zeros(Kfft,obj.L,obj.L);
            gf           = zeros(Kfft,obj.L);
            
            % ----------------------- Compute filters ---------------------        

            % ....... Compute filters using freq. domain formulation ......
            % For each positive freq. bin...
            for f=1:Kfft
                % Compute Hf^pseudo = Hf^H*(Hf*Hf^H+beta*I)^{-1} 
                Hf_ps(:) = ...
                    linsolve(hf(:,:,f)*hf(:,:,f)'+obj.beta*eye(obj.M),...
                             hf(:,:,f),...
                             struct('SYM',true,'POSDEF',true))';
                % Compute Hf^pseudo*Hf
                Hf_ps_x_Hf(f,:,:) = reshape(Hf_ps*hf(:,:,f),1,obj.L,obj.L);
                % Compute optimal freq. filter coefficients    
                gf(f,:)  = Hf_ps*df(:,f);
            end
            % Compute IFFT
            g_fd_full = ifft(gf,Nfft,1,'symmetric');
            % Truncate filters
            g_fd      = g_fd_full(1:obj.Ig,:);

            % ............. Compute and apply correction terms ............
            % Init time and freq. domain residuals
            rt  = g_fd_full;
            rf  = zeros(Nfft,obj.L);
            % Iterate until reach the selected approximation...
            for p=0:obj.P
                % Truncate residual
                rt(1:obj.Ig,:) = 0;
                % Compute its FFT
                rf(:)          = fft(rt,Nfft,1);
                % Correct freq. domain residuals
                rf(1:Kfft,:)   = ...
                          sum(Hf_ps_x_Hf.*reshape(rf(1:Kfft,:),[],1,obj.L),3);
                % Compute time-domain corrected residuals
                rt(:)          = ifft(rf,Nfft,1,'symmetric');
                % Correct filters
                g_fd(:)        = g_fd+rt(1:obj.Ig,:);
            end
            % Store
            g = g_fd;
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g = compFilt_superfast_mex(obj)
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
            % - g: Filters, -> [Ig x L x K/2].  
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['wPMt_Class/compFilt_superfast_mex: ',...
                                 'The object is not setup']);
                             
            % ------------------------- Initialize ------------------------
            % Required FFT size to compute H'*H and H'*d using the FFT
            Nfft           = obj.findFFTsize(obj.Id,11);

            % ----------------------- Compute filters ---------------------        
            % Compute using mex function
            gaux = superfast_LS_real_mex(obj.h,...
                                         obj.d,...
                                         obj.L,...
                                         obj.M,...
                                         obj.Ig,...
                                         obj.beta,...
                                         Nfft,...
                                         obj.P);             
            % Store
            g = reshape(gaux,obj.Ig,[]);    
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
        function FLOPs = getSolverFLOPs(solver,Ig,Ih,L,M,P)
            % -------------------------------------------------------------
            % Method to compute the number of Floating-Point Operations 
            % (FLOPS) required by one of the available solvers to compute 
            % the filters. The FLOPs counts for the different solvers
            % can be found in:
            %
            % V. Moles-Cases, "Filter Optimization for Personal Sound Zones 
            % Systems," PhD thesis, 2022, Available in September 2022 in 
            % UPV repository.
            % --------------------------- Inputs --------------------------
            % - solver:      String idincating the selected solver.
            % - Ig:          Filter length.
            % - Ih:          RIR length.
            % - L:           Number of loudspeakers.
            % - M:           Number of control points.
            % - P:           Approximation order.
            % -------------------------- Outputs --------------------------
            % - FLOPs:       Total number of FLOPs.              
            % -------------------------------------------------------------
            
            if strcmp(solver(1:end-4),'exact')
                % FFT size 
                Nfft     = Ih+Ig-1;
                % 1st order components
                FLOPs1   = ((1/6)*L*Ig)+(2*M*L*(L+3)*Nfft);
                % 2nd order components
                FLOPs2   = (5/2)*(L*Ig)^2;
                % 3rd order components
                FLOPs3   = (1/3)*(L*Ig)^3;
                % Logarithmic components
                FLOPslog = (M*L+M+(1/2)*(L^2)+(3/2)*L)*(5/2)*Nfft*log2(Nfft);
                % Total FLOPS
                FLOPs    = FLOPs1+FLOPs2+FLOPs3+FLOPslog;
            elseif strcmp(solver(1:end-4),'superfast') 
                 % FFT size
                Nfft      = Ig+Ih-1;
                % Logarithmic components
                FLOPslogP = L*P*5*Nfft*log2(Nfft);
                FLOPslog  = (M*L+M+3*L)*(5/2)*Nfft*log2(Nfft);
                % 1st order components
                FLOPs1    = ((2/3)*L^3+5*L^2+(5/6)*L)*Nfft+...
                             L*Ig+(8*(L.^2)+2*L)*P*Nfft+...
                             2*M*L*(5*L+3)*Nfft;
                % Total FLOPS
                FLOPs    = FLOPs0+FLOPs1+FLOPslog+FLOPslogP;
            else
                error('wPMt_Class\getSolverFLOPS: Not supported solver');   
            end
        end
        
    end
end

