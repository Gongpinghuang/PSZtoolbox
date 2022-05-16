classdef gdftFB_Class < handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class with all the required signal processing methods related with
    % Generalized Discrete Fourier Transform (GDFT) filter banks. More 
    % information about GDFT filter banks can be found in:
    %
    % S.Weiss,"On adaptive filtering in oversampled subbands," Ph.D.
    % dissertation, 1998.
    % ---------------------------------------------------------------------
    properties
        
       % ..................... Filter bank properties .....................
       
       K;                % Number of single-sided subbands.
       R;                % Resampling factor.
       Ip;               % Length of the prototype filter.
       hp;               % Prototype filter,  -> [Ip x 1].
       ha;               % Analysis filters,  -> [Ip x K/2].
       hs;               % Synthesis filters, -> [Ip x K/2].
       delay;            % Overrall delay of the filter bank.
       ANp;              % Structure with the required parameters for the 
                         % analysis polyphase filtering.
       SYp;              % Structure with the required parameters for the 
                         % synthesis polyphase filtering.  
       RE_dB;            % Reconstruction Error of the filter bank.                 
       ASR_dB;           % Aliast-To-Signal Ratio of the filter bank.  
       
       % ............................. Flags ..............................
       
       mexFlag = false;  % Flag indicating if mex functions are used with 
                         % the object.
       issetup = false;  % Flag indicating if the object is set up.
    end
    
    methods
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = gdftFB_Class(K,R,Ip,mexFlag)
            % -------------------------------------------------------------                    
            % Constructor with the required parameters. 
            % --------------------------- Inputs --------------------------
            % - K:         Number of single-sided subbands.
            % - R:         Resampling factor.
            % - Ip:        Length of the prototype filter.
            % - mexFlag:   If true, mex functions are used for certain 
            %              operations of the filter bank.
            % -------------------------- Outputs --------------------------
            % - obj:       Handle to object.
            % -------------------------------------------------------------
                  
            % ------------------------ Check Inputs -----------------------
            assert(R<=K,'gdftFB_Class: K must be equal or greater than R');
            assert(mod(K,2)==0,'gdftFB_Class: K must be even');
            
            % ------------------- Filter bank parameters ------------------
            % Set properties   
            obj.K          = K;
            obj.R          = R;
            obj.Ip         = Ip;
            obj.delay      = Ip-1; 
            obj.mexFlag    = mexFlag;

            % ----------------- Generate prototype filter -----------------
            % Path for the prototype filter
            savePath = ['filterBank/database/hp_K_',int2str(K),'_R_',...
                        int2str(R),'_Ip_',int2str(Ip),'.mat'];
            % Generate or load prototype filter
            if exist(savePath,'file')==2
                % If the file exists load the prototype filter
                load(savePath,'hp','RE_dB','ASR_dB');
            else
                % If the file does not exist generate the prototype filter
                [hp,RE_dB,ASR_dB] = gdftFB_Class.comp_PF(Ip,K,R);
                % Save filter 
                save(savePath,'hp','RE_dB','ASR_dB');
            end
            % Store in object properties
            obj.hp     = hp;
            obj.RE_dB  = RE_dB;
            obj.ASR_dB = ASR_dB;
            
            % ---------- Generate analysis and synthesis filters ----------
             % Analysis filters
            obj.ha = hp(:).*exp(j*2*pi*((0:K/2-1)+0.5).*((0:Ip-1)).'/K);
            % Synthesis filters
            obj.hs = conj(flip(obj.ha,1));
            
            % ---------- Setup analysis and synthesis polyphase -----------
            if obj.mexFlag
                obj.setup_polyphase();
            end

            % -------------------------------------------------------------
            % Object is setup
            obj.issetup   = true;
            
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setup_polyphase(obj)
            % -------------------------------------------------------------
            % Method to compute the required parameters to perform the 
            % analysis and synthesis polyphase filtering for the GDFT
            % filter bank. The efficient implementation is based on the 
            % procedure described in section 4.2 of:
            %
            % S.Weiss, "On adaptive filtering in oversampled subbands," 
            % Ph.D. dissertation, 1998. 
            % -------------------------------------------------------------
        
            % -------------------------- General --------------------------
            % Number of branches for the polyphase network
            M   = lcm(2*obj.K,obj.R);
            % Number of groups of branches
            Q   = M/obj.R;
            % Closer multiple of M that is higher than Ip
            IpM = ceil(obj.Ip/M)*M;
            % Upsampled polyphase components
            hp_P            = zeros(((IpM*Q)/M)-(Q-1),M);
            hp_P(1:Q:end,:) = reshape([obj.hp;zeros(IpM-obj.Ip,1)],M,[]).';
            % Shift polyphase components
            Nfft  = size(hp_P,1)+Q-1;
            shift = exp(-j*2*pi*repelem(0:Q-1,obj.R).*(0:Nfft-1).'/Nfft);
            hp_P  = ifft(shift.*fft(hp_P,Nfft,1),Nfft,1,'symmetric');
            % Remove approximation errors produced by the fft
            hp_P(abs(hp_P)<1e-15) = 0;
            % Index of the branches to which we must invert the sign
            m_idx         = find(mod(0:M-1,2*obj.K)>=obj.K);
            % Invert the sign to certain branches
            hp_P(:,m_idx) = -hp_P(:,m_idx);
            % Remove non-needed 0 taps at the end of the responses 
            hp_P          = hp_P(1:max(find(sum(abs(hp_P),2)~=0)),:);
            
            % ------------------------- Analysis --------------------------
            % Store the tap, the coefficient, and the branch for the non-0
            % filter taps in the analysis polyphase network
            [obj.ANp.tap,obj.ANp.branch] = find(abs(hp_P)~=0);
            obj.ANp.coeff                = hp_P(abs(hp_P)~=0);
            obj.ANp.tap                  = obj.ANp.tap-1;
            obj.ANp.branch               = obj.ANp.branch-1;
            
            % ------------------------- Synthesis -------------------------
            % The synthesis polyphase network is a time reversed and 
            % conjugated version of the analysis polyphase network 
            hp_P                         = flip(hp_P,1); 
            % Store the tap, the coefficient, and the branch for the non-0
            % filter taps in the synthesis polyphase network
            [obj.SYp.tap,obj.SYp.branch] = find(abs(hp_P)~=0);
            obj.SYp.coeff                = hp_P(find(abs(hp_P)~=0));
            obj.SYp.tap                  = obj.SYp.tap-1;
            obj.SYp.branch               = obj.SYp.branch-1;
            
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = filter(obj,x,g_k)
            % -------------------------------------------------------------
            % Method used to filter a signal x with the filter bank using a
            % set of L subband filters in each subband. For a PSZ, L is the 
            % number of loudspeakers. Hence, in each subband we filter the 
            % input signal with the subband filters for each loudspeaker, 
            % then, the output of the filter bank is a set of L signals, 
            % one for each loudspeaker. 
            % --------------------------- Inputs --------------------------
            % - x:   Input signal, -> [Ix x 1]
            % - g_k: Subband filters, -> [Ig_k x L x K/2].
            % -------------------------- Outputs --------------------------
            % - y:   Output signal when feeding the filter bank with x,
            %        -> [Iy x L]
            % -------------------------------------------------------------
        
            % Check
            assert(obj.issetup,'gdftFB_Class/filter: Object must be setup');
            assert(any(size(x)==1),...
                         'gdftFB_Class/filter: x must be a single signal');
            assert(size(g_k,3)==obj.K/2,...
                            'gdftFB_Class/filter: Incorrect size for g_k');

            % Analysis filtering
            x_an = obj.analysis_filter(x);
            % Subband filtering 
            x_k  = obj.subband_filter(x_an,g_k);
            % Synthesis filtering 
            y    = obj.synthesis_filter(x_k);
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x_an = analysis_filter(obj,x)
            % -------------------------------------------------------------
            % Method to filter a set of signals with the analysis section 
            % of the filter bank. 
            % --------------------------- Inputs --------------------------
            % - x:    Set of Q signals of length Ix to filter, -> [Ix x Q]
            % -------------------------- Outputs --------------------------
            % - x_an: Analysis signals, -> [ceil((Ix+Ip-1)/R) x Q x K/2]
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['gdftFB_Class/analysis_filter: Object',...
                                'must be setup']);
            % Map higher order dimensions to the second dimension
            x      = x(:,:);
            
            % --------------------------- Filter --------------------------
            % Setup
            if obj.mexFlag
                % .................. Using mex functions .................. 
                % Filter
                x_an  = gdft_fb_analysis_mex(x,...
                                             obj.K,...
                                             obj.R,...
                                             obj.Ip,...
                                             length(obj.ANp.coeff(:)),...
                                             obj.ANp.coeff(:),...
                                             obj.ANp.tap(:),...
                                             obj.ANp.branch(:));  
                % Permute dimensions
                x_an  = permute(x_an(1:obj.K/2,:,:),[2 3 1]);
            else
                % ................. Using matlab functions ................
                % FFT size
                Nfft  = obj.Ip+size(x,1)-1;
                % Filter input signal with analysis filters using FFT
                X_an  = fft(x,Nfft,1).*...
                        reshape(fft(obj.ha,Nfft,1),Nfft,1,[]);
                x_an  = ifft(X_an,Nfft,1);
                % Downsample
                x_an  = x_an(1:obj.R:end,:,:);
            end  
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x_k = subband_filter(obj,x_an,g_k)
            % -------------------------------------------------------------
            % Method to filter the signals at the output of the analysis
            % stage of the filter bank with a set of L filters (for each
            % subband).
            % --------------------------- Inputs --------------------------
            % - x_an:   Analysis signals, -> [Ix_an x 1 x K/2].
            % - g_sb:   Subband filters,  -> [Ig_k x L x K/2].
            % -------------------------- Outputs --------------------------
            % - x_k:   Subband signal,   -> [Ix_k x L x K/2],
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['gdftFB_Class/subband_filter: Object',...
                                'must be setup']);
            assert(size(x_an,2)==1 && size(x_an,3)==obj.K/2 && ...
                   size(g_k,3)==obj.K/2,...
                   'gdftFB_Class/subband_filter: Incorrect input size');

            % Filter using the FFT
            Nfft = size(x_an,1)+size(g_k,1)-1;
            x_k  = ifft(fft(x_an,Nfft,1).*fft(g_k,Nfft,1),Nfft,1);
        end

        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = synthesis_filter(obj,x_k)
            % -------------------------------------------------------------
            % Method to filter L subband signals with the synthesis section
            % of the filter bank.
            % ------------------------- Inputs ----------------------------
            % - x_k: Set of subband signals of length,-> [Ix_sb x L x K/2]
            % ------------------------ Outputs ----------------------------
            % - y:   L output signals -> [Ix_sb*R+Ip-1 x L]
            % -------------------------------------------------------------
           
            % Check
            assert(obj.issetup,['gdftFB_Class/synthesis_filter: Object',...
                                'must be setup']);
                            
            % --------------------------- Filter --------------------------
            % Setup
            if obj.mexFlag
                % .................. Using mex functions ..................
                % Accomodate dimentions and apply hermitian symmetry                         
                x_k    = permute(x_k,[3,1,2]);
                x_k    = cat(1,x_k,conj(x_k(end:-1:1,:,:)));
                % Make sure that the array is complex 
                x_k(1) = complex(real(x_k(1)),imag(x_k(1)));
                % Filter
                y      = gdft_fb_synthesis_mex(x_k,...
                                               obj.K,...
                                               obj.R,...
                                               obj.Ip,...
                                               length(obj.SYp.coeff(:)),...
                                               obj.SYp.coeff(:),...
                                               obj.SYp.tap(:),...
                                               obj.SYp.branch(:));              
            else

                % ................. Using matlab functions ................
                % Upsample subband signals
                x_k_us  = zeros(size(x_k,1)*obj.R,size(x_k,2),obj.K/2);
                x_k_us(1:obj.R:end,:,:) = x_k;
                % Filter using FFTs and apply hermitian symmetry
                Nfft    = size(x_k_us,1)+obj.Ip-1;
                Y       = sum(fft(x_k_us,Nfft,1).*...
                              reshape(fft(obj.hs,Nfft,1),Nfft,1,[]),3);
                y       = 2*real(ifft(Y,Nfft,1));
            end  
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x_k = subband_dec(obj,x)
            % -------------------------------------------------------------
            % Method to compute the subband decomposition of a set of Nx
            % FIR broadband signals of length Ix. The method is based on:
            % 
            % - V. Moles-Cases et al, "Personal Sound Zones by Subband 
            % Filtering and  Time Domain Optimization," in IEEE/ACM 
            % Transactions on Audio, Speech,  and Language Processing, 
            % 2020, doi: 10.1109/TASLP.2020.3023628.
            % 
            % - J. P. Reilly et al "The complex subband decomposition and 
            % its application to the decimation of large adaptive filtering 
            % problems," in IEEE Transactions on Signal Processing, 2002,
            % doi: 10.1109/TSP.2002.804068.
            % --------------------------- Inputs --------------------------
            % - x:   Input signals of length Ix, -> [Ix x Nx].
            % -------------------------- Outputs --------------------------
            % - x_k: Subband component of the Nx input signals in the  
            %        subbands of the positive spectrum of the filter bank, 
            %        -> [Ix_k x Nx x K/2].
            % -------------------------------------------------------------
        
            % ---------------------- Input dimensions ---------------------
            % Length of the FIR broadband signals
            Ix   = size(x,1);
            % Number of different input FIR signals
            Nx   = size(x,2);
            
            % ---------------- Compute analysis components ----------------
            x_an = obj.analysis_filter(x);
            
            % ---------------- Compute subband components -----------------
            if obj.mexFlag
                % .................. Using mex functions ..................
                x_k = gdft_fb_sb_dec_mex(x_an,obj.hp,obj.K,obj.R);
            else
                % ................. Using matlab functions ................
                % Length of the subband components of x
                Ix_k    = obj.getSBlength(Ix);
                % Length of the analysis components of x
                Ix_an   = ceil((obj.Ip+Ix-1)/obj.R);
                % Pseudo-inverse of the convolution matrix for the  
                % downsampled prototype filter
                Wp_pinv = pinv(convmtx(obj.hp(1:obj.R:end),Ix_k));
                % Auxiliary vector
                fk      = (0.5:(obj.K/2)-0.5).';
                % Shift matrices
                F1      = exp(j*2*pi*fk*(0:Ix_an-1)*obj.R/obj.K).';
                F2_inv  = exp(j*2*pi*fk*(0:Ix_k-1)*obj.R/obj.K).';
                % Initialize components for each subband
                x_k = zeros(Ix_k,Nx,obj.K/2);
                % For each subband...
                for k = 1:obj.K/2
                    % Compute subband components for the k-th subband
                    x_k(:,:,k)= F2_inv(:,k).*...
                                 (Wp_pinv*(conj(F1(:,k)).*(x_an(:,:,k))));
                end
            end
            
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Ix_k = getSBlength(obj,Ix)
            % -------------------------------------------------------------
            % Method to compute the length of the subband components of a 
            % Ix-length FIR signal using the current filter bank.
            % ------------------------- Inputs ----------------------------
            % - Ix:   Length of the broadband FIR signal.
            % ------------------------ Outputs ----------------------------
            % - Ix_k: Length of the subband components.
            % -------------------------------------------------------------
            
            % Check
            assert(obj.issetup,['gdftFB_Class/getSBlength: Object must',... 
                                'be setup']);
            % Compute length of the subband components
            Ix_k = obj.getSBlength_static(Ix,obj.Ip,obj.R);
        end
        
    end
    
    methods(Static)
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Ix_k = getSBlength_static(Ix,Ip,R)
            % -------------------------------------------------------------
            % Method to compute the length of the subband components of a 
            % Ix-length FIR signal using a given filter bank.
            % ------------------------- Inputs ----------------------------
            % - Ix:   Length of the broadband FIR signal.
            % ------------------------ Outputs ----------------------------
            % - Ix_k: Length of the subband components.
            % -------------------------------------------------------------

            % Compute length of the subband components
            Ix_k = ceil((Ix+Ip-1)/R)-ceil(Ip/R)+1;
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [hp,RE_dB,ASR_dB] = comp_PF(Ip,K,R)
            % -------------------------------------------------------------
            % Method to design a proptotype filter for a GDFT FB with K 
            % subbands and resampling factor R. The method is a variation
            % of that presented in:
            % 
            % S.Weiss, "On adaptive filtering in oversampled subbands," 
            % Ph.D. dissertation, 1998,
            %
            % in which the prototype filter is obtained by minizing the
            % reconstruction error of the filter bank and the stop band
            % energy of the prototype filters
            % -------------------------- Inputs ---------------------------
            % - Ip:     Prototype filter length (odd).
            % - K:      Number of subbands.
            % - R:      Resampling factor.
            % -------------------------- Outputs --------------------------
            % - hp:     Prototype filter, -> [Ip x 1].
            % - RE_dB:  Reconstruction Error.
            % - ASR_dB: Alias-To-Signal ratio.
            % -------------------------------------------------------------
            
            % ------------------------ Check Inputs -----------------------
            assert(mod(Ip,2)==1,['gdftFB_Class/comp_PF: The length of '...
                                 'the prototype filter must be odd']);
               
            % -------------------- Algorithm parameters -------------------
            % Number of iterations for the algorithm
            Niter    = 20;
            % Define weighting factors to evaluate
            sigma    = 10.^(-2:0.01:2);
            % Define smoothing factor
            w        = 0.5;
            % Stoping threshold
            th       = 10^(-3);
            
            % ----------------------- Initial filter ----------------------   
            if mod(Ip+1,K)==0
                % Square root raised cosine filter
                h0 = (sqrt(R/K))*rcosdesign(R/K,(Ip-1)/K,K).';
            else
                % Find the Square root raised cosine filter with length 
                % longer than Ip
                Ip_aux = ceil(Ip/K)*K-1;
                h0      = (sqrt(R/K))*rcosdesign(R/K,(Ip_aux-1)/K,K).';
                % Clip
                Ip_clip = (Ip_aux-Ip)/2;
                h0      = h0(Ip_clip+1:end-Ip_clip);
            end
            % Keep firt half of initial prototype
            b0 = h0(1:ceil(Ip/2)); 

            % ---------------- Setup optimization elements ----------------
            % Disable warnings
            warning('off','MATLAB:nearlySingularMatrix');
           
            % Number of unique elements of the prototype filter
            Qp       = length(b0);
            % Mapping matrix to convert bp in hp
            PI       = [eye(Qp);flip(eye(Qp-1),1),zeros(Qp-1,1)];
            % FFT size
            Nfft     = 2*Ip-1;
            % FFT matrix 
            F        = dftmtx(Nfft);
            F        = F(:,1:Ip);
            % Fourier transform of matrix PI
            PI_f     = (K/R)*F*PI;
            % Normalized frequency bins
            u        = (0:Nfft)/Nfft;
            % Idx of the frequency bins in the stopband of the prot. filter
            idx_stb  = find(u>(1/(2*R)) & u<(1-(1/(2*R))));  
            % Lags of the correlation of hp that are multiple of K
            corr_idx = [1:K:Ip,Nfft-K+1:-K:Ip+1];
            % Target cross-correlation of hp
            d        = [1;zeros(length(corr_idx)-1,1)];   
            % Auxiliary matrices and vectors
            A        = zeros(Qp,Qp);
            B        = zeros(Nfft,Qp);
            r        = zeros(Qp,1);
            
            % ------------------ Compute prototype filter -----------------
            % Init filters for each weighting factor
            bp       = repmat(b0(:),1,length(sigma));
            % Init auxiliary filter for the iterative algorithm
            bp_i     = zeros(size(b0));
            % Init reference metrics
            RE_dB    = Inf*ones(length(sigma),1);
            ASR_dB   = Inf*ones(length(sigma),1);
            
            % For each weighting factor...
            for s=1:length(sigma)
                % For each iteration...
                for i=1:Niter
                    % Update prototype filter
                    P       = fft(PI*real(bp(:,s)),Nfft,1);
                    B(:)    = ifft(conj(P(:)).*PI_f,Nfft,1);
                    A(:)    = (B(corr_idx,:)'*B(corr_idx,:))+...
                          (sigma(s).^2)*(PI_f(idx_stb,:)'*PI_f(idx_stb,:));
                    r(:)    = B(corr_idx,:)'*d;
                    bp_i(:) = w*bp(:,s)+(1-w)*real(linsolve(A,r));
                    % Compute metrics
                    [RE_dB_i,ASR_dB_i] = ...
                          gdftFB_Class.compute_RE_ASR(PI*real(bp_i),K,R); 
                    % Check stop criterion
                    if RE_dB_i<RE_dB(s) || ASR_dB_i<ASR_dB(s) || ...
                                    (norm(bp(:,s)-bp_i)/norm(bp(:,s)))<th
                        bp(:,s)   = bp_i;
                        RE_dB(s)  = RE_dB_i;
                        ASR_dB(s) = ASR_dB_i;
                    else
                        break;
                    end
                end
            end
            
            % Find weighting factor that leads to more similar RE and ASR
            [~,idx] = min(abs(RE_dB-ASR_dB));
            % Keep prototype filter, RE, and ASR for that weighting factor
            bp      = bp(:,idx);
            RE_dB   = RE_dB(idx);
            ASR_dB  = ASR_dB(idx);
            % Mirror bp to obtain final protoype filter                                                          
            hp      = PI*bp(:);
            % Enable warnings
            warning('on','MATLAB:nearlySingularMatrix');                                                                                                        
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [RE_dB,ASR_dB] = compute_RE_ASR(hp,K,R)
            % -------------------------------------------------------------
            % Method to compute the Reconstruction Error (RE) and the 
            % Alias-To-Signal Ratio (ASR) of the filter bank.
            % -------------------------- Inputs ---------------------------
            % - hp:      Prototype filter, -> [Ip x 1]
            % - K:       Number of subbands.
            % - R:       Resampling factor. 
            % -------------------------- Outputs --------------------------
            % - RE_dB:   Reconstruction error (RE) with respect to a
            %            distortion-less distortion function.
            % - ASR_dB:  Aliasing-To-Signal Ratio (ASR) in the subbands.
            % -------------------------------------------------------------
            
            % Length of the prototype filter
            Ip = length(hp);
            
            % ........................ Compute RE .........................
            % Compute auto-correlation of prototype filter
            Nfft  = 2*Ip;
            Kfft  = (Nfft/2)+1;
            Hp    = fft(hp(:),Nfft,1);
            r     = ifft(Hp.*conj(Hp),Nfft,1,'symmetric');
            % Keep lags that are multiple of K and scale
            r     = (K/R)*[r(1:K:Kfft);r(end-(K-1):-K:Kfft+1)];
            % Compute the reconstruction error
            RE_dB = 10*log10((r(1)-1).^2+sum(r(2:end).^2));
   
            % ........................ Compute ASR ........................
            % Compute convolution of prototype filter and aliasing terms
            Y      = fft(exp(j*2*pi*(1:R-1).*(0:Ip-1).'/R).*hp(:),Nfft,1);
            X      = ifft(Hp.*Y,Nfft,1);
            % Compute aliasing in the subbands normalized by filter energy
            ASR_dB = 10*log10((1/(R-1))*sum(abs(X(:)).^2)/sum(hp.^2));
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function FLOPs = getsbDecFLOPs(Ix,Nx,K,R,Ip)
            % -------------------------------------------------------------
            % Method to compute the number of Floating-Point Operations 
            % (FLOPS) required to compute the subband decomposition of a
            % set of Nx signals of length Ix. We use the subband
            % decomposition presented in:
            % 
            % - V. Moles-Cases et al, "Personal Sound Zones by Subband 
            % Filtering and  Time Domain Optimization," in IEEE/ACM 
            % Transactions on Audio, Speech,  and Language Processing, 
            % 2020, doi: 10.1109/TASLP.2020.3023628.
            %
            % A detail count of the FLOPs for this method can be found in
            % found in the document wPM-S_flops.pdf provided with this
            % toolbox.
            % --------------------------- Inputs --------------------------
            % - Ix:    Length of the FIR signals to decompose.
            % - Nx:    Number of signals of length Ix that we want to
            %          decompose.
            % - K:     Number of subbands.
            % - R:     Resampling factor.
            % - Ip:    Prototype filter length.
            % -------------------------- Outputs --------------------------
            % - FLOPs: Total number of FLOPs required to perform the 
            %          subband decomposition.           
            % -------------------------------------------------------------
            
            % -------------------------------------------------------------
            % Length of analysis components
            Ix_an = ceil((Ix+Ip-1)/R);
            % Length of subband components
            Ix_k  = gdftFB_Class.getSBlength_static(Ix,Ip,R);
            % FFT size
            Nfft  = Ix_an;
            % Compute FLOPs
            FLOPs = (1/3)*Ix_k.^3;
            FLOPs = FLOPs+(2*Nx*K+(1/2))*Ix_k.^2;
            FLOPs = FLOPs+5*(Nx*K+1)*Nfft.*log2(Nfft);
            FLOPs = FLOPs+(3*Nx*K+(1/6))*Ix_k;
            FLOPs = FLOPs+(Nx/R)*(2*Ip+2*K+5*K*log2(K))*Ix;
            FLOPs = FLOPs+6*(Nx*K+1)*Ix_an;
        end

    end
end

