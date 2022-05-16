% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to build the mex files required by the Personal Sound Zones 
% (PSZ) simulator.
% -------------------------------------------------------------------------
function buildmex()
    
    % Clear mex related memory to avoid undesired behaviour 
    clear mex
    % Display
    dispPSZ('Building mex files',1);
    % Check if bin folder exists
    if ~exist('mex/bin','dir')
        % If it does not exist, create it and add it to matlab path
        mkdir('mex/bin'); 
        addpath('mex/bin');
    end
 
    % ------------------------- FFTW3 library setup -----------------------        
    if ispc
        % ......................... Windows OS ............................
        % Path where the fftw3 library is located
        FFTWpath    = 'mex/fftw_3.3.5_win64/';
        % Library name
        FFTWlibfile = 'libfftw3-3.dll';
        % Copy library .dll file to bin folder where the executable mex
        % files will be located (if is not there already)
        if ~exist(['mex/bin/',FFTWlibfile],'file')
            copyfile([FFTWpath,FFTWlibfile],'mex/bin/');
        end
    elseif ismac
        % ........................... MAC OS ..............................
        % Path where the fftw3 library is located
        FFTWpath    = '/opt/local/lib/';
        % Path to fftw headers
        FFTWheader  = '/opt/local/include/';
        % Library name
        FFTWlibfile = 'libfftw3.a';
        % Check if lib file exists
        if ~exist([FFTWpath,FFTWlibfile],'file')
             error(['A valid FFTW library is not found in ',FFTWpath,...
                  '. It can be installed following the steps in: ',...
                  'http://www.fftw.org/install/mac.html']); 
        end
    end  
    
    % ----------------------- Mex functions to build ----------------------
    mex_name = {'superfast_LS_real',...
                'superfast_LS_complex',...
                'exact_LS_real',...
                'exact_LS_complex',...
                'gdft_fb_analysis',...
                'gdft_fb_synthesis',...
                'gdft_fb_sb_dec'};
    % Extension of the mex file for different OS  
    if ispc
        mex_ext = 'mexw64';
    elseif ismac
        mex_ext = 'mexmaci64';
    elseif isunix
        mex_ext = 'mexa64';
    end        
 
    % --------------------------- Build mex filtes ------------------------     
    % For each file to build...
    for i=1:numel(mex_name)
        dispPSZ(mex_name{i},3);
        % Path to the mex file
        mex_path    = ['mex/source/',mex_name{i},'_mex.c'];
        % Get last modifed date of mex file
        file_info   = dir(mex_path);
        date_mex    = file_info.datenum;
        % Path to the source file
        source_path = ['mex/source/',mex_name{i},'.c'];
        % Get last modifed date of source file
        file_info   = dir(source_path);
        date_source = file_info.datenum;
        % Path to the header file
        header_path = ['mex/source/',mex_name{i},'.h'];
        file_info   = dir(header_path);
        date_header = file_info.datenum;
        % Path to the bin file
        bin_path    = ['mex/bin/',mex_name{i},'_mex.',mex_ext];
        % Check if the bin file already exists
        if exist(bin_path,'file')~=0
            file_info = dir(bin_path);
            date_bin  = file_info.datenum;
        else
            date_bin  = -1;
        end
        % Build only if required
        if date_bin<date_mex || date_bin<date_source || date_bin<date_header
            if ispc
                mex('-R2018a',mex_path,'-outdir','mex/bin',source_path,...
                    ['-I',FFTWpath],['-L',FFTWpath],...
                    '-lmwblas','-lmwlapack',['-l',FFTWlibfile(1:end-4)]);
            elseif ismac
                 mex('-R2018a',mex_path,'-outdir','mex/bin',source_path,...
                     [FFTWpath,FFTWlibfile],['-I',FFTWheader],...
                     '-lmwblas','-lmwlapack');
            end
        else
            disp('MEX file already built');
        end
        
    end
end
 



