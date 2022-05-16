% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to initialize the execution of the PSZ simulation. Particularly,  
% this function closes all previous figures, clears the workspace, and adds  
% to the matlab path the required folders for the execution of the 
% simulator.
% -------------------------------------------------------------------------
function init()
    % Clear base work space
    evalin('base','clear all');
    % Close all figures and GUIs
    close all hidden
    % Restore matlab's path to its default value
    restoredefaultpath;
    % Add required paths
    addpath('algorithms',...
            'filterBank',...
            'mex',...
            'utilities');
    if ispc    
        if exist('mex\bin','dir')
            addpath('mex\bin');
        end  
    elseif ismac
        if exist('mex/bin','dir')
            addpath('mex/bin');
        end  
    end
    % Display 
    dispPSZ('PSZ SIMULATOR',0);
end
 

