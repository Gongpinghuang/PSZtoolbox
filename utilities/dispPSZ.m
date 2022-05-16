% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function used to plot messages in the command window with certain style.
% -------------------------------- Inputs ---------------------------------
% - messageIn:   String with the message to show.
% - level:       Hierarchy level of the message to show.
% -------------------------------------------------------------------------
% For message = 'Example', the function outputs:
%
% disp('Example',0):
%
% *************************************************************************
%                                  Example                             
% *************************************************************************
%
% disp('Example',1):
%
% -------------------------------------------------------------------------
%                                  Example                             
% -------------------------------------------------------------------------
%
% disp('Example',2):
%
% ................................ Example ................................
%
% disp('Example',3):
%
%                   ______________ Example ______________     
%
% disp('Example',4):
%
%                               -- Example --
%
% disp('Example',5):
%
%                                  Example

% -------------------------------------------------------------------------
function dispPSZ(messageIn,level)

    % Get the current size of the command window
    commandWindowSize  = get(0, 'CommandWindowSize');
    % Get the current width of the command window
    commandWindowWidth = commandWindowSize(1);
    % Number of empty spaces before message, such that it is centered
    nEmpySpace = ceil((commandWindowWidth/2)-(length(messageIn)/2));

    % For the selected level
    switch level
        case 0
              % String before the message
              str       = repmat(' ',1,nEmpySpace);
              % Upper and lower bounds
              space     = ' ';
              bound     = repmat('*',1,commandWindowWidth);
        case 1
              % String before the message
              str       = repmat(' ',1,nEmpySpace);
              % Upper and lower bounds
              space     = ' ';
              bound     = repmat('-',1,commandWindowWidth);
        case 2
              % String before the message
              str       = [repmat('.',1,nEmpySpace-1),' '];
              % Upper and lower bounds
              space     = [];
              bound     = [];
        case 3
              % String before the message
              str       = [repmat(' ',1,ceil((nEmpySpace-1)/2)),...
                           repmat('_',1,floor((nEmpySpace-1)/2)),...
                           ' '];
              % Upper and lower bounds
              space     = [];
              bound     = [];
        case 4
              % String before the message
              str       = [repmat(' ',1,nEmpySpace-3),'-- '];
              % Upper and lower bounds
              space     = [];
              bound     = [];
         case 5
              % String before the message
              str       = repmat(' ',1,nEmpySpace);
              % Upper and lower bounds
              space     = [];
              bound     = [];
        otherwise
            error('dispPSZ: Not supported level');
    end
    
    % Formatted message to display
    message   = [str,messageIn,flip(str)];
    message   = message(1:commandWindowWidth);
    
    % Display
    disp(space);
    disp(bound);
    disp(message);
    disp(bound);
end





