% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Implemented by Vicent Moles-Cases at GTAC-UPV, 2022         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to resample a set of Room Impulse Responses (RIR).
% -------------------------------- Inputs --------------------------------
% - hin:    RIR between L speakers and M microphones, -> [Ih_in x L x M]  
% - fs_in:  Sampling frequency for the input RIR.
% - fs_out: Sampling frequency for the output RIR.
% -------------------------------- Outputs --------------------------------
% - hout:   Resampled RIR, -> [Ih_out x L x M]
% -------------------------------------------------------------------------
function hout = resampleRIR(h,fs_in,fs_out)

    % Types of RIR
    type = {'ctrl','val'};
    % Resample for all types of RIR
    for i=1:numel(type)
        % Make sure that the input is double
        h.(type{i}) = double(h.(type{i}));
        % Resample (if needed)
        if fs_in~=fs_out
            haux           = resample(h.(type{i})(:,:),1,fs_in/fs_out);    
            hout.(type{i}) = reshape(haux,...
                                     [],...
                                     size(h.(type{i}),2),...
                                     size(h.(type{i}),3));         
        else
            hout.(type{i}) = h.(type{i});
        end
        % Make sure the length of the RIR is even
        if mod(size(hout.(type{i}),1),2)~=0
             hout.(type{i})(end+1,:,:) = 0;
        end
    end
end




