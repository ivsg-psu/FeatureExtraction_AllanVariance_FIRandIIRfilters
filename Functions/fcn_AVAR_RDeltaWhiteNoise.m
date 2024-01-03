function R_delta = fcn_AVAR_RDeltaWhiteNoise(power_spectral_density,...
                   list_of_correlation_intervals,lag,sampling_interval)
%% fcn_AVAR_RDeltaWhiteNoise
%   This function computes auto-correlation of the difference between 
%   successive non-overlapping averages at all the correlation intervals in
%   'list_of_correlation_intervals' for white noise.
%
% FORMAT:
%   R_delta = fcn_AVAR_RDeltaWhiteNoise(power_spectral_density,...
%             list_of_correlation_intervals,lag,sampling_interval)
%
% INPUTS:
%   power_spectral_density: Power spectral density of white noise [unit^2 s].
%   list_of_correlation_intervals: A M x 1 vector containing list of 
%   correlation intervals.
%   lag: Auto-correlation lag.
%   sampling_interval: Sampling interval [s].
%
% OUTPUTS:
%   R_delta: A M x 1 vector containing auto-correlation of the difference 
%   between successive non-overlapping averages.
%
% This function was written on 2022_09_10 by Satya Prasad
% Questions or comments? szm888@psu.edu

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs
    % Are there the right number of inputs?
    if 4~=nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    try
        fcn_AVAR_checkInputsToFunctions(power_spectral_density,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The power_spectral_density input must be a positive number'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'avar interval');
    catch ME
        assert(strcmp(ME.message,...
            'The list_of_correlation_intervals input must be a M x 1 vector of increasing positive integers'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(lag,'non negative integer');
    catch ME
        assert(strcmp(ME.message,...
            'The lag input must be a non-negative integer'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(sampling_interval,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The sampling_interval input must be a positive number'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
end

%% Calculate auto-correlation of the difference between successive 
% non-overlapping averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_delta = NaN(numel(list_of_correlation_intervals),1);
for i = 1:numel(list_of_correlation_intervals)
    correlation_interval = list_of_correlation_intervals(i);
    if lag < correlation_interval
        R_delta(i) = power_spectral_density*(2*correlation_interval - 3*lag)/...
                     (sampling_interval*correlation_interval^2);
    elseif lag < 2*correlation_interval
        R_delta(i) = power_spectral_density*(lag - 2*correlation_interval)/...
                     (sampling_interval*correlation_interval^2);
    else
        R_delta(i) = 0;
    end % NOTE: END IF statement
end % NOTE: END FOR loop 'numel(list_of_correlation_intervals)'

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end