function degrees_of_freedom = fcn_AVAR_empiricalDOFofAVARofNoise(...
                              list_of_correlation_intervals,...
                              number_of_time_steps,noise_type)
%% fcn_AVAR_empiricalDOFofAVARofNoise
%   This function calculates Degrees of Freedom of AVAR of a signal, based 
%   on the noise type.
% 
% FORMAT:
%   degrees_of_freedom = fcn_AVAR_empiricalDOFofAVARofNoise(...
%                        list_of_correlation_intervals,...
%                        number_of_time_steps,noise_type)
% 
% INPUTS:
%   list_of_correlation_intervals: A Mx1 vector containing list of 
%   correlation intervals.
%   number_of_time_steps: Length of the signal.
%   noise_type: 
%              'wn'  -> White Noise
%              'rw'  -> Random Walk
%              'fn'  -> Flicker Noise
% 
% OUTPUTS:
%   degrees_of_freedom: Degrees of Freedom of AVAR of a signal, based 
%   on the noise type.
% 
% Author:  Satya Prasad
% Created: 2022/06/08

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
    if 3~=nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    try
        fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'avar interval');
    catch ME
        assert(strcmp(ME.message,...
            'The list_of_correlation_intervals input must be a M x 1 vector of increasing positive integers'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(number_of_time_steps,'positive integer');
    catch ME
        assert(strcmp(ME.message,...
            'The number_of_time_steps input must be a positive integer'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(noise_type,'string');
    catch ME
        assert(strcmp(ME.message,...
            'The noise_type input must be string'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
end

%% Calculate DOF of AVAR based on noise type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(noise_type,'wn')
    % calculate degrees of freedom of the estimator using empirial equation
    % (white noise)
    degrees_of_freedom = (3*(number_of_time_steps-1)./(2*list_of_correlation_intervals)-...
                         2*(number_of_time_steps-2)/number_of_time_steps).*...
                         (4*list_of_correlation_intervals.^2)./(4*list_of_correlation_intervals.^2+5);
elseif strcmp(noise_type,'rw')
    % calculate degrees of freedom of the estimator using empirial equation
    % (random walk)
    degrees_of_freedom = ((number_of_time_steps-2)./list_of_correlation_intervals).*...
                         (((number_of_time_steps-1)^2-...
                         3*list_of_correlation_intervals*(number_of_time_steps-1)+...
                         4*list_of_correlation_intervals.^2)/...
                         (number_of_time_steps-3)^2);
elseif strcmp(noise_type,'fn')
    % calculate degrees of freedom of the estimator using empirial equation
    % (flicker noise)
    degrees_of_freedom = (5*number_of_time_steps^2)./...
                         (4*list_of_correlation_intervals.*...
                         (number_of_time_steps+3*list_of_correlation_intervals));
    degrees_of_freedom(1==list_of_correlation_intervals) = ...
        2*(number_of_time_steps-2)/(2.3*number_of_time_steps-4.9);
end

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