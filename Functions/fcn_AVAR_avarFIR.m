function allan_variance = ...
    fcn_AVAR_avarFIR(power_spectral_density,random_walk_coefficient,...
    list_of_correlation_intervals,fir_filter_order,fir_filter_num,...
    sampling_interval,noise_model,varargin)
%% fcn_AVAR_avarMA
%   This function calculates AVAR of error in MA filter.
%
% FORMAT:
%   allan_variance = ...
%       fcn_AVAR_avarMA(power_spectral_density,random_walk_coefficient,...
%       list_of_correlation_intervals,window_length,sampling_interval,noise_model)
%
% INPUTS:
%   power_spectral_density: Power spectral density of white noise [unit^2 s].
%   random_walk_coefficient: Noise coefficient for random walk [unit/sqrt(s)].
%   list_of_correlation_intervals: A M x 1 vector containing list of 
%   correlation intervals.
%   window_length: Window length of MA filter.
%   sampling_interval: Sampling interval [s].
%   noise_model:
%       0-> AVAR of output of MA filter
%       1-> Random walk corrupted by white noise
%       2-> Constant corrupted by random walk and white noise
%   varargin: figure number for debugging.
%
% OUTPUTS:
%   allan_variance: A M x 1 vector containing allan variance corresponding 
%   to the correlation intervals.
%
% EXAMPLES:
%   See the script:
%       script_test_fcn_AVAR_avarMA.m for a full test suite.
%
% This script was written on 2023_09_08 by Satya Prasad
% Questions or comments? szm888@psu.edu

flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plot  = 0; % Flag to plot the results
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
    if 7>nargin || 8<nargin
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
        fcn_AVAR_checkInputsToFunctions(random_walk_coefficient,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The random_walk_coefficient input must be a positive number'));
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
        fcn_AVAR_checkInputsToFunctions(fir_filter_order,'positive integer');
    catch ME
        assert(strcmp(ME.message,...
            'The fir_filter_order input must be a positive integer'));
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

% Does the user want to make a plot at the end?
if 8 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Calculate AVAR of MA filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0==noise_model
    % AVAR of output in FIR filter error with Random Walk corrupted by white
    % nosie
    R_delta = fcn_AVAR_RDeltaRandomWalk(random_walk_coefficient,...
              list_of_correlation_intervals,0,sampling_interval) + ...
              fcn_AVAR_RDeltaWhiteNoise(power_spectral_density,...
              list_of_correlation_intervals,0,sampling_interval);
    sum_R_delta = 0.5*sum(fir_filter_num.^2)*R_delta;
    for q = 1:fir_filter_order
        R_delta = fcn_AVAR_RDeltaRandomWalk(random_walk_coefficient,...
                  list_of_correlation_intervals,q,sampling_interval) + ...
                  fcn_AVAR_RDeltaWhiteNoise(power_spectral_density,...
                  list_of_correlation_intervals,q,sampling_interval);
        sum_R_delta = sum_R_delta+...
            sum(fir_filter_num(1:fir_filter_order+1-q).*fir_filter_num(1+q:fir_filter_order+1))*R_delta;
    end % NOTE: END FOR loop 'fir_filter_order'
    allan_variance = sum_R_delta;
elseif 1==noise_model
    % Calculated Allan Variance
    rw_fir_filter_num = fir_filter_num-[1, zeros(1,fir_filter_order)];
    R_delta_WN = fcn_AVAR_RDeltaWhiteNoise(power_spectral_density,...
                 list_of_correlation_intervals,0,sampling_interval);
    R_delta_RW = fcn_AVAR_RDeltaRandomWalk(random_walk_coefficient,...
                 list_of_correlation_intervals,0,sampling_interval);
    sum_R_delta = 0.5*sum(fir_filter_num.^2)*R_delta_WN + ...
                  0.5*sum(rw_fir_filter_num.^2)*R_delta_RW;
    for q = 1:fir_filter_order
        R_delta_WN = fcn_AVAR_RDeltaWhiteNoise(power_spectral_density,...
                     list_of_correlation_intervals,q,sampling_interval);
        R_delta_RW = fcn_AVAR_RDeltaRandomWalk(random_walk_coefficient,...
                     list_of_correlation_intervals,q,sampling_interval);
        sum_R_delta = sum_R_delta+...
            sum(fir_filter_num(1:fir_filter_order+1-q).*fir_filter_num(1+q:fir_filter_order+1))*R_delta_WN + ...
            sum(rw_fir_filter_num(1:fir_filter_order+1-q).*rw_fir_filter_num(1+q:fir_filter_order+1))*R_delta_RW;
    end % NOTE: END FOR loop 'fir_filter_order'
    allan_variance = sum_R_delta;
end % NOTE: END IF statement

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
if flag_do_plot
    figure(fig_num)
    clf
    width = 540; height = 400; right = 100; bottom = 400;
    set(gcf, 'position', [right, bottom, width, height])
    plot(list_of_correlation_intervals,allan_variance,'b','Linewidth',1.2)
    grid on
    set(gca,'XScale','log','YScale','log','FontSize',13)
    ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
    xlabel('Correlation Interval [Number of Samples]','Interpreter','latex','FontSize',18)
end

if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end