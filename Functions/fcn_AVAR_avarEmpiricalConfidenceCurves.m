function [allan_variance,lower_confidence_curve,upper_confidence_curve] = ...
    fcn_AVAR_avarEmpiricalConfidenceCurves(data,list_of_correlation_intervals,...
    confidence_coefficient,noise_type,varargin)
%% fcn_AVAR_avarEmpiricalConfidenceCurves
%   This function computes allan variance of regularly sampled data 'data'
%   for all the correlation intervals in 'list_of_correlation_intervals'.
%   It uses a recursive algorithm, inspired from FFT, over simple averages 
%   along correlation intervals.
%
% FORMAT:
%
%   [allan_variance,lower_confidence_curve,upper_confidence_curve] = ...
%   fcn_AVAR_avarEmpiricalConfidenceCurves(data,list_of_correlation_intervals,...
%   confidence_coefficient,noise_type)
%
% INPUTS:
%
%   data: A Nx1 vector of data points. N should be of form 2^p+1 (p >= 2).
%   list_of_correlation_intervals: A Mx1 vector containing list of 
%   correlation intervals. Each interval must be of the form 2^p (p >= 1).
%   confidence_coefficient: A confidence level between 0 and 1 to estimate
%   confidence and detection curves.
%   noise_type:
%              'wn'  -> White Noise
%              'rw'  -> Random Walk
%              'fn'  -> Flicker Noise
%   varargin: figure number for debugging.
%
% OUTPUTS:
%
%   allan_variance: A Mx1 vector containing allan variance corresponding to 
%   the correlation intervals.
%   lower_confidence_curve: Lower bound on AVAR estimate with a confidence
%   of 'confidence_coefficient'.
%   upper_confidence_curve: Upper bound on AVAR estimate with a confidence
%   of 'confidence_coefficient'.
%
% EXAMPLES:
%
%   See the script:
%       script_test_fcn_AVAR_favarAndConfidenceCurves.m for a full test suite.
%
% This function was written on 2022_02_08 by Satya Prasad
% Questions or comments? szm888@psu.edu
% Updated: 2022/02/15
%% ToDO
% Add DOF for rate random walk and rate ramp

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
    if 4>nargin || 5<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check the input type and domain
    try
        fcn_AVAR_checkInputsToFunctions(data,'favar dataT');
    catch ME
        fprintf(1, '%s\n\n', ME.message)
        assert(strcmp(ME.message,...
            'The data input must be a N x 1 vector of real numbers, where N >= 5'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'favar interval');
    catch ME
        assert(strcmp(ME.message,...
            'The list_of_correlation_intervals input must be a M x 1 vector of increasing numbers of form 2^p (p >= 1)'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(confidence_coefficient,'probability');
    catch ME
        assert(strcmp(ME.message,...
            'The confidence_coefficient input must lie between 0 and 1'));
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

% Does the user want to make a plot at the end?
if 5 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

p = floor(log2(numel(data)-1));
if 2^p+1 ~= numel(data)
    data = data(1:2^p+1);
    list_of_correlation_intervals = 2.^(1:p-1)';
    warning('Data is trimmed to the nearest power of 2 before estimating AVAR using FAVAR.')
end
%% Calculate Allan Variance as well as lower/upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of data points in the INPUT data
number_of_datapoints = numel(data);

% Estimate allan variance
allan_variance = fcn_AVAR_favar(data,list_of_correlation_intervals);

% calculate degrees of freedom of the estimator using empirical equation
dof_empirical = fcn_AVAR_empiricalDOFofAVARofNoise(list_of_correlation_intervals,...
                number_of_datapoints,noise_type);

chi1_squared = icdf('Chisquare',0.5*(1-confidence_coefficient),dof_empirical);
chi2_squared = icdf('Chisquare',0.5*(1+confidence_coefficient),dof_empirical);

% calculate lower and upper confidence surfaces using eq.(31)
lower_confidence_curve = (dof_empirical.*allan_variance)./chi2_squared;
upper_confidence_curve = (dof_empirical.*allan_variance)./chi1_squared;

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
    plot(list_of_correlation_intervals,allan_variance,'r','Linewidth',1.2)
    hold on
    plot(list_of_correlation_intervals,lower_confidence_curve,'g--','Linewidth',1.2)
    plot(list_of_correlation_intervals,upper_confidence_curve,'c-.','Linewidth',1.2)
    grid on
    set(gca,'XScale','log','YScale','log','FontSize',13)
    ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
    xlabel('Correlation Interval [Number of Samples]','Interpreter','latex','FontSize',17)
    legend('Estimated AVAR','Lower Bound','Upper Bound','Interpreter','latex','FontSize',13)
end

if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end