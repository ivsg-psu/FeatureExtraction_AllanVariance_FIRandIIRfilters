function allan_variance = fcn_AVAR_favar(data,list_of_correlation_intervals,...
                          varargin)
%% fcn_AVAR_favar
%   This function computes allan variance of regularly sampled data 'data'
%   for all the correlation intervals in 'list_of_correlation_intervals'.
%   It uses a recursive algorithm, inspired from FFT, over simple averages 
%   along correlation intervals.
%
% FORMAT:
%   allan_variance = fcn_AVAR_favar(data,list_of_correlation_intervals)
%
% INPUTS:
%   data: A Nx1 vector of data points. N should be of form 2^p+1 (p >= 2).
%   list_of_correlation_intervals: A Mx1 vector containing list of 
%   correlation intervals. Each interval must be of the form 2^p (p >= 1).
%   varargin: figure number for debugging.
%
% OUTPUTS:
%   allan_variance: A Mx1 vector containing allan variance corresponding to 
%   the correlation intervals.
%
% EXAMPLES:
%   See the script:
%       script_test_fcn_AVAR_favar.m and 
%       script_compare_fcn_AVAR_favar_with_fcn_AVAR_avar.m 
%       for a full test suite.
%
% This function was written on 2021_05_15 by Satya Prasad
% Questions or comments? szm888@psu.edu
% Updated: 2022/02/15

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
    if 2>nargin || 3<nargin
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
            'The list_of_correlation_intervals input must be a M x 1 vector of increasing numbers of form 2^p (p >= 0)'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
end

% Does the user want to make a plot at the end?
if 3 == nargin
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
    list_of_correlation_intervals = 2.^(0:p-2)';
    warning('Data is trimmed to the nearest power of 2 before estimating AVAR using FAVAR.')
end
%% Calculate Allan Variance using FAVAR
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
number_of_datapoints            = numel(data);
% number of correlation intervals
number_of_correlation_intervals = numel(list_of_correlation_intervals);

vector_of_means = data(1:number_of_datapoints-1); % intialize vector of means with the data
length_of_vector_of_means = numel(vector_of_means); % length of vector_of_means

% initialize variable to store allan variance
allan_variance = nan(number_of_correlation_intervals, 1);
for i = 1:number_of_correlation_intervals  % loop over the list of correlation_intervals
    correlation_interval = list_of_correlation_intervals(i); % correlation_interval
    
    % Recurring step
    if 1~=correlation_interval
        vector_of_means = 0.5*(vector_of_means(1:(length_of_vector_of_means-correlation_interval/2))+...
                               vector_of_means((1+correlation_interval/2):length_of_vector_of_means));
        length_of_vector_of_means = numel(vector_of_means); % length of vector_of_means
    end
    vector_of_means_front = vector_of_means((1+correlation_interval):length_of_vector_of_means);
    vector_of_means_back  = vector_of_means(1:(length_of_vector_of_means-correlation_interval));
    
    allan_variance_sum = sum((vector_of_means_front-vector_of_means_back).^2);
    
    % write Allan Variance to the output
    allan_variance(i) = 0.5*allan_variance_sum/...
                        (number_of_datapoints-2*correlation_interval);
end % END: For loop over correlation_intervals

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
    figure_title = inputname(1);
    figure_title('_' == figure_title) = ' ';
    figure(fig_num)
    clf
    width = 540; height = 400; right = 100; bottom = 400;
    set(gcf, 'position', [right, bottom, width, height])
    plot(list_of_correlation_intervals,allan_variance,'b','Linewidth',1.2)
    grid on
    set(gca,'XScale','log','YScale','log','FontSize',13)
    ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
    xlabel('Correlation Interval [Number of Samples]','Interpreter','latex','FontSize',18)
    title(figure_title,'Interpreter','latex','FontSize',18)    
end

if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end