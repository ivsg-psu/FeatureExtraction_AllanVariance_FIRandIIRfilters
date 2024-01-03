%%%%%%%%%%%%%%%%%%%%%% script_AVAR_effectOfFIROrder.m %%%%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to demonstrate the effect of FIR's order
%   on AVAR of the filter output.
%
% This script was written on 2024_01_02 by Satya Prasad
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
normalized_cutoff_frequency = 0.02;
list_fir_filter_orders      = [(2:5:30), (35:15:100)];
number_of_FIR_filters       = numel(list_fir_filter_orders);

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

%% Synthesize the input signal
white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
               sampling_frequency,number_of_time_steps); % White noise
random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
               sampling_frequency,number_of_time_steps); % Random walk
input_signal = random_walk + white_noise;

%% Estimate AVAR of input to FIR filter
avar_input = fcn_AVAR_favar([input_signal; 0],list_of_correlation_intervals);

%% Plot AVAR of I/O of FIR filter
default_color_map = jet(256);
custom_color_map  = default_color_map(1:floor(256/number_of_FIR_filters):256,:);
legend_cell       = cell(number_of_FIR_filters,1);
figure(12345)
clf
width = 540; height = 448; right = 100; bottom = 100;
axis_position = [0.1354, 70.515/height, 0.7696, 307.32/height];
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(list_of_correlation_intervals,avar_input,'k','Linewidth',1.2)
legend_cell{1} = 'In';
for index_order = 1:number_of_FIR_filters
    fir_filter_order = list_fir_filter_orders(index_order);
    fir_filter_num   = fir1(fir_filter_order,normalized_cutoff_frequency);
    
    %% Synthesize the input signal
    white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
                   sampling_frequency,number_of_time_steps+fir_filter_order); % White noise
    random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
                   sampling_frequency,number_of_time_steps+fir_filter_order); % Random walk
    random_walk  = random_walk - random_walk(fir_filter_order+1);
    input_signal = random_walk + white_noise;
    
    fir_filter_output = filter(fir_filter_num,1,input_signal);
    fir_filter_output = fir_filter_output(fir_filter_order+1:end);
    
    %% Estimate AVAR of output of FIR filter
    avar_output = fcn_AVAR_favar([fir_filter_output; 0],list_of_correlation_intervals);
    plot(list_of_correlation_intervals,avar_output,'Color',custom_color_map(index_order,:),'Linewidth',1.2)
    legend_cell{index_order+1} = ['$p =$ ' num2str(fir_filter_order)];
end % NOTE: END FOR loop 'number_of_iterations'
legend(legend_cell,'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([10^-6.5 1e0])
ax1 = gca;
axes('Position',axis_position,'XAxisLocation','top',...
     'xLim',ax1.XLim*sampling_interval,'XScale','log','ytick',[],...
     'xtick',[1e-1 1e1 1e3],'Color','none','Fontsize',13,'Box','off');
ax2 = gca;
xlabel(ax2,'Correlation Time $[s]$','Interpreter','latex','FontSize',18)
