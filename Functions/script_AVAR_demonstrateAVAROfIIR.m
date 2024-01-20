%%%%%%%%%%%%%%%%%%%% script_AVAR_demonstrateAVAROfIIR.m %%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to demonstrate IIR filter and AVAR of I/O
%   of IIR filter.
%
% This script was written on 2023_12_06 by Satya Prasad
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
number_of_iterations = 500;
iir_filter_order = 3;
normalized_cutoff_frequency = 0.25;

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]
confidence_coefficient  = 0.95;
noise_type              = 'rw';

%% Synthesize the input signal
white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
               sampling_frequency,number_of_time_steps+iir_filter_order); % White noise
random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
               sampling_frequency,number_of_time_steps+iir_filter_order); % Random walk
random_walk  = random_walk - random_walk(iir_filter_order+1);
input_signal = random_walk + white_noise;

% IIR filter coefficients
[iir_filter_num,iir_filter_den] = butter(iir_filter_order,normalized_cutoff_frequency);

%% Demonstrate IIR filter
input_length = 20;
input_data   = input_signal(1:input_length+iir_filter_order);
iir_filter_output = filter(iir_filter_num,iir_filter_den,input_data);
iir_filter_output = iir_filter_output(iir_filter_order+1:end);

% Plot the signal and its IIR filtered output
yLim_min   = min([input_data(iir_filter_order+1:end), iir_filter_output],[],'all');
yLim_max   = max([input_data(iir_filter_order+1:end), iir_filter_output],[],'all');
yLim_Range = 0.1*(yLim_max-yLim_min);
xLim_min   = 0;
xLim_max   = input_length-1;
xLim_Range = 0.1*(xLim_max-xLim_min);
figure(12345)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot((0:input_length-1)',input_data(iir_filter_order+1:end),'k.','Markersize',13)
plot((0:input_length-1)',iir_filter_output,'m*','Markersize',8)
set(gca,'xtick',[3 7 11 15 19],'Fontsize',13)
legend('Input','Output','Location','best','Interpreter','latex','Fontsize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','Fontsize',18)
xlabel('Time Step $[Number \: of \: Samples]$','Interpreter','latex','Fontsize',18)
title(['$M =$ ' num2str(iir_filter_order) ', $\omega_{n} =$ ' num2str(normalized_cutoff_frequency)],...
    'Interpreter','latex','Fontsize',18)
ylim([yLim_min-yLim_Range yLim_max+yLim_Range])
xlim([xLim_min-xLim_Range xLim_max+xLim_Range])

%% Demonstrate AVAR of IIR
%%% Estimate AVAR of input to IIR filter
[avar_input,avar_input_lb,avar_input_ub] = ...
    fcn_AVAR_avarEmpiricalConfidenceCurves([input_signal(1:number_of_time_steps); 0],...
    list_of_correlation_intervals,confidence_coefficient,noise_type);

%%% Estimate AVAR of output of IIR filter with confidence bounds
matrix_AVAR = NaN(p-2,number_of_iterations);
for index_iter = 1:number_of_iterations
    white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
                   sampling_frequency,number_of_time_steps+iir_filter_order); % White noise
    random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
                   sampling_frequency,number_of_time_steps+iir_filter_order); % Random walk
    random_walk  = random_walk - random_walk(iir_filter_order+1);
    input_signal = random_walk + white_noise;
    
    iir_filter_output = filter(iir_filter_num,iir_filter_den,input_signal);
    iir_filter_output = iir_filter_output(iir_filter_order+1:end);
    
    matrix_AVAR(:,index_iter) = fcn_AVAR_favar([iir_filter_output; 0],...
                                list_of_correlation_intervals);
end % NOTE: END FOR loop 'number_of_iterations'
avar_output  = matrix_AVAR(:,1);
% calculate degrees of freedom of the estimator
dof_avar     = 2*(mean(matrix_AVAR,2).^2)./var(matrix_AVAR,1,2);
chi1_squared = icdf('Chisquare',0.5*(1-confidence_coefficient),dof_avar);
chi2_squared = icdf('Chisquare',0.5*(1+confidence_coefficient),dof_avar);
% calculate lower and upper confidence surfaces using eq.(31)
avar_output_lb = (dof_avar.*avar_output)./chi2_squared;
avar_output_ub = (dof_avar.*avar_output)./chi1_squared;

%%% Plot AVAR of I/O of IIR filter
figure(12346)
clf
width = 540; height = 448; right = 100; bottom = 100;
axis_position = [0.1354, 70.515/height, 0.7696, 307.32/height];
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(list_of_correlation_intervals,avar_input,'k','Linewidth',1.2)
plot(list_of_correlation_intervals,avar_input_lb,'k--','Linewidth',1)
plot(list_of_correlation_intervals,avar_input_ub,'k-.','Linewidth',1)
plot(list_of_correlation_intervals,avar_output,'m','Linewidth',1.2)
plot(list_of_correlation_intervals,avar_output_lb,'m--','Linewidth',1)
plot(list_of_correlation_intervals,avar_output_ub,'m-.','Linewidth',1)
legend_cell = cell(6,1);
legend_cell{1} = 'In';
legend_cell{2} = ['In: ' num2str(100*confidence_coefficient) '$\%$ LB'];
legend_cell{3} = ['In: ' num2str(100*confidence_coefficient) '$\%$ UB'];
legend_cell{4} = ['Out: $M =$ ' num2str(iir_filter_order) ', $\omega_{n} =$ ' num2str(normalized_cutoff_frequency)];
legend_cell{5} = ['Out: ' num2str(100*confidence_coefficient) '$\%$ LB'];
legend_cell{6} = ['Out: ' num2str(100*confidence_coefficient) '$\%$ UB'];
legend(legend_cell,'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])
ax1 = gca;
axes('Position',axis_position,'XAxisLocation','top',...
     'xLim',ax1.XLim*sampling_interval,'XScale','log','ytick',[],...
     'xtick',[1e-1 1e1 1e3],'Color','none','Fontsize',13,'Box','off');
ax2 = gca;
xlabel(ax2,'Correlation Time $[s]$','Interpreter','latex','FontSize',18)
