%%%%%%%%%%%%%%%%%% script_AVAR_analyzeRWinFIRfilter_EG.m %%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to create a plot for Area of AVAR + Time 
%   domain demo
%
% Author:  Satya Prasad
% Created: 2024/01/02

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
fir_noise_model = 1;

list_of_fir_filter_orders   = [(1:4), (6:2:16), (24:8:48), 50, (56:8:96), (100:32:1000)]';
number_of_fir_filter_orders = numel(list_of_fir_filter_orders);
list_of_normalized_cutoff_frequencies   = 0.02;
number_of_normalized_cutoff_frequencies = numel(list_of_normalized_cutoff_frequencies);

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals   = 2.^(0:p-3)'; % List of correlation intervals
list_of_correlation_time        = list_of_correlation_intervals*sampling_interval;
number_of_correlation_intervals = numel(list_of_correlation_intervals);

% Noise parameters
power_spectral_density  = 5*0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

%% Estimate MSE and AVAR of error in FIR filter
% Initialize variables
calculated_AVAR = NaN(number_of_correlation_intervals,number_of_fir_filter_orders,number_of_normalized_cutoff_frequencies);

max_fir_filter_order = max(list_of_fir_filter_orders);
time_vector = sampling_interval*(0:number_of_time_steps-1)';
%% Synthesize the input signal
white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
               sampling_frequency,number_of_time_steps+max_fir_filter_order); % White noise
random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
               sampling_frequency,number_of_time_steps+max_fir_filter_order); % Random walk
random_walk  = random_walk - random_walk(max_fir_filter_order+1);
input_signal = random_walk + white_noise;

for i = 1:number_of_fir_filter_orders
    fir_filter_order = list_of_fir_filter_orders(i);
    
    for j = 1:number_of_normalized_cutoff_frequencies
        normalized_cutoff_frequency = list_of_normalized_cutoff_frequencies(j);
        % Design FIR filter and Filter the input signal
        fir_filter_num  = fir1(fir_filter_order,normalized_cutoff_frequency);
        
        % AVAR of error with Random Walk as true input
        calculated_AVAR(:,i,j) = ...
            fcn_AVAR_avarFIR(power_spectral_density,random_walk_coefficient,...
            list_of_correlation_intervals,fir_filter_order,fir_filter_num,...
            sampling_interval,fir_noise_model);
    end % NOTE: END FOR loop 'number_of_normalized_cutoff_frequencies'
end % NOTE: END FOR loop 'number_of_fir_filter_orders'

%% Area of AVAR + Time domain demo
index_freq = 1;
figure(01)
clf
width = 1056.2+10; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
axis_position = [85/width, 0.1567, 415.6/width, 0.7683];
subplot(1,2,1)
hold on
grid on
plot(list_of_fir_filter_orders,...
     0.5*sum(calculated_AVAR(1:end-1,:,index_freq)+calculated_AVAR(2:end,:,1),1),...
     'k','Linewidth',1.2)
xline(3,'--','Color',[0.8500 0.3250 0.0980],'Linewidth',1.2)
xline(132,'b','Linewidth',1.2)
xline(550,'m-.','Linewidth',1.2)
legend(['$\omega_{n} =$ ' num2str(list_of_normalized_cutoff_frequencies(index_freq))],...
       '$p=3$','$p=132$','$p=550$',...
       'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e1 1e2 1e3],'XScale','log','FontSize',13)
ylabel('Area of AVAR $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('FIR Filter Order','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)

normalized_cutoff_frequency = list_of_normalized_cutoff_frequencies(index_freq);
axis_position = [(160+415.584)/width, 0.1567, 415.6/width, 0.7683];
subplot(1,2,2)
hold on
grid on
plot(-1,-1,'Color',[0.7 0.7 0.7],'Linewidth',3)
plot(-1,-1,'.','Color',[0.8500 0.3250 0.0980],'Markersize',10)
plot(-1,-1,'b','Linewidth',1.2)
plot(-1,-1,'m','Linewidth',1.2)

% Design FIR filter and Filter the input signal
fir_filter_num  = fir1(3,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'.','Color',[0.8500 0.3250 0.0980],'Markersize',4)

plot(time_vector,random_walk(end-number_of_time_steps+1:end),'Color',[0.7 0.7 0.7],'Linewidth',3)

% Design FIR filter and Filter the input signal
fir_filter_num  = fir1(132,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'b','Markersize',1.2)

% Design FIR filter and Filter the input signal
fir_filter_num  = fir1(550,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'m','Linewidth',1.2)

legend('Reference Input','$p=3$','$p=132$','$p=550$','NumColumns',2,...
       'Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'FontSize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
xlim([0 100])
