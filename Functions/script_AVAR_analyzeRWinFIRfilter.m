%%%%%%%%%%%%%%%%%%%% script_AVAR_analyzeRWinFIRfilter.m %%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to analyze FIR filter with Random Walk
%   input corrupted by White Noise using Allan Variance and MSE
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
fir_noise_model      = 1;
number_of_iterations = 500;

list_of_fir_filter_orders   = [(1:4), (6:2:16), (24:8:48), 50, (56:8:96), (100:32:1000)]';
number_of_fir_filter_orders = numel(list_of_fir_filter_orders);
list_of_normalized_cutoff_frequencies   = [(0.02:0.02:0.1), (0.2:0.2:0.8)];
number_of_normalized_cutoff_frequencies = numel(list_of_normalized_cutoff_frequencies);

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals   = 2.^(0:p-3)'; % List of correlation intervals
list_of_correlation_time        = list_of_correlation_intervals*sampling_interval;
number_of_correlation_intervals = numel(list_of_correlation_intervals);
target_points = [round(number_of_time_steps/3), round(2*number_of_time_steps/3), ...
                 number_of_time_steps];
number_of_target_points = numel(target_points);

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

%% Estimate MSE and AVAR of error in FIR filter
% Initialize variables
estimated_MSE   = zeros(number_of_time_steps,number_of_fir_filter_orders,number_of_normalized_cutoff_frequencies);
estimated_AVAR  = NaN(number_of_correlation_intervals,number_of_fir_filter_orders,number_of_normalized_cutoff_frequencies);
calculated_AVAR = NaN(number_of_correlation_intervals,number_of_fir_filter_orders,number_of_normalized_cutoff_frequencies);

max_fir_filter_order = max(list_of_fir_filter_orders);
time_vector = sampling_interval*(0:number_of_time_steps-1)';
for index_iter = 1:number_of_iterations
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
        filtered_output = filter(fir_filter_num,1,...
                                 input_signal(end-number_of_time_steps-fir_filter_order+1:end));
        filtered_output = filtered_output(fir_filter_order+1:end);
        actual_error    = random_walk(end-number_of_time_steps+1:end)-filtered_output;
        
        % MSE with Random Walk as true input
        estimated_MSE(:,i,j) = estimated_MSE(:,i,j) + (actual_error.^2);
        if index_iter==1
            % AVAR of error with Random Walk as true input
            estimated_AVAR(:,i,j) = fcn_AVAR_favar([actual_error; 0],list_of_correlation_intervals);
            calculated_AVAR(:,i,j) = ...
                fcn_AVAR_avarFIR(power_spectral_density,random_walk_coefficient,...
                list_of_correlation_intervals,fir_filter_order,fir_filter_num,...
                sampling_interval,fir_noise_model);
        end % NOTE: END IF statement 'index_iter==1'
    end % NOTE: END FOR loop 'number_of_normalized_cutoff_frequencies'
end % NOTE: END FOR loop 'number_of_fir_filter_orders'
end % NOTE: END FOR loop 'number_of_iterations'
estimated_MSE = estimated_MSE/number_of_iterations;

%% MSE and AVAR of error for no filter case
NF_MSE_IN1  = (power_spectral_density/sampling_interval)*ones(number_of_time_steps,1);
NF_AVAR_IN1 = fcn_AVAR_avarWhiteNoise(power_spectral_density,list_of_correlation_time);

%% Plot the results
default_color_map = jet(256);
jf = java.text.DecimalFormat;   % To add comma to numbers
%%% To show MSE is independent of signal length
figure(01)
clf
width = 1469.3; height = 400; right = 20; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
custom_color_map = default_color_map(1:floor(256/number_of_normalized_cutoff_frequencies):256,:);
for index_target = 1:number_of_target_points
    legend_cell   = cell(number_of_normalized_cutoff_frequencies+1,1);
    axis_position = [(75+(index_target-1)*(55+415.584))/width, 0.1567, 415.6/width, 0.7683];
    subplot(1,3,index_target)
    hold on
    grid on
    yline(NF_MSE_IN1(target_points(index_target)),'m--','Linewidth',1.2)
    legend_cell{1} = 'No Filter';
    for i = 1:number_of_normalized_cutoff_frequencies
        temp_var = estimated_MSE(target_points(index_target),:,i);
        plot(list_of_fir_filter_orders,temp_var(:),'Color',custom_color_map(i,:),'Linewidth',1.2)
        legend_cell{i+1} = ['$\omega_{n} =$ ' num2str(list_of_normalized_cutoff_frequencies(i))];
    end % NOTE: END FOR loop 'number_of_normalized_cutoff_frequencies'
    set(gca,'Position',axis_position,'xtick',[1e0 1e1 1e2 1e3],'ytick',[0 0.01 0.02 0.03],'XScale','log','FontSize',13)
    if index_target==1
        legend(legend_cell,'NumColumns',3,'Location','best','Interpreter','latex','FontSize',13)
        ylabel('Mean Squared Error $[Unit^2]$','Interpreter','latex','FontSize',18)
    end % NOTE: END if statement 'index_target==1'
    xlabel('FIR Filter Order','Interpreter','latex','FontSize',18)
    title(['$k = $ ' char(jf.format(target_points(index_target)-1))],'Interpreter','latex','FontSize',18)
    ylim([0 0.03])
end % NOTE: END for loop 'number_of_target_points'

%%% MSE vs Area of AVAR
custom_color_map = default_color_map(1:floor(256/number_of_normalized_cutoff_frequencies):256,:);
legend_cell      = cell(number_of_normalized_cutoff_frequencies,1);
figure(02)
clf
width = 540; height = 770.04; right = 100; bottom = 10;
set(gcf, 'position', [right, bottom, width, height])
axis_position = [0.1354, 432.72/height, 0.7696, 307.32/height];
subplot(2,1,1)
hold on
grid on
yline(NF_MSE_IN1(target_points(end)),'m--','Linewidth',1.2)
for i = 1:number_of_normalized_cutoff_frequencies
    temp_var = estimated_MSE(target_points(end),:,i);
    plot(list_of_fir_filter_orders,temp_var(:),'Color',custom_color_map(i,:),'Linewidth',1.2)
end % NOTE: END FOR loop 'number_of_normalized_cutoff_frequencies'
set(gca,'Position',axis_position,'xtick',[1e0 1e1 1e2 1e3],'ytick',[0 0.01 0.02 0.03],'XScale','log','FontSize',13)
ylabel('Mean Squared Error $[Unit^2]$','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([0 0.03])

axis_position = [0.1354, 62.7/height, 0.7696, 307.32/height];
subplot(2,1,2)
hold on
grid on
yline(0.5*sum(NF_AVAR_IN1(1:end-1)+NF_AVAR_IN1(2:end)),'m--','Linewidth',1.2)
legend_cell{1} = 'No Filter';
for i = 1:number_of_normalized_cutoff_frequencies
    plot(list_of_fir_filter_orders,...
         0.5*sum(estimated_AVAR(1:end-1,:,i)+estimated_AVAR(2:end,:,i),1),...
         'Color',custom_color_map(i,:),'Linewidth',1.2)
    legend_cell{i+1} = ['$\omega_{n} =$ ' num2str(list_of_normalized_cutoff_frequencies(i))];
end % NOTE: END FOR loop 'number_of_normalized_cutoff_frequencies'
legend(legend_cell,'NumColumns',3,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e1 1e2 1e3],'ytick',[0 0.01 0.02 0.03 0.04],'XScale','log','FontSize',13)
ylabel('Area of AVAR $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('FIR Filter Order','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
ylim([0 0.045])

%% Area of AVAR + Time domain demo
index_freq = 1;
figure(03)
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
xline(80,'b','Linewidth',1.2)
xline(550,'m-.','Linewidth',1.2)
legend(['$\omega_{n} =$ ' num2str(list_of_normalized_cutoff_frequencies(index_freq))],...
       '$p=3$','$p=80$','$p=550$',...
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
fir_filter_num  = fir1(80,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'b','Markersize',1.2)

% Design FIR filter and Filter the input signal
fir_filter_num  = fir1(550,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'m','Linewidth',1.2)

legend('Reference Input','$p=3$','$p=80$','$p=550$','NumColumns',2,...
       'Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'FontSize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
xlim([0 100])

%% Area of AVAR + Time domain demo
index_freq = 5;
figure(04)
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
xline(40,'b','Linewidth',1.2)
xline(550,'m-.','Linewidth',1.2)
legend(['$\omega_{n} =$ ' num2str(list_of_normalized_cutoff_frequencies(index_freq))],...
       '$p=3$','$p=40$','$p=550$',...
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

% Design FIR filter and Filter the input signal
fir_filter_num  = fir1(40,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'b','Markersize',0.8)

% Design FIR filter and Filter the input signal
fir_filter_num  = fir1(550,normalized_cutoff_frequency);
filtered_output = filter(fir_filter_num,1,input_signal);
filtered_output = filtered_output(end-number_of_time_steps+1:end);
plot(time_vector,filtered_output,'m','Linewidth',0.8)

plot(time_vector,random_walk(end-number_of_time_steps+1:end),'Color',[0.7 0.7 0.7],'Linewidth',3)

legend('Reference Input','$p=3$','$p=40$','$p=550$','NumColumns',2,...
       'Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'FontSize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
xlim([0 100])
