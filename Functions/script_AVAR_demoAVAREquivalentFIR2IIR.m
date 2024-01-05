%%%%%%%%%%%%%%%%% script_AVAR_demoAVAREquivalentFIR2IIR.m %%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to show the effect of IIR filter order
%   and normalized cutoff frequency on AVAR equivalent FIR filter order.
% 
% Author:  Satya Prasad
% Created: 2024/01/05

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
index_iir_filter = 1;
number_of_iterations = 100;
epsilon = 0.01;

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]
confidence_coefficient  = 0.95;

% IIR filter parameters
if index_iir_filter==1
    list_of_FIR_filter_orders = [(2:4), (15:2:19)];
    % IIR filter parameters
    iir_filter_order = 3; %'n' in the function 'butter' (filter order)
    iir_cutoff_freq  = 0.25; %'Wn = fc/(fs/2)' in the function 'butter' (normalized cutoff frequency)
elseif index_iir_filter==2
    list_of_FIR_filter_orders = [(13:15:43), (103:15:133)];
    % IIR filter parameters
    iir_filter_order = 3; %'n' in the function 'butter' (filter order)
    iir_cutoff_freq  = 0.02; %'Wn = fc/(fs/2)' in the function 'butter' (normalized cutoff frequency)
elseif index_iir_filter==3
    list_of_FIR_filter_orders = [(22:20:62), (162:20:202)];
    % IIR filter parameters
    iir_filter_order = 5; %'n' in the function 'butter' (filter order)
    iir_cutoff_freq  = 0.02; %'Wn = fc/(fs/2)' in the function 'butter' (normalized cutoff frequency)
end % NOTE: END if statement 'index_iir_filter'
number_of_FIR_filters = numel(list_of_FIR_filter_orders);

%% IIR (butterworth) filter coefficients
[A,B,C,D] = butter(iir_filter_order,iir_cutoff_freq);
[b,a]     = butter(iir_filter_order,iir_cutoff_freq);

matrix_AVAR = NaN(p-2,number_of_iterations);
for i = 1:number_of_iterations
    %% Synthesize test signal
    white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
                   sampling_frequency,number_of_time_steps+iir_filter_order); % White noise
    random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
                   sampling_frequency,number_of_time_steps+iir_filter_order); % Random walk
    random_walk  = random_walk - random_walk(iir_filter_order+1);
    input_signal = random_walk + white_noise; % Test signal
    % Filter signal using IIR filter
    filtered_test_signal = filter(b,a,input_signal);
    filtered_test_signal = filtered_test_signal(iir_filter_order+1:end);
    
    % Estimate AVAR
    matrix_AVAR(:,i) = fcn_AVAR_favar([filtered_test_signal; 0],...
                       list_of_correlation_intervals);
end % NOTE: END for loop 'number_of_iterations'
estimated_avar_output = matrix_AVAR(:,1);
% Calculate degrees of freedom of the estimator
degrees_of_freedom    = 2*(mean(matrix_AVAR,2).^2)./var(matrix_AVAR,1,2);
chi1_squared = icdf('Chisquare',0.5*(1-confidence_coefficient),degrees_of_freedom);
chi2_squared = icdf('Chisquare',0.5*(1+confidence_coefficient),degrees_of_freedom);
% Calculate lower and upper confidence surfaces using eq.(31)
lb_estimated_avar = (degrees_of_freedom.*estimated_avar_output)./chi2_squared;
ub_estimated_avar = (degrees_of_freedom.*estimated_avar_output)./chi1_squared;

%% Plot the results
default_color_map = jet(256);
custom_color_map  = default_color_map(1:floor(256/number_of_FIR_filters):256,:);
legend_cell       = cell(number_of_FIR_filters+1,1);
figure(01)
clf
width = 1501.8+42.5; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
%% Compare AVAR of FIR approximation of IIR filter (Varying Order)
subplot(1,3,1)
axis_position = [85/width, 0.1567, 415.6/width, 0.7683];
hold on
grid on
fill([list_of_correlation_intervals; list_of_correlation_intervals(end:-1:1)],...
     [lb_estimated_avar; ub_estimated_avar(end:-1:1)],'m','FaceAlpha',0.5,'EdgeColor','none')
legend_cell{1} = ['IIR: $M =$ ' num2str(iir_filter_order) ', $\omega_{n} =$ ' num2str(iir_cutoff_freq)];
for i = 1:number_of_FIR_filters
    % Calculate FIR filter weights using state matrices and filter order
    filter_order      = list_of_FIR_filter_orders(i);
    filter_weights    = NaN(filter_order+1,1);
    filter_weights(1) = D;
    for j = 2:filter_order+1
        filter_weights(j) = C*(A^(j-2))*B;
    end % NOTE: END for loop 'filter_order+1'
    
    % Calculated Allan Variance
    calculated_avar = fcn_AVAR_avarFIR(power_spectral_density,random_walk_coefficient,...
        list_of_correlation_intervals,filter_order,filter_weights,sampling_interval,0);
    
    % Plot AVAR curves
    plot(list_of_correlation_intervals,calculated_avar,'Color',custom_color_map(i,:),...
         'Linewidth',1.2)
    legend_cell{i+1} = ['FIR: $p =$ ' num2str(filter_order)];
end % NOTE: END for loop 'number_of_FIR_filters'
legend(legend_cell,'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log','YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',17)
title('(a)','Interpreter','latex','FontSize',18)
ylim([10^-8.5 1e0])

%% Compare AVAR of FIR approximation of IIR filter (Varying Order)
list_of_FIR_filter_orders = (1:200);
number_of_FIR_filters     = numel(list_of_FIR_filter_orders);
matrix_AVAR               = NaN(p-2,number_of_FIR_filters);
for index_fir = 1:number_of_FIR_filters
    % Calculate FIR filter weights using state matrices and filter order
    filter_order      = list_of_FIR_filter_orders(index_fir);
    filter_weights    = NaN(filter_order+1,1);
    filter_weights(1) = D;
    for i = 2:filter_order+1
        filter_weights(i) = C*(A^(i-2))*B;
    end % NOTE: END FOR loop 'filter_order+1'
    
    % Calculated Allan Variance
    matrix_AVAR(:,index_fir) = fcn_AVAR_avarFIR(power_spectral_density,random_walk_coefficient,...
        list_of_correlation_intervals,filter_order,filter_weights,sampling_interval,0);
end % NOTE: END FOR loop 'number_of_FIR_filters'

%% Estimate AVAR relative distance
AVAR_relative_dist = sqrt(sum((matrix_AVAR(:,2:number_of_FIR_filters)./matrix_AVAR(:,1:number_of_FIR_filters-1) - 1).^2,1));

axis_position = [(85+85+415.584)/width, 0.1567, 415.6/width, 0.7683];
subplot(1,3,2)
hold on
grid on
plot(list_of_FIR_filter_orders(1:end-1),AVAR_relative_dist,'b','Linewidth',1.2)
set(gca,'Position',axis_position,'YScale','log','FontSize',13)
ylabel('AVAR Relative Distance','Interpreter','latex','FontSize',18)
xlabel('FIR Filter Order','Interpreter','latex','FontSize',18)
title('(b)','Interpreter','latex','FontSize',18)
ylim([1e-6 1e2])

%% Time domain plots for both IIR and AVAR equivalent FIR
%%% Synthesize test signal
time_vector  = (0:number_of_time_steps-1)'/sampling_frequency;
white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
               sampling_frequency,number_of_time_steps+list_of_FIR_filter_orders(end)); % White noise
random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
               sampling_frequency,number_of_time_steps+list_of_FIR_filter_orders(end)); % Random walk
random_walk  = random_walk - random_walk(list_of_FIR_filter_orders(end)+1);
input_signal = random_walk + white_noise; % Test signal
%%% Filter signal using IIR filter
iir_test_signal = filter(b,a,input_signal);
iir_test_signal = iir_test_signal(end-number_of_time_steps+1:end);

%%% Calculate FIR filter weights using state matrices and filter order
fir_filter_order  = find([0, diff(AVAR_relative_dist < epsilon)],1,'last');
filter_weights    = NaN(fir_filter_order+1,1);
filter_weights(1) = D;
for i = 2:fir_filter_order+1
    filter_weights(i) = C*(A^(i-2))*B;
end % NOTE: END FOR loop 'filter_order+1'
%%% Filter signal using AVAR equivalent FIR filter
fir_test_signal = filter(filter_weights,1,input_signal);
fir_test_signal = fir_test_signal(end-number_of_time_steps+1:end);

axis_position = [(85+2*(85+415.584))/width, 0.1567, 415.6/width, 0.7683];
subplot(1,3,3)
hold on
grid on
plot(time_vector,iir_test_signal,'r','Linewidth',1)
plot(time_vector,fir_test_signal,'g--','Linewidth',1)
legend('IIR','Equi FIR','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'FontSize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
title('(c)','Interpreter','latex','FontSize',18)
xlim([0 100])
