%%%%%%%%%%%%%%%%%%%% script_AVAR_findAVAREquiFIR2IIR.m %%%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to find a AVAR equivalent FIR filter of 
%   an IIR filter using relative difference in AVAR.
% 
% Author:  Satya Prasad
% Created: 2023/12/06

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
list_of_FIR_filter_orders = (1:100);

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

% IIR filter parameters
iir_filter_order = 3; %'n' in the function 'butter'
iir_cutoff_freq  = 0.25; %'Wn = fc/(fs/2)' in the function 'butter'

%% IIR (butterworth) filter coefficients
[A,B,C,D] = butter(iir_filter_order,iir_cutoff_freq);

%% Compare AVAR of FIR approximation of IIR filter (Varying Order)
number_of_FIR_filters = numel(list_of_FIR_filter_orders);
matrix_AVAR           = NaN(p-2,number_of_FIR_filters);
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

%% Demonstrate relative distance between AVAR curves
figure(01)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(log2(list_of_correlation_intervals),matrix_AVAR(:,2),'k','Linewidth',1.2)
plot(log2(list_of_correlation_intervals),matrix_AVAR(:,2),'k*','Markersize',8)
plot(log2(list_of_correlation_intervals),matrix_AVAR(:,number_of_FIR_filters),'r','Linewidth',1.2)
plot(log2(list_of_correlation_intervals),matrix_AVAR(:,number_of_FIR_filters),'r*','Markersize',8)
xt = 0:4:p-3;
xticks(xt)
xticklabels(cellstr(num2str(xt(:),'2^{%d}')))
set(gca,'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([10^-4.5 1e0])
xlim([0 p-3])

%% Estimate AVAR relative distance
AVAR_relative_diff = sqrt(sum((matrix_AVAR(:,2:number_of_FIR_filters)./matrix_AVAR(:,1:number_of_FIR_filters-1) - 1).^2,1));

figure(02)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(list_of_FIR_filter_orders(1:end-1),AVAR_relative_diff,'b','Linewidth',1.2)
set(gca,'YScale','log','FontSize',13)
ylabel('AVAR Relative Distance','Interpreter','latex','FontSize',18)
xlabel('FIR Filter Order','Interpreter','latex','FontSize',18)
title(['$M =$ ' num2str(iir_filter_order) ', $\omega_{n} =$ ' num2str(iir_cutoff_freq)],...
      'Interpreter','latex','FontSize',18)
ylim([1e-20 1e5])
