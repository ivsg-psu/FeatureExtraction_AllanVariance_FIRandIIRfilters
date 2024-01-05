%%%%%%%%%%%%%%% script_AVAR_effectOfIIROnAVAREquiFIRfilter.m %%%%%%%%%%%%%%
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
number_of_iterations = 100;

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]
confidence_coefficient  = 0.95;

%% Plot the results
figure(01)
clf
width = 1460; height = 448; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
for index_plot = 1:3
% IIR filter parameters
if index_plot==1
    list_of_FIR_filter_orders = [(2:4), (15:2:19)];
    % IIR filter parameters
    iir_filter_order = 3; %'n' in the function 'butter' (filter order)
    iir_cutoff_freq  = 0.25; %'Wn = fc/(fs/2)' in the function 'butter' (normalized cutoff frequency)
elseif index_plot==2
    list_of_FIR_filter_orders = [(13:15:43), (103:15:133)];
    % IIR filter parameters
    iir_filter_order = 3; %'n' in the function 'butter' (filter order)
    iir_cutoff_freq  = 0.02; %'Wn = fc/(fs/2)' in the function 'butter' (normalized cutoff frequency)
elseif index_plot==3
    list_of_FIR_filter_orders = [(22:20:62), (162:20:202)];
    % IIR filter parameters
    iir_filter_order = 5; %'n' in the function 'butter' (filter order)
    iir_cutoff_freq  = 0.02; %'Wn = fc/(fs/2)' in the function 'butter' (normalized cutoff frequency)
end % NOTE: END if statement 'index_plot'

%% Filter the test signal using IIR (butterworth) filter
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

%% Compare AVAR of FIR approximation of IIR filter (Varying Order)
number_of_FIR_filters = numel(list_of_FIR_filter_orders);
default_color_map = jet(256);
custom_color_map  = default_color_map(1:floor(256/number_of_FIR_filters):256,:);
legend_cell       = cell(number_of_FIR_filters+1,1);
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
    if i==1
        axis_position = [(75+(index_plot-1)*(40+415.584))/width, 70.515/height, 415.6/width, 307.32/height];
        subplot(1,3,index_plot)
        hold on
        grid on
        fill([list_of_correlation_intervals; list_of_correlation_intervals(end:-1:1)],...
             [lb_estimated_avar; ub_estimated_avar(end:-1:1)],'m','FaceAlpha',0.5,'EdgeColor','none')
        legend_cell{i} = ['IIR: $M =$ ' num2str(iir_filter_order) ', $\omega_{n} =$ ' num2str(iir_cutoff_freq)];
    end
    plot(list_of_correlation_intervals,calculated_avar,'Color',custom_color_map(i,:),...
         'Linewidth',1.2)
    legend_cell{i+1} = ['FIR: $p =$ ' num2str(filter_order)];
    if i==number_of_FIR_filters
        legend(legend_cell,'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
        if index_plot==1
            set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log','YScale','log','FontSize',13)
            ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
        else
            set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'YTickLabel',[],'XScale','log','YScale','log','FontSize',13)
        end
        xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',17)
        ylim([10^-8.5 10^0.5])
        ax1 = gca;
        axes('Position',axis_position,'XAxisLocation','top',...
             'xLim',ax1.XLim*sampling_interval,'XScale','log','ytick',[],...
             'xtick',[1e-1 1e1 1e3],'Color','none','Fontsize',13,'Box','off');
        ax2 = gca;
        xlabel(ax2,'Correlation Time $[s]$','Interpreter','latex','FontSize',18)
        if index_plot==1
            title('(a)','Interpreter','latex','FontSize',18)
        elseif index_plot==2
            title('(b)','Interpreter','latex','FontSize',18)
        elseif index_plot==3
            title('(c)','Interpreter','latex','FontSize',18)
        end
    end
end % NOTE: END for loop 'number_of_FIR_filters'
end
