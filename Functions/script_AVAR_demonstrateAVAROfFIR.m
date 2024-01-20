%%%%%%%%%%%%%%%%%%% script_AVAR_demonstrateAVAROfFIR.m %%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to demonstrate FIR filter and AVAR of
%   I/O of FIR filter.
%
% This script was written on 2023_09_19 by Satya Prasad
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
number_of_iterations = 100;
% fir_noise_model -> 0: AVAR of output of FIR filter
% fir_noise_model -> 1: Random walk corrupted by white noise
% fir_noise_model -> 2: Constant corrupted by random walk and white noise
fir_noise_model = 0;

% Design an FIR filter
fir_filter_order            = 50;
normalized_cutoff_frequency = 0.25;
fir_filter_num              = fir1(fir_filter_order,normalized_cutoff_frequency);

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

%%% Loop over 'number_of_iterations'
matrix_fir_avar = NaN(p-2,number_of_iterations);
for index_iter = 1:number_of_iterations
    %% Synthesize the input signal
    white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
                   sampling_frequency,number_of_time_steps+fir_filter_order); % White noise
    random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
                   sampling_frequency,number_of_time_steps+fir_filter_order); % Random walk
    random_walk  = random_walk - random_walk(fir_filter_order+1);
    input_signal = random_walk + white_noise;
    
    %% Estimate AVAR of FIR output
    fir_output_signal = filter(fir_filter_num,1,input_signal);
    fir_output_signal = fir_output_signal(fir_filter_order+1:end);
    matrix_fir_avar(:,index_iter) = ...
        fcn_AVAR_favar([fir_output_signal; 0],list_of_correlation_intervals);
    
    if index_iter==number_of_iterations
        %% Demonstrate FIR filter
        input_length = 20;
        input_data   = input_signal(1:fir_filter_order+input_length);
        fir_output_signal = filter(fir_filter_num,1,input_data);
        fir_output_signal = fir_output_signal(fir_filter_order+1:end);
        
        % Plot the input and its FIR output
        yLim_min   = min([input_data(fir_filter_order+1:end), fir_output_signal],[],'all');
        yLim_max   = max([input_data(fir_filter_order+1:end), fir_output_signal],[],'all');
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
        plot((0:input_length-1)',input_data(fir_filter_order+1:end),'k.','Markersize',13)
        plot((0:input_length-1)',fir_output_signal,'m*','Markersize',8)
        set(gca,'xtick',[3 7 11 15 19],'Fontsize',13)
        legend('Input','Output','Location','best','Interpreter','latex','Fontsize',13)
        ylabel('Amplitude $[Unit]$','Interpreter','latex','Fontsize',18)
        xlabel('Time Step $[Number \: of \; Samples]$','Interpreter','latex','Fontsize',18)
        title(['$p =$ ' num2str(fir_filter_order) ', $\omega_{n} =$ ' num2str(normalized_cutoff_frequency)],...
              'Interpreter','latex','Fontsize',18)
        ylim([yLim_min-yLim_Range yLim_max+yLim_Range])
        xlim([xLim_min-xLim_Range xLim_max+xLim_Range])
        
        %% Demonstrate AVAR of FIR filter
        %%% Estimate AVAR of input to MA filter
        [avar_input,avar_input_lb,avar_input_ub] = ...
            fcn_AVAR_avarEmpiricalConfidenceCurves([input_signal(1:number_of_time_steps); 0],...
            list_of_correlation_intervals,confidence_coefficient,noise_type);
        
        %%% Estimate AVAR of output of MA filter with confidence bounds
        avar_output = matrix_fir_avar(:,1);
        % calculate degrees of freedom of the estimator using empirical equation
        dof_avar     = 2*(mean(matrix_fir_avar,2).^2)./var(matrix_fir_avar,1,2);
        chi1_squared = icdf('Chisquare',0.5*(1-confidence_coefficient),dof_avar);
        chi2_squared = icdf('Chisquare',0.5*(1+confidence_coefficient),dof_avar);
        % calculate lower and upper confidence surfaces using eq.(31)
        avar_output_lb = (dof_avar.*avar_output)./chi2_squared;
        avar_output_ub = (dof_avar.*avar_output)./chi1_squared;
        
        %%% Plot AVAR of I/O of FIR filter
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
        legend_cell{4} = ['Out: $p =$ ' num2str(fir_filter_order) ', $\omega_{n} =$ ' num2str(normalized_cutoff_frequency)];
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

        %% Calculated Allan Variance of FIR output
        calculated_avar = fcn_AVAR_avarFIR(power_spectral_density,random_walk_coefficient,...
            list_of_correlation_intervals,fir_filter_order,fir_filter_num,...
            sampling_interval,fir_noise_model);
        
        %%% Plot calculated AVAR against estimated AVAR
        figure(12347)
        clf
        width = 540; height = 448; right = 100; bottom = 100;
        axis_position = [0.1354, 70.515/height, 0.7696, 307.32/height];
        set(gcf, 'position', [right, bottom, width, height])
        hold on
        grid on
        plot(list_of_correlation_intervals,avar_output,'m','Linewidth',1.2)
        plot(list_of_correlation_intervals,avar_output_lb,'m--','Linewidth',1)
        plot(list_of_correlation_intervals,avar_output_ub,'m-.','Linewidth',1)
        plot(list_of_correlation_intervals,calculated_avar,'k','Linewidth',1.2)
        legend('Estimated',[num2str(100*confidence_coefficient) '$\%$ LB'],...
               [num2str(100*confidence_coefficient) '$\%$ UB'],'Calculated',...
               'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
        set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log','YScale','log','FontSize',13)
        ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
        xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
        ylim([1e-4 1e0])
        ax1 = gca;
        axes('Position',axis_position,'XAxisLocation','top',...
             'xLim',ax1.XLim*sampling_interval,'XScale','log','ytick',[],...
             'xtick',[1e-1 1e1 1e3],'Color','none','Fontsize',13,'Box','off');
        ax2 = gca;
        xlabel(ax2,'Correlation Time $[s]$','Interpreter','latex','FontSize',18)
    end % NOTE: END IF statement 'index_iter==number_of_iterations'
end % NOTE: END FOR loop 'number_of_iterations'
