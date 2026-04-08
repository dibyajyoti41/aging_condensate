% =============================================================================
% fitting2anda0e.m
%
% Post-analysis of the pairwise overlap function F(tau) computed by
% droplet_overlap.f90 and shellwise_overlap.f90.
%
% The overlap function for each aging time window t_w is fit to a
% bi-exponential decay with a long-time plateau:
%
%   F(tau) = A0 + A1*exp(-tau/tau1) + A2*exp(-tau/tau2)
%
% with the constraint A0 + A1 + A2 = 1, i.e. A2 = 1 - A0 - A1.
%
%   A0   -- non-ergodic parameter: the long-time plateau of F(tau)
%   A1   -- amplitude of the fast (beta) relaxation process
%   A2   -- amplitude of the slow (alpha) relaxation process
%   tau1 -- fast relaxation timescale (cage-rattling / local rearrangements)
%   tau2 -- slow relaxation timescale (cooperative/structural relaxation)
%
% The fits are performed window-by-window as a function of the waiting
% time t_w (age of the system after quench), allowing the evolution of
% relaxation timescales with system age to be tracked.
%
% Inputs:  bs1, bs2, ..., bs29  -- workspace variables, each containing
%          the F(tau) vs tau data for one waiting-time window.
% Output:  res2f matrix (rows = windows, columns = fit parameters + errors),
%          subplot panel showing all fits, and a final plot of tau1 vs t_w.
% =============================================================================

% --- Load overlap function data for all waiting-time windows ---
% Each bsN is a matrix with columns [tau_index, F(tau), ...] as written
% by droplet_overlap.f90 / shellwise_overlap.f90. Windows correspond to
% successive non-overlapping segments of the trajectory.
datasets = {bs1,bs2,bs3,bs4, bs6, bs7, bs8, bs9, bs10,bs11,bs12,bs13,bs14,bs15,bs16,bs17,bs18,bs19,bs20,bs21,bs22,bs23,bs24,bs25,bs26,bs27,bs28,bs29};  %time windows data of overlap functions saved as cell array bs1,bs2,..,bs29

% x1: waiting times t_w corresponding to each window.
% The snapshot indices are spaced 6592 frames apart starting from frame 19000.
% Each snapshot is converted to simulation time: snapshot * 0.05 * 10000.
x1=19000:6592:205000; 
x1=x1*0.05*10000;

% res2f: results matrix. Each row stores fit parameters for one window.
% Columns: [A0, err_A0, A1, err_A1, A2, tau1, err_tau1, tau2, err_tau2, R^2]

res2f = zeros(length(datasets), 10);

chs=[13 30 30 30 39 40 51 49 49 49 51 51 51 49 51 51 53 53 57 57 59 59 60 60 60 60 62 62];%  droplet size in chains window-wise
%     --- Define the bi-exponential fit function ---
% F(tau) = A0 + A1*exp(-tau/tau1) + (1-A0-A1)*exp(-tau/tau2)
% Parameters: [A0, A1, tau1, tau2]

    fitFunc = @(params, x) params(1) + params(2)*exp(-x/params(3)) + (1-params(1)-params(2))*exp(-x/params(4));
    
    % Initial guesses for parameters: [A0, A1, tc1, tc2]
    initial_params = [0.4, 0.4, 5000, 500000];
    cuto=0.95; % minimum R^2 threshold: fits with R^2 < cuto are flagged as failed
    % Bounds for parameters
    lb = [0,0,0,0];
    ub = [Inf, Inf, Inf, Inf];
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

for i = 1:length(datasets)
    
    data=[];
    fc=[];
    data = datasets{i};
    
     % Prepend the exact initial condition F(tau=0) = 1,pMSD(tau=0)=0
    % The row format matches the data file: [tau_index, F, 0, 1, 0, 1, 0].
    fc=[0 1 0 1 0 1 0]; 
    data=[fc;data];
    
   xdata = data(:, 1)*0.05*10000;
     ydata = data(:, 2); 
    options = optimset('Display', 'off');
    
     % Nonlinear least-squares fit; J is the Jacobian at the solution,
    % needed for confidence interval calculation via nlparci.
    [params_opt, ~, residuals_nonlinear, ~, ~, ~, J] = lsqcurvefit(fitFunc, initial_params, xdata, ydata, lb, ub, options);
    
    A0 = params_opt(1);
    A1 = params_opt(2);
    A2 = 1 - params_opt(1) - params_opt(2);
    tc1 = params_opt(3);
    tc2 = params_opt(4);
    
        % R^2 goodness-of-fit: 1 indicates perfect fit; values below cuto = 0.95
    % typically indicate the bi-exponential model is insufficient or the
    % data are too noisy for reliable parameter extraction.
    ss_res_nonlinear = sum(residuals_nonlinear.^2);
    ss_tot = sum((ydata - mean(ydata)).^2);
    R_squared_nonlinear = 1 - (ss_res_nonlinear / ss_tot);
    
    if R_squared_nonlinear > cuto
        
        ci = nlparci(params_opt, residuals_nonlinear, 'Jacobian', J);
        
        param_errors = (ci(:,2) - ci(:,1)) / 2;
        err_A0 = param_errors(1);
        err_A1 = param_errors(2);
        err_tc1 = param_errors(3);
        err_tc2 = param_errors(4);

        subplot(5, 6, i);
        hold on;
        plot(xdata, ydata, 'ko', 'MarkerSize', 4);
        plot(xdata, fitFunc(params_opt, xdata), '-', 'LineWidth', 2);
        
        xlabel('\tau'); 
        ylabel('F(\tau)');
        
        title_str = sprintf('w=%d A0=%.2f\n\\tau_1=%.1f \\pm %.1f\n\\tau_2=%.1f \\pm %.1f',...
                            i+3, A0, tc1, err_tc1, tc2, err_tc2);
        title(title_str, 'FontSize', 8);
        
	 % Store fit results: [A0, err_A0, A1, err_A1, A2,
        %                     tau1, err_tau1, tau2, err_tau2, R^2]
        
        res2f(i,:)=[A0 err_A0 A1 err_A1 A2 tc1 err_tc1 tc2 err_tc2 R_squared_nonlinear];
    else
        subplot(5, 6, i);
        plot(xdata, ydata, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', [0.7 0.7 0.7]);
        title(sprintf('Fit Failed (R^2 = %.2f)', R_squared_nonlinear), 'FontSize', 8);
        xlabel('\tau'); 
        ylabel('F(\tau)');
    end
end

% =============================================================================
% Plot: fast relaxation timescale tau1 as a function of waiting time t_w
%
% A systematic increase in tau1 with t_w is the hallmark of physical aging:
% the system's dynamics slow down as it evolves toward deeper energy minima
% after the quench.
% =============================================================================

% Colour scheme for line and shaded-error-band plots

line_colors = [
    0.9059, 0.5608, 0.5608;
    0.7255, 0.7843, 0.9059;
    0.7176, 0.8745, 0.6863;
];

patch_colors = [
    0.9059, 0.5608, 0.5608;
    0.7255, 0.7843, 0.9059;
    0.7176, 0.8745, 0.6863;
      
];

figure('Color', 'w'); 
hold on;

% Add the zero line
yline(0, 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'HandleVisibility', 'off');



h_lines = gobjects(1, 1); 


patch_alpha = 0.3; 
% j=2 selects the color. fast timescale tau1 (column 6 of res2f) and its error (column 7).
j=2;
 y_mean = res2f(:,6);
    y_std = res2f(:,7); 
    
    valid_idx = 1:size(res2f,1);
    x_curr = x1(valid_idx);
    x_curr=x_curr';
    y_curr = y_mean(valid_idx);
    err_curr = y_std(valid_idx);
    
    curve1 = y_curr + err_curr;      
    curve2 = y_curr - err_curr;      
    
    x_poly = [x_curr; flipud(x_curr)];
    y_poly = [curve1; flipud(curve2)];
    
    fill(x_poly(:), y_poly(:), patch_colors(j,:), ...
        'FaceAlpha', patch_alpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off'); % <--- Vital: Tells legend to ignore this
        
    % --- CREATE LINE ---
    % We save the handle 'h' into our array 'h_lines'
    h_lines(j) = plot(x_curr(:), y_curr(:), ...
        'Color', line_colors(j,:), ...
        'LineWidth', 3);
%%%%%%%%%%%%%%%%%%

