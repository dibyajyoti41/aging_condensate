% pmsdanalysis.m
%
% Three-step analysis of window-resolved pairwise MSD (pMSD), producing
% the condensate aging exponent alpha (Eq. 8 in the paper):
%   <dR^2(tau)>|_tw = D(tw) * tau^kappa ,   D(tw) = D0 * tw^(-alpha)
%
% INPUT (must exist in the workspace before running):
%   datasets : 1xN cell array, one entry per aging window, in age order.
%              Each entry is an M-by-3+ matrix with columns matching the
%              output of droplet_overlap.f90:
%                  col 1 = lag_frames
%                  col 3 = msd_whole      (this is the pMSD used below)
%   t_obs    : 1xN vector, the physical observation time (window age) for
%              each entry of `datasets`, same order/length as `datasets`
%        
%
% STEP 1: for each window, fit log10(pMSD) vs log10(tau) as a straight
%         line through the origin (i.e. pMSD ~ tau^k) to get the local
%         subdiffusion exponent k for that window. Windows whose log-log 
%	  fit falls below R^2 = 0.96 are dropped from all later steps.
%         
% STEP 2: compute the effective diffusion prefactor D = pMSD / tau^kappa
%         at every point in each (valid) window, then average over tau
%         to get one D value per window. kappa is FIXED at 0.35 here
%         (which is near average of per-window kappa in step 1) - this matches the fixed
%         subdiffusive exponent used throughout the paper (Eq. 8-10).
% STEP 3: fit D(t_obs) ~ t_obs^p (again log-log, through the origin)
%         across all valid windows; alpha = -p is the aging exponent.

n_windows = length(datasets);
%tau_trim  = @(x) x(300:end-20);   % drop short-time and long-time ends

% =====================================================================
% STEP 1: local subdiffusion exponent k per window, log-log linear fit
% =====================================================================
k_fit      = zeros(n_windows, 1);
R2_fit     = zeros(n_windows, 1);
rmse_fit   = zeros(n_windows, 1);
intercept  = zeros(n_windows, 1);
valid_idx  = [];

R2_min = 0.96;   % quality cutoff: windows below this are excluded

figure; hold on;
colors = turbo(n_windows);

for i = 1:n_windows
    data = datasets{i};
    %tau  = tau_trim(data(:, 1));
    tau = data(:, 1);
    tau = tau*0.05*10000; %conversion from snapshot to md timesteps
    pmsd = tau_trim(data(:, 3));       % pMSD = column 3 (msd_whole)

    log_tau  = log10(tau);
    log_pmsd = log10(pmsd);

    pp = polyfit(log_tau, log_pmsd, 1);   % pp(1) = local k, pp(2) = intercept
    fit_vals = polyval(pp, log_tau);
    R2 = 1 - sum((log_pmsd - fit_vals).^2) / sum((log_pmsd - mean(log_pmsd)).^2);
    rmse = sqrt(mean((log_pmsd - fit_vals).^2));

    k_fit(i)     = pp(1);
    intercept(i) = pp(2);
    R2_fit(i)    = R2;
    rmse_fit(i)  = rmse;

    if R2 > R2_min
        valid_idx = [valid_idx; i]; %#ok<AGROW>
        plot(tau, pmsd, '-', 'Color', colors(i, :));
    end
end
xlabel('\tau'); ylabel('pMSD(\tau)');
title(sprintf('Step 1: pMSD vs. \\tau (only R^2 > %.2f windows shown)', R2_min));
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

fprintf('Step 1: kept %d of %d windows (R^2 > %.2f)\n', numel(valid_idx), n_windows, R2_min);

% =====================================================================
% STEP 2: effective diffusion coefficient per window, kappa fixed at 0.35
% =====================================================================
kappa = 0.35;
D_mean = zeros(n_windows, 1);
D_std  = zeros(n_windows, 1);

for idx = 1:numel(valid_idx)
    i    = valid_idx(idx);
    data = datasets{i};
    tau  = tau_trim(data(:, 1));
    pmsd = tau_trim(data(:, 3));

    D_points  = pmsd ./ (tau .^ kappa);
    D_mean(i) = mean(D_points);
    D_std(i)  = std(D_points);
end

% =====================================================================
% STEP 3: fit D(t_obs) ~ t_obs^p across valid windows; alpha = -p
% =====================================================================
log_tobs = log10(t_obs(valid_idx));
log_D    = log10(D_mean(valid_idx));

pp3 = polyfit(log_tobs, log_D, 1);
D_fit_vals = polyval(pp3, log_tobs);
R2_p  = 1 - sum((log_D(:) - D_fit_vals(:)).^2) / sum((log_D(:) - mean(log_D)).^2);
p_fit = pp3(1);
alpha = -p_fit;

fprintf('Step 3: D ~ t_obs^p,  p = %.4f,  alpha = %.4f,  R^2 = %.4f\n', p_fit, alpha, R2_p);

figure('Color', 'w'); hold on;
errorbar(t_obs(valid_idx), D_mean(valid_idx), D_std(valid_idx), 'bo', 'MarkerFaceColor', 'b');
t_line = logspace(log10(min(t_obs(valid_idx))), log10(max(t_obs(valid_idx))), 100);
plot(t_line, 10.^polyval(pp3, log10(t_line)), 'r-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('t_{obs}'); ylabel('D(t_{obs})');
title(sprintf('Step 3: D(t_{obs}) fit,  \\alpha = %.3f,  R^2 = %.3f', alpha, R2_p));
legend('D per window', sprintf('fit: D \\propto t_{obs}^{%.3f}', p_fit), 'Location', 'southwest');
hold off;
