%% Monte Carlo Simulation: SDID vs TWFE
% Verifying unbiasedness and coverage of SDID under "Divergent Controls"

clear; clc;
rng(12345);

% --- Simulation Settings ---
nSim = 500;                  % Number of repetitions
N = 200;                    % Total units
T_len = 10;                 % Time periods
N_treated = 50;
cohortTimes = [6];
TrueEffect  = 1.0;

% Divergent Controls Setup
% Treated Trend: +0.5
% Good Controls: +0.5 (Parallel)
% Bad Controls : -0.2 (Divergent)
% TWFE compares Treated (+0.5) to Mix (+0.15) -> Bias ~ +0.35 per period?
% SDID should match Good Controls -> Unbiased.

% --- Storage ---
est_sdid = NaN(nSim, 1);
se_sdid  = NaN(nSim, 1);
cov_sdid = NaN(nSim, 1);

est_twfe = NaN(nSim, 1);
se_twfe  = NaN(nSim, 1);
cov_twfe = NaN(nSim, 1);

% Prepare parallel pool if needed
if isempty(gcp('nocreate')), parpool('Processes'); end

fprintf('Starting MC Simulation (%d reps)...\n', nSim);

tic;
parfor r = 1:nSim
    if mod(r, 10) == 0, fprintf('Rep %d/%d...\n', r, nSim); end

    % 1. Generate Data
    data = did.genDIDdata(T_len, N, N_treated/N, ...
        "ATT", TrueEffect, ...
        "treatType", "constantTime", ...
        "cohortTimes", cohortTimes, ...
        "meanError", 0, "errorStd", 1.0, ...
        "preTrendType", "divergentControls", ...
        "preTrendMeanTreated", 0.5, ...
        "preTrendMeanControl", -0.2, ...
        "preTrendSd", 0.5,...
        "preTrendDivergentShare", 0.8); 

    ds = did.Dataset.fromTable(data, ...
        "idVar", "id", "timeVar", "time", "dVar", "D", "yVar", "y", ...
        "Display", false);

    % 2. Fit SDID
    % Use Placebo inference for validity
    try
        res_sdid = did.fit("sdid", ds, ...
            "Display", false, ...
            "SEMethod", "Placebo", ...
            "B", 50, ...             % Low B for speed in MC
            "UseParallel", false);   % Avoid nested parfor

        est_sdid(r) = res_sdid.coef.Estimate;
        se_sdid(r)  = res_sdid.coef.SE;

        ci_lo = est_sdid(r) - 1.96 * se_sdid(r);
        ci_hi = est_sdid(r) + 1.96 * se_sdid(r);
        cov_sdid(r) = (TrueEffect >= ci_lo) && (TrueEffect <= ci_hi);
    catch ME
        warning('SDID failed in rep %d: %s', r, ME.message);
    end

    % 3. Fit TWFE
    try
        res_twfe = did.fit("twfe", ds, "Display", false);
        est_twfe(r) = res_twfe.coef.Estimate;
        se_twfe(r)  = res_twfe.coef.SE;

        ci_lo = est_twfe(r) - 1.96 * se_twfe(r);
        ci_hi = est_twfe(r) + 1.96 * se_twfe(r);
        cov_twfe(r) = (TrueEffect >= ci_lo) && (TrueEffect <= ci_hi);
    catch
    end
end
simTime = toc;

%% Results Analysis
fprintf('\n=== Monte Carlo Results (R=%d) ===\n', nSim);
fprintf('Simulation Time: %.1f seconds\n', simTime);

% Helper Stats
calc_bias = @(est) mean(est, 'omitnan') - TrueEffect;
calc_rmse = @(est) sqrt(mean((est - TrueEffect).^2, 'omitnan'));
calc_cov  = @(cv)  mean(cv, 'omitnan') * 100;

% TWFE
bias_twfe = calc_bias(est_twfe);
rmse_twfe = calc_rmse(est_twfe);
covg_twfe = calc_cov(cov_twfe);
mean_se_twfe = mean(se_twfe, 'omitnan');
sd_est_twfe  = std(est_twfe, 'omitnan');

fprintf('\n[TWFE]\n');
fprintf('  Bias:     %.4f\n', bias_twfe);
fprintf('  RMSE:     %.4f\n', rmse_twfe);
fprintf('  Coverage: %.1f%%\n', covg_twfe);
fprintf('  Mean SE:  %.4f  (vs SD of Est: %.4f)\n', mean_se_twfe, sd_est_twfe);

% SDID
bias_sdid = calc_bias(est_sdid);
rmse_sdid = calc_rmse(est_sdid);
covg_sdid = calc_cov(cov_sdid);
mean_se_sdid = mean(se_sdid, 'omitnan');
sd_est_sdid  = std(est_sdid, 'omitnan');

fprintf('\n[SDID]\n');
fprintf('  Bias:     %.4f\n', bias_sdid);
fprintf('  RMSE:     %.4f\n', rmse_sdid);
fprintf('  Coverage: %.1f%%\n', covg_sdid);
fprintf('  Mean SE:  %.4f  (vs SD of Est: %.4f)\n', mean_se_sdid, sd_est_sdid);

% Bias Reduction
reduc = (abs(bias_twfe) - abs(bias_sdid)) / abs(bias_twfe) * 100;
fprintf('\nBias Reduction: %.1f%%\n', reduc);

% Plot Distribution
figure('Name', 'MC Distribution: SDID vs TWFE', 'Color', 'w');
histogram(est_twfe, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'TWFE');
hold on;
histogram(est_sdid, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'SDID');
xline(TrueEffect, 'k--', 'True Effect');
legend('Location', 'best');
title('Estimator Distribution (Divergent Trends)');
xlabel('Estimate'); ylabel('Frequency');
grid on;
