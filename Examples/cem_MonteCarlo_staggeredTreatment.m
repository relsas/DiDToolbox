% VALIDATE_TOOLBOX_CEM_MIXED_WOOLDRIDGE
% Monte Carlo Simulation: Staggered Adoption with Heterogeneous Effects
% Testing the "Mixed" Strategy and Cohort-level Bias.

clear classes;
clear functions;
rehash toolbox;

% 0. Setup Path
if isempty(which('did.fit'))
    toolboxPath = 'd:\Temp\MCP\DID\did-toolbox';
    if exist(toolboxPath, 'dir')
        addpath(toolboxPath);
    end
end

% 1. Simulation Parameters
nReps = 50;
N_units = 800;
T_periods = 12;

TargetBaseATT = 2.0;
CohortInc = 0.5;

fprintf('Starting Mixed Simulation (Wool vs BJS) with Cohort Analysis (M=%d)...\n', nReps);

bias_Naive       = zeros(nReps, 1);
bias_CEM         = zeros(nReps, 1);
bias_Wool_Raw    = zeros(nReps, 1);
bias_Wool_Matched= zeros(nReps, 1);
bias_Wool_Matched= zeros(nReps, 1);
bias_BJS         = zeros(nReps, 1);
bias_BJS_Cov     = zeros(nReps, 1);

% Cohort Bias Storage: Struct of arrays?
% Simplify: We know cohorts will be stable in this DGP configuration?
% Usually random assignment, but over 50 reps, we'll see all cohorts.
% Let's store bias in a cell array or table accumulator.
cohort_bias_acc = struct('Wool', [], 'BJS', [], 'BJS_Cov', []);

tic;
for i = 1:nReps
    if mod(i, 10) == 0, fprintf('.'); end

    rng(3000 + i);

    % Generate Data
    T_sim = did.genDIDdata(T_periods, N_units, 0.5, ...
        'xNum', 1, ...
        'SelectionBias', 1.0, ...
        'TrendEffectX', 0.5, ...
        'treatType',    "constantTime", ...
        'startPeriod',  4, ...
        'endPeriod',    9, ...
        'ATT',          TargetBaseATT, ...
        'CohortIncrease', CohortInc, ...
        'xTimeStd', 0.1, ...
        'xShockStd', 0.5, ...
        'errorStd', 0.5);

    % True ATT (Global)
    true_tau_global = mean(T_sim.ATT(T_sim.D==1));

    % True ATT per Cohort
    % T_sim.g contains cohort (0 for control).
    treat_G = T_sim(T_sim.D==1, :);
    true_att_by_g = groupsummary(treat_G, "g", "mean", "ATT");
    % Key: g, mean_ATT

    % 1. Naive TWFE
    try
        res_naive = did.fit("twfe", T_sim, 'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', 'Display',false, 'vcov','none');
        bias_Naive(i) = res_naive.ATT - true_tau_global;
    catch, bias_Naive(i) = NaN; end

    % 2. CEM (Standard)
    % Details=true needed for Model return?? We patched CEM.m to silence print.
    res_cem = did.fit("CEM", T_sim, ...
        'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', ...
        'Covariates', ["x1"], ...
        'nBins', 5, ...
        'Robust', true, ...
        'Details', true, 'Display', false);
    bias_CEM(i) = res_cem.tau - true_tau_global;

    % 3. Wooldridge (Raw)
    try
        res_wool = did.wooldridge_TB(T_sim, ...
            'idVar','id', 'timeVar','time', 'yVar','y', 'gVar','g', ...
            'Display', false);
        bias_Wool_Raw(i) = res_wool.overall.Estimate - true_tau_global;
    catch, bias_Wool_Raw(i) = NaN; end

    % Prepare Matched Data
    if isfield(res_cem, 'Model')
        if isa(res_cem.Model, 'LinearModel')
            T_matched = res_cem.Model.Variables;
        elseif isa(res_cem.Model, 'did.Model')
            T_matched = res_cem.Model.ds.T;
        else
            T_matched = table();
        end

        idx_keep = T_matched.Weight > 0;
        T_subset = T_matched(idx_keep, :);
        w_subset = T_matched.Weight(idx_keep);

        % 4. Wooldridge (Matched Sample)
        try
            res_wool_m = did.wooldridge_TB(T_subset, ...
                'idVar','id', 'timeVar','time', 'yVar','y', 'gVar','g', ...
                'Weights', w_subset, ...
                'Display', false);

            bias_Wool_Matched(i) = res_wool_m.overall.Estimate - true_tau_global;

            % Cohort Bias (Wool)
            est_G = res_wool_m.summaryTable; % Cohort, Estimate
            % Compare with true_att_by_g
            % Join on Cohort/g
            T_merged = outerjoin(true_att_by_g, est_G, 'Type','inner', 'LeftKeys','g', 'RightKeys','Cohort');
            % Bias = Estimate - mean_ATT
            bias_vals = T_merged.Estimate - T_merged.mean_ATT;

            % Store: [Rep, Cohort, Bias]
            block = [repmat(i, height(T_merged), 1), T_merged.g, bias_vals];
            cohort_bias_acc.Wool = [cohort_bias_acc.Wool; block];

        catch
            bias_Wool_Matched(i) = NaN;
        end

        % 5. Weighted BJS
        try
            w_bjs = T_matched.Weight; % Include zeros (BJS handles) or subset?
            % BJS expects full table structure usually, let's pass T_matched (full)
            % and Weights vector.

            res_bjs = did.estimators.bjs_imputation(T_matched, ...
                'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', ...
                'Weights', w_bjs, ...
                'Display', false, 'SEMethod', 'None');

            bias_BJS(i) = res_bjs.ATT_overall.ATT - true_tau_global;

            % Cohort Bias (BJS)
            est_G_bjs = res_bjs.ATT_cohort_obs; % cohort, ATT_obs
            T_merged_b = outerjoin(true_att_by_g, est_G_bjs, 'Type','inner', 'LeftKeys','g', 'RightKeys','cohort');
            bias_vals_b = T_merged_b.ATT_obs - T_merged_b.mean_ATT;

            block_b = [repmat(i, height(T_merged_b), 1), T_merged_b.g, bias_vals_b];
            cohort_bias_acc.BJS = [cohort_bias_acc.BJS; block_b];

        catch
            bias_BJS(i) = NaN;
        end

        % 6. Weighted BJS + Covariates (Trend)
        try
            % Construct Trend Covariate (Match what Robust CEM does: x * time)
            % T_matched needs numeric time?
            % T_sim.time is usually categorical or string?
            % Check genDIDdata: time is 1..T (double).

            T_matched.x1_trend = T_matched.x1 .* T_matched.time;

            res_bjs_c = did.estimators.bjs_imputation(T_matched, ...
                'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', ...
                'Weights', w_bjs, ...
                'Covariates', "x1_trend", ...
                'Display', false, 'SEMethod', 'None');

            bias_BJS_Cov(i) = res_bjs_c.ATT_overall.ATT - true_tau_global;

            % Cohort Bias (BJS+Cov)
            est_G_c = res_bjs_c.ATT_cohort_obs;
            T_merged_c = outerjoin(true_att_by_g, est_G_c, 'Type','inner', 'LeftKeys','g', 'RightKeys','cohort');
            bias_vals_c = T_merged_c.ATT_obs - T_merged_c.mean_ATT;

            block_c = [repmat(i, height(T_merged_c), 1), T_merged_c.g, bias_vals_c];
            cohort_bias_acc.BJS_Cov = [cohort_bias_acc.BJS_Cov; block_c];

        catch ME
            if i==1, warning('BJS+Cov failed: %s', ME.message); end
            bias_BJS_Cov(i) = NaN;
        end
    else
        bias_Wool_Matched(i) = NaN;
        bias_BJS(i) = NaN;
        bias_BJS_Cov(i) = NaN;
    end
end
elapsed = toc;
fprintf('\nSimulation Complete in %.2f seconds.\n', elapsed);

% Report Overall
fprintf('\n--- Combined Results (Global Bias) ---\n');
fprintf('Naive TWFE:       %.4f\n', mean(bias_Naive, 'omitnan'));
fprintf('CEM (Standard):   %.4f\n', mean(bias_CEM, 'omitnan'));
fprintf('Wooldridge (Raw): %.4f\n', mean(bias_Wool_Raw, 'omitnan'));
fprintf('Wool on CEM Data: %.4f\n', mean(bias_Wool_Matched, 'omitnan'));
fprintf('Weighted BJS:     %.4f\n', mean(bias_BJS, 'omitnan'));
fprintf('W-BJS + Trend:    %.4f  << TARGET\n', mean(bias_BJS_Cov, 'omitnan'));

% Plot Overall
local_plot_results(bias_Naive, bias_CEM, bias_Wool_Raw, bias_Wool_Matched, bias_BJS, bias_BJS_Cov);

% Analyze & Plot Cohort Bias
local_plot_cohort_bias(cohort_bias_acc);


% -------------------------------------------------------------------------
function local_plot_results(b_naive, b_cem, b_wool, b_mixed, b_bjs, b_bjs_c)
fig = figure('Name', 'Mixed Strategy Bias', 'Color', 'w');
hold on;
    function plot_kde(data, color, name, lw)
        data = data(~isnan(data));
        if isempty(data), return; end
        [f, x] = ksdensity(data);
        plot(x, f, color, 'LineWidth', lw, 'DisplayName', name);
    end
plot_kde(b_naive, 'r', 'Naive', 1.5);
plot_kde(b_cem,   'b', 'CEM', 1.5);
plot_kde(b_wool,  'g', 'Wooldridge', 1.5);
plot_kde(b_mixed, 'k', 'Wooldridge (W)', 2.0);
plot_kde(b_bjs,   'm', 'BJS (W)', 2.5);
plot_kde(b_bjs_c, 'c', 'BJS (W+Trend)', 3.0);

xline(0, '--k', 'Unbiased');
legend('Location','best');
title('Global Estimator Performance');
xlabel('Bias'); ylabel('Density');
grid on;
saveas(fig, 'd:\Temp\MCP\CEM_Analysis\Mixed_Sim_Coefficients.png');
end

function local_plot_cohort_bias(acc)
% acc.Wool: [Rep, Cohort, Bias]

W = array2table(acc.Wool, 'VariableNames',{'Rep','Cohort','Bias'});
B = array2table(acc.BJS,  'VariableNames',{'Rep','Cohort','Bias'});
C = array2table(acc.BJS_Cov, 'VariableNames',{'Rep','Cohort','Bias'});

statsW = groupsummary(W, "Cohort", {"mean","std"}, "Bias");
statsB = groupsummary(B, "Cohort", {"mean","std"}, "Bias");
statsC = groupsummary(C, "Cohort", {"mean","std"}, "Bias");

fig = figure('Name','Cohort Bias Analysis', 'Color','w');
hold on;

% Wooldridge
errorbar(statsW.Cohort - 0.15, statsW.mean_Bias, statsW.std_Bias, ...
    'o-', 'Color','k', 'LineWidth', 1.5, 'DisplayName', 'Wooldridge (W)');

% BJS (Weights Only)
errorbar(statsB.Cohort, statsB.mean_Bias, statsB.std_Bias, ...
    's-', 'Color','m', 'LineWidth', 1.5, 'DisplayName', 'BJS (W)');

% BJS (Weights + Trend)
errorbar(statsC.Cohort + 0.15, statsC.mean_Bias, statsC.std_Bias, ...
    '^-', 'Color','c', 'LineWidth', 1.5, 'DisplayName', 'BJS (W+Trend)');

yline(0, '--k');
xlabel('Cohort (Treatment Start Time)');
ylabel('Bias (Estimate - True)');
legend('Location','best');
title('Bias by Cohort: Weighted Estimators');
grid on;
saveas(fig, 'd:\Temp\MCP\CEM_Analysis\Cohort_Bias_Analysis.png');
end
