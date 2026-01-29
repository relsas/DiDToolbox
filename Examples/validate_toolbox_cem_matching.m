% VALIDATE_TOOLBOX_CEM (Simulation Mode)
% Conducts a Monte Carlo Simulation using the integrated DID Toolbox estimators.
% Compares Naive TWFE, CEM, and CEM+Robust (Trends).
% Generates plots for Coefficient Distribution, Covariate Balance, and Event Trends.


% 0. Setup Path
if isempty(which('did.fit'))
    toolboxPath = 'd:\Temp\MCP\DID\did-toolbox';
    if exist(toolboxPath, 'dir')
        addpath(toolboxPath);
        fprintf('Added DID Toolbox to path: %s\n', toolboxPath);
    else
        warning('DID Toolbox not found at %s.', toolboxPath);
    end
end

% 1. Simulation Parameters
nReps = 50;         % Number of Repetitions (Reduced for speed, set to 100+ for smooth plots)
N_units = 500;
T_periods = 10;
TargetATT = 2.0;

fprintf('Starting Toolbox CEM Simulation (M=%d)...\n', nReps);

% Storage
est_Naive     = zeros(nReps, 1);
est_CEM       = zeros(nReps, 1);
est_CEM_Rob   = zeros(nReps, 1);

% Keep last iteration data for diagnostics
T_last = [];
res_cem_last = [];

tic;
for i = 1:nReps
    if mod(i, 10) == 0, fprintf('.'); end

    % A. Generate Data (Relevant Scenario: Selection on X + Trend on X)
    % Using Standard did.genDIDdata with REDUCED NOISE for clearer trends
    rng(1000 + i); % Reproducible seeds
    T_sim = did.genDIDdata(10, N_units, 0.4, ...
        'xNum', 1, ...
        'SelectionBias', 1.0, ...
        'TrendEffectX', 0.5, ...
        'ATT', TargetATT, ...
        'xTimeStd', 0.1, ...   % Reduced from 1.0 to 0.1
        'xShockStd', 0.5, ...  % Reduced
        'errorStd', 0.5);      % Reduced outcome noise

    % B. Naive TWFE
    % We can use did.fit("twfe") or a quick manual regression for speed if did.fit is heavy.
    % Let's use did.fit to fully validate the toolbox.
    % Suppress output typically.
    try
        % Fast fit if possible.
        res_naive = did.fit("twfe", T_sim, 'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', ...
            'Display',false, 'vcov','none');
        est_Naive(i) = res_naive.ATT;
    catch ME
        warning('Naive fit failed at rep %d: %s', i, ME.message);
        est_Naive(i) = NaN;
    end

    % C. CEM (Base)
    % Robust = false
    res_cem = did.fit("CEM", T_sim, ...
        'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', ...
        'Covariates', ["x1"], ...
        'nBins', 5, ...
        'Robust', false, ...
        'Details', false, ...
        'Display', false); % Add Details flag to suppress print output in CEM.m
    est_CEM(i) = res_cem.tau;

    % D. CEM (Robust)
    % Robust = true
    res_cem_rob = did.fit("CEM", T_sim, ...
        'idVar','id', 'timeVar','time', 'yVar','y', 'dVar','D', ...
        'Covariates', ["x1"], ...
        'nBins', 5, ...
        'Robust', true, ...
        'Details', false, ...
        'Display', false);
    est_CEM_Rob(i) = res_cem_rob.tau;

    if i == nReps
        T_last = T_sim;
        res_cem_last = res_cem_rob;
    end
end
elapsed = toc;
fprintf('\nSimulation Complete in %.2f seconds.\n', elapsed);

% 2. Analysis
bias_Naive   = mean(est_Naive - TargetATT, 'omitnan');
bias_CEM     = mean(est_CEM - TargetATT, 'omitnan');
bias_CEM_Rob = mean(est_CEM_Rob - TargetATT, 'omitnan');

fprintf('\n--- Results (Target ATT = %.2f) ---\n', TargetATT);
fprintf('Naive TWFE Bias:   %.4f\n', bias_Naive);
fprintf('CEM (Base) Bias:   %.4f\n', bias_CEM);
fprintf('CEM (Robust) Bias: %.4f\n', bias_CEM_Rob);

% 3. Plotting

fprintf('Generating plots...\n');

% Extract data from last run for diagnostic plots
if isfield(res_cem_last, 'Model')
    mdl = res_cem_last.Model;
    T_matched = mdl.Variables;

    % Weights
    if ismember('Weight', T_matched.Properties.VariableNames)
        w_matched = T_matched.Weight;
    else
        w_matched = ones(height(T_matched), 1);
    end

    local_plot_event_study(T_last, T_matched, w_matched);
    local_plot_covariate_balance(T_last, T_matched, w_matched);
end

local_plot_coefficients(est_Naive, est_CEM, est_CEM_Rob, TargetATT);



% -------------------------------------------------------------------------
% LOCAL PLOTTING FUNCTIONS
% -------------------------------------------------------------------------

function local_plot_coefficients(b_naive, b_cem, b_ctrl, true_att)
fig = figure('Name', 'Toolbox Simulation: Coefficients', 'Color', 'w');
hold on;

[f1, x1] = ksdensity(b_naive);
[f2, x2] = ksdensity(b_cem);
[f3, x3] = ksdensity(b_ctrl);

plot(x1, f1, 'r', 'LineWidth', 2, 'DisplayName', 'Naive TWFE');
plot(x2, f2, 'b', 'LineWidth', 2, 'DisplayName', 'CEM');
plot(x3, f3, 'g', 'LineWidth', 2, 'DisplayName', 'CEM + Trend');

xline(true_att, '--k', 'True ATT', 'LineWidth', 1.5, 'HandleVisibility','off');

legend('Location', 'best');
title('Methods Comparison: Coefficient Distribution');
xlabel('Estimated ATT'); ylabel('Density');
grid on; hold off;

end

function local_plot_event_study(T_full, T_match, ~)
fig = figure('Name', 'Toolbox Simulation: Event Trends', 'Color', 'w');

% 1. Unmatched
% Use 'eventTime' instead of 'time'
isTreated = (T_full.g > 0);
[G, ~, tr_full] = findgroups(T_full.eventTime, isTreated);
y_mean_full = splitapply(@mean, T_full.y, G);

% Get unique event times and unique group indicators
% G maps to unique rows of [eventTime, isTreated]
unique_keys = unique([T_full.eventTime, double(isTreated)], 'rows');
times = unique_keys(:,1);
is_tr = unique_keys(:,2) == 1;

% Separate Treated and Control
% Ensure alignment
u_times = unique(times);
y_tr_full = NaN(size(u_times));
y_co_full = NaN(size(u_times));

% Map back
for k = 1:numel(u_times)
    t = u_times(k);
    % Treated
    idx_tr = (times == t) & is_tr;
    if any(idx_tr), y_tr_full(k) = y_mean_full(idx_tr); end
    % Control
    idx_co = (times == t) & ~is_tr;
    if any(idx_co), y_co_full(k) = y_mean_full(idx_co); end
end

subplot(1, 2, 1);
plot(u_times, y_tr_full, '-or', 'LineWidth', 2, 'DisplayName', 'Treated');
hold on;
plot(u_times, y_co_full, '-ob', 'LineWidth', 2, 'DisplayName', 'Control');
xline(-0.5, '--k', 'Treatment', 'HandleVisibility','off');
title('Unmatched Data (Event Time)');
xlabel('Event Time'); ylabel('Mean Y');
legend('Location', 'best'); grid on;

% 2. Matched (Weighted)
isTreatedM = (T_match.g > 0);
w = T_match.Weight;

[G_m, ~, tr_match] = findgroups(T_match.eventTime, isTreatedM);
w_mean = @(y, w) sum(y.*w)/sum(w);
y_means_w = splitapply(w_mean, T_match.y, w, G_m);

% Re-map for matched
unique_keys_m = unique([T_match.eventTime, double(isTreatedM)], 'rows');
times_m = unique_keys_m(:,1);
is_tr_m = unique_keys_m(:,2) == 1;

u_times_m = unique(times_m);
y_tr_w = NaN(size(u_times_m));
y_co_w = NaN(size(u_times_m));

for k = 1:numel(u_times_m)
    t = u_times_m(k);
    % Treated
    idx_tr = (times_m == t) & is_tr_m;
    if any(idx_tr), y_tr_w(k) = y_means_w(idx_tr); end
    % Control
    idx_co = (times_m == t) & ~is_tr_m;
    if any(idx_co), y_co_w(k) = y_means_w(idx_co); end
end

subplot(1, 2, 2);
plot(u_times_m, y_tr_w, '-or', 'LineWidth', 2, 'DisplayName', 'Treated (Matched)');
hold on;
plot(u_times_m, y_co_w, '-om', 'LineWidth', 2, 'DisplayName', 'Control (Weighted)');
xline(-0.5, '--k', 'Treatment', 'HandleVisibility','off');
title('Matched Data (Event Time)');
xlabel('Event Time'); ylabel('Mean Y');
legend('Location', 'best'); grid on;

end

function local_plot_covariate_balance(T_full, T_match, ~)
fig = figure('Name', 'Toolbox Simulation: Balance', 'Color', 'w');

% Collapse to unit
[~, idx] = unique(T_full.id);
T_units_full = T_full(idx, :);

% Matched units
if ismember('id', T_match.Properties.VariableNames)
    uid = T_match.id;
else
    uid = double(string(T_match.id_cat));
end
[~, idxM] = unique(uid);
T_units_match = T_match(idxM, :);

subplot(1, 2, 1);
histogram(T_units_full.x1(T_units_full.g > 0), 'Normalization','pdf','FaceColor','r','FaceAlpha',0.5); hold on;
histogram(T_units_full.x1(T_units_full.g == 0), 'Normalization','pdf','FaceColor','b','FaceAlpha',0.5);
title('Unmatched Balance'); legend('Treated','Control');

subplot(1, 2, 2);
histogram(T_units_match.x1(T_units_match.g > 0), 'Normalization','pdf','FaceColor','r','FaceAlpha',0.5); hold on;
histogram(T_units_match.x1(T_units_match.g == 0), 'Normalization','pdf','FaceColor','m','FaceAlpha',0.5);
title('Matched Balance'); legend('Treated','Control');


end
