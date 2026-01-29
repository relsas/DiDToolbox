function verify_monte_carlo()
% VERIFY_MONTE_CARLO  Backtest all toolbox estimators against known truth.
%
%   Scenario:
%     - N=200, T=10
%     - Treatment: Staggered, constant over time but increasing by cohort.
%       (Scenario "Primary Use Case": constantTime + CohortIncrease)
%     - M=50 Simulations (adjustable)
%
%   Estimators:
%     1. TWFE (Static) - Expect Bias due to heterogeneity
%     2. Wooldridge (Mundlak)
%     3. CS (Callaway & Sant'Anna)
%     4. IW (Sun & Abraham)
%     5. BJS (Borusyak et al.)
%     6. CH (de Chaisemartin)
%
% -------------------------------------------------------------------------

% Force clear and path refresh
clear; clc;
rng(42); % reproducibility

% --- Configuration ---
M = 100;             % Fast check
N = 200;            % Units
T_periods = 10;     % Time periods

% "Primary Use Case": Cohort-specific constant effects
genOpts = {
    'treatType',      "constantTime", ...
    'ATT',            1.0, ...
    'CohortIncrease', 0.2, ...
    'startPeriod',    3, ...    % Allow early treatment
    'endPeriod',      10, ...
    'numCohorts',     4, ...    % Staggered adoption
    'errorStd',       1.0 ...
    };

estimators = {
    "TWFE",         @() did.estimators.TWFE("details",false, "Cluster",true);
    "Wooldridge",   @() did.estimators.Wooldridge("Display",false);
    "CS",           @() did.estimators.CS("Display",false, "B",99);
    "IW",           @() did.estimators.IW("Display",false);
    "BJS",          @() did.estimators.BJS("Display",false, "useParallel",24); % Serial!
    "CH",           @() did.estimators.CH("Display",false, "B",99);
    };

nEst = size(estimators, 1);

% Storage
% --- Storage (Matrix form for parfor) ---
% est_mat: M x nEst
est_mat   = zeros(M, nEst);
se_mat    = zeros(M, nEst);
bias_mat  = zeros(M, nEst);
cover_mat = false(M, nEst);
true_att_store = zeros(M,1);

% --- Parallelization Config ---
% "SimLevel" (Recommended): Outer parfor (M workers), Inner Serial. Best for M=100.
% "EstLevel": Outer for (Serial), Inner Parallel. Best for debugging BJS logic with M=1.
% "Serial":   All serial. Debugging.
ParallelMode = "SimLevel";

estimators = {
    "TWFE",         @() did.estimators.TWFE("details",false, "Cluster",true);
    "Wooldridge",   @() did.estimators.Wooldridge("Display",false);
    "CS",           @() did.estimators.CS("Display",false, "B",99);
    "IW",           @() did.estimators.IW("Display",false);
    "BJS",          @() did.estimators.BJS("Display",false, "useParallel",0); % Default
    "CH",           @() did.estimators.CH("Display",false, "B",99);
    };

nEst = size(estimators, 1);

% Configure Based on Mode
switch ParallelMode
    case "SimLevel"
        % Enable Parfor for M, force BJS serial (to avoid nested overhead)
        forceBJSSerial = true;
        RunParforM     = Inf; % Max workers
    case "EstLevel"
        % Serial M, Parallel BJS (enable explicitly)
        forceBJSSerial = false;
        RunParforM     = 0;   % Serial
        % Update BJS to use parallel
        for i=1:nEst
            if estimators{i,1} == "BJS"
                estimators{i,2} = @() did.estimators.BJS("Display",false, "useParallel",24);
            end
        end
    case "Serial"
        forceBJSSerial = true;
        RunParforM     = 0;   % Serial
    otherwise
        error("Unknown ParallelMode");
end

if forceBJSSerial
    for i=1:nEst
        if estimators{i,1} == "BJS"
            estimators{i,2} = @() did.estimators.BJS("Display",false, "useParallel",0);
        end
    end
end

% Storage
% est_mat: M x nEst
est_mat   = zeros(M, nEst);
se_mat    = zeros(M, nEst);
bias_mat  = zeros(M, nEst);
cover_mat = false(M, nEst);
true_att_store = zeros(M,1);

if RunParforM
    fprintf('Starting Monte Carlo (M=%d) with PARFOR (SimLevel)...\n', M);
else
    fprintf('Starting Monte Carlo (M=%d) with FOR (Serial/EstLevel)...\n', M);
end

parfor (m = 1:M, RunParforM) % Conditional parfor
    % No printing in parfor

    % 1. Data Generation
    try
        T = did.genDIDdata(T_periods, N, 0.5, genOpts{:});

        % Dynamic Truth
        isTreated = (T.D == 1);
        if ~any(isTreated)
            currentTrueATT = NaN;
        else
            currentTrueATT = mean(T.ATT(isTreated));
        end
        true_att_store(m) = currentTrueATT;

        if isnan(currentTrueATT)
            continue;
        end

        % Create Dataset (Silence!)
        ds = did.Dataset.fromTable(T, "idVar","id", "timeVar","time", ...
            "yVar","y", "dVar","D", "Display",false, "describe",false);

        % 2. Run Estimators
        m_est   = zeros(1, nEst);
        m_se    = zeros(1, nEst);
        m_bias  = zeros(1, nEst);
        m_cover = false(1, nEst);

        for i = 1:nEst
            name = estimators{i,1};
            factory = estimators{i,2};

            try
                estObj = factory();

                % Fit (Silence output via argument, evalc not allowed in parfor)
                res = estObj.fit(ds);

                [b, se] = extract_result(res, name);

                % --- TWFE Special Handling for SE ---
                if name == "TWFE" && (isnan(se) || se==0)
                    % Calculate Cluster SE manually using diagnostics
                    se = calc_twfe_se(res, ds);
                end
                % ------------------------------------

                cover = false;
                if ~isnan(se) && se~=0
                    ci_lower = b - 1.96*se;
                    ci_upper = b + 1.96*se;
                    cover = (currentTrueATT >= ci_lower && currentTrueATT <= ci_upper);
                end

                m_est(i)  = b;
                m_se(i)   = se;
                m_bias(i) = b - currentTrueATT;
                m_cover(i)= cover;

            catch ME
               
                m_est(i)  = NaN;
                m_se(i)   = NaN;
                m_bias(i) = NaN;
            end
        end

        est_mat(m,:)   = m_est;
        se_mat(m,:)    = m_se;
        bias_mat(m,:)  = m_bias;
        cover_mat(m,:) = m_cover;

    catch ME
        fprintf('Critical Data Gen Error: %s\n', ME.message);
    end
end

% --- Unpack to results struct for reporting ---
for i = 1:nEst
    name = estimators{i,1};
    results.(name).est   = est_mat(:,i);
    results.(name).se    = se_mat(:,i);
    results.(name).bias  = bias_mat(:,i);
    results.(name).cover = cover_mat(:,i);
end

% --- Reporting ---
fprintf('\n%s\n', repmat('-',1,85));
fprintf('%-12s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
    'Estimator','Bias','Emp. SD','Avg SE','RMSE','Coverage');
fprintf('%s\n', repmat('-',1,85));

mean_truth = mean(true_att_store, 'omitnan');
fprintf('TRUE ATT (Mean over Sim): %.4f\n', mean_truth);
fprintf('%s\n', repmat('-',1,85));

for i = 1:nEst
    name = estimators{i,1};
    vals = results.(name);

    valid = ~isnan(vals.est);
    nValid = sum(valid);

    if nValid < 5
        fprintf('%-12s | (Failed or Insufficient Data)\n', name);
        continue;
    end

    bias   = mean(vals.bias(valid));
    emp_sd = std(vals.est(valid));
    avg_se = mean(vals.se(valid));
    rmse   = sqrt(mean(vals.bias(valid).^2));
    cov_rate = mean(vals.cover(valid));

    fprintf('%-12s | %10.4f | %10.4f | %10.4f | %10.4f | %10.1f%%\n', ...
        name, bias, emp_sd, avg_se, rmse, cov_rate*100);
end
fprintf('%s\n', repmat('-',1,85));
end

function [b, se] = extract_result(res, method)
b = NaN; se = NaN;

% Attempt to find 'Estimate' and 'SE' in .overall or .beta or standard tables

% 0. BJS specific: ATT_overall
if isfield(res, 'ATT_overall') && (isstruct(res.ATT_overall) || istable(res.ATT_overall))
    try
        b  = res.ATT_overall.ATT;
        if ismember('SE', res.ATT_overall.Properties.VariableNames)
            se = res.ATT_overall.SE;
        else
            se = NaN;
        end

        if ~isscalar(b), b=b(1); end
        if ~isscalar(se), se=se(1); end
        return;
    catch
    end
end

% 0b. CS specific: Aggregates.overall
if isfield(res, 'Aggregates') && isfield(res.Aggregates, 'overall')
    try
        r = res.Aggregates.overall;
        if istable(r)
            b = r.Estimate(1);
            if ismember('SE', r.Properties.VariableNames), se = r.SE(1); end
        elseif isstruct(r)
            b = r.Estimate;
            se = r.SE;
        end
        return;
    catch
    end
end

% 1. Explicit .overall property (e.g., Wooldridge wrapper, CS wrapper new)
if isfield(res, 'overall') && istable(res.overall)
    % usually table with Estimate, SE
    try
        b  = res.overall.Estimate(1);
        se = res.overall.SE(1);
        return;
    catch
    end
elseif isfield(res, 'overall') && isstruct(res.overall)
    try
        b  = res.overall.Estimate;
        se = res.overall.SE;
        return;
    catch
    end
end

% 2. TWFE result (usually .coef table, first param is ATT)
% Or a field named 'ATT'
if isfield(res, 'ATT') && ~isempty(res.ATT) && isnumeric(res.ATT) && isscalar(res.ATT)
    b = res.ATT;
    if isfield(res, 'SE'), se = res.SE; end
    return;
end

% 3. Check summaryTable for a row named "ATT", "Treatment", or "D"
if isfield(res, 'summaryTable')
    T = res.summaryTable;
    % Try finding row
    idx = find(contains(string(T.Name), ["ATT","D","Treatment","Overall"],'IgnoreCase',true), 1);
    if ~isempty(idx)
        b  = T.Estimate(idx);
        if ismember('SE', T.Properties.VariableNames)
            se = T.SE(idx);
        end
        return;
    end
    % Fallback: take first row if only one or if logic dictates
    if height(T)==1
        b = T.Estimate(1);
        if ismember('SE', T.Properties.VariableNames)
            se = T.SE(1);
        else
            % Maybe in vcov?
        end
    end
end

% 4. Coefficients table fallback
if isnan(b) && isfield(res,'coef')
    b = res.coef.Estimate(1);
    if ismember('SE', res.coef.Properties.VariableNames)
        se = res.coef.SE(1);
    end
end
end

function se = calc_twfe_se(res, ds)
se = NaN;
try
    if ~isfield(res, 'Diagnostics') || ~isfield(res.Diagnostics, 'design'), return; end

    des = res.Diagnostics.design;
    X = des.X;
    y = des.y;
    b = res.coef.Estimate;

    if numel(b) ~= size(X,2), return; end

    resid = y - X * b;
    % Access T from Dataset (public property)
    ids = ds.T.id; % Hardcoded idVar='id' as per script usage

    if numel(ids) ~= numel(resid), return; end

    [g_id, ~] = findgroups(ids);
    K = size(X,2);
    meat = zeros(K,K);
    u_g = unique(g_id);
    for k = 1:numel(u_g)
        idx = (g_id == u_g(k));
        X_g = X(idx, :);
        e_g = resid(idx);
        meat = meat + (X_g' * (e_g * e_g') * X_g);
    end
    XTX_inv = (X'*X) \ eye(K);
    V = XTX_inv * meat * XTX_inv;
    G = numel(u_g);
    adj = G / (G - 1);
    V = V * adj;
    se = sqrt(V(1,1));
catch
end
end

