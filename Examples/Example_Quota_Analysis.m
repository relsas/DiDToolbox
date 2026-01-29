%% Example: Women's Quota and Maternal Mortality (Bhalotra et al., 2023)
%
% This script analyzes the effect of gender quotas in parliament on
% maternal mortality capability.
% Data source: "quota_stata.xlsx" (Balanced Panel)
%
% Estimators compared:
% 1. Standard Two-Way Fixed Effects (TWFE) - Static
% 2. Synthetic Difference-in-Differences (SDID) - Static
% 3. Wooldridge (2021) Event Study - For Pre-trend Diagnostics

clear; clc; close all;
% Ensure we use the local toolbox version
toolboxPath = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(toolboxPath));

% --- 1. Load Data ---
dataPath = fullfile(fileparts(mfilename('fullpath')), '..', 'Data', 'quota_stata.xlsx');
fprintf('Loading data from: %s\n', dataPath);
T = readtable(dataPath);

% Display basic info
head(T)
summary(T)

% Drop missing values (BUT allow quotaYear to be NaN for controls)
% Logic: If quotaYear is NaN, it implies "Never Treated" (so quota should be 0).
if ismember('quotaYear', T.Properties.VariableNames)
    idxControl = isnan(T.quotaYear);
    % Check if truly never treated
    if any(T.quota(idxControl) == 1)
        warning('Some units with NaN quotaYear have quota=1. Dropping them.');
        T(idxControl & T.quota==1, :) = [];
    end

    % For clean controls, set quotaYear = Inf (or a large number)
    % This ensures did.Dataset handles them as never-treated if supported,
    % or we can rely on gVar being > max(time).
    T.quotaYear(isnan(T.quotaYear)) = 9999;
end

% Check remaining vars for random missingness
% User Note: Stata drops cases where lngdp is NaN (N=115 groups for lnmmrt analysis).
varsToCheck = {'lnmmrt', 'year', 'quota', 'lngdp'};
for i = 1:length(varsToCheck)
    vn = varsToCheck{i};
    if ismember(vn, T.Properties.VariableNames)
        nMiss = sum(isnan(T.(vn)));
        if nMiss > 0
            fprintf('Dropping %d rows with missing %s.\n', nMiss, vn);
            T = T(~isnan(T.(vn)), :);
        end
    end
end


% --- 2. Setup Dataset ---
% Strict Balanced Panel Enforcement (to match Stata's N=115, T=26)
[G, TID] = findgroups(T.country);
num_obs = splitapply(@numel, T.year, G);
max_T = max(num_obs);
% Count valid obs post-filtering
num_obs_post = splitapply(@numel, T.year, G);
valid_groups = num_obs_post == max_T;

if any(~valid_groups)
    idx_drop = ismember(G, find(~valid_groups));
    fprintf('Dropping %d unbalanced units (T < %d).\n', sum(~valid_groups), max_T);
    T = T(~idx_drop, :);
end

% Recalculate QuotaYear (Cohort) just in case
T.quotaYear = zeros(height(T), 1);
[G_clean, ID_clean] = findgroups(T.country);
for i = 1:max(G_clean)
    idx = (G_clean == i);
    df_i = T(idx, :);
    if sum(df_i.quota) > 0
        % Cohort = First year of treatment
        first_yr = min(df_i.year(df_i.quota==1));
        T.quotaYear(idx) = first_yr;
    else
        T.quotaYear(idx) = 9999;
    end
end

% idVar:   country
% timeVar: year
% yVar:    lnmmrt
% dVar:    quota
% gVar:    quotaYear

% Use strings to avoid char concatenation issues in sortrows
% Debug: check surviving units and their treatment status
[gUnits, gId] = findgroups(T.country);
gQuota = splitapply(@(x) mean(x, 'omitnan'), T.quotaYear, gUnits);
gTreated = splitapply(@(x) max(x, [], 'omitnan'), T.quota, gUnits);
fprintf('\nSurviving Units (%d):\n', numel(gId));
% Tdebug = table(gId, gQuota, gTreated, 'VariableNames', {'Country', 'QuotaYear', 'EverTreated'});
% disp(Tdebug);

% Check if we have controls
if all(gTreated == 1)
    warning('NO CONTROL UNITS available in the balanced/cleaned sample! SDID cannot run.');
end

ds = did.Dataset.fromTable(T, ...
    "idVar", "country", ...
    "timeVar", "year", ...
    "yVar", "lnmmrt", ...
    "dVar", "quota", ...
    "gVar", "quotaYear"); % Needed for Wooldridge

fprintf('\nDataset Summary:\n');
ds.describe();

% --- 3. Pre-trend Diagnostics (Event Study) ---
% We use Wooldridge's method to estimate dynamic effects and test for pre-trends.
fprintf('\n=========================================\n');
fprintf('   Pre-trend Diagnostic (Wooldridge ES)  \n');
fprintf('=========================================\n');

% Fit Event Study (includes leads and lags)
resES = did.fit("wooldridge", ds, "EventStudy", true, "Display", false);

% Check Pre-trends using the diagnostics module
% PreCut=-1 because e=-1 is reference? No, usually e=-1 is reference (0 coeff).
% We check e < -1 or e < 0 depending on implementation.
% Wooldridge drops -1 by default.
diag = did.diagnostics.pretrend_event(resES, 'PreCut', 0);

fprintf('\n[Diagnostic Result]\n');
fprintf('Slope of Pre-trends: %.4f (t=%.2f)\n', diag.slope, diag.slope_t_approx);
fprintf('Joint Test p-value:  %.4f\n', diag.pJoint_approx);
if diag.pJoint_approx < 0.05
    fprintf('VALIDITY WARNING: Significant pre-trends detected!\n');
else
    fprintf('VALIDITY CHECK: No significant pre-trends detected.\n');
end

% Plot Event Study
% figure('Name', 'Pre-trend Event Study');
% resES.plot(); % Wooldridge struct does not support .plot() method yet
% title('Event Study (Wooldridge)');

% Check if SEs are valid
fprintf('Event Study Table Head:\n');
disp(resES.summaryTable(1:min(10,height(resES.summaryTable)), :));

if any(isnan(resES.summaryTable.SE))
    warning('Wooldridge SEs contain NaNs.');
end

% --- 4. Standard DiD (TWFE) ---
fprintf('\n=========================================\n');
fprintf('      Standard DiD (TWFE Static)         \n');
fprintf('=========================================\n');

resTWFE = did.fit("twfe", ds, "Display", true);

% --- 5. Synthetic DiD (SDID) ---
fprintf('\n=========================================\n');
fprintf('      Synthetic DiD (SDID Static)        \n');
fprintf('=========================================\n');

% Using Placebo inference
try
    % Seed 1213 to match Stata example
    resSDID = did.fit("sdid", ds, "SEMethod", "Placebo", "B", 50, "Seed", 1213, "Display", true);
catch ME
    warning('SDID failed: %s', ME.message);
    resSDID = struct('tau', NaN, 'se', NaN);
end

% --- 6. Comparison ---
fprintf('\n=========================================\n');
fprintf('             COMPARISON                  \n');
fprintf('=========================================\n');

% Fix TWFE SE extraction:
twfe_se = NaN; twfe_ci = [NaN, NaN];
try
    if isfield(resTWFE, 'summaryTable') && height(resTWFE.summaryTable) > 0
        twfe_se = resTWFE.summaryTable.SE(1);
        twfe_ci = [resTWFE.summaryTable.Estimate(1) - 1.96*twfe_se, resTWFE.summaryTable.Estimate(1) + 1.96*twfe_se];
    end
catch
end

fprintf('%-15s %-10s %-10s %-20s\n', 'Method', 'Estimate', 'SE', '95% CI');
fprintf('%-15s %-10.4f %-10.4f [%.4f, %.4f]\n', ...
    "TWFE", resTWFE.ATT, twfe_se, twfe_ci(1), twfe_ci(2));

if ~isnan(resSDID.tau)
    fprintf('%-15s %-10.4f %-10.4f [%.4f, %.4f]\n', ...
        "SDID", resSDID.tau, resSDID.se, resSDID.tau - 1.96*resSDID.se, resSDID.tau + 1.96*resSDID.se);

    fprintf('\nInterpretation:\n');
    diff = abs(resTWFE.ATT - resSDID.tau);
    fprintf('Difference: %.4f\n', diff);

    % Plot SDID Curves
    try
        fprintf('Plotting SDID Trajectories...\n');
        if isfield(resSDID, 'Curves')
            fprintf('Curves Data: T=%d, Treated_Avg=%.4f\n', ...
                length(resSDID.Curves.Time), mean(resSDID.Curves.Treated));
            if all(isnan(resSDID.Curves.Treated))
                warning('Treated Curve is all NaN!');
            end
        else
            warning('resSDID.Curves field is missing.');
        end

        figure('Name', 'SDID Trajectories');
        did.estimators.SDID.plot(resSDID);
    catch ME
        warning('Plot failed: %s', ME.message);
    end
else
    fprintf('%-15s %-10s %-10s %-20s\n', "SDID", "N/A", "N/A", "No clean controls");
end
