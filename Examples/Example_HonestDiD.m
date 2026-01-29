%% Honest DiD with Real Data (Medicaid Expansion)
% This script demonstrates Honest DiD sensitivity analysis using the
% Callaway & Sant'Anna (2021) estimator on the Medicaid Expansion dataset.
%
% We will:
% 1. Estimate dynamic effects (Event Study) using CS.
% 2. Check sensitivity of the results to parallel trend violations.
% 3. Illustrate both "Relative Magnitude" (RM) and "Smoothness" (SD) constraints.
% 4. Filter and interpret valid sensitivity bounds.

clear; close all;
% Force clear to ensure fresh code


%% 1. Load and Prepare Data
fprintf('--- Loading Medicaid Data ---\n');
load('.\Data\Medicaid.mat'); % Variables: outcome=crude_rate_20_64

% Create Dataset
ds = did.Dataset.fromTable(data, ...
    'idVar', 'county_code', ...
    'timeVar','year_code', ...
    'yVar',   'crude_rate_20_64', ...
    'dVar',   'treatPost');

%% 2. Estimate Initial Event Study (Callaway & Sant'Anna)
fprintf('\n--- Running CS Estimator (may take a moment) ---\n');
% We use CS to get Group-Time Average Treatment Effects aggregated into an Event Study
% User request: Drop first pre-period (2009) as base, NOT 2013.
resCS = did.fit('CS', ds, 'Display', true, 'PreTrendBase', 'first');

%% 3. Honest DiD Analysis
fprintf('\n--- Initializing Honest DiD ---\n');
hd = did.diagnostics.HonestDiD.fromFit(resCS);

% Define broad sensitivity grid (M)
M_grid = 0:0.05:2;

%% 4. Sensitivity Analysis 1: Relative Magnitude (RM)
% Constraint: Max post-period bias <= M * Max pre-period bias
fprintf('\n--- Analysis 1: Relative Magnitude (RM) ---\n');
resRM = hd.fit('M', M_grid, ...
    'DeltaType', 'RM', ...
    'Target', 'last', ... % Effect at last event time
    'Alpha', 0.05);

% Display feasible results (HonestDiD automatically filters NaNs)
fprintf('Results table (M, CI_Lo, CI_Hi):\n');

nShow = min(10, numel(resRM.M));
if nShow > 0
    disp(table(resRM.M(1:nShow), resRM.CI_Lo(1:nShow), resRM.CI_Hi(1:nShow), 'VariableNames',{'M','LB','UB'}));
else
    fprintf('No feasible M values found for RM.\n');
end

% Plot
hd.plot(resRM);
title('Honest DiD (Relative Magnitude)');

%% 5. Sensitivity Analysis 2: Smoothness (SD)
% Constraint: Bias cannot change more than M between periods.
fprintf('\n--- Analysis 2: Smoothness (SD) ---\n');
M_grid_SD = 0:0.01:0.5;

resSDRM = hd.fit('M', M_grid_SD, ...
    'DeltaType', 'SDRM', ...
    'Target', 'last', ...
    'Alpha', 0.05);

fprintf('Results table (SDRM) (Partial):\n');
nShow = min(5, numel(resSDRM.M));
if nShow > 0
    disp(table(resSDRM.M(1:nShow), resSDRM.CI_Lo(1:nShow), resSDRM.CI_Hi(1:nShow), 'VariableNames',{'M','LB','UB'}));
else
    fprintf('No feasible M values found for SD.\n');
end




%% Honest SDRM
resRM = hd.fit('M', M_grid, ...
    'DeltaType', 'SDRM', ...
    'Target', 'last', ... % Effect at last event time
    'Alpha', 0.05);




% Comparison
fprintf('Valid M count: RM=%d, SDRM=%d\n', numel(resRM.M), numel(resSDRM.M));
