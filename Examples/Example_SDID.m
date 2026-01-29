%% Example: Synthetic Difference-in-Differences (SDID)
% This script demonstrates how to use the SDID estimator with simulated data.
% We generate a dataset where Parallel Trends might be violated or simply
% noisy, and apply SDID to find an optimal weighting.

clear; clc;
rng(42); % Reproducibility

% Ensure package is on path
scriptDir = fileparts(mfilename('fullpath'));
toolboxDir = fileparts(scriptDir);
addpath(toolboxDir);

%% 1. Simulate Data (Block Design)
% N = 200 units, T = 30 periods.
% Treatment starts at T=16 for the treated group.
N = 200;
T_len = 10;
N_treated = 50;

% Use genDIDdata to create a basic structure
% TreatType="constantTime" -> Block design (all treated units start at same time)
cohortTimes = [3,5,7];

% We introduce a violation of Parallel Trends:
% Treated units have a linear upward trend (slope=0.5) distinct from control (slope=0).
% TWFE should be biased. SDID should reweight to zero out this trend.
data = did.genDIDdata(T_len, N, N_treated/N, ...
    "treatType", "constantTime", ...
    "cohortTimes", cohortTimes, ...
    "meanError", 0, "errorStd", 0.5, ...
    "preTrendType", "divergentControls", ...
    "preTrendSd", 0.2, ...            % Add variance (the key to the "distributed" look)
    "preTrendMeanTreated", 0.5, ...   % Treated units (and "Good" controls)
    "preTrendMeanControl", -0.2);     % "Bad" divergent controls

% Convert to did.Dataset
ds = did.Dataset.fromTable(data, ...
    "idVar", "id", "timeVar", "time", "dVar", "D", "yVar", "y", ...
    "Display", true);

% Extract Simulation Truth (Latent Slopes) for Validation Plots
% We need one value per unit (using the first observation per unit)
[~, uniqueIdx] = unique(data.id);
% Sort by ID to match did.Dataset internal ordering if needed, but unique returns in order of appearance usually.
% Let's be safe: sort uniqueIdx
% The simulation generates IDs 1..N sequentially, so unique is safe.
true_slopes = data.latent_slope(uniqueIdx);
is_treated  = data.everTreated(uniqueIdx) == 1;

fprintf('\nData generated. Treatment starts at t=%d.\n', cohortTimes(1));

%% 2. Fit SDID Estimator
% SDID reweights control units to match the trends of treated units.
fprintf('\n--- Running SDID ---\n');

res_sdid = did.fit("sdid", ds, ...
    "Display", true, ...
    "SEMethod", "Placebo", ...   % Options: "Placebo", "Clustered"
    "B", 50);                    % Low bootstrap reps for example speed

%% 3. Fit Standard TWFE for Comparison
fprintf('\n--- Running TWFE (Comparison) ---\n');
res_twfe = did.fit("twfe", ds, "Display", true);

%% 4. Display Summary Comparison
fprintf('\n=== Results Comparison ===\n');
fprintf('True Effect (approx): 2.0 (from genDIDdata default)\n');
fprintf('SDID Estimate: %.4f (SE: %.4f)\n', res_sdid.coef.Estimate, res_sdid.coef.SE);
fprintf('TWFE Estimate: %.4f (SE: %.4f)\n', res_twfe.coef.Estimate, res_twfe.coef.SE);

%% 5. Plotting
% The SDID estimator provides a static plot method
fprintf('\nGenerating SDID plots...\n');
did.estimators.SDID.plot(res_sdid);


%%
% -------------------------------------------------------------------------
% SDID provides unit weights (omega).
% Hypothesis: SDID should assign HIGH weights to Control units that, by chance,
% have steep slopes (similar to Treated), and LOW weights to flat controls.

% Extract Weights for the first (and only) cohort
% Since we have simultaneous adoption here for illustration, weights are simple.
% SDID.fit for staggered aggregates, but we can look at the first cohort details.
% Extract Weights for the first (and only) cohort
% Since we have simultaneous adoption here for illustration, weights are simple.
% SDID.fit for staggered aggregates, but we can look at the first cohort details.
cDetails = res_sdid.Details.Cohorts(1);
fprintf('\nVisualizing details for Cohort starting at G=%d\n', cDetails.G);
fprintf('(Note: In staggered designs, weights differ by cohort. Showing the first one.)\n');

% Get Control indices and their weights
co_idx = cDetails.Controls;  % Indices in the Y matrix
weights = cDetails.Omega;    % Vector of weights summing to 1

% Get the TRUE slopes for these control units
% Note: co_idx contains the indices of the Control units in the original Y matrix.
% We must ensure 'true_slopes' is indexed similarly.
% did.Dataset typically preserves the order of IDs if they are numeric sortable.
co_slopes = true_slopes(co_idx);

% Plot: Weight vs True Slope
figure('Name','SDID Unit Weights Analysis','Color','w','Position',[100 100 1200 500]);

subplot(1, 2, 1);
scatter(co_slopes, weights, 50, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('True Latent Slope (Control Units)');
ylabel('SDID Assigned Weight (\omega)');
title(['Weight vs. Latent Trend (Cohort G=' num2str(cDetails.G) ')']);
grid on;
xline(mean(true_slopes(is_treated)), '--r', 'Avg Treated Slope', 'LineWidth', 2);
subtitle('SDID picks controls with slopes matching the Treated group');
annotation('textbox', [0.15 0.7 0.3 0.1], 'String', 'Bimodal weights expected:\nHigh for Good Match, Low for Bad', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'BackgroundColor', 'w');

% -------------------------------------------------------------------------
% Interpretation of Time Weights (Lambda)
% -------------------------------------------------------------------------
% Time weights (\lambda) emphasize pre-treatment periods that best predict
% the post-treatment average.
% In a pure linear trend case, this might be balanced.

lambda = cDetails.Lambda;
% Lambda corresponds to the pre-treatment periods used in estimation.
% Length should be G-1?
% Let's make x-axis dynamic based on Lambda length.
pre_times = (1:length(lambda))';

subplot(1, 2, 2);
bar(pre_times, lambda, 'FaceColor', [0.2 0.6 0.6]);
xlabel('Pre-Treatment Time Period');
ylabel('Time Weight (\lambda)');
title(['SDID Time Weights (Cohort G=' num2str(cDetails.G) ')']);
grid on;
xline(cDetails.G, '-', ['Treatment ' num2str(cDetails.G)]);

sgtitle('Anatomy of Synthetic DiD (Staggered Example)');