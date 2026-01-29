
%% Example_CEM_Showcase.m
% -------------------------------------------------------------------------
% Illustration of Selection Bias and CEM Correction
% -------------------------------------------------------------------------
% This script demonstrates:
% 1. Simulating data where treatment probability depends on a covariate X.
% 2. Showing how naive TWFE is biased because Treated/Control differ on X.
% 3. Using Coarsened Exact Matching (CEM) to re-weight the control group.
% 4. Validating the reduction in bias.
% -------------------------------------------------------------------------

clear; clc;
rng(42);

fprintf('=== CEM Simulation Showcase ===\n\n');

%% 1. Generate Data with Selection Bias
%
% Attributes:
% - SelectionBias=2.5: Strong dependency (High X -> High Prob(Treat))
% - betaX=2: X strongly affects Outcome Y
% - True ATT = 2
% - Naive comparison will be biased upwards because Treated have high X -> high Y.

fprintf('1. Generating Data...\n');
N_units = 1000;
T_periods = 10;
TRUE_ATT = 2

T = did.genDIDdata(T_periods, N_units, 0.3, ...
    "treatType",     "constant", ...
    "startPeriod",   6, ...
    "ATT",           TRUE_ATT, ...
    "xNum",          1, ...   % One covariate
    "xUnitStd",      1.5, ... % Variation in X across units
    "SelectionBias", 2.5, ... % Bias: Treatment depends on X
    "betaX",         2, ...   % Confounding: Outcome depends on X
    "TrendEffectX",  0.5, ... % Differential Trends: X affects growth
    "Seed",          101);

% Convert to Dataset
ds = did.Dataset.fromTable(T, "idVar","id", "timeVar","time", "yVar","y", "dVar","D");

% Verify Imbalance (Using ever-treated status at t=1 base period)
% Note: T.D is time-varying treatment. We want to compare groups based on assignment.
is_treated = (T.everTreated == 1);
% Take one slice (t=1)
X_tr = T.x1(is_treated & T.time==1);
X_co = T.x1(~is_treated & T.time==1);

fprintf('   - Treated Units: %d\n', numel(unique(T.id(T.everTreated==1))));
fprintf('   - Control Units: %d\n', numel(unique(T.id(T.everTreated==0))));
fprintf('   - Mean X (Treated): %.3f\n', mean(X_tr));
fprintf('   - Mean X (Control): %.3f\n', mean(X_co));
fprintf('   -> Selection Bias generates imbalance in X!\n\n');


%% 2. Naive Estimation (Standard TWFE)
% Since X affects Y and X varies between groups, parallel trends may hold
% conditionally but the levels are very different, and if effects are heterogeneous
% or if we just want to control for X, TWFE without X might be biased
% (or at least inefficient).
% Actually in this simulation, since X is time-invariant unit-fixed effect,
% strictly speaking standard TWFE handles the intercept difference.
%
% BUT: If trend depends on X (TrendEffectX), then parallel trends is violated.
% Let's assume for this example we just want to show balancing.
% Or let's imply a trend bias.
% Re-generate with TrendBias for a stronger case?
% Let's stick to the current one. If betaX is large, Unit Fixed Effects
% handle the level difference. Bias in TWFE usually comes from
% differential trends or heterogeneous effects.
%
% However, CEM is often used to ensure "Common Support".
% Let's proceed and see. If standard TWFE works perfectly, we need to add TrendEffectX.

fprintf('2. Estimating Naive TWFE...\n');
res_twfe = did.fit("TWFE", ds, "Display", false);
% Extract TWFE estimate (usually first row if simple)
if isfield(res_twfe, 'summaryTable')
    % Find "D" or treatment var
    idx = contains(res_twfe.summaryTable.Name, "D");
    if any(idx)
        est = res_twfe.summaryTable.Estimate(idx);
        bias_twfe = est(1) - TRUE_ATT;
        fprintf('   - Naive TWFE Estimate: %.4f (Bias: %.4f)\n', est(1), bias_twfe);
    else
        fprintf('   - Naive TWFE: Could not find treatment coeff.\n');
    end
else
    % Fallback if older structure
    fprintf('   - Naive TWFE: Structure unseen.\n');
end


%% 3. CEM Estimation
% We explicitly match on "x1".
% - nBins: Number of bins to coarsen X into.
% - Robust: If true, runs TWFE on matched sample (Weighted TWFE).

fprintf('\n3a. Estimating CEM (Standard - Weights Only)...\n');
res_cem_base = did.fit("CEM", ds, "Covariates", "x1", "nBins", 10, "Robust", false, "Display", false);
bias_cem_base = res_cem_base.tau - TRUE_ATT;
fprintf('   - CEM (Weights Only):  %.4f (Bias: %.4f)\n', res_cem_base.tau, bias_cem_base);

fprintf('\n3b. Estimating CEM (Robust - Weights + Control)...\n');
res_cem = did.fit("CEM", ds, ...
    "Covariates", "x1", ...
    "nBins",      10, ...
    "Robust",     true, ...
    "Display",    true);
bias_cem = res_cem.tau - TRUE_ATT;
fprintf('   - CEM (Robust):        %.4f (Bias: %.4f)\n', res_cem.tau, bias_cem);


%% 4. Analyze Weights
% CEM assigns weights to Control units to match the Treated distribution.
% Treated units get weight 1 (if matched). Unmatched get 0.

% Visualize Weights
% Check for 'WeightsVector' (standard output) or 'weights'
weights_found = [];
if isfield(res_cem, 'WeightsVector')
    weights_found = res_cem.WeightsVector;
elseif isfield(res_cem, 'weights')
    weights_found = res_cem.weights;
end

if ~isempty(weights_found)
    fprintf('\nVisualizing weights...\n');

    % Map weights back to full T (N x 1)
    if isfield(res_cem, 'MatchedIDs')
        w_ids = res_cem.MatchedIDs;
        w_vec = res_cem.WeightsVector;

        % Create map: ID -> Weight
        % Since MatchedIDs corresponds to observations in the matched sample (inner join),
        % and weights are unit-specific (constant for unit).
        % We extract unique unit weights.

        [uIDs, idx] = unique(w_ids);
        uWeights = w_vec(idx);

        % Find indices in T to get Treatment Status (from everTreated)
        [Lia, Locb] = ismember(uIDs, T.id);
        if ~all(Lia)
            warning('Some matched IDs not found in T?');
        end

        uTreated = T.everTreated(Locb);

        % Extract X for scatter plot
        uX       = T.x1(Locb);
        uX_tr    = uX(uTreated==1);
        uX_co    = uX(uTreated==0);

        w_tr     = uWeights(uTreated==1);
        w_co     = uWeights(uTreated==0);

        figure('Name', 'CEM Weight Distribution', 'Color', 'w');

        subplot(2,1,1);
        histogram(X_tr, 20, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'Normalization', 'pdf');
        hold on;
        histogram(X_co, 20, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'Normalization', 'pdf');
        title('Original Covariate Balance (Unweighted)');
        legend('Treated', 'Control');
        ylabel('Density');
        xlabel('Covariate X1');

        subplot(2,1,2);
        % Scatter Plot: Weight vs X
        scatter(uX_co, w_co, 25, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
        hold on;
        scatter(uX_tr, w_tr, 25, 'r', 'filled', 'MarkerFaceAlpha', 0.6);

        title('Weight vs Covariate (Individual Units)');
        ylabel('Weight (w)');
        xlabel('Covariate X1');
        legend('Control', 'Treated (w=1)', 'Location', 'best');
        grid on;

    else
        fprintf('Warning: MatchedIDs not found in result. Skipping detailed plot.\n');
    end

end

fprintf('\nDone.\n');
