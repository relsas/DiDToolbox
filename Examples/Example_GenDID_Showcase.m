% Example_GenDID_Showcase.m
% -------------------------------------------------------------------------
% Showcase of Synthetic Difference-in-Differences Data Generation
% This script illustrates all major options of did.genDIDdata.
% -------------------------------------------------------------------------
clear; clc;
fprintf('=== DID Toolbox Data Generation Showcase ===\n\n');

%% Scenario 1: Basic Staggered Adoption
% The standard case:
% - Units adopt treatment at different times.
% - Treatment effect is constant (ATT = 2).
% - Parallel trends hold.
fprintf('---------------------------------------------------\n');
fprintf('Scenario 1: Basic Single Cohort Treatment\n');
fprintf('Attributes: Random start times, Constant ATT=2\n');

T1 = did.genDIDdata(10, 200, 0.4, ...
    "treatType",   "constant", ...
    "startPeriod", 5, ...
    "ATT",         2, ...
    "Seed",        1);

% Quick check: Mean outcome for treated vs controls at end
ds1 = did.Dataset.fromTable(T1, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
did.did_plot(ds1, "trends");
title('Scenario 1: Basic Staggered');
pause(1); % Pause to let user see plot

%% Scenario 2: Dynamic Treatment Effects
% Effects grow over time ("event time").
% - "timeIncrease" adds `dynEffect * (t - g)` to the base ATT.
fprintf('\nScenario 2: Dynamic Treatment Effects\n');
fprintf('Attributes: Effect grows by 0.5 per period\n');

T2 = did.genDIDdata(10, 200, 0.4, ...
    "treatType",   "timeIncrease", ...
    "dynEffect",   0.5, ...
    "ATT",         1, ...  % Intercept at t=g
    "startPeriod", 5, ...
    "Seed",        2);

ds2 = did.Dataset.fromTable(T2, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
% Plot "event" view to see the slope
did.fit("TWFE", ds2, "Display", true, "preTrendMode", "event");
title('Scenario 2: Dynamic Effects (Event Study)');
pause(1);

%% Scenario 3: On/Off Treatment (Reversible)
% Units turn on treatment, then turn it off later.
fprintf('\nScenario 3: On/Off Treatment\n');
fprintf('Attributes: Treatment ends at period 15\n');

T3 = did.genDIDdata(10, 200, 0.4, ...
    "treatType",   "onOff", ...
    "startPeriod", 5, ...
    "endPeriod",   7, ...
    "ATT",         3, ...
    "Seed",        3);

ds3 = did.Dataset.fromTable(T3, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
did.did_plot(ds3, "trends");
title('Scenario 3: On/Off Treatment');
pause(1);

%% Scenario 4: Explicit Cohorts (Block Design)
% Instead of random staggered times, we specify exactly 3 start dates.
% Useful for replicating block-assignment designs.
fprintf('\nScenario 4: Explicit Cohorts (Block Design)\n');
fprintf('Attributes: Cohorts start at t=8, t=12, t=16\n');

T4 = did.genDIDdata(10, 300, 0.5, ...
    "treatType",   "constant", ...
    "cohortTimes", [8, 12, 16], ...
    "ATT",         2, ...
    "Seed",        4);

ds4 = did.Dataset.fromTable(T4, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
did.did_plot(ds4, "cohort");
title('Scenario 4: Explicit Cohorts');
pause(1);

%% Scenario 5: Divergent Controls (Violated Trends)
% A fraction of the control group is "bad" (divergent trend).
% Standard TWFE will fail here; SDID should perform better.
fprintf('\nScenario 5: Divergent Controls (Violated Trends)\n');
fprintf('Attributes: 50%% of controls have a trend slope +0.3\n');

T5 = did.genDIDdata(10, 500, 0.3, ...
    "preTrendType",           "divergentControls", ...
    "preTrendDivergentShare", 0.5, ...
    "preTrendMeanControl",    0.3, ... % Bad controls drift up
    "Seed",                   5);

ds5 = did.Dataset.fromTable(T5, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
did.did_plot(ds5, "trends");
title('Scenario 5: Divergent Controls');
pause(1);

%% Scenario 6: Unit-Specific Linear Trends
% Every unit has its own random time trend.
% Allows testing estimators robust to unit-trends.
fprintf('\nScenario 6: Unit-Specific Linear Trends\n');
fprintf('Attributes: Random slope per unit (SD=0.2)\n');

T6 = did.genDIDdata(10, 200, 0.4, ...
    "preTrendType", "unitLinear", ...
    "preTrendSd",   0.2, ...
    "Seed",         6);

ds6 = did.Dataset.fromTable(T6, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
did.did_plot(ds6, "trends");
title('Scenario 6: Unit Linear Trends');
pause(1);

%% Scenario 7: Selection Bias & Covariates
% Treatment probability depends on a covariate X.
% Outcome also depends on X.
% Naive comparison is biased.
fprintf('\nScenario 7: Selection Bias & Covariates\n');
fprintf('Attributes: Prob(Treat) depends on X; Y depends on X\n');

T7 = did.genDIDdata(10, 1000, 0.2, ...
    "xNum",          1, ... % 1 Covariate
    "SelectionBias", 2, ... % Strong selection on X
    "betaX",         1, ... % X affects Y
    "Seed",          7);

ds7 = did.Dataset.fromTable(T7, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");
fprintf('  - Treated Mean X: %.2f\n', mean(T7.x1(T7.D==1)));
fprintf('  - Control Mean X: %.2f\n', mean(T7.x1(T7.D==0)));
% Plot wouldn't show much diff visually, but stats confirm bias.

fprintf('\n=== Showcase Completed ===\n');
