% Example_Trends.m
% Demonstrates how to visualize raw outcome trends for multiple cohorts
% using the "trends" option in did_plot.

clear; rng(2025);

%% 1. Simulate Data
% We generate data with:
% - 3 treated cohorts (plus never-treated)
% - "constantTime" treatment effect (same effect at same calendar time)
% - Parallel trends hold (no divergent controls)
fprintf('Simulating data...\n');
opts = struct();
opts.numCohorts = 3;           % 3 distinct cohorts
opts.treatType  = "constantTime";
opts.ATT        = 2;           % True effect size
opts.startPeriod = 10;
opts.numPeriods = 30;
opts.numIds     = 500;

% Generate table
T = did.genDIDdata(opts.numPeriods, opts.numIds, 0.4, ...
    "numCohorts", opts.numCohorts, ...
    "treatType", opts.treatType, ...
    "ATT", opts.ATT, ...
    "startPeriod", opts.startPeriod);

%% 2. Wrap in did.Dataset
% The "trends" plot requires a did.Dataset object to handle grouping and labeling automatically.
fprintf('Creating Dataset object...\n');
ds = did.Dataset.fromTable(T, "idVar", "id", "timeVar", "time", "yVar", "y", "dVar", "D");

%% 3. Plot Trends
% Plot the raw outcome trends for each cohort (and never-treated).
% The function automatically handles styling and legend.
fprintf('Plotting trends...\n');
figure('Color','w', 'Name', 'DiD Trends Example');
did.did_plot(ds, "trends");

fprintf('Done.\n');
