%% MC Simulation: TWFE Event Pre-trends - Size vs Power
%
% This script assesses the performance of the TWFE Event Study estimator
% in detecting pre-trend violations.
%
% Scenarios:
% 1. H0 (Null): Parallel Trends hold. Measures "Size" (False Positive Rate).
% 2. H1 (Alternative): Treated units have a linear pre-trend. Measures "Power" (True Positive Rate).

clear; clc; close all;
rng(12345);

%% Simulation Parameters
M = 100;      % Number of repetitions (Higher = more precise size/power)
N = 500;      % Number of units
T_len = 20;   % Total periods
cohortTimes = [8 14]; % Staggered adoption

% H1 Violation Parameter:
% Treated units have a linear trend slope difference of 0.1 per period.
slope_violation = 0.1;

% Storage matrices [Iteration x CoeffIndex]
% We'll map coefficient names dynamically later.
resultsH0 = [];
resultsH1 = [];

fprintf('Running MC Simulation (M=%d)...\n', M);

for m = 1:M
    if mod(m, 10) == 0, fprintf('.'); end

    %% Scenario 1: H0 (Parallel Trends)
    T0 = did.genDIDdata(T_len, N, 0.5, ...
        "treatType", "constantTime", ...
        "cohortTimes", cohortTimes, ...
        "CohortIncrease", 1, ...
        "meanError", 0, "errorStd", 1.5, ...
        "preTrendType", "none"); % <--- Valid Parallel Trends

    ds0 = did.Dataset.fromTable(T0, "idVar","id", "timeVar","time", "dVar","D", "yVar","y", "Display",false);

    try
        res0 = did.fit("twfe", ds0, "preTrend",true, "preTrendMode","event", ...
            "baseline",-1, "details",false, "Display",false);

        % Store relevant stats (Name, Estimate, PValue)
        % We flatten the table to specific "Event k" rows later
        resultsH0{m} = res0.summaryTable;
    catch
        resultsH0{m} = [];
    end

    %% Scenario 2: H1 (Linear Trend Violation)
    T1 = did.genDIDdata(T_len, N, 0.5, ...
        "treatType", "constantTime", ...
        "cohortTimes", cohortTimes, ...
        "CohortIncrease", 1, ...
        "meanError", 0, "errorStd", 1.5, ...
        "preTrendType", "unitLinear", ...         % <--- Violation
        "preTrendMeanTreated", slope_violation, ...
        "preTrendMeanControl", 0);

    ds1 = did.Dataset.fromTable(T1, "idVar","id", "timeVar","time", "dVar","D", "yVar","y", "Display",false);

    try
        res1 = did.fit("twfe", ds1, "preTrend",true, "preTrendMode","event", ...
            "baseline",-1, "details",false, "Display",false);
        resultsH1{m} = res1.summaryTable;
    catch
        resultsH1{m} = [];
    end
end
fprintf('\nSimulation Complete.\n');

%% Analysis Function
function S = analyze_results(resultsCell, M)
% Extract all coefficient names
all_names = {};
for m = 1:M
    if ~isempty(resultsCell{m})
        all_names = [all_names; resultsCell{m}.Name];
    end
end
unique_names = unique(all_names);

% Filter for "Event k=..."
valid_names = strings(0);
k_vals = [];
for i = 1:numel(unique_names)
    s = string(unique_names{i});
    if contains(s, "Event k=")
        val = str2double(extractAfter(s, "Event k="));
        if ~isnan(val) && val < 0 % Pre-trends only
            k_vals(end+1) = val;
            valid_names(end+1) = s;
        end
    end
end
[k_vals, sIdx] = sort(k_vals);
valid_names = valid_names(sIdx);

% Build stats
count = zeros(numel(valid_names), 1);
reject_count = zeros(numel(valid_names), 1);
sum_est = zeros(numel(valid_names), 1);

for m = 1:M
    Tbl = resultsCell{m};
    if isempty(Tbl), continue; end

    for j = 1:numel(valid_names)
        idx = find(string(Tbl.Name) == valid_names(j), 1);
        if ~isempty(idx)
            est = Tbl.Estimate(idx);
            pval = Tbl.pValue(idx);

            sum_est(j) = sum_est(j) + est;
            count(j) = count(j) + 1;

            if pval < 0.05
                reject_count(j) = reject_count(j) + 1;
            end
        end
    end
end

AvgEst = sum_est ./ count;
RejRate = reject_count ./ count;

S = table(k_vals(:), AvgEst, RejRate, count, ...
    'VariableNames', {'k', 'Avg_Estimate', 'RejectionRate', 'N_Sims'});
end

statsH0 = analyze_results(resultsH0, M);
statsH1 = analyze_results(resultsH1, M);

%% Joint Reporting
% Merge tables on k
finalTbl = table(statsH0.k, statsH0.Avg_Estimate, statsH0.RejectionRate, ...
    statsH1.Avg_Estimate, statsH1.RejectionRate, ...
    'VariableNames', {'k', 'H0_AvgEst', 'Size_Rate', 'H1_AvgEst', 'Power_Rate'});

fprintf('\n=== MC Results: Size vs Power ===\n');
fprintf('Nominal Size: 0.05\n');
fprintf('H1 Slope Violation: Base + %.2f * t (for Treated)\n', slope_violation);
disp(finalTbl);

%% Visualization
figure('Name', 'MC Pre-trend Power Analysis', 'Color','w');

subplot(1,2,1);
plot(statsH0.k, statsH0.RejectionRate, '-o', 'LineWidth', 2, 'DisplayName', 'H0 (Size)');
yline(0.05, '--k', 'Nominal 5%');
ylim([0 1]);
xlabel('Event Time (k)');
ylabel('Rejection Rate (Pr p < 0.05)');
title('Size: False Positive Rate');
grid on;

subplot(1,2,2);
plot(statsH1.k, statsH1.RejectionRate, '-rs', 'LineWidth', 2, 'DisplayName', 'H1 (Power)');
yline(0.05, '--k', 'Nominal 5%');
ylim([0 1]);
xlabel('Event Time (k)');
ylabel('Rejection Rate (Pr p < 0.05)');
title('Power: True Positive Rate');
legend('Location','best');
grid on;

sgtitle('Testing for Pre-trends: Size vs Power');
