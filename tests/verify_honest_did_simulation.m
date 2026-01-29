function verify_honest_did_simulation()
% VERIFY_HONEST_DID_SIMULATION
%
%   Compares two scenarios to prove HonestDiD works for both:
%     1. Single Cohort (Standard DiD)
%     2. Staggered Adoption (CS Aggregated Event Study)
%
% -------------------------------------------------------------------------
clear; close all; clc;

%% Scenario A: Single Cohort
fprintf('\n=================================================\n');
fprintf(' SCENARIO A: Single Cohort (Standard DiD)\n');
fprintf('=================================================\n');
run_scenario(1);

%% Scenario B: Staggered
fprintf('\n=================================================\n');
fprintf(' SCENARIO B: Staggered Adoption (3 Cohorts)\n');
fprintf('=================================================\n');
run_scenario(3);

end

function run_scenario(nCohorts)
N = 1000;
T_periods = 8;
trueATT = 1.0;

% Generate data
if nCohorts > 1
    tType = 'constantTime'; % Supports multiple cohorts
else
    tType = 'constant';     % Single cohort
end

T = did.genDIDdata(T_periods, N, 0.5, ...
    'treatType', tType, ...
    'numCohorts', nCohorts, ...
    'startPeriod', 3, ...
    'endPeriod', 8, ...
    'ATT', trueATT, ...
    'errorStd', 1.0, ...
    'Seed', 42 + nCohorts); % Diff seed for variety in staggered

% Calculate realized true ATT (weighted average on treated)
isTreated = (T.D == 1);
if any(isTreated)
    realizedATT = mean(T.ATT(isTreated));
else
    realizedATT = NaN;
end

ds = did.Dataset.fromTable(T, "yVar","y", "idVar","id", "timeVar","time", "dVar","D");
fprintf('Data: N=%d, Periods=%d, Cohorts=%d. Target ATT=%.2f, Realized ATT=%.3f\n', ...
    N, T_periods, nCohorts, trueATT, realizedATT);

% Update truth variable for checks
trueATT = realizedATT;

% 1. Run CS Estimator
fprintf('Running CS Estimator...\n');
% Suppress display for cleaner output
resCS = did.fit('CS', ds, 'Display', false, 'PreTrendBase', 'universal');

% 2. Honest DiD
hd = did.diagnostics.HonestDiD.fromFit(resCS);
fprintf('HonestDiD: NumPre=%d, NumPost=%d\n', hd.NumPre, hd.NumPost);

% 3. Check RM (M=0)
resRM = hd.fit('M', 0:0.1:2, 'DeltaType', 'RM', 'Alpha', 0.05);

if any(~isnan(resRM.CI_Lo))
    % Check M=0
    idx0 = find(resRM.M == 0, 1);
    if ~isempty(idx0)
        lb = resRM.CI_Lo(idx0);
        ub = resRM.CI_Hi(idx0);
        hit = (trueATT >= lb && trueATT <= ub);
        fprintf('[RM Result] M=0 Interval [%.3f, %.3f] (Width=%.3f). Contains Truth? %d\n', ...
            lb, ub, ub-lb, hit);
    else
        fprintf('[RM Result] M=0 not feasible.\n');
    end
else
    fprintf('[RM Result] No feasible bounds found.\n');
end

% 4. Check SDRM
resSDRM = hd.fit('M', 0:0.1:2, 'DeltaType', 'SDRM', 'Alpha', 0.05);
if any(~isnan(resSDRM.CI_Lo))
    idx0 = find(resSDRM.M == 0, 1);
    if ~isempty(idx0)
        fprintf('[SDRM Result] M=0 Interval [%.3f, %.3f].\n', resSDRM.CI_Lo(idx0), resSDRM.CI_Hi(idx0));
        hd.plot(resSDRM)
    end
else
    fprintf('[SDRM Result] No feasible bounds found.\n');
end

end
