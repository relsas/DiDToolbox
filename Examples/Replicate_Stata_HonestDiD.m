%% Replicate Stata HonestDiD Example (Full: Baseline + Staggered + Wooldridge)
%
% This script replicates both the Baseline (2014-only) and Staggered
% analyses from the Stata HonestDiD README.
%
% Dataset: \Data\medic_mixtape.xlsx
%
% -------------------------------------------------------------------------
clear; try clear classes; catch; end; close all; clc; rehash toolbox;
addpath(genpath('../'));
dataFile = 'd:\Temp\MCP\DID\did-toolbox\Data\medic_mixtape.xlsx';

%% PART 1: Baseline DiD (Only 2014 vs Never Treated)
fprintf('\n=========================================\n');
fprintf(' PART 1: Baseline DiD (2014 vs Never)\n');
fprintf('=========================================\n');

% 1. Load & Prepare Data (Baseline Filter)
data = readtable(dataFile);
if iscell(data.year), data.year = str2double(data.year); elseif isstring(data.year), data.year = double(data.year); end
if ismember('yexp2', data.Properties.VariableNames)
    if iscell(data.yexp2), data.yexp2 = str2double(data.yexp2); elseif isstring(data.yexp2), data.yexp2 = double(data.yexp2); end
    data.g = data.yexp2; data.g(isnan(data.g)) = 0;
else
    error('Variable "yexp2" not found.');
end

% Filter: Year < 2016, Drop treated in 2015
idx_year = (data.year <= 2015);
idx_g = (data.g ~= 2015);
data_base = data(idx_year & idx_g, :);
data_base.g(data_base.g > 2015) = 0;
data_base.D = (data_base.year >= data_base.g) & (data_base.g > 0);

ds_base = did.Dataset.fromTable(data_base, "yVar", "dins", "idVar", "stfips", "timeVar", "year", "dVar", "D", "gVar", "g");

% 2. Estimate CS (Baseline)
% PreTrendBase='universal' (2013 reference for 2014 cohort)
resBase = did.fit('CS', ds_base, 'Display', true, 'PreTrendBase', 'universal', ...
    'SEMethod', 'clustered', 'ClusterVar', 'stfips');

% Identify time variable (t)
tVar = 't';
if ~ismember(tVar, resBase.summaryTable.Properties.VariableNames)
    if ismember('year', resBase.summaryTable.Properties.VariableNames), tVar = 'year';
    elseif ismember('Year', resBase.summaryTable.Properties.VariableNames), tVar = 'Year';
    end
end
fprintf('Baseline Estimate (2014): \n');
try disp(resBase.summaryTable(resBase.summaryTable.g==2014 & resBase.summaryTable.(tVar)==2014, :)); catch; end

% 3. HonestDiD (RM & SD) for Baseline
hdBase = did.diagnostics.HonestDiD.fromFit(resBase);

% TRIM NaNs
valid_idx = ~isnan(hdBase.Betas);
if any(~valid_idx)
    hdBase.Betas = hdBase.Betas(valid_idx);
    hdBase.Sigma = hdBase.Sigma(valid_idx, valid_idx);
    hdBase.EventTimes = hdBase.EventTimes(valid_idx);
    hdBase.NumPre = sum(hdBase.EventTimes < 0);
    hdBase.NumPost = sum(hdBase.EventTimes >= 0);
end

target_idx = find(hdBase.EventTimes(hdBase.NumPre+1:end) == 0, 1) + hdBase.NumPre;
if isempty(target_idx), target_idx = hdBase.NumPre + 1; end
l_vec_base = zeros(hdBase.NumPost, 1);
l_vec_base(target_idx - hdBase.NumPre) = 1;

fprintf('\n--- Baseline Sensitivity (RM) ---\n');
outRM_base = hdBase.fit('M', [0:0.5:2], 'DeltaType', 'RM', 'Target', l_vec_base, 'Alpha', 0.05);
disp(outRM_base.summaryTable);

fprintf('\n--- Baseline Sensitivity (SD) ---\n');
outSD_base = hdBase.fit('M', [0:0.01:0.05], 'DeltaType', 'SD', 'Target', l_vec_base, 'Alpha', 0.05);
disp(outSD_base.summaryTable);

figure('Name', 'Baseline RM'); hdBase.plot(outRM_base, gca); title('Baseline (2014): RM Sensitivity');
figure('Name', 'Baseline SD'); hdBase.plot(outSD_base, gca); title('Baseline (2014): SD Sensitivity');


%% PART 2: Staggered DiD (Full Sample)
fprintf('\n=========================================\n');
fprintf(' PART 2: Staggered DiD (Full Sample)\n');
fprintf('=========================================\n');

% 1. Load Full Data (No filtering besides cleaning)
data = readtable(dataFile);
if iscell(data.year), data.year = str2double(data.year); elseif isstring(data.year), data.year = double(data.year); end
if ismember('yexp2', data.Properties.VariableNames)
    if iscell(data.yexp2), data.yexp2 = str2double(data.yexp2); elseif isstring(data.yexp2), data.yexp2 = double(data.yexp2); end
    data.g = data.yexp2; data.g(isnan(data.g)) = 0;
end
data.D = (data.year >= data.g) & (data.g > 0);

ds_stag = did.Dataset.fromTable(data, "yVar", "dins", "idVar", "stfips", "timeVar", "year", "dVar", "D", "gVar", "g");

% 2. Estimate CS (Staggered)
resStag = did.fit('CS', ds_stag, 'Display', true, 'PreTrendBase', 'varying', ...
    'SEMethod', 'clustered', 'ClusterVar', 'stfips');

% 3. HonestDiD (Staggered)
hdStag = did.diagnostics.HonestDiD.fromFit(resStag);

% TRIM NaNs
valid_idx_stag = ~isnan(hdStag.Betas);
if any(~valid_idx_stag)
    hdStag.Betas = hdStag.Betas(valid_idx_stag);
    hdStag.Sigma = hdStag.Sigma(valid_idx_stag, valid_idx_stag);
    hdStag.EventTimes = hdStag.EventTimes(valid_idx_stag);
    hdStag.NumPre = sum(hdStag.EventTimes < 0);
    hdStag.NumPost = sum(hdStag.EventTimes >= 0);
end

l_vec_stag = zeros(hdStag.NumPost, 1);
idx0_stag = find(hdStag.EventTimes(hdStag.NumPre+1:end) == 0, 1);
if ~isempty(idx0_stag), l_vec_stag(idx0_stag) = 1; else, l_vec_stag(1) = 1; end

fprintf('\n--- Staggered Sensitivity (RM) ---\n');
outRM_stag = hdStag.fit('M', [0:0.5:2], 'DeltaType', 'RM', 'Target', l_vec_stag, 'Alpha', 0.05);
disp(outRM_stag.summaryTable);

fprintf('\n--- Staggered Sensitivity (SD) ---\n');
outSD_stag = hdStag.fit('M', [0:0.01:0.05], 'DeltaType', 'SD', 'Target', l_vec_stag, 'Alpha', 0.05);
disp(outSD_stag.summaryTable);

figure('Name', 'Staggered RM'); hdStag.plot(outRM_stag, gca); title('Staggered: RM Sensitivity (e=0)');
figure('Name', 'Staggered SD'); hdStag.plot(outSD_stag, gca); title('Staggered: SD Sensitivity (e=0)');


%% PART 3: Wooldridge (Mundlak) Estimator
fprintf('\n=========================================\n');
fprintf(' PART 3: Wooldridge Estimator\n');
fprintf('=========================================\n');

% Estimate Wooldridge (Enable EventStudy for Leads/Lags)
resW = did.fit('Wooldridge', ds_stag, 'Display', true, 'EventStudy', true);

% Wooldridge output has 'details.attByEventTime' (EventStudy table)
% HonestDiD.fromFit can read a table with EventTime, Estimate, SE
esW = resW.details.attByEventTime;

fprintf('\nWooldridge Event Study (First 5 rows):\n');
try disp(esW(1:min(5,height(esW)), :)); catch; end

hdW = did.diagnostics.HonestDiD.fromFit(esW);

% TRIM NaNs (should be rare in Wooldridge OLS unless singular)
valid_idx_W = ~isnan(hdW.Betas);
if any(~valid_idx_W)
    hdW.Betas = hdW.Betas(valid_idx_W);
    hdW.Sigma = hdW.Sigma(valid_idx_W, valid_idx_W);
    hdW.EventTimes = hdW.EventTimes(valid_idx_W);
    hdW.NumPre = sum(hdW.EventTimes < 0);
    hdW.NumPost = sum(hdW.EventTimes >= 0);
end

l_vec_W = zeros(hdW.NumPost, 1);
idx0_W = find(hdW.EventTimes(hdW.NumPre+1:end) == 0, 1);
if ~isempty(idx0_W), l_vec_W(idx0_W) = 1; else, l_vec_W(1) = 1; end

fprintf('\n--- Wooldridge Sensitivity (RM) ---\n');
outRM_W = hdW.fit('M', [0:0.5:2], 'DeltaType', 'RM', 'Target', l_vec_W, 'Alpha', 0.05);
disp(outRM_W.summaryTable);

figure('Name', 'Wooldridge RM'); hdW.plot(outRM_W, gca); title('Wooldridge: RM Sensitivity (e=0)');
figure('Name', 'Wooldridge SD');
outSD_W = hdW.fit('M', [0:0.01:0.05], 'DeltaType', 'SD', 'Target', l_vec_W, 'Alpha', 0.05);
hdW.plot(outSD_W, gca); title('Wooldridge: SD Sensitivity (e=0)');

fprintf('Done.\n');


%% =========================================
% PART 4: DID_M (Chaisemartin & d'Haultfoeuille)
% =========================================
fprintf('\n=========================================\n');
fprintf(' PART 4: DID_M Estimator (Dynamic)\n');
fprintf('=========================================\n');

% Estimate DID_M with Horizons and Placebos
% We ask for lags 0..3 and 2 placebos (-1,-2)
% Note: DID_M is computationally intensive; B=2 for quick check
resM = did.fit('CH', ds_stag, 'Display',true, 'Horizons',[0 1 2 3], 'Placebos',2, 'B',2);

fprintf('\nDID_M Event Study:\n');
try disp(resM.EventStudy); catch; end

% HonestDiD Analysis
% We need to check if we have enough pre-periods. Placebos=2 gives e=-1,-2.
% HonestDiD requires ref period e=-1 to be 0 (or dropped).
% Our DID_M code generates e=-1 estimate.
% HonestDiD might need anchoring.
% Let's see if fromFit handles it.
% If e=-1 is estimated, HonestDiD.fit might complain "Reference period required".
% Unless we set PreTrendBase logic. But HonestDiD logic anchors at e=-1 if it finds a zero there.
% If it doesn't find a zero, it might error.
% We might need to manually insert a reference row if DID_M estimates everything relative to something else.
% But DID_M Placebos ARE deviations from trend.
% Let's run and see.

fprintf('\n--- DID_M Sensitivity Analysis ---\n');
try
    hdM = did.diagnostics.HonestDiD.fromFit(resM);

    % Requires at least 2 pre-periods for RM (trend) and usually 3 points for SD (curvature).
    % The 'medic' dataset with this DID_M setup only yields e=-1.
    if hdM.NumPre >= 2
        fprintf('Running RM and SD analysis...\n');

        % RM
        outRM_M = hdM.fit('M', [0:0.5:2], 'DeltaType', 'RM', 'Alpha', 0.05);
        disp(outRM_M.summaryTable);
        figure('Name','DID_M RM'); hdM.plot(outRM_M, gca); title('DID_M RM');

        % SD
        outSD_M = hdM.fit('M', [0:0.01:0.05], 'DeltaType', 'SD', 'Alpha', 0.05);
        disp(outSD_M.summaryTable);
        figure('Name','DID_M SD'); hdM.plot(outSD_M, gca); title('DID_M SD');

    else
        fprintf('NOTE: Skipping HonestDiD for DID_M.\n');
        fprintf('      Reason: Insufficient pre-treatment coefficients (NumPre=%d).\n', hdM.NumPre);
        fprintf('      Sensitivity analysis (RM/SD) requires at least 2 pre-periods to estimate trends.\n');
    end

catch ME
    fprintf('Error in HonestDiD DID_M: %s\n', ME.message);
end

fprintf('Done with PART 4.\n');
