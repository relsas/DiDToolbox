function out = wooldridge_TB(tbl, args)
% WOOLDRIDGE_TB  Staggered DiD (Wooldridge two-way Mundlak) via pooled OLS.
%
%   out = wooldridge_TB(tbl, yvar="y", idvar="id", tvar="time", gvar="g", ...
%                       Cluster=true, Covariates=[], DropTimeBase=1, details=false,Print=true)
%
% Coefficients on I{cohort=r}×I{t=j} (for j>=r) are the cohort/time ATTs.
%
% NOTE (Toolbox integration):
% - Extra NV-args 'vcov' is accepted for did.fit compatibility but ignored;
%   standard errors are governed by 'Cluster' (true = clustered) and 'clusters'.
% - 'clusters' (name or vector) IS USED when Cluster=true; defaults to unit id.
% 
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/30/2025  
% ------------------------------------------------------------------------
arguments
    tbl table
    args.yVar (1,1) string = "y"
    args.idVar (1,1) string = "id"
    args.timeVar (1,1) string = "time"
    args.gVar (1,1) string = "g"
    args.Covariates (1,:) string = string.empty
    args.DropTimeBase (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(args.DropTimeBase,1)} = 1
    args.details (1,1) logical = false % prints further details 
    args.Print logical=true  % prints just standard results
    % --- Minimal additions for toolbox integration ---
    args.vcov (1,1) string = "clustered"
    args.clusters = []                 % name or vector aligned with rows (any type)
end

% ----- Columns -----
y    = tbl.(args.yVar);
id   = categorical(tbl.(args.idVar));                       % grouping-friendly categoricals
tRaw = tbl.(args.timeVar);
gRaw = tbl.(args.gVar);
n    = height(tbl);

% Resolve cluster variable (name or vector); default to id
if isempty(args.clusters)
    clust = id;  % default: cluster by unit id
else
        cn = string(args.clusters);
        if all(ismember(cn, tbl.Properties.VariableNames))
            clust = args.clusters ;
        elseif any(ismember(cn, tbl.Properties.VariableNames))
                clust = args.clusters(ismember(cn, tbl.Properties.VariableNames));
        else
            error('wooldridge_TB:clusterVarNotFound', ...
                'clusters "%s" not found in table.', cn);
        end
end

% ----- Time → stable integer indices 1..T -----
[timeVals, ~, tIdx] = unique(tRaw, 'stable');               % preserve observed order
T = numel(timeVals);
baseTime = min(max(1, args.DropTimeBase), T);

% ----- Cohort per unit: first-treatment index; 0 = never-treated -----
[gid, ~] = findgroups(id);                                  % unit ids (groups)
gUnit = splitapply(@(x)x(1), gRaw, gid);                    % one g per unit
[tf, loc] = ismember(gUnit, timeVals);
cohortIdxUnit = zeros(numel(gUnit),1,'uint32');             % 0 => NT
cohortIdxUnit(tf) = uint32(loc(tf));
cohortIdx = cohortIdxUnit(gid);

treatedCoh = unique(cohortIdxUnit(cohortIdxUnit>0),'stable');
hasNT      = any(cohortIdxUnit==0);

% ----- Build design matrix -----
names = string([]); cols = {};

% Intercept
cols{end+1} = ones(n,1); names(end+1) = "const";

% Cohort dummies (treated only). If NT exists, NT is baseline; else drop last treated cohort.
if hasNT
    dropCoh = uint32(0);
elseif isempty(treatedCoh)
    dropCoh = uint32(0);
else
    dropCoh = treatedCoh(end);
end
for r = reshape(treatedCoh,1,[])
    if r==dropCoh, continue; end
    v = double(cohortIdx==r);
    key = "Coh_"+string(timeVals(double(r)));
    cols{end+1} = v; names(end+1) = key;
end

% Time FE (drop baseline time)
for j = 1:T
    if j==baseTime, continue; end
    v = double(tIdx==j);
    key = "F_"+string(timeVals(j));
    cols{end+1} = v; names(end+1) = key;
end

% Optional level covariates
for k = 1:numel(args.Covariates)
    vn = args.Covariates(k);
    if ~ismember(vn, tbl.Properties.VariableNames)
        warning('wooldridge_TB:covarMissing','Covariates "%s" not found; skipping.', vn);
        continue
    end
    cols{end+1} = tbl.(vn);
    names(end+1) = vn;
end

% Cohort×time ATTs: 1{cohort=r}*1{t=j} for j>=r
attKeys = strings(0,1);
attRT   = zeros(0,2,'uint32');    % rows: [r j]
for r = reshape(treatedCoh,1,[])
    for j = r:T
        v = double((cohortIdx==r) & (tIdx==j));
        key = "ATT_"+string(timeVals(double(r)))+"_x_"+string(timeVals(j));
        cols{end+1} = v; names(end+1) = key;
        attKeys(end+1,1) = key;
        attRT(end+1,:)   = [r j];
    end
end

% Stack X
p = numel(cols);
X = zeros(n,p);
for j = 1:p, X(:,j) = cols{j}; end

% ----- OLS -----
b    = X \ y;
u    = y - X*b;
XTXi = (X.'*X) \ eye(p);     % numerically stable, no explicit inv

% ----- SEs (cluster uses clusters; default = id) -----
    % Cluster-robust (CR0) using chosen cluster variable
    [cgid, ~] = findgroups(clust);
    G = max(cgid);
    
    S = zeros(p,p);
    for g = 1:G
        sel = find(cgid==g);
        Xg = X(sel,:); ug = u(sel);
        S = S + (Xg' * (ug*ug') * Xg);
    end
    V = XTXi * S * XTXi;
    df = max(1, G - 1);                          % clusters-1

se   = real(sqrt(diag(V)));
se(se<1e-6)=NaN;     % avoid divide-by-zero in t
t    = b ./ se;
pval = 2*tcdf(-abs(t), df);

coefTbl = table(names(:), b, se, t, pval, ...
    'VariableNames', {'Name','Estimate','SE','tStat','pValue'});

% ----- Extract ATTs -----
% Build a fast map from name -> position
name2pos = containers.Map(cellstr(names), num2cell(1:numel(names)));
% Positions of ATT coefficients in the same order as attKeys/attRT were constructed
posAtt = zeros(numel(attKeys),1);
for k = 1:numel(attKeys), posAtt(k) = name2pos(char(attKeys(k))); end

isATT = ismember(coefTbl.Name, attKeys);
ATT = coefTbl(isATT, :);
CohLab = string(timeVals(double(attRT(:,1))));
TimLab = string(timeVals(attRT(:,2)));
ATT.Cohort = CohLab; ATT.Time = TimLab;
ATT = movevars(ATT, {'Cohort','Time'}, 'Before', 1);
[~,ord] = sortrows(attRT, [1 2]);
ATT = ATT(ord,:);
posAtt = posAtt(ord,:);           % keep posAtt aligned with ATT rows & attRT

% ----- Details (always computed & stored; optionally displayed) -----
% 1) Units per cohort (including NT)
[Gcoh, cohvals] = findgroups(cohortIdxUnit);                % groups over units (cohort labels incl. NT=0)
nUnitsPerCoh = splitapply(@numel, cohortIdxUnit, Gcoh);     % count units
cohLabel = strings(numel(cohvals),1);
for ii = 1:numel(cohvals)
    r = cohvals(ii);
    if r==0
        cohLabel(ii) = "NT";
    else
        cohLabel(ii) = string(timeVals(double(r)));
    end
end
unitsPerCohort = table(cohLabel, nUnitsPerCoh, ...
    'VariableNames', {'Cohort','Units'});

% 2) Support by event time k using raw data (treated obs only)
validObs = (cohortIdx>0);                         % include all treated obs
k_all = double(tIdx(validObs)) - double(cohortIdx(validObs));
[Gk, kuniq] = findgroups(k_all);
nObs_k = splitapply(@numel, k_all, Gk);
nCoh_k = splitapply(@(g) numel(unique(g)), double(cohortIdx(validObs)), Gk);
supportK = table(kuniq, nObs_k, nCoh_k, ...
    'VariableNames', {'k','nObs','nCohorts'});

% 3) ATT-by-event time (averaging ATT coefficients across cohorts)
k_hat = double(attRT(:,2)) - double(attRT(:,1));            % event time per ATT coefficient (aligned with ATT rows)
[Gkh, khuniq] = findgroups(k_hat);
att_mean = splitapply(@mean, ATT.Estimate, Gkh);
att_med  = splitapply(@median, ATT.Estimate, Gkh);
nCells   = splitapply(@numel, ATT.Estimate, Gkh);
attByK = table(khuniq, att_mean, att_med, nCells, ...
    'VariableNames', {'k','ATT_hat_mean','ATT_hat_median','nCells'});

% 4) Cohort-level averages WITH SEs (delta method) + Overall (cohort-share weighted)
% Build cohort-wise equal-weight means over their ATT cells
treatedMask = (cohvals>0);
treatedCohAll = cohvals(treatedMask);                          % numeric cohort indices (time indices), treated only
treatedUnits = nUnitsPerCoh(treatedMask);                      % counts per treated cohort
sumUnits = sum(treatedUnits);
if sumUnits<=0, sumUnits = 1; end
wCoh = treatedUnits / sumUnits;                                % cohort-share weights for overall

K_att = numel(posAtt);
V_att = V(posAtt, posAtt);                                     % vcov over ATT cells only
b_att = b(posAtt);

% For each cohort r: indices of its ATT cells in ATT/V_att
cohortEffect_rows = [];   % to collect rows
for rr = 1:numel(treatedCohAll)
    r = treatedCohAll(rr);
    I_r = find(attRT(ord,1) == r);                             % positions in ATT (and in posAtt/V_att)
    m_r = numel(I_r);
    if m_r==0, continue; end
    v_r = (1/m_r) * ones(m_r,1);                               % equal weights across that cohort's post cells
    est_r = v_r' * b_att(I_r);
    var_r = v_r' * V_att(I_r, I_r) * v_r;
    se_r  = sqrt(max(var_r,0));
    t_r   = est_r / se_r;
    p_r   = 2*tcdf(-abs(t_r), df);
    cohortEffect_rows = [cohortEffect_rows; [double(r), m_r, est_r, se_r, t_r, p_r]]; %#ok<AGROW>
end

cohortEffect = array2table(cohortEffect_rows, ...
    'VariableNames', ["Cohort","#TimePeriods","ATT(k)","SE","tStat","pValue"]);
cohortEffect = sortrows(cohortEffect,"Cohort","ascend");

% Overall ATT: weighted average across cohorts (weights = cohort unit shares)
% Construct a single linear combo over ATT cells: for cohort r, every cell gets weight wCoh(r) * 1/m_r
a_over = zeros(K_att,1);
for rr = 1:numel(treatedCohAll)
    r = treatedCohAll(rr);
    I_r = find(attRT(ord,1) == r);
    m_r = numel(I_r);
    if m_r==0, continue; end
    a_over(I_r) = a_over(I_r) + (wCoh(rr) / m_r);
end
overall_est = a_over' * b_att;
overall_var = a_over' * V_att * a_over;
overall_se  = sqrt(max(overall_var,0));
overall_t   = overall_est / overall_se;
overall_p   = 2*tcdf(-abs(overall_t), df);
overallTbl  = table(overall_est, overall_se, overall_t, overall_p, ...
    'VariableNames', {'Estimate','SE','tStat','pValue'});

% 5) Metadata & dropped dummies
if dropCoh==0
    droppedCohLabel = "NT"; droppedCohName = "Coh_NT";
else
    droppedCohLabel = string(timeVals(double(dropCoh)));
    droppedCohName  = "Coh_"+droppedCohLabel;
end
cohortDummyNames = "Coh_"+string(timeVals(double(treatedCoh)));
cohortDummyNames = setdiff(cohortDummyNames, droppedCohName, 'stable');

details = struct( ...
    'unitsPerCohort', unitsPerCohort, ...
    'supportByEventTime', supportK, ...
    'attByEventTime', attByK, ...
    'ATTbyCohort', cohortEffect, ...             % now includes SE/t/p
    'Overall', overallTbl, ...                   % overall ATT with SE/t/p
    'TimeValues', string(timeVals), ...
    'DropTimeBase', string(timeVals(baseTime)), ...
    'HasNeverTreated', hasNT, ...
    'DroppedCohort', droppedCohLabel, ...
    'CohortDummies', cohortDummyNames, ...
    'TimeDummies', "F_"+string(timeVals(setdiff(1:T, baseTime))), ...
    'ATT_Names', attKeys);



% ----- Optional display -----
if args.Print
    fprintf('[wooldridge] ATT by cohort (mean of cohorts coefficients):\n'); disp(cohortEffect);
    fprintf('[wooldridge] Overall ATT (cohort-share weighted):\n'); disp(overallTbl);
    fprintf('[wooldridge] ATT by event time (mean/median across cohorts):\n'); disp(attByK);
end

if args.details
    % Standard Output

    fprintf('[wooldridge_TB] ATT cell estimates:\n');
    disp(ATT(:, {'Cohort','Time','Estimate','SE','tStat','pValue'}));
end

% ----- Pack output -----
out = struct();
out.coef    = coefTbl;
out.att     = ATT(:, {'Cohort','Time','Estimate','SE','tStat','pValue'});
out.details = details;
out.summaryTable = cohortEffect;     % keep cohort summary as main table
out.overall = overallTbl;            % expose overall at top level as well
end
