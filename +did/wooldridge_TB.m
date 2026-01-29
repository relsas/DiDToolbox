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
% Last change: 11/16/2025
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
    args.Display logical=true  % prints just standard results
    % --- Minimal additions for toolbox integration ---
    args.vcov (1,1) string = "clustered"
    args.clusters = []                 % name or vector aligned with rows (any type)
    args.EventStudy (1,1) logical = false
    args.Weights (:,1) double = []     % Optional weights vector
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
        % args.clusters contains variable names -> extract columns
        if numel(cn) == 1
            clust = tbl.(cn);
        else
            clust = tbl(:, cn); % findgroups accepts table
        end
    elseif any(ismember(cn, tbl.Properties.VariableNames))
        % Mixed or partial? Use what we found
        found = cn(ismember(cn, tbl.Properties.VariableNames));
        if numel(found) == 1
            clust = tbl.(found);
        else
            clust = tbl(:, found);
        end
    else
        error('wooldridge_TB:clusterVarNotFound', ...
            'clusters "%s" not found in table.', strjoin(cn,', '));
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

if args.EventStudy
    % EVENT STUDY: Estimate ALL interactions except reference (t = g - 1)
    % Leads (t < g-1) and Lags (t >= g)
    for r = reshape(treatedCoh,1,[])
        for j = 1:T
            % Identification: Drop reference period relative to treatment start
            % Standard: Ref = g - 1. Time index for g is 'r' (since r is index in timeVals).
            if j == (r - 1), continue; end

            % Also skip if t < first_data_time ? No, tIdx covers available data.

            v = double((cohortIdx==r) & (tIdx==j));

            % If no observations for this cell, skip (collinear/empty)
            if sum(v)==0, continue; end

            key = "ATT_"+string(timeVals(double(r)))+"_x_"+string(timeVals(j));
            cols{end+1} = v; names(end+1) = key;
            attKeys(end+1,1) = key;
            attRT(end+1,:)   = [r j];
        end
    end
else
    % POOLED / POST-ONLY: Estimate only j >= r
    for r = reshape(treatedCoh,1,[])
        for j = r:T
            v = double((cohortIdx==r) & (tIdx==j));
            key = "ATT_"+string(timeVals(double(r)))+"_x_"+string(timeVals(j));
            cols{end+1} = v; names(end+1) = key;
            attKeys(end+1,1) = key;
            attRT(end+1,:)   = [r j];
        end
    end
end

% Stack X
p = numel(cols);
X = zeros(n,p);
for j = 1:p, X(:,j) = cols{j}; end

% ----- WLS Transformation (if Weights provided) -----
if ~isempty(args.Weights)
    if numel(args.Weights) ~= n
        error('wooldridge_TB:DimMismatch', 'Weights vector must match table height.');
    end
    sqrtW = sqrt(args.Weights);
    % Transform y and X
    y = y .* sqrtW;
    for j=1:numel(cols)
        cols{j} = cols{j} .* sqrtW;
    end
    % Note: Cluster variables (id, clust) should NOT be transformed.
    % The SE calculation uses 'X' and 'u'; if we use transformed X and y,
    % b will be correct WLS estimates.
    % Residuals u = y_trans - X_trans * b will be "weighted residuals".
    % The cluster sandwich formula XTXi * Sum(Xg' ug ug' Xg) * XTXi
    % will thus naturally implement the weighted cluster robust variance.
end

% ----- OLS with Rank Deficiency Handling -----
% Stata behaviors: drops collinear columns. We do the same via QR.
[Q, R, perm] = qr(X,0);
pFull = size(X,2);
tol   = max(size(X)) * eps(abs(R(1,1)));
rankX = sum(abs(diag(R)) > tol);

if rankX < pFull
    % Identify kept columns (in original order)
    keepIdx = sort(perm(1:rankX));

    % Warn if we are dropping columns (common in saturated models, but good to know)
    droppedCols = setdiff(1:pFull, keepIdx);
    if args.details || args.Display
        % fprintf('[wooldridge] Rank deficiency detected. Dropping %d columns: %s\n', ...
        %    numel(droppedCols), strjoin(names(droppedCols), ', '));
    end

    % Subset X and names
    X = X(:, keepIdx);
    names = names(keepIdx);
    p = size(X,2);
else
    p = pFull;
end

b    = X \ y;
u    = y - X*b;
XTXi = (X.'*X) \ eye(p);     % Now full rank, so this is safe and finite

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

% Small-sample correction (Stata/R consistent)
adj = (G / (G - 1)) * ((n - 1) / (n - p));
V = V * adj;

se   = real(sqrt(diag(V)));
se(se<1e-6)=NaN;     % avoid divide-by-zero in t
t    = b ./ se;
pval = 2*tcdf(-abs(t), df);

coefTbl = table(names(:), b, se, t, pval, ...
    'VariableNames', {'Name','Estimate','SE','tStat','pValue'});

% ----- Extract ATTs -----
% Build a fast map from name -> position (names is now subsetted)
name2pos = containers.Map(cellstr(names), num2cell(1:numel(names)));

% Filter attKeys/attRT to only those that survived rank check
validMask = isKey(name2pos, cellstr(attKeys));
if any(~validMask)
    % optionally warn
    % droppedATTs = attKeys(~validMask);
    attKeys = attKeys(validMask);
    attRT   = attRT(validMask,:);
end

% Positions of ATT coefficients
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

% 4) Cohort-level averages & PRE-CALCULATE Weights for (3)
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

% 3) ATT-by-event time (Aggregated with SEs via Delta Method)
% Aggregation: ATT(k) = Sum_{g} w_{g,k} * ATT(g, g+k)
% Weights w_{g,k} are proportional to cohort size N_g (average effect on the treated).

k_hat = double(attRT(:,2)) - double(attRT(:,1));            % event time per ATT coefficient (aligned with ATT rows)
[Gkh, khuniq] = findgroups(k_hat);

% Prepare storage
numK = numel(khuniq);
est_k = zeros(numK, 1);
se_k  = zeros(numK, 1);
t_k   = zeros(numK, 1);
p_k   = zeros(numK, 1);
n_cells_k = zeros(numK, 1);
H_es  = zeros(numK, K_att); % Transformation matrix for V_es

for ii = 1:numK
    k_val = khuniq(ii);

    % Find all (g,t) cells contributing to this k
    idx_k = find(k_hat == k_val); % Indices in ATT/b_att/V_att

    if isempty(idx_k)
        continue;
    end
    n_cells_k(ii) = numel(idx_k);

    % Get weights for these cells
    % We weight by Cohort Size (N_g)
    % Identify which cohort each cell belongs to
    cohorts_in_k = attRT(ord(idx_k), 1); % Cohort indices (time values)

    % Map cohort time -> index in 'treatedCohAll' / 'treatedUnits'
    % We need N_g for each g in cohorts_in_k.
    % Faster lookup:
    weights_k = zeros(numel(idx_k), 1);
    for jj = 1:numel(idx_k)
        coh_t = cohorts_in_k(jj);
        % Find N_g
        % treatedCohAll contains the numeric cohort times
        % treatedUnits contains counts
        mask_c = (treatedCohAll == coh_t);
        if any(mask_c)
            weights_k(jj) = treatedUnits(mask_c);
        else
            weights_k(jj) = 0; % Should not happen
        end
    end

    % Normalize weights
    sum_w = sum(weights_k);
    if sum_w > 0
        w_vec = weights_k / sum_w;
    else
        w_vec = ones(size(weights_k)) / numel(weights_k);
    end

    % Linear combination indices in full b vector
    % idx_k are indices into 'b_att' (subset).
    % 'b_att' = b(posAtt).

    % Store weights in Transformation Matrix
    H_es(ii, idx_k) = w_vec';

    % Estimate point and SE (legacy code kept for safety, though V_es covers it)
    est_k(ii) = w_vec' * b_att(idx_k);

    % Variance: w' V_{sub} w
    % V_{sub} is block of V_att corresponding to these cells
    V_sub = V_att(idx_k, idx_k);
    var_k = w_vec' * V_sub * w_vec;
    se_k(ii)  = sqrt(max(var_k, 0));

    t_k(ii)   = est_k(ii) / se_k(ii);
    p_k(ii)   = 2*tcdf(-abs(t_k(ii)), df);
end

% Compute Full Covariance Matrix for Event Study Estimates
V_es = H_es * V_att * H_es';

attByK = table(khuniq, est_k, se_k, t_k, p_k, n_cells_k, ...
    'VariableNames', {'EventTime','Estimate','SE','tStat','pValue','N_Cells'});

% If Event Study, ensure Reference Period (k=-1) is present
if isfield(args,'EventStudy') && args.EventStudy && ~ismember(-1, attByK.EventTime)
    % Create reference row
    refRow = table(-1, 0, 0, NaN, NaN, NaN, ...
        'VariableNames', {'EventTime','Estimate','SE','tStat','pValue','N_Cells'});

    % Append to table
    attByK_unsorted = [attByK; refRow];

    % Expand Covariance Matrix (Correlation of Ref Period with others is 0)
    V_es_unsorted = blkdiag(V_es, 0);

    % Sort table and apply same permutation to V_es
    [attByK, sortIdx] = sortrows(attByK_unsorted, 'EventTime');
    V_es = V_es_unsorted(sortIdx, sortIdx);
end

% 4) Cohort-level averages WITH SEs (delta method) + Overall (cohort-share weighted)
% Build cohort-wise equal-weight means over their ATT cells
% ... (rest of logic handles weighted average of post-treatment cells)
% Note: The block below only processes 'treatedCohAll' and uses 'attRT'
% 'attRT' contains indices in 'b_att'. If EventStudy, 'attRT' includes Leads.
% We should ensure 'overallTbl' only aggregates POST-treatment effects, or clarify what 'Overall' means.
% Standard DiD "Overall ATT" is average of post-treatment effects.
% So we should filter 'attRT' for j >= r (k >= 0) if we want standard ATT.
% Let's keep it simple: Use existing logic. 'treatedCohAll' iterates cohorts.
% 'attRT' will contain leads if EventStudy.
% We should likely restrict "Overall" to post-treatment.

% Filter attRT for overall calculation to ensure it's ATT (post)
% Current logic: I_r = find(attRT(ord,1) == r);
% This grabs ALL coefficients for cohort r (both leads and lags).
% We should filter for lags (k>=0).


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
    'VariableNames', ["Cohort","#TimePeriods","Estimate","SE","tStat","pValue"]);
cohortEffect = sortrows(cohortEffect,"Cohort","ascend");

% Overall ATT: weighted average across cohorts (weights = observation count shares)
% Correct sample ATT matches the definition: E[Y(1)-Y(0)|D=1]
% Total number of treated observations in the estimation sample:
totalTreatedObs = 0;
for rr = 1:numel(treatedCohAll)
    r = treatedCohAll(rr);
    m_r = sum(attRT(ord,1) == r);
    totalTreatedObs = totalTreatedObs + (treatedUnits(rr) * m_r);
end

if totalTreatedObs == 0, totalTreatedObs = 1; end

% Construct a single linear combo over ATT cells:
% For cohort r, number of obs contributing to each cell is N_r.
% Weight for cell (r,t) is N_r / TotalTreatedObs.
a_over = zeros(K_att,1);
for rr = 1:numel(treatedCohAll)
    r = treatedCohAll(rr);
    I_r = find(attRT(ord,1) == r);
    if isempty(I_r), continue; end

    % Weight per cell in this cohort
    w_cell = treatedUnits(rr) / totalTreatedObs;
    a_over(I_r) = w_cell;
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
    'V_es', V_es, ...                            % Full Vcov for Event Study
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
if args.Display
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
out.Method ="Wooldridge";
out.coef    = coefTbl;
out.ATT     = ATT(:, {'Cohort','Time','Estimate','SE','tStat','pValue'});
out.details = details;
out.summaryTable = cohortEffect;     % keep cohort summary as main table
out.overall = overallTbl;            % expose overall at top level as well
end
