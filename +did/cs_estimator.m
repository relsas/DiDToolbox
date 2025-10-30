function out = cs_estimator(T, opts)
% DID.CS_ESTIMATOR  Callaway & Sant'Anna (2021) ATT(g,t) with inference
%
% Accepts both:
%   out = did.cs_estimator(T, idVar="id", timeVar="time", ...)
%   out = did.cs_estimator(T, struct('idVar',"id", 'timeVar',"time", ...))
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/13/2025  (patched: Weighting="cohortShare"|"treatedObs")
% ------------------------------------------------------------------------

arguments
    T table

    % Required names in T
    opts.idVar       (1,1) string = "id"
    opts.timeVar     (1,1) string = "time"
    opts.yVar        (1,1) string = "y"
    opts.dVar        (1,1) string = "D"

    % Optional columns
    opts.WeightVar          string = string.empty
    opts.Covariates         string = string.empty

    % Core controls
    opts.Approach    (1,1) string {mustBeMember(opts.Approach,["unconditional","or","ipw","dr"])} = "unconditional"
    opts.Comparison  (1,1) string {mustBeMember(opts.Comparison,["never","notyet"])} = "never"
    opts.Delta       (1,1) double {mustBeInteger, mustBeNonnegative} = 0

    % NEW: Aggregation weighting toggle
    % "cohortShare" (current default behavior) vs "treatedObs" (post-only treated counts per (g,t))
    opts.Weighting   (1,1) string {mustBeMember(opts.Weighting,["cohortShare","treatedObs"])} = "treatedObs"

    % Inference controls
    opts.SEMethod    (1,1) string {mustBeMember(opts.SEMethod,["multiplier","clustered","clustered2"])} = "multiplier"
    opts.B           (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.Seed        (1,1) double = randi([1,1000],1,1)
    opts.Multiplier  (1,1) string {mustBeMember(opts.Multiplier,["rademacher","mammen"])} = "rademacher"
    opts.Studentize  (1,1) logical = true
    opts.ClusterVar          string = string.empty
    opts.ClusterVar2         string = string.empty
    opts.SmallSample (1,1) logical = true
    opts.UseParallel (1,1) logical = false
    opts.Details     (1,1) logical = false
    opts.Print       (1,1) logical = true

    % Cross-fitting controls
    opts.CrossFit    (1,1) logical = false
    opts.Kfolds      (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.Kfolds,2)} = 5
    opts.StratifyFoldsBy (1,1) string {mustBeMember(opts.StratifyFoldsBy,["none","cohort"])} = "none"
end

% ----- Guards / normalization -----
Approach   = lower(string(opts.Approach));
Comparison = lower(string(opts.Comparison));
SEMethod   = lower(string(opts.SEMethod));

if Approach=="dr" && isempty(opts.Covariates)
    error('did:cs:drNeedsX','DR requires Covariates.');
end

% ---------- Extract variables ----------
id = T.(opts.idVar);
tt = T.(opts.timeVar);
y  = T.(opts.yVar);
D  = T.(opts.dVar);
if islogical(D), D = double(D); end
if ~isfloat(D),  D = double(D); end

% Optional weights
if ~isempty(opts.WeightVar) && any(strcmp(opts.WeightVar, T.Properties.VariableNames))
    w = double(T.(opts.WeightVar));
else
    w = ones(height(T),1);
end

% ---------- Indexing ----------
[~,~,idIdx] = unique(id,'stable');        % unit index 1..N
N = max(idIdx);

% Map time -> 1..Tnum robustly
Tvals = unique(tt,'stable');
if iscellstr(Tvals) || isstring(Tvals) || iscategorical(Tvals)
    keys  = cellstr(string(Tvals));
    t2num = containers.Map(keys, num2cell(1:numel(keys)));
    tnum  = arrayfun(@(x) t2num(char(string(x))), tt);
elseif isnumeric(Tvals)
    keys  = num2cell(Tvals);
    t2num = containers.Map(keys, num2cell(1:numel(keys)));
    tnum  = arrayfun(@(x) t2num(x), tt);
else
    error('did:cs:timeType','Time must be numeric, string, or categorical.');
end
Tnums = unique(tnum,'stable')';

% ---------- Group assignment: first treat time G; never-treated C ----------
G = inf(N,1);
for i=1:N
    ii   = (idIdx==i);
    Dt_i = D(ii);
    tt_i = tnum(ii);
    first = find(Dt_i>0,1,'first');
    if ~isempty(first), G(i) = tt_i(first); end
end
Cunit = isinf(G);
Gs = unique(G(~isinf(G)))';
if isempty(Gs), error('did:cs:noTreated','No treated units found.'); end

% ---------- Precompute period means per unit ----------
Yt  = NaN(N, numel(Tnums));
Wt  = NaN(N, numel(Tnums));
for ti = 1:numel(Tnums)
    t   = Tnums(ti);
    rows = (tnum==t);
    Yt(:,ti) = accumarray(idIdx(rows), y(rows), [N,1], @mean, NaN);
    Wt(:,ti) = accumarray(idIdx(rows), w(rows), [N,1], @mean, NaN);
end

% ---------- Covariates (only if needed) ----------
needsX = ismember(Approach, ["or","ipw","dr"]);
if ~needsX && ~isempty(opts.Covariates)
    warning('did:cs:CovariatesIgnored', ...
        'Covariates provided but Approach="%s" ignores them; continuing without X.', Approach);
end
hasX = ~isempty(opts.Covariates) && needsX;

Xt = [];   % each cell ti holds N x K means at time t
if hasX
    K = numel(opts.Covariates);
    % row-level matrix
    Xfull = zeros(height(T), K);
    for k=1:K
        vn = opts.Covariates(k);
        if ~any(strcmp(vn, T.Properties.VariableNames))
            error('did:cs:covarMissing','Covariate "%s" missing in table.', vn);
        end
        xk = T.(vn);
        if islogical(xk), xk = double(xk); end
        Xfull(:,k) = double(xk(:));
    end
    % collapse per column
    Xt = cell(1, numel(Tnums));
    for ti = 1:numel(Tnums)
        rows = (tnum==Tnums(ti));
        ii   = idIdx(rows);
        Xi   = Xfull(rows,:);
        Xmean = NaN(N, K);
        for k=1:K
            Xmean(:,k) = accumarray(ii, Xi(:,k), [N,1], @mean, NaN);
        end
        Xt{ti} = Xmean;   % N x K
    end
end

% ---------- (Build unit-level folds for cross-fitting ----------
if opts.CrossFit
    if opts.StratifyFoldsBy=="cohort"
        unitFold = unit_kfold_stratified_cohort_(G, opts.Kfolds, opts.Seed);  % N x 1
    else
        unitFold = unit_kfold_(N, opts.Kfolds, opts.Seed);                     % N x 1
    end
else
    unitFold = [];   % signal: no cross-fit
end

% ---------- Main loop: ATT(g,t) + per-unit ψ_i(g,t) ----------
names   = strings(0,1);
g_list  = []; t_list = [];
ATT     = [];
psiCols = {};                   % each is N x 1 (ψ_i for this cell)
cell_wt = [];                   % NEW: per-(g,t) treated-post weight (Wtreat)

for g = Gs
    refT  = g - opts.Delta - 1;
    refIx = find(Tnums == refT, 1);
    if isempty(refIx), continue; end

    for ti = 1:numel(Tnums)
        t = Tnums(ti);
        if t < g, continue; end   % post

        % Long difference ΔY = Y_t - Y_ref
        dY   = Yt(:,ti) - Yt(:,refIx);
        wt_u = Wt(:,ti);                    % period-t weights
        valid = ~isnan(dY);

        % Groups
        treat_u = (G==g) & valid;
        switch Comparison
            case "never"
                comp_u = (Cunit==1) & valid;
            case "notyet"
                comp_u = (G > (t + opts.Delta)) & valid;   % includes never-treated since inf > ...
        end
        eligible = (treat_u | comp_u);
        if ~any(treat_u) || ~any(comp_u), continue; end

        W_all  = nansum(wt_u(eligible));
        Wtreat = nansum(wt_u(treat_u));      % ← NEW: store this for treatedObs weighting
        Wcomp  = nansum(wt_u(comp_u));
        if W_all<=0 || Wtreat<=0 || Wcomp<=0, continue; end

        p_treat = Wtreat / W_all;
        p_comp  = Wcomp  / W_all;

        % Treated mean
        mu_treat = nansum(wt_u(treat_u).*dY(treat_u))/Wtreat;

        % Branch by Approach
        psi = zeros(N,1);   % ψ_i contributions so that θ̂ ≈ Σ_i ψ_i (O(1) total)
        switch Approach
            case "unconditional"
                mu_comp = nansum(wt_u(comp_u).*dY(comp_u))/Wcomp;
                att = mu_treat - mu_comp;

                psi(treat_u) = (wt_u(treat_u)/W_all) .* ((dY(treat_u)-mu_treat)/p_treat);
                psi(comp_u)  = psi(comp_u) - (wt_u(comp_u)/W_all) .* ((dY(comp_u)-mu_comp)/p_comp);

            case "or"
                if ~hasX, error('did:cs:orNeedsX','OR requires Covariates.'); end
                if ~opts.CrossFit
                    % ---- original OR (no CF) ----
                    Xc = Xt{ti}(comp_u & eligible, :);
                    yc = dY(comp_u & eligible);
                    wc = wt_u(comp_u & eligible);
                    Xc1 = [ones(size(Xc,1),1), Xc];
                    beta = wls_(Xc1, yc, wc);

                    Xt_t = Xt{ti}(treat_u, :);
                    mt    = [ones(sum(treat_u),1), Xt_t] * beta;
                    w_t   = wt_u(treat_u);
                    mu_reg_t = nansum(w_t .* mt) / nansum(w_t);

                    att = mu_treat - mu_reg_t;

                    psi_t = (w_t / W_all) .* ( ((dY(treat_u) - mu_treat) - (mt - mu_reg_t)) / p_treat );
                    psi(treat_u) = psi_t;

                else
                    % ---- cross-fitted OR (CF-OR) ----
                    [mt_all, mu_reg_t] = cf_or_predict_(Xt{ti}, dY, wt_u, treat_u, comp_u, eligible, idIdx, unitFold);
                    att = mu_treat - mu_reg_t;

                    % ψ remains treated-only with CF predictions (orthogonal)
                    w_t  = wt_u(treat_u);
                    psi_t = (w_t / W_all) .* ( ((dY(treat_u) - mu_treat) - (mt_all(treat_u) - mu_reg_t)) / p_treat );
                    psi(treat_u) = psi_t;
                end

            case "ipw"
                if ~hasX, error('did:cs:ipwNeedsX','IPW requires Covariates.'); end
                ge = (G(eligible)==g);
                Xe = Xt{ti}(eligible, :);
                we = wt_u(eligible);

                if ~opts.CrossFit
                    pi_e = fit_logit_prob_(Xe, ge, we);
                else
                    rowFold = unitFold(idIdx(eligible));
                    pi_e = fit_logit_prob_cf_(Xe, ge, we, rowFold, max(unitFold));
                end

                pi_bar     = nansum(we .* pi_e)/nansum(we);
                one_pi_bar = 1 - pi_bar;
                omega_e = zeros(sum(eligible),1);
                omega_e(~ge) = (pi_e(~ge) ./ (1 - pi_e(~ge))) * (one_pi_bar / pi_bar);

                omega = zeros(N,1);
                idx_el = find(eligible);
                omega(idx_el) = omega_e;

                w_eff = wt_u .* omega;
                mu_c_ipw = nansum(w_eff(comp_u).*dY(comp_u)) / nansum(w_eff(comp_u));
                att = mu_treat - mu_c_ipw;

                psi(treat_u) = (wt_u(treat_u)/W_all) .* ((dY(treat_u)-mu_treat)/p_treat);
                denom = nansum(w_eff(comp_u));
                if denom<=0, continue; end
                psi_c = zeros(N,1);
                psi_c(comp_u) = (w_eff(comp_u)/denom) .* (dY(comp_u) - mu_c_ipw);
                psi = psi - (denom/W_all) * (psi_c / p_comp);

            case "dr"
                if ~hasX, error('did:cs:drNeedsX','DR requires Covariates.'); end

                if ~opts.CrossFit
                    % OR on controls
                    Xc = Xt{ti}(comp_u & eligible, :);
                    yc = dY(comp_u & eligible);
                    wc = wt_u(comp_u & eligible);
                    Xc1 = [ones(size(Xc,1),1), Xc];
                    beta = wls_(Xc1, yc, wc);

                    % propensity on eligible
                    ge = (G(eligible)==g);
                    Xe = Xt{ti}(eligible, :);
                    we = wt_u(eligible);
                    pi_e = fit_logit_prob_(Xe, ge, we);

                    pi_bar     = nansum(we .* pi_e)/nansum(we);
                    one_pi_bar = 1 - pi_bar;
                    omega_e = zeros(sum(eligible),1);
                    omega_e(~ge) = (pi_e(~ge) ./ (1 - pi_e(~ge))) * (one_pi_bar / pi_bar);

                    omega = zeros(N,1);
                    idx_el = find(eligible);
                    omega(idx_el) = omega_e;

                    Xt_t = Xt{ti}(treat_u & eligible, :);
                    mt   = [ones(sum(treat_u & eligible),1), Xt_t] * beta;
                    w_t  = wt_u(treat_u & eligible);
                    mu_reg_t = nansum(w_t .* mt) / nansum(w_t);
                    t_part   = mu_treat - mu_reg_t;

                    Xc_full = Xt{ti}(comp_u & eligible, :);
                    m0_c    = [ones(sum(comp_u & eligible),1), Xc_full] * beta;
                    resid_c = dY(comp_u & eligible) - m0_c;
                    w_eff   = wt_u .* omega;
                    denom   = nansum(w_eff(comp_u & eligible));
                    if denom<=0, continue; end
                    c_part  = nansum(w_eff(comp_u & eligible) .* resid_c) / denom;

                    att = t_part - c_part;

                    psi = zeros(N,1);
                    mt_all = zeros(N,1);
                    mt_all(treat_u & eligible) = mt;
                    treated_res = dY - mt_all;
                    mu_tr_res  = (mu_treat - mu_reg_t);
                    psi(treat_u) = (wt_u(treat_u)/W_all) .* ((treated_res(treat_u) - mu_tr_res)/p_treat);
                    psi_c = zeros(N,1);
                    resid_all = zeros(N,1);
                    resid_all(comp_u & eligible) = resid_c;
                    psi_c(comp_u & eligible) = (w_eff(comp_u & eligible)/denom) .* (resid_all(comp_u & eligible) - c_part);
                    psi = psi - (denom/W_all) * (psi_c / p_comp);

                else
                    % ---- cross-fitted DR (CF-DR) ----
                    % cross-fitted outcome model on controls
                    [mt_all, mu_reg_t] = cf_or_predict_(Xt{ti}, dY, wt_u, treat_u, comp_u, eligible, idIdx, unitFold);
                    % cross-fitted propensity on eligible
                    ge = (G(eligible)==g);
                    Xe = Xt{ti}(eligible, :);
                    we = wt_u(eligible);
                    rowFold = unitFold(idIdx(eligible));
                    pi_e = fit_logit_prob_cf_(Xe, ge, we, rowFold, max(unitFold));

                    pi_bar     = nansum(we .* pi_e)/nansum(we);
                    one_pi_bar = 1 - pi_bar;
                    omega_e = zeros(sum(eligible),1);
                    omega_e(~ge) = (pi_e(~ge) ./ (1 - pi_e(~ge))) * (one_pi_bar / pi_bar);
                    omega = zeros(N,1);
                    idx_el = find(eligible);
                    omega(idx_el) = omega_e;

                    % treated part
                    w_t = wt_u(treat_u);
                    t_part = mu_treat - mu_reg_t;

                    % control residuals with CF m0(X)
                    m0_cf = cf_or_predict_controls_(Xt{ti}, dY, wt_u, comp_u, eligible, idIdx, unitFold); % N x 1
                    resid_c = zeros(N,1);
                    ce = (comp_u & eligible);
                    resid_c(ce) = dY(ce) - m0_cf(ce);

                    w_eff   = wt_u .* omega;
                    denom   = nansum(w_eff(ce));
                    if denom<=0, continue; end
                    c_part  = nansum(w_eff(ce) .* resid_c(ce)) / denom;

                    att = t_part - c_part;

                    % ψ (orthogonal form)
                    psi = zeros(N,1);
                    treated_res = dY - mt_all;     % only nonzero on treated slice
                    mu_tr_res  = (mu_treat - mu_reg_t);
                    psi(treat_u) = (wt_u(treat_u)/W_all) .* ((treated_res(treat_u) - mu_tr_res)/p_treat);
                    psi_c = zeros(N,1);
                    psi_c(ce) = (w_eff(ce)/denom) .* (resid_c(ce) - c_part);
                    psi = psi - (denom/W_all) * (psi_c / p_comp);
                end
        end

        names(end+1,1) = "ATT(g="+string(g)+", t="+string(t)+")"; %#ok<AGROW>
        g_list(end+1,1) = g;                                       %#ok<AGROW>
        t_list(end+1,1) = t;                                       %#ok<AGROW>
        ATT(end+1,1)    = att;                                     %#ok<AGROW>
        psiCols{end+1}  = psi;                                     %#ok<AGROW>
        cell_wt(end+1,1)= Wtreat;                                  %#ok<AGROW> % NEW
    end
end

R = numel(ATT);
if R==0, error('did:cs:noCells','No valid (g,t) cells.'); end

% ----- Influence matrices -----
Psi     = zeros(N,R);
for r = 1:R, Psi(:,r) = psiCols{r}; end
Psi_unc = Psi;                      % uncentered ψ for clustered SE
Psi_ctr = Psi - mean(Psi,1);        % centered ψ for multiplier bootstrap
Psi_unc(~isfinite(Psi_unc)) = 0;

% ---- Diagnostics for Vcov engines (CS Wild Cluster, etc.) ----
diagCS = struct( ...
    'Psi_unc', Psi_unc, ...
    'idIdx',   idIdx, ...
    'ATT',     ATT, ...
    'names',   names, ...
    'g_list',  g_list, ...
    't_list',  t_list );

% ---------- Inference ----------
SE   = NaN(R,1); tStat = NaN(R,1); pValue = NaN(R,1);
Bands = struct('alpha',0.95,'crit',NaN,'lower',NaN(R,1),'upper',NaN(R,1));
att_b = [];   % R x B draws (for aggregates if multiplier is used)
S_for_aggs = []; % struct for clustered aggregates

switch SEMethod
    case "multiplier"
        if opts.B > 0
            sd = double(opts.Seed);
            if isnan(sd) ||isempty(sd) || ~isfinite(sd) || sd < 0 || sd >= 2^32, sd = randi([1,1000],1,1); end
            rng(floor(sd));

            XI = draw_multipliers_(N, max(0,floor(double(opts.B))), lower(string(opts.Multiplier)));
            % θ* = θ̂ + Σ_i (ψ_i − \barψ) ξ_i
            att_b = ATT + (Psi_ctr.' * XI);          % R x B
            att_b_centered = att_b - mean(att_b,2);
            SE = std(att_b_centered, 0, 2);

            zeroSE = (SE==0 | ~isfinite(SE));
            tStat(~zeroSE) = ATT(~zeroSE) ./ SE(~zeroSE);

            if opts.Studentize
                Z = abs((att_b - ATT) ./ max(SE, eps));     % R x B
                thr = abs(tStat);
                pValue(~zeroSE) = mean(Z(~zeroSE,:) >= thr(~zeroSE), 2);
                pValue(zeroSE)  = NaN;
            else
                pValue(~zeroSE) = 2 * (1 - normcdf(abs(tStat(~zeroSE))));
                pValue(zeroSE)  = NaN;
            end

            % Simultaneous max-t bands (95%)
            Z = abs((att_b - ATT) ./ max(SE, eps));
            maxZ = max(Z, [], 1);
            crit = did.quantile(maxZ, 0.95);
            Bands.crit  = crit;
            Bands.lower = ATT - crit .* SE;
            Bands.upper = ATT + crit .* SE;
        else
            warning('did:cs:noBootstrap','SEMethod="multiplier" but B=0 — SE/t/p will be NaN.');
        end

    case "clustered"
        % one-way (unit or custom unit-constant var)
        [cIdx1, C1] = cluster_index_per_unit_(T, idIdx, N, opts.ClusterVar, isempty(opts.ClusterVar));
        S1 = cluster_sums_(Psi_unc, cIdx1);
        [SE, tStat, pValue, Bands, S1] = crv1_from_S_(ATT, S1);
        S_for_aggs = struct('kind',"clustered",'S1',S1,'C1',C1);

    case "clustered2"
        % two-way CGM: unit x second unit-constant cluster (e.g., industry)
        [cIdx1, C1] = cluster_index_per_unit_(T, idIdx, N, opts.ClusterVar, true);   % default: unit
        if isempty(opts.ClusterVar2)
            error('did:cs:clustered2','ClusterVar2 must be provided for SEMethod="clustered2".');
        end
        [cIdx2, C2] = cluster_index_per_unit_(T, idIdx, N, opts.ClusterVar2, false);
        % Intersection clusters (pair labels)
        c12 = double(categorical( string(cIdx1) + "|" + string(cIdx2) ));
        [~,~,cIdx12] = unique(c12,'stable'); C12 = max(cIdx12);

        % Cluster sums for each margin and joint
        S1  = cluster_sums_(Psi_unc, cIdx1);    % R x C1
        S2  = cluster_sums_(Psi_unc, cIdx2);    % R x C2
        S12 = cluster_sums_(Psi_unc, cIdx12);   % R x C12

        % CRV1 margin/joint (center across clusters)
        S1c  = S1  - mean(S1, 2);
        S2c  = S2  - mean(S2, 2);
        S12c = S12 - mean(S12,2);

        V1   = (C1 /(C1 -1)) * sum(S1c .^2, 2);
        V2   = (C2 /(C2 -1)) * sum(S2c .^2, 2);
        V12  = (C12/(C12-1)) * sum(S12c.^2, 2);

        Var_hat = V1 + V2 - V12;                % CGM
        SE = sqrt(Var_hat);

        zeroSE = (SE==0 | ~isfinite(SE));
        tStat = NaN(numel(ATT),1); pValue = NaN(numel(ATT),1);
        tStat(~zeroSE)  = ATT(~zeroSE) ./ SE(~zeroSE);
        pValue(~zeroSE) = 2*(1-normcdf(abs(tStat(~zeroSE))));
        crit = 1.96;
        Bands = struct('alpha',0.95,'crit',crit,'lower',ATT-crit.*SE,'upper',ATT+crit.*SE);

        S_for_aggs = struct('kind',"clustered2",'S1',S1,'C1',C1,'S2',S2,'C2',C2,'S12',S12,'C12',C12);
end

% ---------- Summary (cells) ----------
summaryTable = table(names, repmat("ATT(g,t)",R,1), ATT, g_list, t_list, ...
    SE, tStat, pValue, ...
    'VariableNames', {'Name','Effect','Estimate','g','t','SE','tStat','pValue'});

% ---------- Aggregations + year labels (one line) ----------
[summaryTable, Agg] = make_aggregates_( ...
    summaryTable, ATT, g_list, t_list, G, ...
    att_b, S_for_aggs, ...
    'TimeValues', Tvals, ...   % ordered unique labels of your timeVar
    'Delta', opts.Delta, ...   % anticipation window used to form refYear
    'CellWeights', cell_wt, ...% NEW: per-(g,t) treated-post weights
    'Weighting',  opts.Weighting);  % NEW: toggle

% ---------- Output ----------
out = struct();
out.Estimator    = "CS2021";
if isfield(opts,'CrossFit') && logical(opts.CrossFit)
    cfTag = " CF";
else
    cfTag = "";
end

out.Method = "CS2021 (" + Approach + "; " + Comparison + ...
             "; δ=" + string(opts.Delta) + ", SE=" + SEMethod + cfTag + ...
             "; W=" + string(opts.Weighting) + ")";

out.Params       = opts;
out.summaryTable = summaryTable;
out.SimBands     = Bands;
out.Aggregates   = Agg;
if ~isfield(out,'Diagnostics') || ~isstruct(out.Diagnostics)
    out.Diagnostics = struct();
end
out.Diagnostics.cs = diagCS;

if opts.Print
    fprintf('[CS] R=%d cells; Approach=%s; Comp=%s; δ=%d; SE=%s%s; W=%s. ', ...
        R, Approach, Comparison, opts.Delta, SEMethod, cfTag, string(opts.Weighting));
    if SEMethod=="multiplier" && opts.B>0
        fprintf('Bootstrap B=%d, crit(95%%)=%.3f\n', opts.B, Bands.crit);
    else
        fprintf('\n');
    end

    fprintf('[CS] Overall ATT \n');
    if isstruct(out.Aggregates.overall)
        display(struct2table(out.Aggregates.overall));
    else
        display(out.Aggregates.overall);
    end

    fprintf('[CS] Estimates \n');
    display(out.summaryTable)

    fprintf('[CS] By Cohort \n');
    display(out.Aggregates.byCohort);
end
end

% ======================= Helpers =======================

function beta = wls_(X, y, w)
if isempty(w), w = ones(size(y)); end
W = sqrt(max(w,0));
Xw = X .* W; yw = y .* W;
beta = Xw \ yw;
end

function p = fit_logit_prob_(X, y, w)
if exist('glmfit','file')==2
    if isempty(w), w = ones(size(y)); end
    ps = glmfit(X, double(y), 'binomial','link','logit','weights',w);
    p  = 1./(1+exp(-[ones(size(X,1),1), X]*ps));
else
    X1 = [ones(size(X,1),1), X];
    beta = zeros(size(X1,2),1);
    for it=1:25
        z = X1*beta;
        mu = 1./(1+exp(-z));
        g  = X1'*(double(y)-mu);
        W  = diag(mu.*(1-mu) + 1e-6);
        H  = X1'*W*X1;
        step = H \ g;
        beta = beta + step;
        if norm(step) < 1e-6, break; end
    end
    p = 1./(1+exp(-X1*beta));
end
p = min(max(p, 1e-6), 1-1e-6);
end

% ---- cross-fitted propensity (by unit folds) ----
function p = fit_logit_prob_cf_(X, y, w, rowFold, K)
if isempty(rowFold), p = fit_logit_prob_(X,y,w); return; end
p = NaN(size(y));
for k=1:K
    test = (rowFold==k);
    train = ~test;
    if ~any(test) || ~any(train)
        p(test) = fit_logit_prob_(X(test,:), y(test), w(test));
    else
        if exist('glmfit','file')==2
            ps = glmfit(X(train,:), double(y(train)), 'binomial','link','logit','weights',w(train));
            p(test) = 1./(1+exp(-[ones(sum(test),1), X(test,:)]*ps));
        else
            X1 = [ones(sum(train),1), X(train,:)];
            beta = zeros(size(X1,2),1);
            for it=1:25
                z = X1*beta; mu = 1./(1+exp(-z));
                g  = X1'*(double(y(train))-mu);
                W  = diag(mu.*(1-mu) + 1e-6);
                H  = X1'*W*X1;
                step = H \ g;
                beta = beta + step;
                if norm(step) < 1e-6, break; end
            end
            p(test) = 1./(1+exp(-[ones(sum(test),1), X(test,:)]*beta));
        end
    end
end
p = min(max(p, 1e-6), 1-1e-6);
end

function XI = draw_multipliers_(n, B, kind)
if B<=0, XI = zeros(n,0); return; end
switch kind
    case "mammen"
        a = (1 - sqrt(5))/2; b = (1 + sqrt(5))/2;
        p = (sqrt(5)+1)/(2*sqrt(5));
        U = rand(n,B);
        XI = b*(U>p) + a*(U<=p);
    otherwise
        XI = 2*(rand(n,B)>0.5) - 1;
end
end

function [cIdx, C, report] = cluster_index_per_unit_(T, idIdx, N, varname, defaultUnit)
% Map row-level cluster var -> unit-constant cluster index
    report = struct('var',string(varname),'nUnits',N,'nClusters',NaN, ...
                    'nConflictedUnits',0,'conflictedUnitExamples',[]);

    useUnit = defaultUnit || isempty(varname) || ~any(strcmp(varname, T.Properties.VariableNames));
    if useUnit
        cIdx = (1:N).';   % each unit its own cluster
        C    = N;
        report.nClusters = C;
        return;
    end

    v = T.(varname);
    if iscellstr(v) || isstring(v) || iscategorical(v)
        v = string(v); isStr = true;
    else
        v = double(v); isStr = false;
    end

    cPerUnit = NaN(N,1); nConf = 0; ex = [];
    for i = 1:N
        rows = (idIdx == i);
        vi   = v(rows);
        if isStr, vi = vi(~ismissing(vi)); else, vi = vi(~isnan(vi)); end
        if isempty(vi), continue; end
        if isStr
            uv = unique(vi);
            if numel(uv) > 1, nConf = nConf + 1; if numel(ex) < 5, ex = [ex; i]; end, end %#ok<AGROW>
            [~,idx] = max(histcounts(categorical(vi), categorical(uv)));
            cPerUnit(i) = double(categorical(uv(idx)));
        else
            uv = unique(vi);
            if numel(uv) > 1, nConf = nConf + 1; if numel(ex) < 5, ex = [ex; i]; end, end %#ok<AGROW>
            [vals,~,ic] = unique(vi);
            cnts = accumarray(ic, 1);
            [~,k] = max(cnts);
            cPerUnit(i) = vals(k);
        end
    end

    if any(isnan(cPerUnit))
        nanMask = isnan(cPerUnit);
        cPerUnit(nanMask) = - (1:sum(nanMask));
    end

    [~, ~, idx] = unique(cPerUnit, 'stable');
    cIdx = idx; C = max(cIdx);

    report.nClusters = C;
    report.nConflictedUnits = nConf;
    report.conflictedUnitExamples = ex;

    if nConf > 0
        warning('did:cs:clusterVarNotUnitConstant', ...
            'ClusterVar "%s" is not unit-constant for %d units (e.g., unit(s) %s). Using modal value per unit.', ...
            string(varname), nConf, mat2str(ex(:).'));
    end
end

function S = cluster_sums_(Psi_unc, cIdx)
if size(Psi_unc,1) < size(Psi_unc,2), Psi_unc = Psi_unc.'; end
N = size(Psi_unc,1);
if numel(cIdx) ~= N
    error('did:cs:clusterIndexSize','cIdx length (%d) must match number of units N (%d).', numel(cIdx), N);
end
R = size(Psi_unc,2); C = max(cIdx);
S = zeros(R, C);
for c = 1:C
    I = (cIdx == c);
    if any(I), S(:,c) = sum(Psi_unc(I, :), 1).'; end
end
end

function [SE, tStat, pValue, Bands, S_out] = crv1_from_S_(ATT, S)
R = numel(ATT);
C = size(S,2);
if C < 2
    SE = NaN(R,1); tStat = NaN(R,1); pValue = NaN(R,1);
    Bands = struct('alpha',0.95,'crit',NaN,'lower',NaN(R,1),'upper',NaN(R,1));
    S_out = S; return;
end
S_centered = S - mean(S,2);
ssq = sum(S_centered.^2, 2);
Var_hat = (C/(C-1)) * ssq;
SE = sqrt(Var_hat);
zeroSE = (SE==0 | ~isfinite(SE));
tStat = NaN(R,1); pValue = NaN(R,1);
tStat(~zeroSE) = ATT(~zeroSE)./SE(~zeroSE);
pValue(~zeroSE) = 2*(1-normcdf(abs(tStat(~zeroSE))));
crit = 1.96;
Bands = struct('alpha',0.95,'crit',crit,'lower',ATT-crit.*SE,'upper',ATT+crit.*SE);
S_out = S;
end

function [summaryTable, Agg] = make_aggregates_( ...
        summaryTable_in, ...
        ATT, g_list, t_list, ...
        G, att_b, Sinfo, varargin)

ip = inputParser;
addParameter(ip,'TimeValues',[]);
addParameter(ip,'Delta',0);
% NEW params:
addParameter(ip,'CellWeights',[]);         % length R: Wtreat per (g,t)
addParameter(ip,'Weighting',"cohortShare");
parse(ip,varargin{:});
TimeValues = ip.Results.TimeValues;
Delta      = ip.Results.Delta;
CellWts    = ip.Results.CellWeights(:);
Weighting  = lower(string(ip.Results.Weighting));

labelFrom = @(idx) map_labels_(idx, TimeValues);

summaryTable = summaryTable_in;
summaryTable.gYear   = labelFrom(summaryTable.g);
summaryTable.tYear   = labelFrom(summaryTable.t);
summaryTable.refYear = labelFrom(summaryTable.g - Delta - 1);

% cohort-share weights (existing behavior)
Gfinite = G(~isinf(G));
[uG,~,ic] = unique(Gfinite,'stable');
countG = accumarray(ic, 1, [numel(uG),1], @sum, 0);
wG = countG / max(sum(countG),1);
wGmap = containers.Map(num2cell(uG), num2cell(wG));
wg = @(gv) arrayfun(@(x) wGmap(x), gv);

% --- NEW: select a weight accessor w_fun(I) that returns a vector of weights per subset of cells I
% * cohortShare: w_fun(I) = wg(g_list(I))
% * treatedObs : w_fun(I) = CellWts(I)
w_fun = @(I) wg(g_list(I));   % default
if Weighting == "treatedobs"
    if isempty(CellWts) || numel(CellWts) ~= numel(ATT)
        warning('did:cs:treatedObsWeightsMissing','CellWeights missing or wrong length; falling back to cohortShare.');
    else
        w_fun = @(I) CellWts(I);
    end
end

% event-time
e_list   = t_list - g_list;
Evals    = unique(e_list(e_list>=0));
theta_es = NaN(numel(Evals),1);
for k=1:numel(Evals)
    ek = Evals(k); I  = find(e_list==ek);
    if ~isempty(I)
        w = w_fun(I);
        theta_es(k) = nansum( w .* ATT(I) ) / max(nansum(w), eps);
    end
end

% calendar-time
Tuniq   = unique(t_list);
theta_c = NaN(numel(Tuniq),1);
for k=1:numel(Tuniq)
    tk = Tuniq(k);
    I  = find(t_list==tk & g_list<=tk);
    if ~isempty(I)
        w = w_fun(I);
        theta_c(k) = nansum( w .* ATT(I) ) / max(nansum(w), eps);
    end
end

% by-cohort (keep your original simple mean; can be swapped to weighted if desired)
Guniq    = unique(g_list);
theta_g  = NaN(numel(Guniq),1);
for k=1:numel(Guniq)
    gk = Guniq(k);
    I  = find(g_list==gk & t_list>=gk);
    theta_g(k) = mean(ATT(I),'omitnan');
end

% overall
pick_post   = (t_list>=g_list);
I_all       = find(pick_post);
w_all_cells = w_fun(I_all);
theta_OW = nansum( w_all_cells .* ATT(I_all) ) / max(nansum(w_all_cells),eps);

% inference for aggregates
if ~isempty(att_b)
    [es_tab, es_crit]   = agg_infer_(theta_es, build_draws_(Evals,  e_list, g_list, t_list, att_b, w_fun, 'event'));
    [cal_tab, cal_crit] = agg_infer_(theta_c,  build_draws_(Tuniq,  t_list, g_list, t_list, att_b, w_fun, 'calendar'));
    [sel_tab, sel_crit] = agg_infer_(theta_g,  build_draws_(Guniq,  g_list, t_list, t_list, att_b, w_fun, 'cohort'));
    [ow_tab,  ow_crit]  = agg_infer_(theta_OW, build_draws_(NaN,    [],     g_list, t_list, att_b, w_fun, 'overall', pick_post));
    Agg.es       = addvars(table(Evals,'VariableNames',{'e'}), es_tab.Estimate, es_tab.SE, es_tab.tStat, es_tab.pValue, es_tab.LB, es_tab.UB, ...
        'NewVariableNames',{'Estimate','SE','tStat','pValue','LB','UB'});
    Agg.calendar = addvars(table(Tuniq,'VariableNames',{'t'}), cal_tab.Estimate, cal_tab.SE, cal_tab.tStat, cal_tab.pValue, cal_tab.LB, cal_tab.UB, ...
        'NewVariableNames',{'Estimate','SE','tStat','pValue','LB','UB'});
    Agg.byCohort = addvars(table(Guniq,'VariableNames',{'g'}), sel_tab.Estimate, sel_tab.SE, sel_tab.tStat, sel_tab.pValue, sel_tab.LB, sel_tab.UB, ...
        'NewVariableNames',{'Estimate','SE','tStat','pValue','LB','UB'});
    Agg.overall  = ow_tab;
    Agg.Crit     = struct('es',es_crit,'calendar',cal_crit,'byCohort',sel_crit,'overall',ow_crit);
else
    if ~isempty(Sinfo)
       W_es   = weight_matrix_(Evals,  e_list, g_list, t_list, w_fun, 'event');
       W_cal  = weight_matrix_(Tuniq,  t_list, g_list, t_list, w_fun, 'calendar');
       W_sel  = weight_matrix_(Guniq,  g_list, t_list, t_list, w_fun, 'cohort');
       W_over = weight_matrix_(NaN,    [],     g_list, t_list, w_fun, 'overall', pick_post);

       [SE_es, t_es, p_es, LB_es, UB_es] = var_from_S_lincombo_(theta_es, W_es,  Sinfo);
       [SE_c,  t_c,  p_c,  LB_c,  UB_c ] = var_from_S_lincombo_(theta_c,  W_cal, Sinfo);
       [SE_g,  t_g,  p_g,  LB_g,  UB_g ] = var_from_S_lincombo_(theta_g,  W_sel, Sinfo);
       [SE_ow, t_ow, p_ow, LB_ow, UB_ow] = var_from_S_lincombo_(theta_OW, W_over, Sinfo);

        Agg.es       = table(Evals, theta_es, SE_es, t_es, p_es, LB_es, UB_es, ...
            'VariableNames',{'e','Estimate','SE','tStat','pValue','LB','UB'});
        Agg.calendar = table(Tuniq, theta_c,  SE_c,  t_c,  p_c,  LB_c,  UB_c, ...
            'VariableNames',{'t','Estimate','SE','tStat','pValue','LB','UB'});
        Agg.byCohort = table(Guniq, theta_g,  SE_g,  t_g,  p_g,  LB_g,  UB_g, ...
            'VariableNames',{'g','Estimate','SE','tStat','pValue','LB','UB'});
        Agg.overall  = struct('Estimate',theta_OW,'SE',SE_ow,'tStat',t_ow,'pValue',p_ow,'LB',LB_ow,'UB',UB_ow,'crit',1.96);
        Agg.Crit     = struct('es',1.96,'calendar',1.96,'byCohort',1.96,'overall',1.96);
    else
        Agg.es       = table(Evals, theta_es, 'VariableNames',{'e','Estimate'});
        Agg.calendar = table(Tuniq, theta_c,  'VariableNames',{'t','Estimate'});
        Agg.byCohort = table(Guniq, theta_g,  'VariableNames',{'g','Estimate'});
        Agg.overall  = struct('Estimate',theta_OW,'SE',NaN,'tStat',NaN,'pValue',NaN,'LB',NaN,'UB',NaN,'crit',NaN);
        Agg.Crit     = struct('es',NaN,'calendar',NaN,'byCohort',NaN,'overall',NaN);
    end
end

% labels on aggregates
if istable(Agg.byCohort) && ismember('g', Agg.byCohort.Properties.VariableNames)
    Agg.byCohort.gYear = labelFrom(Agg.byCohort.g);
end
if istable(Agg.calendar) && ismember('t', Agg.calendar.Properties.VariableNames)
    Agg.calendar.tYear = labelFrom(Agg.calendar.t);
end
end

% ====== nested helpers (scoped) ======
function outL = map_labels_(idx, L)
    if isempty(L), outL = idx; return; end
    L = L(:);
    outL = repmat(missingLike_(L), size(idx));
    valid = idx >= 1 & idx <= numel(L) & isfinite(idx);
    if any(valid)
        if isdatetime(L)
            tmp = NaT(size(idx)); tmp(valid) = L(idx(valid)); outL = tmp;
        elseif isstring(L) || iscellstr(L)
            tmp = strings(size(idx)); tmp(valid) = string(L(idx(valid))); outL = tmp;
        elseif iscategorical(L)
            tmp = categorical(missing(size(idx))); %#ok<CTPCT>
            tmp(valid) = L(idx(valid)); outL = tmp;
        else
            tmp = NaN(size(idx)); tmp(valid) = double(L(idx(valid))); outL = tmp;
        end
    end
end

function m = missingLike_(L)
    if isdatetime(L),      m = NaT;
    elseif isstring(L),    m = string(missing);
    elseif iscategorical(L), m = categorical(missing);
    else,                  m = NaN;
    end
end


function A = weight_matrix_(axisVals, a_list, b_list, t_list, w_fun, kind, pick_post)
% Builds linear-combo matrix A so that theta = A * ATT (for variance via S)
% a_list / b_list correspond to the first / second index that define cells.
R = numel(b_list);  %#ok<NASGU> % (only used implicitly through find on conditions)
if nargin==7 && strcmpi(string(kind),'overall')
    % overall: single row
    A = zeros(1, numel(t_list));
    I = find(pick_post);
    if ~isempty(I)
        w = w_fun(I);
        A(1, I) = w / max(sum(w), eps);
    end
    return
end

knd = "event"; if nargin>=6 && ~isempty(kind), knd = lower(string(kind)); end
Rcells = numel(t_list);
A = zeros(0, Rcells); % will set below

switch knd
    case "calendar"
        Tuniq = axisVals;
        A = zeros(numel(Tuniq), Rcells);
        for k=1:numel(Tuniq)
            tk = Tuniq(k);
            I  = find(t_list==tk & b_list<=tk);
            if isempty(I), continue; end
            w = w_fun(I);
            A(k,I) = w / max(sum(w), eps);
        end
    case "cohort"
        Guniq = axisVals;
        A = zeros(numel(Guniq), Rcells);
        for k=1:numel(Guniq)
            gk = Guniq(k);
            I  = find(b_list==gk & t_list>=gk);
            if isempty(I), continue; end
            % keep simple average for "byCohort" as in original
            A(k,I) = 1/numel(I);
        end
    otherwise % "event"
        Evals  = axisVals;
        e_list = t_list - b_list;
        A = zeros(numel(Evals), Rcells);
        for k=1:numel(Evals)
            ek = Evals(k);
            I  = find(e_list==ek);
            if isempty(I), continue; end
            w = w_fun(I);
            A(k,I) = w / max(sum(w), eps);
        end
end
end

function draws = build_draws_(axisVals, a_list, b_list, t_list, att_b, w_fun, kind, pick_post)
B = size(att_b,2);

% overall (single-number) case
if nargin>=8 && strcmpi(string(kind),'overall')
    draws = NaN(1,B);
    I = find(pick_post);
    if isempty(I), return; end
    w = w_fun(I);
    for b=1:B
        draws(1,b) = nansum(w .* att_b(I,b)) / max(sum(w), eps);
    end
    return
end

knd = "event";
if nargin>=7 && ~isempty(kind), knd = lower(string(kind)); end

switch knd
    case "calendar"
        Tuniq = axisVals;
        draws = NaN(numel(Tuniq), B);
        for k=1:numel(Tuniq)
            tk = Tuniq(k);
            I  = find(t_list==tk & b_list<=tk);
            if isempty(I), continue; end
            w = w_fun(I);
            for b=1:B
                draws(k,b) = nansum(w .* att_b(I,b)) / max(sum(w), eps);
            end
        end

    case "cohort"
        Guniq = axisVals;
        draws = NaN(numel(Guniq), B);
        for k=1:numel(Guniq)
            gk = Guniq(k);
            I  = find(b_list==gk & t_list>=gk);
            if isempty(I), continue; end
            for b=1:B
                draws(k,b) = mean(att_b(I,b), 'omitnan');  % keep simple average (as original)
            end
        end

    otherwise  % "event"
        Evals  = axisVals;
        e_list = t_list - b_list;   % event-time offsets
        draws  = NaN(numel(Evals), B);
        for k=1:numel(Evals)
            ek = Evals(k);
            I  = find(e_list==ek);
            if isempty(I), continue; end
            w = w_fun(I);
            for b=1:B
                draws(k,b) = nansum(w .* att_b(I,b)) / max(sum(w), eps);
            end
        end
end
end

function [SEv, tSv, pV, LB, UB] = var_from_S_lincombo_(theta, A, Sinfo)
K = size(A,1); SEv = NaN(K,1); tSv = NaN(K,1); pV = NaN(K,1); LB = NaN(K,1); UB = NaN(K,1);
if isempty(Sinfo), return; end
crit = 1.96;
if Sinfo.kind=="clustered"
    S1 = Sinfo.S1; C1 = Sinfo.C1;
    S1c = S1 - mean(S1,2);
    for j=1:K
        a = A(j,:).';
        s = (a.' * S1c);
        Var = (C1/(C1-1)) * sum(s.^2);
        se = sqrt(Var); SEv(j)=se;
        if isfinite(se) && se>0
            tSv(j) = theta(j)/se;
            pV(j)  = 2*(1-normcdf(abs(tSv(j))));
            LB(j)  = theta(j) - crit*se;
            UB(j)  = theta(j) + crit*se;
        end
    end
else
    S1 = Sinfo.S1; C1 = Sinfo.C1;
    S2 = Sinfo.S2; C2 = Sinfo.C2;
    S12= Sinfo.S12; C12= Sinfo.C12;
    S1c  = S1  - mean(S1, 2);
    S2c  = S2  - mean(S2, 2);
    S12c = S12 - mean(S12,2);
    for j=1:K
        a = A(j,:).';
        s1  = (a.' * S1c);   V1  = (C1 /(C1 -1)) * sum(s1.^2);
        s2  = (a.' * S2c);   V2  = (C2 /(C2 -1)) * sum(s2.^2);
        s12 = (a.' * S12c);  V12 = (C12/(C12-1)) * sum(s12.^2);
        Var = V1 + V2 - V12;
        se = sqrt(Var); SEv(j)=se;
        if isfinite(se) && se>0
            tSv(j) = theta(j)/se;
            pV(j)  = 2*(1-normcdf(abs(tSv(j))));
            LB(j)  = theta(j) - crit*se;
            UB(j)  = theta(j) + crit*se;
        end
    end
end
end

function [tab, crit] = agg_infer_(theta, draws)
K = numel(theta);
if isempty(draws)
    tab = table(theta, NaN(K,1), NaN(K,1), NaN(K,1), NaN(K,1), NaN(K,1), ...
        'VariableNames', {'Estimate','SE','tStat','pValue','LB','UB'});
    crit = NaN; return;
end
dc = draws - mean(draws,2);
SE = std(dc,0,2);
tStat = theta ./ max(SE, eps);
Z = abs((draws - theta) ./ max(SE, eps));
thr = abs(tStat);
pValue = mean(Z >= thr, 2);
maxZ = max(Z, [], 1);
crit = did.quantile(maxZ, 0.95);
LB = theta - crit .* SE;
UB = theta + crit .* SE;
tab = table(theta, SE, tStat, pValue, LB, UB, ...
    'VariableNames', {'Estimate','SE','tStat','pValue','LB','UB'});
end

% ---------- Cross-fitting utilities ----------

function unitFold = unit_kfold_(N, K, seed)
% Assigns each unit 1..N to one of K folds (balanced, deterministic given seed)
if nargin<3 || isnan(seed)||isempty(seed) || ~isfinite(seed), seed = randi([1,1e7],1,1); end
rng(floor(double(seed)));
perm = randperm(N);
foldSizes = repmat(floor(N/K), K, 1);
rem = N - sum(foldSizes);
foldSizes(1:rem) = foldSizes(1:rem) + 1;
unitFold = zeros(N,1);
idx = 1;
for k=1:K
    take = foldSizes(k);
    unitFold(perm(idx:idx+take-1)) = k;
    idx = idx + take;
end
end

function [mt_all, mu_reg_t] = cf_or_predict_(Xt_ti, dY, wt_u, treat_u, comp_u, eligible, idIdx, unitFold)
% Cross-fitted outcome regression on controls, predictions for treated
K = max(unitFold);
mt_all = zeros(size(dY));   % N×1
w_t_all = 0; s_num = 0;
for k = 1:K
    test_treat = treat_u & (unitFold == k);
    if ~any(test_treat), continue; end
    train_ctrl = comp_u & eligible & (unitFold ~= k);
    if ~any(train_ctrl), train_ctrl = comp_u & eligible; if ~any(train_ctrl), continue; end, end

    Xc  = Xt_ti(train_ctrl,:); yc = dY(train_ctrl); wc = wt_u(train_ctrl);
    beta = wls_([ones(sum(train_ctrl),1), Xc], yc, wc);

    Xt_k = Xt_ti(test_treat,:);
    mt_k = [ones(sum(test_treat),1), Xt_k]*beta;

    mt_all(test_treat) = mt_k;

    wt_k = wt_u(test_treat);
    w_t_all = w_t_all + nansum(wt_k);
    s_num   = s_num   + nansum(wt_k .* mt_k);
end
if w_t_all > 0
    mu_reg_t = s_num / w_t_all;
else
    mu_reg_t = 0;
end
end

function m0_cf = cf_or_predict_controls_(Xt_ti, dY, wt_u, comp_u, eligible, idIdx, unitFold)
% Cross-fitted m0(X) for controls (unit-level)
K = max(unitFold);
m0_cf = zeros(size(dY));  % N×1
for k = 1:K
    test_ctrl  = comp_u & eligible & (unitFold == k);
    if ~any(test_ctrl), continue; end
    train_ctrl = comp_u & eligible & (unitFold ~= k);
    if ~any(train_ctrl), train_ctrl = comp_u & eligible; end

    Xc  = Xt_ti(train_ctrl,:); yc = dY(train_ctrl); wc = wt_u(train_ctrl);
    beta = wls_([ones(sum(train_ctrl),1), Xc], yc, wc);

    Xc_te = Xt_ti(test_ctrl,:);
    m0_cf(test_ctrl) = [ones(sum(test_ctrl),1), Xc_te]*beta;
end
end

function unitFold = unit_kfold_stratified_cohort_(G, K, seed)
% Cohort-stratified K-fold assignment at the *unit* level.
if nargin<3 || isempty(seed) ||isnan(seed) || ~isfinite(seed), seed = randi([1,1e7],1,1), end
rng(floor(double(seed)));

N = numel(G);
unitFold = zeros(N,1);

% Map Inf (never-treated) to a dedicated label
Glab = G;
if any(isinf(Glab))
    maxFinite = max(Glab(~isinf(Glab)));
    if isempty(maxFinite) || ~isfinite(maxFinite), maxFinite = 0; end
    Glab(isinf(Glab)) = maxFinite + 1;   % NT becomes its own highest label
end

% Stable unique cohort labels (including NT label)
[uCoh, ~, ic] = unique(Glab, 'stable');  %#ok<ASGLU>
S = numel(uCoh);

for s = 1:S
    I = find(ic==s);
    ns = numel(I);
    if ns==0, continue; end
    perm = I(randperm(ns));
    folds_s = repmat((1:K).', ceil(ns/K), 1);
    folds_s = folds_s(1:ns);
    unitFold(perm) = folds_s;
end

unitFold = max(1, min(K, unitFold));
end
