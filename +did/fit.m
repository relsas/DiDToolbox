function res = fit(method, T_or_ds, opts)
%DID.FIT  Unified interface across estimators (orchestration).
%
% Accepts either:
%   did.fit(method, T,   idVar=..., timeVar=..., yVar=..., dVar=..., ...)
% or
%   did.fit(method, ds,  ...)     % where ds is did.Dataset
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 11/18/2025
% ------------------------------------------------------------------------

arguments
    % ---- Positional ----
    method (1,1) string
    T_or_ds

    % ---- Name–value in 'opts' ----
    % Required only if T_or_ds is a TABLE (legacy path)
    opts.idVar (1,1) string = ""
    opts.timeVar (1,1) string = ""
    opts.yVar (1,1) string = ""
    opts.dVar (1,1) string = ""

    % Common VCOV options
    opts.vcov (1,1) string = "clustered"
    opts.clusters string = []

    % Wild bootstrap options (for methods that use it)
    opts.clusterWild string = string.empty
    opts.B (1,1) double {mustBeInteger,mustBePositive} = 100
    opts.multiplier (1,1) string = "mammen"
    opts.studentize (1,1) logical = true
    opts.Seed double = randi([1,1e7],1,1)

    % ---- BJS passthrough options ----
    opts.Covariates string = string.empty(1,0)
    opts.Horizons double = []
    opts.Placebos double = []
    opts.Balanced (1,1) logical = false
    % EXTEND SEMethod admissible values to include CS:
    opts.SEMethod (1,1) string {mustBeMember(opts.SEMethod, {'LOO','BootstrapUnit','None','multiplier','clustered','clustered2','Placebo'})} = "clustered"
    opts.BootReps (1,1) double {mustBeInteger, mustBeNonnegative} = 199
    opts.Display (1,1) logical = true
    opts.CovarSample (1,1) string {mustBeMember(opts.CovarSample,["D0","never","all"])} = "D0"
    opts.useParallel double = 0
    opts.details = false % Wooldridge, TWFE

    % ---- clustering options for CS ----
    opts.ClusterVar string = string.empty
    opts.ClusterVar2 string = string.empty

    % ---- CS controls (so users can pass them via did.fit) ----
    opts.Approach (1,1) string = "unconditional"
    opts.Comparison (1,1) string = "never"
    opts.Delta (1,1) double = 0
    opts.CrossFit (1,1) logical=false
    opts.Kfolds (1,1) double=5
    opts.StratifyFoldsBy (1,1) string = "none";
    opts.Weighting string {mustBeMember(opts.Weighting,["cohortShare","treatedObs"])} = "treatedObs"
    opts.PreTrendBase string {mustBeMember(opts.PreTrendBase, ["universal","varying","first","last"])} = "universal"
    opts.describe logical=true

    % ---- SDID / SC options ----
    opts.Weights (1,1) string = "SDID"
    opts.Regularization (1,1) double = NaN % Auto
    opts.Solver (1,1) string = "auto"

    % ---- CEM Options ----
    opts.nBins (1,1) double = 5
    opts.Robust (1,1) logical = false
    opts.Estimator (1,1) string = "TWFE"

    % ---- Absorption for TWFE with large individual FE
    opts.absorbMethod (1,1) string = "halperin"
    opts.absorbUseGPU (1,1) logical = false
    opts.absorbTol   (1,1) double  = 1e-8
    opts.absorbMaxIter (1,1) double = 1000

    % --- Diagnostics
    opts.preTrend (1,1) logical = false
    opts.preTrendMode (1,1) string = "pooled_preMinG"

    opts.EventStudy (1,1) logical = false
    opts.TimeEffects (1,1) logical = true

    % ---- Generic Options ----
    opts.baseline double = NaN
end

% ---- Normalize method & resolve dataset/table ----
m = lower(string(method));

if isa(T_or_ds, 'did.Dataset')
    % Dataset path (no repeated var names needed)
    ds = T_or_ds;
else
    % Legacy table path
    T = T_or_ds;
    idVar   = opts.idVar;
    timeVar = opts.timeVar;
    yVar    = opts.yVar;
    dVar    = opts.dVar;

    if strlength(idVar)==0 || strlength(timeVar)==0 || strlength(yVar)==0 || strlength(dVar)==0
        error('did:fit:MissingVarNames', ...
            'When passing a table, idVar/timeVar/yVar/dVar must be provided.');
    end
    % Sort and build Dataset (reuses your new class)
    T = sortrows(T,[idVar,timeVar]);

    ds = did.Dataset.fromTable(T, opts);
end

% ---- SPECIAL: Wooldridge staggered DiD ----
if ismember(m, ["wooldridge","twm"])
    TT = ds.materialize("gVar");     % guarantees a column with the canonical name
    gVar = ds.cohortVarName;         % usually "gVar"

    out = did.wooldridge_TB(TT, ...
        'yVar',ds.yVar, 'idVar',ds.idVar, 'timeVar',ds.timeVar, 'gVar',gVar, 'Covariates',opts.Covariates,...
        'vcov',opts.vcov, 'clusters',opts.clusters, 'details',opts.details,'Display',opts.Display, ...
        'EventStudy', opts.EventStudy);

    res = out;
    if isfield(out,'vcov') && ~isfield(out,'Vcov'), res.Vcov = out.vcov; end
    return
end

% ---- Common pipeline for everything else (TWFE, CS, CH, BJS, IW, …) ----
% Build estimator from opts
est = did.factories.getEstimator(method, opts);

% ---- NEW: toggle TWFE pretrend mode when requested ----------------------
if opts.preTrend
    % We currently only implement preTrend for TWFE
    if ~isa(est, 'did.estimators.TWFE')
        error('did:fit:preTrendOnlyTWFE', ...
            'The preTrend option is only implemented for the TWFE estimator (method="twfe").');
    end
    if isprop(est, 'preTrend')
        est.preTrend = true;
    else
        warning('did:fit:preTrendNoSupport', ...
            'TWFE estimator has no Pretrend property; preTrend option will be ignored.');
    end
end
% ------------------------------------------------------------------------

% Attach VCOV engine if available
try
    est.VcovEngine = did.factories.getVcov(opts);
catch
    % keep estimator default; not all estimators use a VCOV engine object
end

% Fit via Model(ds, est) -> estimator.fit(ds)
mdl = did.Model(ds, est);
res = mdl.fit();

% Normalize field name if your estimators store it in .vcov
if isfield(res,'vcov') && ~isfield(res,'Vcov')
    res.Vcov = res.vcov;
end

% Pretty print for TWFE and TWFE_pretrend (respect Display/Print/describe)
if isfield(res,'Method') && res.Method=="TWFE"
    if opts.Display
        fprintf('[TWFE] DepVar: %s\n', ds.yVar);
        fprintf('[TWFE] Obs: %d | Units: %d | Time: %d\n', ...
            res.Details.Nobs, res.Details.Nunits, res.Details.Ntime);
        if isfield(res.Details, 'R2')
            fprintf('[TWFE] Within-R2: %.4f\n', res.Details.R2);
        end
        if ismissing(res.Details.DroppedTimeFE)
            fe_str = "Time";
        else
            fe_str = "None (" + res.Details.DroppedTimeFE + ")";
        end
        fprintf('[TWFE] Fixed Effects: Unit, %s\n', fe_str);
        if isfield(res,'Vcov') && isfield(res.Vcov,'clusters')
            fprintf('[TWFE] Standard errors clustered by %s \n', join(res.Vcov.clusters, ", "));
        end
        display(res.summaryTable);
    end

elseif isfield(res,'Method') && res.Method=="TWFE_pretrend"
    if opts.Display

        % ---- Mode + window info ----
        modeStr = "pooled_preMinG";
        if isfield(res,'Details') && isfield(res.Details,'Mode') ...
                && ~isempty(res.Details.Mode)
            modeStr = string(res.Details.Mode);
        end
        fprintf('[TWFE pretrend] Mode: %s\n', modeStr);

        % Optional: print pre window and baseline if available (mainly pooled modes)
        if isfield(res,'Details')
            dDet = res.Details;
            if isfield(dDet,'PreStart') && isfield(dDet,'PreEnd') ...
                    && isfinite(dDet.PreStart) && isfinite(dDet.PreEnd)
                fprintf('  Pre window: [%g, %g]\n', dDet.PreStart, dDet.PreEnd);
            end
            if isfield(dDet,'Baseline') && isfinite(dDet.Baseline)
                fprintf('  Baseline period: %g\n', dDet.Baseline);
            end
        end

        % Clustering info (if VCOV engine produced it)
        if isfield(res,'Vcov') && isfield(res.Vcov,'clusters') ...
                && ~isempty(res.Vcov.clusters)
            fprintf('[TWFE pretrend] Standard errors clustered by %s \n', ...
                join(res.Vcov.clusters, ", "));
        end

        % Header depends on mode
        if opts.details
            % Header depends on mode
            if modeStr == "cohortwise_notYetTreated"
                fprintf('[TWFE pretrend] Pretrend coefficients (cohort × pre_t):\n');
            else
                fprintf('[TWFE pretrend] Pretrend coefficients (isTreated × pre_t):\n');
            end

            % Full table of pretrend coefficients
            display(res.summaryTable);
        end

        % ---- Joint tests ----
        if isfield(res,'PretrendDiag')
            d = res.PretrendDiag;

            % Cohortwise mode: per-cohort joint tests
            if isfield(res,'Details') && isfield(res.Details,'Mode') ...
                    && res.Details.Mode == "cohortwise_notYetTreated" ...
                    && isfield(d,'type') && d.type == "twfe_pretrend_cohortwise"

                fprintf('[TWFE pretrend] Cohortwise joint tests H0: all pretrend coefficients = 0 (per cohort)\n');
                fprintf('  Cohort    Nobs   nGamma      F(df1,df2)        p-value\n');
                for j = 1:numel(d.cohort)
                    fprintf('  %7g  %6d  %7d   F(%d,%d)=%.3f   p=%.4f\n', ...
                        d.cohort(j), d.Nobs(j), d.nGamma(j), ...
                        d.df1(j), d.df2(j), d.F(j), d.pValue(j));
                end

            else
                % Pooled pretrend (preMinG or notYetTreated): single joint test
                if all(isfield(d, ["df1","df2","F","pValue"]))
                    fprintf('[TWFE pretrend] Joint test H0: all pretrend coefficients = 0\n');
                    fprintf('  F(%d,%d) = %.3f,  p = %.4f\n', d.df1, d.df2, d.F, d.pValue);
                end
            end
        end

        % ---- Optional graphical diagnostics (light scheme) ----
        if opts.details
            try
                if modeStr == "pooled_notYetTreated"
                    localPlotTWFEPretrendPooled(res);
                elseif modeStr == "cohortwise_notYetTreated"
                    localPlotTWFEPretrendCohortwise(res);
                elseif modeStr == "event"
                    localPlotTWFEPretrendEvent(res);
                end
            catch ME
                % Do not crash fit(); just warn if plotting fails
                warning('did:fit:PretrendPlotFailed', ...
                    'Pretrend plot failed: %s', ME.message);
            end
        end

    end

elseif isfield(res,'Method') && res.Method=="SDID"
    if opts.Display
        fprintf('[SDID] Weights: %s\n', res.WeightsType);
        t_stat = res.tau / res.se;
        p_val  = 2 * (1 - normcdf(abs(t_stat)));

        fprintf('[SDID] Aggregate Estimate: %.4f (SE: %.4f, t: %.2f, p: %.4f)\n', ...
            res.tau, res.se, t_stat, p_val);

        if isfield(res, 'CohortTable')
            fprintf('\n[SDID] Cohort Estimates:\n');
            disp(res.CohortTable);
        end
    end
end





end % function fit


function localPlotTWFEPretrendPooled(res)
% Plot pooled pretrend coefficients: ever-treated × pre_t or × t (D=0)

tab = res.summaryTable;
needed = ["Name","Estimate","SE","pValue"];
if ~all(ismember(needed, tab.Properties.VariableNames))
    return;
end

nameStr = string(tab.Name);
T = NaN(height(tab),1);

for i = 1:height(tab)
    % Pattern 1: "... pre_t=2"
    tok = regexp(nameStr(i), 'pre_t=([\-0-9\.]+)', 'tokens', 'once');

    % Pattern 2: "... t=2 (D=0)"
    if isempty(tok)
        tok = regexp(nameStr(i), 't=([\-0-9\.]+)', 'tokens', 'once');
    end

    if ~isempty(tok)
        T(i) = str2double(tok{1});
    end
end

ok = isfinite(T) & isfinite(tab.Estimate) & isfinite(tab.SE);
if ~any(ok), return; end

T   = T(ok);
est = tab.Estimate(ok);
se  = tab.SE(ok);
p   = tab.pValue(ok);

ciLo = est - 1.96*se;
ciHi = est + 1.96*se;

% Sort by time for a nicer plot
[T, idx] = sort(T);
est  = est(idx);
se   = se(idx);
p    = p(idx);
ciLo = ciLo(idx);
ciHi = ciHi(idx);

fig = figure('Name','TWFE pretrend (pooled)', 'Color','w');
% Enforce light theme where available

set(fig, 'Theme', 'light');


hold on;

% Whiskers (light grey)
for i = 1:numel(T)
    plot([T(i) T(i)], [ciLo(i) ciHi(i)], '-', ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Significance: all markers filled, light gray vs black
sig = p <= 0.1 & ~isnan(p);

% Non-significant: light gray filled circles
plot(T(~sig), est(~sig), 'o', ...
    'MarkerFaceColor',[0.5 0.5 0.5], ...
    'MarkerEdgeColor',[0 0 0], ...
    'LineStyle','none', ...
    'MarkerSize',6);

% Significant: black filled circles
plot(T(sig),  est(sig),  'o', ...
    'MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0 0 0], ...
    'LineStyle','none', ...
    'MarkerSize',6);

yline(0, '--', 'Color',[0.4 0.4 0.4]);

xlabel('Pre-period t');
ylabel('Pretrend coefficient');
title('TWFE pretrend (pooled)');
grid on;
box on;
hold off;
end

% --------------------------------------------------------------
function localPlotTWFEPretrendCohortwise(res)
% Plot cohortwise pretrend coefficients: cohort × pre_t

tab = res.summaryTable;
needed = ["Name","Estimate","SE","pValue"];
if ~all(ismember(needed, tab.Properties.VariableNames))
    return;
end

nameStr = string(tab.Name);
n = height(tab);

cohort = NaN(n,1);
T      = NaN(n,1);

for i = 1:n
    tokC = regexp(nameStr(i), 'cohort=([\-0-9\.]+)', 'tokens', 'once');
    if ~isempty(tokC)
        cohort(i) = str2double(tokC{1});
    end
    tokT = regexp(nameStr(i), 'pre_t=([\-0-9\.]+)', 'tokens', 'once');
    if ~isempty(tokT)
        T(i) = str2double(tokT{1});
    end
end

ok = isfinite(cohort) & isfinite(T) & ...
    isfinite(tab.Estimate) & isfinite(tab.SE);
if ~any(ok), return; end

cohort = cohort(ok);
T      = T(ok);
est    = tab.Estimate(ok);
se     = tab.SE(ok);
p      = tab.pValue(ok);

ciLo   = est - 1.96*se;
ciHi   = est + 1.96*se;

uC = unique(cohort);
nC = numel(uC);

% Modest layout: up to 3 columns
nCols = min(3, nC);
nRows = ceil(nC / nCols);

fig = figure('Name','TWFE pretrend (cohortwise)', 'Color','w');

set(fig,'Theme','light');


tiledlayout(nRows, nCols, 'Padding','compact', 'TileSpacing','compact');

for j = 1:nC
    c = uC(j);
    idx = (cohort == c);
    if ~any(idx), continue; end

    Tj    = T(idx);
    estj  = est(idx);
    sej   = se(idx);
    pj    = p(idx);
    ciLoj = ciLo(idx);
    ciHij = ciHi(idx);

    [Tj, order] = sort(Tj);
    estj  = estj(order);
    sej   = sej(order);
    pj    = pj(order);
    ciLoj = ciLoj(order);
    ciHij = ciHij(order);

    nexttile;
    hold on;

    % Whiskers (light grey)
    for i = 1:numel(Tj)
        plot([Tj(i) Tj(i)], [ciLoj(i) ciHij(i)], '-', ...
            'Color',[0.7 0.7 0.7], 'LineWidth',1);
    end

    % Significance markers: both filled, light gray vs black
    sig = pj <= 0.1 & ~isnan(pj);

    % Non-significant: light gray filled
    plot(Tj(~sig), estj(~sig), 'o', ...
        'MarkerFaceColor',[0.5 0.5 0.5], ...
        'MarkerEdgeColor',[0 0 0], ...
        'LineStyle','none', ...
        'MarkerSize',6);

    % Significant: black filled
    plot(Tj(sig),  estj(sig),  'o', ...
        'MarkerFaceColor',[0 0 0], ...
        'MarkerEdgeColor',[0 0 0], ...
        'LineStyle','none', ...
        'MarkerSize',6);

    yline(0, '--', 'Color',[0.4 0.4 0.4]);

    title(sprintf('Cohort g = %g', c));
    grid on;
    box on;
    hold off;

    if j > (nRows-1)*nCols
        xlabel('Pre-period t');
    end
    ylabel('\gamma_g(t)');
end
end
% --------------------------------------------------------------
function localPlotTWFEPretrendEvent(res)
% Plot event-time coefficients: Event k=...

tab = res.summaryTable;
needed = ["Name","Estimate","SE","pValue"];
if ~all(ismember(needed, tab.Properties.VariableNames))
    return;
end

nameStr = string(tab.Name);
n = height(tab);

k = NaN(n,1);

for i = 1:n
    % Pattern: "Event k=-3"
    tok = regexp(nameStr(i), 'Event k=([\-0-9\.]+)', 'tokens', 'once');
    if ~isempty(tok)
        k(i) = str2double(tok{1});
    end
end

ok = isfinite(k) & isfinite(tab.Estimate) & isfinite(tab.SE);
if ~any(ok), return; end

k   = k(ok);
est = tab.Estimate(ok);
se  = tab.SE(ok);
p   = tab.pValue(ok);

ciLo = est - 1.96*se;
ciHi = est + 1.96*se;

[k, idx] = sort(k);
est  = est(idx);
se   = se(idx);
p    = p(idx);
ciLo = ciLo(idx);
ciHi = ciHi(idx);

fig = figure('Name','TWFE pretrend (Event)', 'Color','w');
set(fig, 'Theme', 'light');

hold on;

% Whiskers (light grey)
for i = 1:numel(k)
    plot([k(i) k(i)], [ciLo(i) ciHi(i)], '-', ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Significance: all markers filled
sig = p <= 0.1 & ~isnan(p);

% Non-significant: light gray filled circles
plot(k(~sig), est(~sig), 'o', ...
    'MarkerFaceColor',[0.5 0.5 0.5], ...
    'MarkerEdgeColor',[0 0 0], ...
    'LineStyle','none', ...
    'MarkerSize',6);

% Significant: black filled circles
plot(k(sig),  est(sig),  'o', ...
    'MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0 0 0], ...
    'LineStyle','none', ...
    'MarkerSize',6);

yline(0, '--', 'Color',[0.4 0.4 0.4]);

xlabel('Event Time k');
ylabel('Pretrend Coefficient \gamma_k');
title('TWFE pretrend (Event Mode)');
grid on;
box on;
hold off;
end
