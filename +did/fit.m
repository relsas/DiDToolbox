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
% Last change: 10/09/2025
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
    opts.cluster string = string.empty
    opts.B (1,1) double {mustBeInteger,mustBePositive} = 100
    opts.multiplier (1,1) string = "mammen"
    opts.studentize (1,1) logical = true
    opts.Seed double = randi([1,1e7],1,1)

    % ---- BJS passthrough options ----
    opts.Covariates string = string.empty(1,0)
    opts.Horizons double = []
    opts.Balanced (1,1) logical = false
    % EXTEND SEMethod admissible values to include CS:
    opts.SEMethod (1,1) string {mustBeMember(opts.SEMethod,["LOO","BootstrapUnit","None","multiplier","clustered","clustered2"])} = "clustered"
    opts.BootReps (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.Display (1,1) logical = true
    opts.CovarSample (1,1) string {mustBeMember(opts.CovarSample,["D0","never","all"])} = "D0"
    opts.useParallel double = 1
    opts.details = false % wooldridge

    % ---- NEW: clustering options for CS ----
    opts.ClusterVar string = string.empty
    opts.ClusterVar2 string = string.empty

    % ---- NEW: CS controls (so users can pass them via did.fit) ----
    opts.Approach (1,1) string = "unconditional"
    opts.Comparison (1,1) string = "never"
    opts.Delta (1,1) double = 0
    opts.CrossFit (1,1) logical=false
    opts.Kfolds (1,1) double=5
    opts.StratifyFoldsBy (1,1) string = "none";
    opts.Print (1,1) logical=true
    opts.Weighting string {mustBeMember(opts.Weighting,["cohortShare","treatedObs"])} = "treatedObs"
    opts.describe logical=true
end

% ---- Normalize method & resolve dataset/table ----
m = lower(string(method));

if isa(T_or_ds, 'did.Dataset')
    % Dataset path (no repeated var names needed)
    ds = T_or_ds;
    % T  = ds.T;
    % idVar   = ds.idVar;
    % timeVar = ds.timeVar;
    % yVar    = ds.yVar;
    % dVar    = ds.dVar;
    % gVar    = ds.gVar;
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
        'vcov',opts.vcov, 'clusters',opts.clusters, 'details',opts.details,'Print',opts.Print);

    res = out; if isfield(out,'vcov') && ~isfield(out,'Vcov'), res.Vcov = out.vcov; end
    return
end


% ---- Common pipeline for everything else (TWFE, CS, CH, BJS, IW, …) ----
% Build estimator from opts
est = did.factories.getEstimator(method, opts);

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

% Pretty print for TWFE (unchanged)
if isfield(res,'Method') && res.Method=="TWFE"
    fprintf('[TWFE] Standard errors clustered by %s \n',join(res.Vcov.clusters,", "));
    fprintf('[TWFE]\n'); display(res.summaryTable);
end
end
