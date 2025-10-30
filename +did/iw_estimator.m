function out = iw_estimator(T, opts)
% DID.IW_ESTIMATOR  Sun & Abraham (2021) Interaction-Weighted (IW) event-study.
%
%  out = did.iw_estimator(T, ...
%        idVar="id", timeVar="time", yVar="y", dVar="D", WeightVar="", ...
%        Delta=0, Comparison="notyet", ...
%        SEMethod="multiplier", B=999, Seed=42, Multiplier="rademacher", ...
%        Studentize=true, ClusterVar="", ClusterVar2="", ...
%        Weighting="cohortShare", ...
%        Print=true)
%
% Notes:
% - Uses the CS engine to estimate ATT(g,t) and then relies on its event-time
%   aggregation (θ(e)) under the chosen Weighting. For SA IW, use "cohortShare".
% - Approach is fixed to "unconditional" (no X-modeling) to match IW.
% - Inference (SE/bands) comes directly from the CS call you choose.
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/13/2025  (Dataset-aware passthrough; Weighting forward; arguments NV-pairs)
% ------------------------------------------------------------------------

arguments
    T table

    % Data columns
    opts.idVar       (1,1) string = "id"
    opts.timeVar     (1,1) string = "time"
    opts.yVar        (1,1) string = "y"
    opts.dVar        (1,1) string = "D"
    opts.WeightVar         string = string.empty

    % IW/CS controls
    opts.Delta       (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.Comparison  (1,1) string {mustBeMember(opts.Comparison,["notyet","never"])} = "notyet"
    opts.SEMethod    (1,1) string {mustBeMember(opts.SEMethod,["multiplier","clustered","clustered2"])} = "multiplier"
    opts.B           (1,1) double {mustBeInteger, mustBeNonnegative} = 999
    opts.Seed        (1,1) double = NaN
    opts.Multiplier  (1,1) string {mustBeMember(opts.Multiplier,["rademacher","mammen"])} = "rademacher"
    opts.Studentize  (1,1) logical = true
    opts.ClusterVar        string = string.empty
    opts.ClusterVar2       string = string.empty

    % Aggregation weighting (forwarded to CS):
    %   "cohortShare" → classic SA IW (weights constant within cohort across post t)
    %   "treatedObs"  → Deb–Norton–Wooldridge–Zabel post-only treated counts per (g,t)
    opts.Weighting   (1,1) string {mustBeMember(opts.Weighting,["cohortShare","treatedObs"])} = "cohortShare"

    % I/O
    opts.Print       (1,1) logical = true
end

% ---------- 1) Run CS engine to get ATT(g,t) + inference carriers ----------
% IW uses unconditional CS ATT(g,t); forward relevant options 1:1.
nv = {
    'idVar',opts.idVar, 'timeVar',opts.timeVar, 'yVar',opts.yVar, 'dVar',opts.dVar, ...
    'WeightVar',opts.WeightVar, ...
    'Approach',"unconditional", ...
    'Comparison',lower(string(opts.Comparison)), ...
    'Delta',opts.Delta, ...
    'SEMethod',lower(string(opts.SEMethod)), ...
    'B',opts.B, 'Seed',opts.Seed, 'Multiplier',lower(string(opts.Multiplier)), ...
    'Studentize',opts.Studentize, ...
    'ClusterVar',opts.ClusterVar, 'ClusterVar2',opts.ClusterVar2, ...
    'Weighting',opts.Weighting, ...   % ← aggregation scheme (SA or DNWZ)
    'Print',false ...
    };

cs = did.cs_estimator(T, nv{:});

% ---------- 2) Package as Sun–Abraham IW ----------
% IW estimand: θ(e) = Σ_g w_g ATT(g,g+e); CS Aggregates.es provides θ(e)
out = struct();
out.Estimator = "IW";
out.Method    = "IW ("+lower(string(opts.Comparison))+"; δ="+string(opts.Delta)+ ...
                ", SE="+lower(string(opts.SEMethod))+"; W="+string(opts.Weighting)+")";
out.Params    = opts;

% Relay CS aggregates into a familiar block
out.Aggregates = cs.Aggregates;

% Also expose top-level fields for adapters/plotters
if isfield(cs,'Aggregates') && istable(cs.Aggregates.es)
    out.es = cs.Aggregates.es;
end
if isfield(cs,'Aggregates') && istable(cs.Aggregates.byCohort)
    out.byCohort = cs.Aggregates.byCohort;
end
if isfield(cs,'Aggregates') && istable(cs.Aggregates.calendar)
    out.calendar = cs.Aggregates.calendar;
end

% Forward overall (if present)
if isfield(cs,'Aggregates') && isfield(cs.Aggregates,'overall')
    out.overall = cs.Aggregates.overall;  % struct with Estimate, SE, LB, UB, crit
end

% Diagnostics passthrough (optional)
out.Diagnostics = struct();
if isfield(cs,'Diagnostics'), out.Diagnostics.cs = cs.Diagnostics; end

% ---------- Build a summaryTable (for Model printing) ----------
sumRows = strings(0,1);
effLab  = strings(0,1);
est     = [];
se      = [];
tstat   = [];
pval    = [];

% Event-time rows θ(e)
if isfield(out,'es') && istable(out.es) && all(ismember({'e','Estimate'}, out.es.Properties.VariableNames))
    Ne = height(out.es);
    sumRows = [sumRows; "ES e="+string(out.es.e)];
    effLab  = [effLab;  repmat("θ(e)", Ne, 1)];
    est     = [est;     out.es.Estimate];
    if ismember('SE', out.es.Properties.VariableNames)
        se = [se; out.es.SE];
    else
        se = [se; NaN(Ne,1)];
    end
    if all(ismember({'tStat','pValue'}, out.es.Properties.VariableNames))
        tstat = [tstat; out.es.tStat];
        pval  = [pval;  out.es.pValue];
    else
        ts = NaN(Ne,1); pv = NaN(Ne,1);
        nz = isfinite(se) & se>0;
        ts(nz) = est(nz)./se(nz);
        pv(nz) = 2*(1-normcdf(abs(ts(nz))));
        tstat = [tstat; ts];
        pval  = [pval;  pv];
    end
end

% Optional overall row
if isfield(out,'overall') && isstruct(out.overall) && isfield(out.overall,'Estimate')
    sumRows = [sumRows; "Overall"];
    effLab  = [effLab;  "θ(overall)"];
    est     = [est;     out.overall.Estimate];
    if isfield(out.overall,'SE')
        se = [se; out.overall.SE];
        if isfinite(out.overall.SE) && out.overall.SE>0
            ts = out.overall.Estimate / out.overall.SE;
            pv = 2*(1-normcdf(abs(ts)));
        else
            ts = NaN; pv = NaN;
        end
    else
        se = [se; NaN];
        ts = NaN; pv = NaN;
    end
    tstat = [tstat; ts];
    pval  = [pval;  pv];
end

out.summaryTable = table(sumRows, effLab, est, se, tstat, pval, ...
    'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue'});

% Pretty print (compact)
if opts.Print
    fprintf('[IW] %s\n', out.Method);
    if isfield(out,'es') && istable(out.es)
        K = height(out.es);
        fprintf('[IW] K=%d event-time points. ', K);
        if isfield(out,'overall') && ~isempty(out.overall) && isfield(out.overall,'Estimate')
            estO = out.overall.Estimate;
            if isfield(out.overall,'SE') && isfinite(out.overall.SE)
                seO  = out.overall.SE;
                fprintf('Overall=%.4f (SE=%.4f)\n', estO, seO);
            else
                fprintf('Overall=%.4f\n', estO);
            end
        else
            fprintf('\n');
        end
    end
    if isfield(out,'es') && istable(out.es)
        toShow = intersect({'e','Estimate','SE','tStat','pValue','LB','UB'}, out.es.Properties.VariableNames,'stable');
        disp(out.es(:,toShow));
    end
end
end
