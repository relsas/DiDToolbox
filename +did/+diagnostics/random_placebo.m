function plc = random_placebo(ds, method, varargin)
% RANDOM_PLACEBO  Randomization (placebo) test for DiD estimators.
%
%   plc = did.diagnostics.random_placebo(ds, method, Name=Value, ...)
%
% INPUTS
%   ds      : did.Dataset object (original data)
%   method  : string, e.g. "twfe", "BJS", "CS", ...
%
% NAME–VALUE PAIRS
%   Placebo-specific (handled here, NOT passed to did.fit):
%     B               : # placebo replications (default 199)
%     Seed            : RNG seed (default NaN = do not reset)
%     scheme          : "permuteEverTreated" (default) or "permuteAll"
%     treatClusterVar : name of treatment cluster variable in ds.T
%                       (default: ds.idVar)
%     Verbose         : logical (default true)
%
%   All other Name=Value pairs are passed directly to did.fit,
%   e.g.: Covariates, vcov, clusters, Display, Print, details, ...
%
% OUTPUT
%   plc : struct with fields
%       .type            : "random_placebo"
%       .method          : method string
%       .B               : # placebo replications
%       .scheme          : scheme used
%       .treatClusterVar : treatment cluster variable actually used
%       .statObs         : observed statistic on original data
%       .statPlacebo     : B×1 vector of placebo statistics
%       .okMask          : logical mask of successful replications
%       .pTwoSided       : two-sided randomization p-value
%       .pGreater        : P(statPlacebo >= statObs)
%       .pLess           : P(statPlacebo <= statObs)
%       .quantiles       : [0.025 0.5 0.975] quantiles of placebo stats
% -----------------------------------------------------------------------------------
%  Dr. Ralf Elsas-Nicolle, LMU, Germany / Last change: 11/17/2025
% -----------------------------------------------------------------------------------
arguments
    ds (1,1) did.Dataset
    method (1,1) string
end

arguments (Repeating)
    varargin
end

% ---- Default placebo-specific options -----------------------------------
B_placebo    = 199;
Seed_placebo = NaN;
scheme       = "permuteEverTreated";
treatClVar   = "";        % default: will become ds.idVar
Verbose      = true;

% ---- Split varargin into placebo-specific vs fit-related ----------------
if mod(numel(varargin),2) ~= 0
    error('did:random_placebo:NVpairs', ...
        'Name–value pairs must come in key/value pairs.');
end

nv_fit = {};   % forwarded to did.fit

placeboNames = ["B","Seed","scheme","treatClusterVar","Verbose"];

for k = 1:2:numel(varargin)
    name = string(varargin{k});
    val  = varargin{k+1};

    lowerName = lower(name);

    if any(strcmpi(name, placeboNames))
        % Handle placebo-specific options here
        switch lowerName
            case "b"
                B_placebo = val;
            case "seed"
                Seed_placebo = val;
            case "scheme"
                scheme = string(val);
            case "treatclustervar"
                treatClVar = string(val);
            case "verbose"
                Verbose = logical(val);
        end
    else
        % Forward everything else to did.fit unmodified
        nv_fit(end+1:end+2) = {char(name), val}; %#ok<AGROW>
    end
end

% Validate scheme
if ~ismember(scheme, ["permuteEverTreated","permuteAll"])
    error('did:random_placebo:BadScheme', ...
        'Unknown scheme "%s". Use "permuteEverTreated" or "permuteAll".', scheme);
end

% ---- Ensure simulation runs quietly by default --------------------------
% (respect user overrides if present)
fitNames = string(nv_fit(1:2:end));

if ~any(strcmpi(fitNames, "Display"))
    nv_fit(end+1:end+2) = {'Display', false}; %#ok<AGROW>
end
if ~any(strcmpi(fitNames, "describe"))
    nv_fit(end+1:end+2) = {'describe', false}; %#ok<AGROW>
end
if ~any(strcmpi(fitNames, "details"))
    nv_fit(end+1:end+2) = {'details', false}; %#ok<AGROW>
end

% -------------------------------------------------------------------------
%    From here on, nv_fit is exactly what we pass to did.fit
% -------------------------------------------------------------------------

% ---- Extract core data from Dataset -------------------------------------
T       = ds.T;
idVar   = ds.idVar;
timeVar = ds.timeVar;
dVar    = ds.dVar;

id = T.(idVar);
tt = T.(timeVar);
D  = T.(dVar);

% Canonical helpers (cohort g; 0 for never-treated)
g = ds.get("g");

% ---- Determine treatment "cluster" level --------------------------------
if strlength(treatClVar) > 0
    if ~ismember(treatClVar, T.Properties.VariableNames)
        error('did:random_placebo:NoTreatClusterVar', ...
            'treatClusterVar "%s" not found in ds.T.', treatClVar);
    end
    cl = T.(treatClVar);
else
    % Default: each unit is its own treatment cluster
    treatClVar = idVar;
    cl = id;
end

[cl_u, ~, cl_idx] = unique(cl);    % cl_idx(i) in {1..G}
G = numel(cl_u);

% Cluster-level first-treatment time: gC (0 = never-treated cluster)
g_vec = g(:);
gC = accumarray(cl_idx, g_vec, [G,1], @minPosOrZero);

% Treated vs never-treated clusters
isTreatedCluster = (gC > 0);
idxTreatedCl     = find(isTreatedCluster);
idxNeverCl       = find(~isTreatedCluster);

nTreatCl = numel(idxTreatedCl);
if nTreatCl == 0
    error('did:random_placebo:NoTreatedClusters', ...
        'Dataset has no treated clusters; randomized placebo not meaningful.');
end

if Verbose
    fprintf('[random_placebo] %d treatment clusters (%d ever-treated, %d never-treated).\n', ...
        G, nTreatCl, numel(idxNeverCl));
    fprintf('[random_placebo] Scheme: %s, Cluster level: %s\n', ...
        scheme, treatClVar);
end

% ---- Fit the estimator on the original data (get statObs) --------------
if Verbose
    fprintf('[random_placebo] Fitting estimator "%s" on original data...\n', method);
end

resObs  = did.fit(method, ds, nv_fit{:});
statObs = extract_stat_(resObs);
if ~isfinite(statObs)
    warning('did:random_placebo:StatObsNaN', ...
        'Observed statistic (statObs) is NaN or Inf for method "%s".', method);
end

% ---- RNG handling -------------------------------------------------------
if ~isnan(Seed_placebo) && isfinite(Seed_placebo)
    rngOld = rng;
    rng(Seed_placebo,'twister');
else
    rngOld = [];
end

% ---- Placebo loop -------------------------------------------------------
B = B_placebo;
statPlacebo = NaN(B,1);

if Verbose
    fprintf('[random_placebo] Running %d placebo replications...\n', B);
end

for b = 1:B
    % 1) Generate placebo cluster-level gC^(b)
    gC_b = gC;

    switch scheme
        case "permuteEverTreated"
            perm = idxTreatedCl(randperm(nTreatCl));
            gC_b(idxTreatedCl) = gC(perm);

        case "permuteAll"
            perm = randperm(G);
            gC_b = gC(perm);
    end

    % 2) Map cluster-level gC_b back to observation-level g_b
    g_b = g;
    g_b(:) = gC_b(cl_idx);

    % 3) Rebuild placebo D^(b) as 1{t >= g_b > 0}
    D_b = zeros(size(D));
    D_b(g_b > 0 & tt >= g_b) = 1;

    % 4) Build placebo table and Dataset
    T_b = T;
    T_b.(dVar) = D_b;

    ds_b = did.Dataset.fromTable(T_b, ...
        idVar=ds.idVar, timeVar=ds.timeVar, yVar=ds.yVar, dVar=ds.dVar, ...
        weightVar=ds.weightVar, describe=false, Print=false);

    % 5) Fit estimator on placebo data and extract statistic
    try
        res_b = did.fit(method, ds_b, nv_fit{:});
        statPlacebo(b) = extract_stat_(res_b);
    catch ME
        warning('did:random_placebo:FitFail', ...
            'Placebo replication %d failed: %s', b, ME.message);
        statPlacebo(b) = NaN;
    end

    if Verbose && (mod(b, max(1, floor(B/10))) == 0)
        fprintf('  [random_placebo] Completed %d / %d replications.\n', b, B);
    end
end

% Restore RNG
if ~isempty(rngOld)
    rng(rngOld);
end

% Drop NaN stats
ok = isfinite(statPlacebo);
statPlacebo_ok = statPlacebo(ok);

if ~any(ok)
    warning('did:random_placebo:AllNaN', ...
        'All placebo statistics are NaN; check estimator/inputs.');
end

% ---- Randomization p-values --------------------------------------------
if any(ok)
    absObs = abs(statObs);
    absPlc = abs(statPlacebo_ok);

    B_eff  = numel(statPlacebo_ok);
    pTwo   = (1 + sum(absPlc >= absObs)) / (B_eff + 1);
    pGt    = (1 + sum(statPlacebo_ok >= statObs)) / (B_eff + 1);
    pLt    = (1 + sum(statPlacebo_ok <= statObs)) / (B_eff + 1);

    qs     = quantile(statPlacebo_ok, [0.025 0.5 0.975]);
else
    pTwo = NaN; pGt = NaN; pLt = NaN;
    qs   = [NaN NaN NaN];
end

% ---- Pack output struct -----------------------------------------------
plc = struct();
plc.type            = "random_placebo";
plc.method          = string(method);
plc.B               = B;
plc.scheme          = scheme;
plc.treatClusterVar = treatClVar;
plc.statObs         = statObs;
plc.statPlacebo     = statPlacebo;
plc.okMask          = ok;
plc.pTwoSided       = pTwo;
plc.pGreater        = pGt;
plc.pLess           = pLt;
plc.quantiles       = qs;
if ~isnan(Seed_placebo) && isfinite(Seed_placebo)
    plc.Seed = Seed_placebo;
end

if Verbose
    fprintf('[random_placebo] Done. Two-sided p = %.4f (based on %d valid replications).\n', ...
        pTwo, nnz(ok));
    if pTwo<=0.1 & scheme=="permuteEverTreated";
        display("[random_placebo] Since p<=10%, amongst all assignments of the same adoption pattern across treated clusters, "...
          + newline + "the actual assignment produces by far the strongest alignment between treatment and outcomes.");
    elseif pTwo>0.1 & scheme=="permuteEverTreated";
          display("[random_placebo] Since p>10%, a random assignment of treatment times for treated clusters" ...
              + newline +"would have generated similar ATTs.");
    end
end

end

% ===================================================================
% Local helper: extract scalar statistic from estimator result
% ===================================================================
function s = extract_stat_(res)

s = did.utils.getATT(res);
end

% ===================================================================
% Local helper: min positive value, or 0 if none exist
% ===================================================================
function m = minPosOrZero(x)
x = x(x > 0);
if isempty(x)
    m = 0;
else
    m = min(x);
end
end
