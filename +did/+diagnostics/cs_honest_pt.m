function out = cs_honest_pt(csRes, opts)
% DID.DIAGNOSTICS.CS_HONEST_PT
%
% Wrapper: Callaway–Sant'Anna (2021) -> Honest PT (Rambachan–Roth ID step)
% using the dynamic event-study table csRes.Aggregates.es.
%
%   out = did.diagnostics.cs_honest_pt(csRes, Name=Value, ...)
%
% OPTIONS
%   DeltaType : "SD", "RM", "SDRM"   (currently only "SD" used here)
%   M         : vector of tuning parameters (e.g., 0:0.5:2) or [] to
%               use honest_pt's default fine grid.
%   L_vec     : custom weight vector over post e>=0 (else determined by Target)
%   Target    : "avgPost" (default) or "lastPost"
%   PreBand   : scalar c ≥ 0 for pre-intervals [β_pre ± c*SE_pre] (default 1.96)
%   Verbose   : logical (default true)
%
% OUTPUT
%   out : struct from honest_pt (plus cs meta)
% -------------------------------------------------------------------------

arguments
    csRes struct

    opts.DeltaType (1,1) string {mustBeMember(opts.DeltaType,["SD","RM","SDRM"])} = "SD"
    opts.M double = []

    opts.L_vec double = []
    opts.Target (1,1) string {mustBeMember(opts.Target,["avgPost","average","lastPost","firstPost"])} = "avgPost"

    opts.PreBand (1,1) double = 1.96

    opts.Verbose (1,1) logical = true
end

Verbose = opts.Verbose;

% ---- Basic checks on csRes ----------------------------------------------
if ~isfield(csRes, 'Aggregates') || ~isfield(csRes.Aggregates, 'es')
    error('did:cs_honest_pt:NoES', ...
        'csRes.Aggregates.es not found. Run cs_estimator with standard aggregation.');
end

esTab = csRes.Aggregates.es;

if ~ismember('e', esTab.Properties.VariableNames) || ...
        ~ismember('Estimate', esTab.Properties.VariableNames)
    error('did:cs_honest_pt:BadESTable', ...
        'Aggregates.es must contain variables "e" and "Estimate".');
end
if ~ismember('SE', esTab.Properties.VariableNames)
    warning('did:cs_honest_pt:NoSE', ...
        'Aggregates.es.SE not found; pre-standard-error bands cannot be used.');
end

% ---- Extract full event-study from CS (pre + post) ----------------------
e_all    = double(esTab.e(:));         % event times (can be negative)
beta_all = double(esTab.Estimate(:));  % θ̂(e)
if ismember('SE', esTab.Properties.VariableNames)
    se_all = double(esTab.SE(:));      % SE(e)
else
    se_all = NaN(size(beta_all));
end

% Split into pre (e<0) and post (e>=0)
maskPre  = (e_all < 0);
maskPost = (e_all >= 0);

if ~any(maskPost)
    error('did:cs_honest_pt:NoPost', ...
        'No post-period event times (e>=0) found in Aggregates.es.');
end

% ---- Pre-period coefficients and SEs ------------------------------------
if any(maskPre)
    hasPre   = true;
    preTimes = e_all(maskPre);
    beta_pre = beta_all(maskPre);
    se_pre   = se_all(maskPre);

    % sort by event time (ascending)
    [preTimes, idxPreSort] = sort(preTimes);
    beta_pre = beta_pre(idxPreSort);
    se_pre   = se_pre(idxPreSort);
    numPre   = numel(beta_pre);
else
    hasPre   = false;
    preTimes = -1;
    beta_pre = 0;
    se_pre   = NaN;
    numPre   = 1;
end

% ---- Post-period coefficients and SEs -----------------------------------
postTimes = e_all(maskPost);
beta_post = beta_all(maskPost);
se_post   = se_all(maskPost);

[postTimes, idxPostSort] = sort(postTimes);
beta_post = beta_post(idxPostSort);
se_post   = se_post(idxPostSort);
numPost   = numel(beta_post);

% ---- Choose actual DeltaType (downgrade if no pre) ----------------------
DeltaTypeUsed = opts.DeltaType;
if ~hasPre && DeltaTypeUsed ~= "SD"
    if Verbose
        warning('did:cs_honest_pt:NoPreForDelta', ...
            ['No pre-period event-study (e<0) found in Aggregates.es. ', ...
            'DeltaType="%s" requires pre-information, switching to DeltaType="SD".'], ...
            string(DeltaTypeUsed));
    end
    DeltaTypeUsed = "SD";
end

% ---- Weight vector over post periods (l_vec) ----------------------------
if ~isempty(opts.L_vec)
    l_vec = opts.L_vec(:);
    if numel(l_vec) ~= numPost
        error('did:cs_honest_pt:BadLvec', ...
            'L_vec must have length equal to #post periods (%d).', numPost);
    end
else
    switch lower(string(opts.Target))
        case "lastpost"
            l_vec = zeros(numPost,1);
            l_vec(end) = 1;
        case "firstpost"
            l_vec = zeros(numPost,1);
            l_vec(1) = 1;
        otherwise  % "avgPost" or "average"
            l_vec = ones(numPost,1) ./ numPost;
    end
end

% ---- Approximate SE(θ̂) for this target ---------------------------------
% θ̂ = l' β_post,  Var(θ̂) ≈ Σ l_j^2 Var(β_post_j)  (ignoring covariances)
if all(~isnan(se_post)) && all(isfinite(se_post))
    theta_se = sqrt(sum((l_vec.^2) .* (se_post.^2)));
else
    theta_se = NaN;
end

% ---- Logging ------------------------------------------------------------
if Verbose
    nPreDisplay = numPre - (~hasPre);
    fprintf('[cs_honest_pt] Using CS event-study with %d pre and %d post periods.\n', ...
        nPreDisplay, numPost);
    if hasPre
        fprintf('  Pre e in [%g, ..., %g]; Post e in [%g, ..., %g].\n', ...
            min(preTimes), max(preTimes), min(postTimes), max(postTimes));
    else
        fprintf('  No genuine pre-period e<0 in Aggregates.es.\n');
    end
    fprintf('  Target = %s, DeltaType = %s, PreBand=%.3f.\n', ...
        string(opts.Target), string(DeltaTypeUsed), opts.PreBand);
end

% ---- Call Honest PT ID routine ------------------------------------------
betahat       = [beta_pre; beta_post];
eventTimesAll = [preTimes; postTimes];

idRes = did.diagnostics.honest_pt( ...
    betahat, numPre, numPost, ...
    DeltaType = DeltaTypeUsed, ...
    M         = opts.M, ...
    l_vec     = l_vec, ...
    eventTimes= eventTimesAll, ...
    preSE     = se_pre, ...
    PreBand   = opts.PreBand, ...
    ThetaSE   = theta_se, ...
    Display   = Verbose, ...
    StoreDelta= true);

% ---- Pack output with CS meta -------------------------------------------
out = idRes;
out.baseEstimator  = "CS2021";
if isfield(csRes, 'Method')
    out.csMethod = csRes.Method;
else
    out.csMethod = "";
end
out.DeltaType_used = DeltaTypeUsed;
out.Target         = opts.Target;
out.hasPre         = hasPre;

if isfield(csRes, 'Diagnostics') && isfield(csRes.Diagnostics, 'cs')
    out.csDiagnostics = csRes.Diagnostics.cs;
else
    out.csDiagnostics = struct();
end

end
