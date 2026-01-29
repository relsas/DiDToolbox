function out = honest_pt(betahat, numPre, numPost, opts)
% HONEST_PT  Honest parallel trends (Rambachan & Roth style) for a scalar target.
%
%   out = did.diagnostics.honest_pt(betahat, numPre, numPost, Name=Value,...)
%
% INPUTS
%   betahat : K x 1 vector of event-study coefficients β(e), ordered as
%             [β_pre; β_post], with first numPre entries "pre" and last
%             numPost entries "post".
%   numPre  : number of pre-period event times (can be 0 or a dummy>0).
%   numPost : number of post-period event times (must be >= 1).
%
% NAME–VALUE OPTIONS (opts)
%   DeltaType : "SD", "RM", or "SDRM"
%               SD   = smooth-drift: |Δ^2 δ| <= M  (second differences)
%               RM   = relative magnitudes: post first diffs bounded
%                      by M times the maximum pre first diff (calibrated
%                      from pre-period intervals).
%               SDRM = intersection of SD and RM.
%
%   M         : 1 x J vector of curvature / magnitude bounds M
%               (default: 0.01:0.01:3).
%
%   l_vec     : numPost x 1 weight vector over post e>=0 such that
%               θ = sum_j l_vec(j) β_post(j). (Required, or at least used.)
%
%   eventTimes: K x 1 vector of event times e corresponding to betahat.
%               (If empty, defaults to [-numPre:-1, 0:(numPost-1)].')
%
%   preSE     : numPre x 1 vector of standard errors for pre-period β_pre(e).
%
%   PreBand   : scalar c ≥ 0 for pre-intervals [β_pre ± c * preSE].
%               If preSE missing or c<=0, we fix δ_pre = β_pre.
%
%   ThetaSE   : scalar SE(θ̂) for the target (used for 95% CI check).
%
%   Display   : logical; if true, prints diagnostics.
%
%   StoreDelta: logical; if true, stores δ paths and plots extremal paths.
%
% OUTPUT (struct)
%   out.theta_hat  : scalar θ̂ (target, e.g. avgPost).
%   out.ThetaSE    : SE(θ̂) as passed in.
%   out.ci95       : [ciLo, ciHi] usual 95% CI for θ̂.
%   out.Mvec       : grid of M values used.
%   out.theta_lb   : J x 1 vector of lower bounds θ_min(M_j).
%   out.theta_ub   : J x 1 vector of upper bounds θ_max(M_j).
%   out.M_feasible : smallest feasible M (LP succeeds).
%                    (Minimum allowed violation required to fit the data)
%   out.M_sign     : largest M where honest interval excludes 0 and has
%                    the same sign as θ̂.
%   out.delta_min  : J x K matrix, δ paths achieving θ_min(M_j) (if stored).
%   out.delta_max  : J x K matrix, δ paths achieving θ_max(M_j) (if stored).
%   out.DeltaType  : "SD", "RM" or "SDRM".
%   out.PreBand    : PreBand value used.
%   out.eventTimes : event times used.
%   out.betahat    : β(e) in same order as eventTimes.
%
% -------------------------------------------------------------------------
% Honest PT implementation (SD + RM + SDRM), last change: 2025-11-18
% -------------------------------------------------------------------------

arguments
    betahat    (:,1) double
    numPre     (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(numPre,0)}
    numPost    (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(numPost,1)}

    opts.DeltaType (1,1) string {mustBeMember(opts.DeltaType,["SD","RM","SDRM"])} = "SD"
    opts.M         (1,:) double                             = []
    opts.l_vec     (:,1) double                             = double.empty(0,1)
    opts.eventTimes(:,1) double                             = double.empty(0,1)
    opts.preSE     (:,1) double                             = double.empty(0,1)
    opts.PreBand   (1,1) double                             = 1.96
    opts.ThetaSE   (1,1) double                             = NaN
    opts.Display   (1,1) logical                            = true
    opts.StoreDelta(1,1) logical                            = true
end

Display   = opts.Display;
DeltaType = string(opts.DeltaType);
PreBand   = opts.PreBand;

% --- M grid --------------------------------------------------------------
Mvec = opts.M(:).';
if isempty(Mvec)
    % default fine grid
    Mvec = 0.01:0.01:3;
end
ThetaSE   = opts.ThetaSE;

K = numel(betahat);
if numPre + numPost ~= K
    error('honest_pt:dimMismatch', ...
        'numPre + numPost (%d) must equal length(betahat) (%d).', ...
        numPre + numPost, K);
end

% ---------- Event times --------------------------------------------------
if ~isempty(opts.eventTimes)
    e = opts.eventTimes(:);
    if numel(e) ~= K
        error('honest_pt:dimMismatch', ...
            'Length of eventTimes (%d) must match betahat (%d).', ...
            numel(e), K);
    end
else
    % Default: -numPre:-1, 0:(numPost-1)
    % Convention: e=0 is first post period (treatment time).
    e = [(-numPre:-1) 0:(numPost-1)].';
end


% ---------- L vector (target weights) ------------------------------------
if isempty(opts.l_vec)
    % Default target: average of all post-periods?
    l_vec = ones(numPost, 1) ./ numPost;
end

l_vec = opts.l_vec(:);
if numel(l_vec) ~= numPost
    error('honest_pt:dimMismatch', ...
        'Length of l_vec (%d) must match numPost (%d).', ...
        numel(l_vec), numPost);
end

% We define theta = l_vec' * beta_post.
% Vector L (size K) such that theta = L' * beta.
% pre-period weights are 0.
L = [zeros(numPre, 1); l_vec];


% ---------- Solver setup -------------------------------------------------
% We solve min/max L'(beta - delta) subject to constraints on delta.
% This is equivalent to min/max (-L' delta)  (since L'beta is constant).
%
% Variables in LP: delta (K x 1).
%
% Constraints:
% 1. SD:  | (D*delta)_t | <= M
%    D is second-difference matrix.
% 2. RM:  | (delta_{t+1} - delta_t) | <= M * max_pre | beta_{t+1} - beta_t |
%    (calibrated from pre-period).

% Build Difference Matrix D (second differences)
% D has size (K-2) x K.
% Row r corresponds to delta_{r+2} - 2*delta_{r+1} + delta_r.
if K < 3
    % Cannot do second diffs?
    D = zeros(0, K);
else
    D = spdiags(ones(K-2,1)*[1 -2 1], 0:2, K-2, K);
    D = full(D);
end

Nd = size(D, 1);

% Define A, b for LP: A*delta <= b
A = [];
b = [];

% Anchor delta to 0 at reference period (e = -1) to prevent unboundedness (level shift)
% If -1 is not in eventTimes, we might be in trouble (unbounded).
% Ideally, we should look for an option 'ReferencePeriod'. Defaulting to -1.
% Anchor delta to 0 at reference period
% Strategy: Find event time where pre-period SE is 0 (or beta is exactly 0).
refIdx = [];
if ~isempty(opts.preSE)
    % preSE corresponds to first numPre elements of e
    zeroSE = (opts.preSE == 0 | isnan(opts.preSE));
    if any(zeroSE)
        % Reference is typically close to treatment? Or first?
        % We look for the index in 1:numPre
        idxInPre = find(zeroSE);
        % If multiple zeros, which one?
        % Standard practice: The one closest to treatment? Or the base period.
        % Actually, if multiple, anchoring all of them might be valid if they are all refs.
        % But usually just one. Let's pick the last one found (closest to 0 if sorted?)
        % No, 'e' is sorted [-numPre ... -1].
        % If 'varying' base, multiple might be low SE? No, only one ref.
        refIdx = idxInPre(end);
    end
end

% Fallback: Check for exact zero in betahat (pre-period)
if isempty(refIdx)
    preBeta = betahat(1:numPre);
    zeroBeta = (abs(preBeta) < 1e-12);
    if any(zeroBeta)
        refIdx = find(zeroBeta, 1, 'last');
    end
end

% Fallback: Standard -1 if nothing else found
if isempty(refIdx)
    refIdx = find(e == -1);
end

% We implement anchor via LB/UB constraint (0 <= delta <= 0) if found.
lb = -inf(K,1);
ub =  inf(K,1);

if ~isempty(refIdx)
    lb(refIdx) = 0;
    ub(refIdx) = 0;
    if Display
        fprintf('  [honest_pt] Anchoring delta(%d) = 0 (Reference Period, e=%g)\n', refIdx, e(refIdx));
    end

    % For SD (Smoothness), we need a second anchor to prevent unbounded linear trend.
    % 2nd difference = 0 implies linear trend is allowed.
    % To identify the trend (slope), we usually pin delta at (ref-1) to 0 as well,
    % assuming bias is 0 in pre-period (parallel trends holds exactly there).
    if (DeltaType == "SD" || DeltaType == "SDRM") && (refIdx > 1)
        lb(refIdx-1) = 0;
        ub(refIdx-1) = 0;
        if Display
            fprintf('  [honest_pt] Anchoring delta(%d) = 0 (Pre-Reference, e=%g) to fix SD trend.\n', ...
                refIdx-1, e(refIdx-1));
        end
    end
end

% Equality constraints?
Aeq = [];
beq = [];


% ---------- Pre-period calibration (for RM) ------------------------------
Bmax = 0;
BmaxIdx = NaN;
Bmax_e_start = NaN;
Bmax_e_end   = NaN;

if DeltaType == "RM" || DeltaType == "SDRM"
    if numPre < 2
        error('honest_pt:NotEnoughPre','RM requires at least 2 pre-periods.');
    end

    preBeta = betahat(1:numPre);

    if isempty(opts.preSE)
        warning('honest_pt:NoPreSE','No pre-period SEs provided. Assuming 0 noise for RM calibration.');
        preSE = zeros(numPre, 1);
    else
        preSE = opts.preSE(:);
    end

    % Loop over pre-period intervals
    for r = 1:(numPre-1)
        diffVal = preBeta(r+1) - preBeta(r);
        % SE of difference?
        % Assuming independent estimation errors:
        diffSE  = sqrt(preSE(r+1)^2 + preSE(r)^2);

        % CI for the difference
        lb1 = diffVal - PreBand * diffSE;
        ub1 = diffVal + PreBand * diffSE;

        diffMin = lb1;
        diffMax = ub1;
        supAbs  = max(abs(diffMin), abs(diffMax));
        if supAbs > Bmax
            Bmax         = supAbs;
            BmaxIdx      = r;
            Bmax_e_start = e(r);
            Bmax_e_end   = e(r+1);
        end
    end

    if Display
        if ~isnan(BmaxIdx)
            fprintf(['  RM calibration: max pre |δ_{t+1}-δ_t| ≤ Bmax = %.4f ', ...
                '(attained on interval [%g, %g]).\n'], ...
                Bmax, Bmax_e_start, Bmax_e_end);
        else
            fprintf('  RM calibration: no genuine pre differences found for RM.\n');
        end
    end

    if Bmax == 0
        warning('honest_pt:RM_BmaxZero', ...
            ['RM calibration gives Bmax = 0. Pre-period first differences are ', ...
            'pinned to zero by the constraints, so the RM restriction implies ', ...
            'all post-period first differences of δ must be zero. ', ...
            'In this case, M has no effect on the RM part of the honest region.']);
    end
end


% ---------- Prepare grid & storage --------------------------------------
Mvec      = Mvec(:).';
nM        = numel(Mvec);
theta_lb  = NaN(nM,1);
theta_ub  = NaN(nM,1);
delta_min = NaN(nM,K);
delta_max = NaN(nM,K);

idxFeas         = [];   % first feasible M
lastSignSafeIdx = [];   % last M with sign-robust honest interval (no zero)

signTheta = sign(betahat(:)' * L);  % +1 or -1 (target est sign)
theta_hat = betahat(:)' * L;

lpOpts = optimoptions('linprog', ...
    'Algorithm','dual-simplex', ...
    'Display','off');

if Display
    fprintf('[honest_pt] Scanning curvature grid with %d values of M.\n', nM);
end

% ---------- Loop over M --------------------------------------------------

for mIdx = 1:nM
    M = Mvec(mIdx);

    % --- SD part: |D δ| <= M
    if (DeltaType == "SD" || DeltaType == "SDRM") && Nd > 0
        A_sd = [ D; -D ];
        b_sd = M * ones(2*Nd,1);
    else
        A_sd = [];
        b_sd = [];
    end

    % --- RM part: |(δ_{t+1}-δ_t)| <= M * Bmax for post t >= 0
    if (DeltaType == "RM" || DeltaType == "SDRM") && (Bmax > 0)
        F = zeros(K-1, K);
        for r = 1:(K-1)
            F(r, r)   = -1;
            F(r, r+1) =  1;
        end
        % We want to constrain diffs where the *destination* is post-treatment (or 0).
        % i.e., e(r+1) >= 0. This captures the delta(0)-delta(-1) jump.
        % FIXED (2025): If reference is earlier (e.g. "first"), we must constrain the PATH
        % from the reference to the post-period to prevent unboundedness.
        % So we constrain all diffs where the interval is after the reference period.

        idxMap = 1:(K-1);
        if ~isempty(refIdx)
            isRelevant = (idxMap >= refIdx); % Constrain jumps starting at or after Ref
        else
            isRelevant = e(2:K) >= 0; % Fallback (should not happen with anchor)
        end

        F_post     = F(isRelevant, :);
        A_rm       = [ F_post; -F_post ];
        b_rm       = M * Bmax * ones(2*nnz(isRelevant), 1);
    elseif (DeltaType == "RM" || DeltaType == "SDRM") && Bmax == 0
        % If Bmax==0, RM would force all post differences to zero.
        F = zeros(K-1, K);
        for r = 1:(K-1)
            F(r, r)   = -1;
            F(r, r+1) =  1;
        end

        idxMap = 1:(K-1);
        if ~isempty(refIdx)
            isRelevant = (idxMap >= refIdx);
        else
            isRelevant = e(2:K) >= 0;
        end

        F_post     = F(isRelevant, :);
        A_rm       = [ F_post; -F_post ];
        b_rm       = zeros(2*nnz(isRelevant), 1);
    else
        A_rm = [];
        b_rm = [];
    end

    % Combine all inequality constraints
    A = [A_sd; A_rm];
    b = [b_sd; b_rm];

    % 1) θ_min(M): minimize θ = L'*(β - δ) = const - L'δ
    f_minTheta = -L;
    [delta_for_min, ~, flag_min] = linprog(f_minTheta, A, b, Aeq, beq, lb, ub, lpOpts);

    % 2) θ_max(M): maximize θ = L'*(β - δ) → minimize -θ = L'δ - const
    f_maxTheta =  L;
    [delta_for_max, ~, flag_max] = linprog(f_maxTheta, A, b, Aeq, beq, lb, ub, lpOpts);

    % Check constraints
    feas_min = (flag_min > 0);
    feas_max = (flag_max > 0);

    if ~(feas_min && feas_max)
        if Display
            fprintf('[honest_pt] linprog infeasible or failed for M=%.4f (flags max=%d, min=%d).\n', ...
                M, flag_max, flag_min);
        end
        continue;  % no bounds for this M
    end

    % First feasible M
    if isempty(idxFeas)
        idxFeas = mIdx;
    end

    % Compute bounds on θ(M)
    theta_lb(mIdx) = sum( L .* (betahat - delta_for_min) ); % lower bound (θ_min)
    theta_ub(mIdx) = sum( L .* (betahat - delta_for_max) ); % upper bound (θ_max)

    if opts.StoreDelta
        delta_min(mIdx,:) = delta_for_min(:).';
        delta_max(mIdx,:) = delta_for_max(:).';
    end

    % Print line for this M
    if Display
        fprintf('  M=%.3f: θ ∈ [%.4f, %.4f]\n', M, theta_lb(mIdx), theta_ub(mIdx));
    end

    % Check if interval contains zero
    containsZero = (theta_lb(mIdx) <= 0) && (theta_ub(mIdx) >= 0);

    % Check if interval has same sign as θ̂ and does NOT contain zero
    sameSign = false;
    if signTheta > 0 && theta_lb(mIdx) > 0
        sameSign = true;
    elseif signTheta < 0 && theta_ub(mIdx) < 0
        sameSign = true;
    end

    if sameSign
        lastSignSafeIdx = mIdx;  % largest M with sign-robust honest interval so far
    end

    % If zero enters the interval, stop scanning the grid
    if containsZero
        if Display
            fprintf('[honest_pt] Zero enters honest interval at M=%.4f. Stopping grid here.\n', M);
        end
        break;
    end
end

% ---------- Summarise critical M's --------------------------------------
if ~isempty(idxFeas)
    M_feas = Mvec(idxFeas);
else
    M_feas = NaN;
end

if ~isempty(lastSignSafeIdx)
    M_sign = Mvec(lastSignSafeIdx);
else
    M_sign = NaN;
end

if Display
    if ~isnan(M_feas)
        fprintf('[honest_pt] Smallest feasible M (pre-trends compatible): M_feasible = %.4f, θ ∈ [%.4f, %.4f]\n', ...
            M_feas, theta_lb(idxFeas), theta_ub(idxFeas));
        fprintf('             (M_feasible > 0 implies data contradicts strict parallel trends/restrictions)\n');
    else
        fprintf('[honest_pt] No feasible M in grid (LP always infeasible).\n');
    end

    if ~isnan(M_sign)
        fprintf('[honest_pt] Largest M with sign(θ) preserved (honest interval excludes 0): M_sign = %.4f, θ ∈ [%.4f, %.4f]\n', ...
            M_sign, theta_lb(lastSignSafeIdx), theta_ub(lastSignSafeIdx));
    else
        fprintf('[honest_pt] Even at the smallest feasible M, the honest interval includes 0; sign is not robust on this grid.\n');
    end
end
% ---------- Build summary diagnostics table ------------------------------
if ~isnan(M_feas)
    theta_lb_feas = theta_lb(idxFeas);
    theta_ub_feas = theta_ub(idxFeas);
else
    theta_lb_feas = NaN;
    theta_ub_feas = NaN;
end

if ~isnan(M_sign)
    theta_lb_sign = theta_lb(lastSignSafeIdx);
    theta_ub_sign = theta_ub(lastSignSafeIdx);
else
    theta_lb_sign = NaN;
    theta_ub_sign = NaN;
end

if DeltaType == "RM" || DeltaType == "SDRM"
    Bmax_val   = Bmax;
    Bmax_start = Bmax_e_start;
    Bmax_end   = Bmax_e_end;
else
    Bmax_val   = NaN;
    Bmax_start = NaN;
    Bmax_end   = NaN;
end

summaryTbl = table( ...
    theta_hat, ThetaSE, [NaN NaN], [NaN NaN], ... % Placeholder for ciLo (can't calc here without alpha? No we need ThetaSE)
    M_feas, theta_lb_feas, theta_ub_feas, ...
    M_sign, theta_lb_sign, theta_ub_sign, ...
    Bmax_val, Bmax_start, Bmax_end, ...
    'VariableNames', { ...
    'ThetaHat', ...
    'ThetaSE', ...
    'CI_Lower', ... % Note: we don't compute FLCI here, just honest interval stats
    'CI_Upper', ...
    'M_feasible', ...
    'Theta_Lower_at_M_feasible', ...
    'Theta_Upper_at_M_feasible', ...
    'M_sign', ...
    'Theta_Lower_at_M_sign', ...
    'Theta_Upper_at_M_sign', ...
    'Bmax', ...
    'Bmax_e_start', ...
    'Bmax_e_end' ...
    });
% Note: ci95 is simpler:
ciLo = theta_hat - 1.96 * ThetaSE;
ciHi = theta_hat + 1.96 * ThetaSE;
summaryTbl.CI_Lower = ciLo;
summaryTbl.CI_Upper = ciHi;

summaryTbl.Properties.RowNames = {'Summary'};

if Display
    fprintf('\n[honest_pt] Summary diagnostics table:\n');
    disp(summaryTbl);
end

% ---------- Prepare output struct ---------------------------------------
out = struct();
out.theta_hat  = theta_hat;
out.ThetaSE    = ThetaSE;
out.ci95       = [ciLo, ciHi];
out.Mvec       = Mvec;
out.theta_lb   = theta_lb;
out.theta_ub   = theta_ub;
out.M_feasible = M_feas;
out.M_sign     = M_sign;
out.delta_min  = delta_min;
out.delta_max  = delta_max;
out.DeltaType  = DeltaType;
out.PreBand    = PreBand;
out.eventTimes = e;
out.betahat    = betahat;

% RM-specific diagnostics
if DeltaType == "RM" || DeltaType == "SDRM"
    out.Bmax         = Bmax;
    out.BmaxIdx      = BmaxIdx;
    out.BmaxInterval = [Bmax_e_start, Bmax_e_end];
else
    out.Bmax         = NaN;
    out.BmaxIdx      = NaN;
    out.BmaxInterval = [NaN, NaN];
end

% Summary table
out.summaryTable = summaryTbl;

% ---------- Plot extremal treatment paths for largest sign-robust M -----
if Display && opts.StoreDelta && ~isnan(M_sign) && ~isempty(lastSignSafeIdx)
    plotIdx    = lastSignSafeIdx;
    Msel       = Mvec(plotIdx);
    delta_max_ = delta_max(plotIdx,:);   % δ giving θ_max
    delta_min_ = delta_min(plotIdx,:);   % δ giving θ_min
    beta_all   = betahat(:)';

    % Implied treatment paths τ(e) = β(e) − δ(e)
    tau_min_ = beta_all - delta_min_;    % treatment path associated with θ_min
    tau_max_ = beta_all - delta_max_;    % treatment path associated with θ_max

    f  = figure('Color','w');
    ax = axes(f);
    hold(ax,'on');

    set(ax, 'Color','w', ...
        'XColor','k','YColor','k', ...
        'GridColor',[0.7 0.7 0.7], ...
        'GridAlpha',0.5);

    plot(ax, e, beta_all, '--k', 'LineWidth', 1);
    plot(ax, e, tau_max_, 'LineWidth', 1.5);
    plot(ax, e, tau_min_, 'LineWidth', 1.5);

    if any(e == 0)
        xline(ax, 0, ':', 'Pre/Post split');
    end

    xlabel(ax, 'Event time e');
    ylabel(ax, '\tau_e (treatment effect)');
    ttl = title(ax, sprintf('Extremal treatment paths for largest sign-robust M = %.3f', Msel));
    ttl.Color = 'k';

    legend(ax, {'Estimated \beta(e)', '\tau^{max}(e)', '\tau^{min}(e)'}, 'Location','Best');
    grid(ax,'on');
    hold(ax,'off');
end

end
