function diagOut = pretrend_event(res, varargin)
% PRETREND_EVENT  Diagnostics for parallel trends using event-study ATTs.
%
%   diag = did.diagnostics.pretrend_event(res)
%   diag = did.diagnostics.pretrend_event(res, 'Name',Value,...)
%
% INPUT
%   res   : estimator result (TWFE, BJS, CS, IW, etc.), as returned by did.fit.
%
% NAME–VALUE PAIRS
%   'Alpha'        : Confidence level used in did_standardize (default: 0.95).
%   'PreCut'       : Threshold for pre-treatment event times (default: 0).
%                    Pre periods are those with e < PreCut.
%   'MinPrePeriods': Minimum number of pre periods required to run diagnostics
%                    (default: 2). If fewer, diag.nPre < MinPrePeriods and
%                    tests are set to NaN.
%   'Verbose'      : If true (default), print a one-line summary to the console.
%
% OUTPUT (struct diag)
%   diag.type           : 'event_pretrend'
%   diag.Alpha          : Alpha used
%   diag.PreCut         : PreCut used
%   diag.nPre           : number of pre-treatment event times used
%   diag.ePre           : vector of pre-treatment event times (column)
%   diag.betaPre        : ATT(e) estimates for e < PreCut
%   diag.sePre          : standard errors for ATT(e) if available, else NaN
%   diag.tPre           : t-stats betaPre ./ sePre if available, else NaN
%   diag.meanBetaPre    : mean(betaPre)
%   diag.sdBetaPre      : std(betaPre)
%   diag.maxAbsT        : max(abs(tPre))
%   diag.meanAbsT       : mean(abs(tPre))
%   diag.slope          : approximate linear pre-trend slope (WLS of betaPre on ePre)
%   diag.slope_t_approx : approximate t-stat for slope (ignoring covariance across e)
%   diag.pJoint_approx  : approximate p-value of joint test H0: betaPre = 0
%                         (chi-square using sum(tPre.^2), independence assumption)
%   diag.note           : explanatory note on approximations used
%
% IMPORTANT:
%   - This diagnostic uses ONLY the event-study ATTs and their SEs as returned
%     (after did_standardize). It DOES NOT use the full covariance matrix across
%     event times. Hence the joint test and slope t-stat are approximate and
%     may overstate significance if covariances are large.
%   - For a more formal test, use estimator-specific pre-trend diagnostics
%     (e.g. a TWFE regression with ever-treated × time dummies) in addition.

ip = inputParser;
addParameter(ip, 'Alpha',        0.95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(ip, 'PreCut',       0,    @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'MinPrePeriods',2,    @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(ip, 'Verbose',      true, @(x) islogical(x) && isscalar(x));
parse(ip, varargin{:});
P = ip.Results;

Alpha        = P.Alpha;
PreCut       = P.PreCut;
MinPre       = P.MinPrePeriods;
Verbose      = P.Verbose;

% -------------------------------------------------------------------------
% 1) Standardize estimator output to unified schema (S.es, etc.)
% -------------------------------------------------------------------------
S = did.did_standardize(res, 'Alpha', Alpha);

if ~isfield(S,'es') || isempty(S.es) || ~istable(S.es)
    warning('did.diagnostics.pretrend_event:NoEventStudy', ...
        'No event-study table (S.es) available for this estimator.');
    diagOut = make_empty_diag_(Alpha, PreCut);
    diagOut.note = sprintf('%s | No event-study info found (S.es is empty).', diagOut.note);
    return;
end

T = S.es;
% Expect at least variables: e, Estimate, and possibly SE or (LB,UB)
if ~ismember('e', T.Properties.VariableNames)
    warning('did.diagnostics.pretrend_event:MissingE', ...
        'Event-time column ''e'' not found in S.es. Cannot run pretrend diagnostics.');
    diagOut = make_empty_diag_(Alpha, PreCut);
    diagOut.note = sprintf('%s | Column ''e'' missing in S.es.', diagOut.note);
    return;
end
if ~ismember('Estimate', T.Properties.VariableNames)
    warning('did.diagnostics.pretrend_event:MissingEstimate', ...
        'Column ''Estimate'' not found in S.es. Cannot run pretrend diagnostics.');
    diagOut = make_empty_diag_(Alpha, PreCut);
    diagOut.note = sprintf('%s | Column ''Estimate'' missing in S.es.', diagOut.note);
    return;
end

% -------------------------------------------------------------------------
% 2) Subset to pre-treatment event times: e < PreCut
% -------------------------------------------------------------------------
idxPre = T.e < PreCut;
Tpre   = T(idxPre, :);
nPre   = height(Tpre);

% Initialize diag struct with defaults
diagOut = make_empty_diag_(Alpha, PreCut);
diagOut.nPre = nPre;

if nPre == 0
    % No pre periods at all
    diagOut.note = sprintf('%s | No pre-treatment event times (no e < PreCut).', diagOut.note);
    if Verbose
        fprintf('Pretrend (event): no pre-treatment event times (no e < %.3g).\n', PreCut);
    end
    return;
end

if nPre < MinPre
    % Too few pre periods for meaningful diagnostics
    diagOut.ePre        = Tpre.e(:);
    diagOut.betaPre     = Tpre.Estimate(:);
    diagOut.note        = sprintf('%s | Fewer than MinPrePeriods (%d) pre-treatment periods.', ...
        diagOut.note, MinPre);
    if Verbose
        fprintf('Pretrend (event): only %d pre periods (< MinPrePeriods = %d). Limited diagnostics.\n', ...
            nPre, MinPre);
    end
    return;
end

% -------------------------------------------------------------------------
% 3) Extract estimates and standard errors (if possible)
% -------------------------------------------------------------------------
ePre    = Tpre.e(:);
betaPre = Tpre.Estimate(:);

% Try to get SE:
sePre = NaN(size(betaPre));

if ismember('SE', Tpre.Properties.VariableNames)
    sePre = Tpre.SE(:);
elseif all(ismember({'LB','UB'}, Tpre.Properties.VariableNames))
    % If we have confidence bands, back out SE from LB/UB assuming normal approx
    crit = norminv(0.5 + Alpha/2);
    sePre = (Tpre.UB(:) - Tpre.LB(:)) / (2*crit);
end

% t-stats (approx) if SE available
tPre = NaN(size(betaPre));
hasSE = all(isfinite(sePre));
if hasSE
    % Avoid division by zero
    sePre_safe = sePre;
    sePre_safe(sePre_safe == 0) = NaN;
    tPre = betaPre ./ sePre_safe;
end

% 4) Summary statistics & approximate tests
% -------------------------------------------------------------------------
diagOut.ePre        = ePre;
diagOut.betaPre     = betaPre;
diagOut.sePre       = sePre;
diagOut.tPre        = tPre;
diagOut.meanBetaPre = mean(betaPre, 'omitnan');
diagOut.sdBetaPre   = std(betaPre,  'omitnan');

if hasSE
    diagOut.maxAbsT  = max(abs(tPre), [], 'omitnan');
    diagOut.meanAbsT = mean(abs(tPre), 'omitnan');
else
    diagOut.maxAbsT  = NaN;
    diagOut.meanAbsT = NaN;
end

% ---- Approximate joint test H0: betaPre = 0 using chi-square on t-squared ----
% Improved: Use Full Covariance if available
pJoint = NaN;
if hasSE


    if isfield(S, 'V_es') && ~isempty(S.V_es) && all(size(S.V_es) == [height(T), height(T)])
        % We have full covariance!
        % Subset to pre-periods
        V_pre = S.V_es(idxPre, idxPre);

        % Drop NaNs (e.g. baseline periods with 0 variance)
        validIdx = isfinite(betaPre) & (diag(V_pre) > 1e-12);

        if any(validIdx)
            b_valid = betaPre(validIdx);
            V_valid = V_pre(validIdx, validIdx);

            % Wald Statistic: b' V^{-1} b
            % Use pinv for rank-deficient matrices
            [U, Sig, ~] = svd(V_valid);
            tol   = max(size(V_valid)) * eps(max(diag(Sig)));
            rankV = sum(diag(Sig) > tol);

            invV = pinv(V_valid, tol);
            chi2 = b_valid' * invV * b_valid;

            if exist('chi2cdf','file') == 2
                pJoint = 1 - chi2cdf(chi2, rankV);
            end
            diagOut.note = sprintf('%s | Uses Full Covariance Matrix for Joint Test.', diagOut.note);
        else
            pJoint = NaN;
        end
    else
        % Fallback: Independence Assumption
        k    = numel(betaPre);
        chi2 = sum(tPre.^2, 'omitnan');
        if exist('chi2cdf','file') == 2
            pJoint = 1 - chi2cdf(chi2, k);
        else
            % Statistics Toolbox not available
            pJoint = NaN;
        end
        diagOut.note = sprintf(['%s | Uses event-study ATTs only (no full covariance across e). ', ...
            'Joint test and slope t-stat are approximate (independence assumption).'], ...
            diagOut.note);
    end
end
diagOut.pJoint_approx = pJoint;

% ---- Approximate linear pre-trend slope via WLS ----
slope      = NaN;
slope_t    = NaN;

if hasSE && nPre >= 2
    % Weighted least squares of betaPre on [1, ePre], weights = 1./sePre^2
    X   = [ones(nPre,1), ePre];
    w   = 1 ./ (sePre.^2);
    w(~isfinite(w) | w <= 0) = NaN;
    W   = diag(w);
    % Drop any rows with NaN weights
    valid = all(isfinite(X),2) & isfinite(betaPre) & isfinite(diag(W));
    if nnz(valid) >= 2
        Xv = X(valid,:);
        yv = betaPre(valid);
        wv = w(valid);
        Wv = diag(wv);
        % WLS coefficient: (X'WX)^{-1} X'W y
        XtW = Xv' * Wv;
        betaHat = (XtW * Xv) \ (XtW * yv);
        slope   = betaHat(2);

        % Approximate variance of betaHat ignoring covariance across betas
        % Var(betaHat) ~ (X'WX)^{-1}
        Vbeta   = inv(XtW * Xv);
        seSlope = sqrt(Vbeta(2,2));
        slope_t = slope / seSlope;
    end
end

diagOut.slope          = slope;
diagOut.slope_t_approx = slope_t;

% -------------------------------------------------------------------------
% 5) Note on approximations & optional verbose output
% -------------------------------------------------------------------------
% (Logic moved to respective blocks)

if Verbose
    fprintf('Pretrend (event): n_pre=%d, mean(ATT_pre)=%.4g, ', nPre, diagOut.meanBetaPre);
    if hasSE
        fprintf('max|t_pre|=%.2f, approx p_joint=%.3f\n', diagOut.maxAbsT, diagOut.pJoint_approx);
    else
        fprintf('SEs not available, only point summaries.\n');
    end
end

end

% =========================================================================
% Helper: empty diagnostic skeleton
% =========================================================================
function diag = make_empty_diag_(Alpha, PreCut)
diag = struct( ...
    'type',           'event_pretrend', ...
    'Alpha',          Alpha, ...
    'PreCut',         PreCut, ...
    'nPre',           0, ...
    'ePre',           [], ...
    'betaPre',        [], ...
    'sePre',          [], ...
    'tPre',           [], ...
    'meanBetaPre',    NaN, ...
    'sdBetaPre',      NaN, ...
    'maxAbsT',        NaN, ...
    'meanAbsT',       NaN, ...
    'slope',          NaN, ...
    'slope_t_approx', NaN, ...
    'pJoint_approx',  NaN, ...
    'note',           'Event-study-based pretrend diagnostic.' ...
    );
end
