function T = genDIDdata(numPeriods,numIds,percTreated, opts)
% DID.GENDIDDATA  Simulate panel data for DiD with staggered / on-off treatment,
%                 optional covariates, pre-trend drift, and cohort controls.
%
% Usage:
%   T = did.genDIDdata(numPeriods, numIds, percTreated, Name, Value, ...)
%
% Arguments:
%   numPeriods  (int)    Number of time periods (default 11)
%   numIds      (int)    Number of cross-sectional units (default 100)
%   percTreated (double) Fraction of units to be ever-treated (0..1)
%
% Name-Value Options:
%   -- Treatment Process --
%   ATT           (double) True Average Treatment Effect on Treated (default 2)
%   startPeriod   (int)    Period when treatment starts (default 5)
%   endPeriod     (int)    Period when treatment ends (for "onOff") (default 8)
%   treatType     (string) "constant" (default), "constantTime", "timeIncrease", "onOff"
%   dynEffect     (double) Slope of dynamic effect in "timeIncrease" (default 0.2)
%   cohortTimes   (vector) Specific start times for cohorts (e.g. [5, 10])
%
%   -- Errors / Fixed Effects --
%   meanError     (double) Mean of error term (default 0.5)
%   errorStd      (double) Std dev of error term (default 1.5)
%   Seed          (int)    Random seed (default random)
%
%   -- Pre-Trends --
%   preTrendType  (string) "none" (default), "unitLinear", "preOnly", "divergentControls"
%   preTrendSd    (double) Std dev of random slopes (default 0)
%   preTrendMeanTreated (double) Mean slope added to treated units (default 0)
%   preTrendMeanControl (double) Mean slope added to standard/divergent controls (default 0)
%   preTrendDivergentShare (double) Fraction of controls that are divergent (0..1, default 0.5)
%                                   Used only with "divergentControls".
%
%   -- Covariates & Selection --
%   xNum, xUnitStd, xTimeStd... (See code for details)
%   SelectionBias (double) Coeff of X on treatment probability (default 0)
%   TrendEffectX  (double) Coeff of X on linear time trend (default 0)
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 01/26/2026 (Added divergentControls)
% ------------------------------------------------------------------------

arguments
    numPeriods  (1,1) double {mustBeInteger, mustBePositive} = 11
    numIds      (1,1) double {mustBeInteger, mustBePositive} = 100
    percTreated (1,1) double {mustBeInRange(percTreated,0,1)} = 0.4

    % Treatment process
    opts.ATT (1,1) double = 2
    opts.startPeriod (1,1) double {mustBeInteger, mustBePositive} = 5
    opts.endPeriod  double {mustBeInteger, mustBePositive} = 8
    opts.treatType (1,1) string {mustBeMember(opts.treatType,["constant","constantTime","timeIncrease","onOff"])} = "constant"
    opts.dynEffect (1,1) double = 0.2

    % Errors / FE
    opts.meanError (1,1) double = 0.5
    opts.errorStd  (1,1) double = 1.5
    opts.Seed double = randi([1,1e7],1,1)

    % Covariates X(t,i) = unit FE + time FE + shock
    opts.xNum (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.xUnitStd (1,1) double {mustBeNonnegative} = 1
    opts.xTimeStd (1,1) double {mustBeNonnegative} = 1
    opts.xShockStd (1,1) double {mustBeNonnegative} = 1
    % betaX: scalar applied to all X or vector length xNum
    opts.betaX double = 0

    % Pre-trend drift
    opts.preTrendType (1,1) string {mustBeMember(opts.preTrendType,["none","unitLinear","preOnly","divergentControls"])} = "none"
    opts.preTrendSd (1,1) double {mustBeNonnegative} = 0
    opts.preTrendMeanTreated (1,1) double = 0
    opts.preTrendMeanControl (1,1) double = 0
    opts.preCenterTime double = []

    % Control share for divergentControls mode (0..1 fraction of controls that are divergent/bad)
    opts.preTrendDivergentShare (1,1) double {mustBeInRange(opts.preTrendDivergentShare,0,1)} = 0.5

    % Cohort controls (staggered only)
    opts.cohortTimes double = []   % explicit adoption times to use
    opts.numCohorts double = NaN   % if set, auto-pick this many evenly spaced times
    opts.cohortProbs double = []   % length == numel(cohortTimes) probabilities
    opts.cohortSize  double = []   % length == numel(cohortTimes) exact counts; sum == #treated
    opts.treatedNum double = NaN   % optional: force exact number treated

    % Cohort-level ATT controls (for constantTime)
    opts.CohortIncrease (1,1) double = 0.5
    opts.CohortLevels double = []

    % CEM / Selection Bias extensions
    opts.SelectionBias (1,1) double = 0   % Coefficient of X1 on Treatment Probability
    opts.TrendEffectX  (1,1) double = 0   % Coefficient of X1 on Linear Time Trend
end

% RNG
if ~isnan(opts.Seed), rng(opts.Seed); end

% Guards
if opts.startPeriod > numPeriods
    error('did:genDIDdata:startPeriod','startPeriod must be <= numPeriods.');
end
if opts.endPeriod > numPeriods
    error('did:genDIDdata:endPeriod','endPeriod must be <= numPeriods.');
end

N  = numIds * numPeriods;
Id = (1:numIds)';
Tm = (1:numPeriods)';

% Unit & time FEs
i_FE = randn(numIds,1) * opts.errorStd;
t_FE = opts.meanError/2 + randn(numPeriods,1) * (opts.errorStd/2);

% Covariates X (xNum of them)
% MOVED UP to allow Selection Bias based on X
xK = opts.xNum;
X  = [];
betaX = [];
if xK > 0
    if isscalar(opts.betaX)
        betaX = repmat(opts.betaX, 1, xK);
    else
        betaX = opts.betaX(:)'; if numel(betaX) ~= xK, error('did:genDIDdata:betaXlen','betaX must be scalar or length xNum.'); end
    end
    X = zeros(N, xK);
    % We need unit-level X for selection.
    % The legacy code generated X as (unit + time + shock).
    % For selection, we likely care about the Unit FE component of X.
    % We'll generate the components here but construct the full panel X later or now?
    % X is N x xK (Long).
    % Let's generate Unit-level X first.

    % Store unit-level components for selection
    u_i_mat = zeros(numIds, xK);

    for k = 1:xK
        u_i  = randn(numIds,1)    * opts.xUnitStd;   % unit FE component
        u_i_mat(:,k) = u_i;

        v_t  = randn(numPeriods,1)* opts.xTimeStd;   % time FE component
        e_it = randn(N,1)         * opts.xShockStd;  % idio shock

        X(:,k) = repelem(u_i,numPeriods) + repmat(v_t,numIds,1) + e_it;
    end
else
    % If SelectionBias is requested but xNum=0, forcing xNum=1 would break "output columns" expectation?
    % Better to throw error if SelectionBias != 0 and xNum == 0, OR create a latent X that isn't output?
    % The user plan said "Add argument... if SelectionBias is active...".
    if opts.SelectionBias ~= 0
        warning('did:genDIDdata:BiasWithoutX', 'SelectionBias requested but xNum=0. Bias will be based on a latent unobserved variable.');
        % Generate latent X for selection
        u_i_mat = randn(numIds, 1);
    else
        u_i_mat = [];
    end
end

% Treatment assignment at unit level (everTreated)
if ~isnan(opts.treatedNum)
    % exact count requested by user
    nTreated = round(opts.treatedNum);
    if nTreated < 0 || nTreated > numIds
        error('did:genDIDdata:treatedNum','treatedNum must be in [0, %d].', numIds);
    end
    everTreated = false(numIds,1);

    if opts.SelectionBias ~= 0 && ~isempty(u_i_mat)
        % Prob score p ~ sigmoid(bias * X)
        z = -1 + opts.SelectionBias * u_i_mat(:,1); % Use first X
        % Pick units with highest propensity? Or weighted sampling?
        % Standard practice for "Exact Number" with bias:
        % exact treated are those with highest latent index? Or weighted sample.
        % Let's do weighted sample without replacement.
        prob = 1 ./ (1 + exp(-z));
        everTreated(randsample(numIds, nTreated, true, prob)) = true;
    else
        everTreated(randperm(numIds, nTreated)) = true;
    end
else
    if opts.SelectionBias ~= 0 && ~isempty(u_i_mat)
        % Logistic selection
        % Calibrate intercept? No, user provided percTreated is ignored?
        % Or we assume the bias modifies the base rate.
        % z = logit(percTreated) + bias * X
        base_z = log(percTreated / (1 - percTreated));
        z = base_z + opts.SelectionBias * u_i_mat(:,1);
        prob = 1 ./ (1 + exp(-z));
        everTreated = rand(numIds,1) < prob;
    else
        everTreated = rand(numIds,1) < percTreated;   % Binomial draw
    end
    nTreated = sum(everTreated);
end

% If exact cohort sizes are provided, override treated count to match them
if ~isempty(opts.cohortSize)
    sizes   = round(opts.cohortSize(:)');
    nTarget = sum(sizes);
    if nTarget > numIds
        error('did:genDIDdata:cohortSizeTooBig','sum(cohortSize)=%d exceeds numIds=%d.', nTarget, numIds);
    end
    % Force exactly nTarget treated units
    everTreated = false(numIds,1);
    everTreated(randperm(numIds, nTarget)) = true;
    nTreated = nTarget;
end

% First treatment time per unit (g_unit: 0 for never)
switch opts.treatType

    case "onOff"
        % Support staggered STARTS & per-cohort ENDS if user provided cohort info;
        % otherwise fall back to legacy single-cohort behavior.

        % Determine candidate START times (cohorts)
        userWantsStaggered = (~isempty(opts.cohortTimes)) || (~isnan(opts.numCohorts) && opts.numCohorts>=1) ...
            || ~isempty(opts.cohortSize) || ~isempty(opts.cohortProbs);

        if userWantsStaggered
            if ~isempty(opts.cohortTimes)
                cand = unique(round(opts.cohortTimes(:)'));
                % User overrides startPeriod. Only bound by 1..numPeriods
                cand = cand(cand >= 1 & cand <= numPeriods);
            elseif ~isnan(opts.numCohorts) && opts.numCohorts >= 1
                K  = min(numPeriods - opts.startPeriod + 1, round(opts.numCohorts));
                cand = unique(round(linspace(opts.startPeriod, numPeriods, K)));
            else
                cand = opts.startPeriod:numPeriods;
            end
            if isempty(cand)
                error('did:genDIDdata:onOffCohorts','No valid cohort times derived for onOff.');
            end
            K = numel(cand);

            % Allocate treated units to cohorts (exact sizes or probabilistic)
            g_unit = zeros(numIds,1);  % 0 for never by default
            treatedIds = find(everTreated);

            if ~isempty(opts.cohortSize)
                sizes = round(opts.cohortSize(:)');
                if numel(sizes) ~= K
                    error('did:genDIDdata:cohortSizeLen','cohortSize length must equal number of cohorts (%d).', K);
                end
                if sum(sizes) ~= nTreated
                    error('did:genDIDdata:cohortSizeSum','sum(cohortSize) must equal number of treated units (%d).', nTreated);
                end
                bag = repelem(cand, sizes);
                bag = bag(randperm(numel(bag)));
                g_unit(treatedIds) = bag;
            else
                % probabilistic by cohortProbs (or equal)
                if isempty(opts.cohortProbs)
                    pr = ones(1,K) / K;
                else
                    pr = opts.cohortProbs(:)'; pr = pr / sum(pr);
                    if numel(pr) ~= K
                        error('did:genDIDdata:cohortProbsLen','cohortProbs length must equal number of cohorts (%d).', K);
                    end
                end
                u = rand(nTreated,1);
                edges = [0, cumsum(pr)]; edges(end) = 1;
                bins = discretize(u, edges);         % 1..K
                g_unit(treatedIds) = cand(bins);
            end

        else
            % Legacy single cohort at startPeriod
            g_unit = opts.startPeriod * double(everTreated);
            K = 1; cand = opts.startPeriod;
        end

        onOffMode = true;

        % Build cohort-specific END periods vector (length K)
        % Allow scalar endPeriod (legacy), or vector length K (per-cohort).
        if isscalar(opts.endPeriod)
            endPerCohort = repmat(round(opts.endPeriod), 1, K);
        else
            endVec = round(opts.endPeriod(:)');     % <-- keep duplicates!
            if numel(endVec) == K
                endPerCohort = endVec;
            else
                error('did:genDIDdata:onOffEndVector',...
                    'endPeriod must be scalar or length equal to number of cohorts (%d).', K);
            end
        end

        % Guard: end must be within [start,numPeriods+1] (numPeriods+1 means "never off")
        endPerCohort = min(max(endPerCohort, cand+1), numPeriods+1);

        % Map each treated unit's cohort start → its cohort end
        h_unit = (numPeriods+1) * ones(numIds,1); % default: never off
        for k = 1:K
            inCoh = (g_unit == cand(k));
            h_unit(inCoh) = endPerCohort(k);
        end
    case {"constant"}
        g_unit = opts.startPeriod * double(everTreated);
        K = 1; cand = opts.startPeriod;


    case {"constantTime","timeIncrease"}  % staggered starts allowed
        % Determine candidate cohort times
        cand = [];
        if ~isempty(opts.cohortTimes)
            cand = unique(round(opts.cohortTimes(:)'));
            % User overrides startPeriod. Only bound by 1..numPeriods
            cand = cand(cand >= 1 & cand <= numPeriods);
            if isempty(cand)
                error('did:genDIDdata:cohortTimes','cohortTimes must contain at least one integer in [startPeriod, numPeriods].');
            end
        elseif ~isnan(opts.numCohorts) && opts.numCohorts >= 2
            K  = min(numPeriods - opts.startPeriod + 1, round(opts.numCohorts));
            cand = unique(round(linspace(opts.startPeriod, numPeriods, K)));
        elseif ~isnan(opts.numCohorts) && opts.numCohorts== 1
            cand =opts.startPeriod;
        else
            cand = opts.startPeriod:numPeriods;   % default: any period in window
        end
        K = numel(cand);


        % Allocate treated units to cohorts
        g_unit = zeros(numIds,1);
        treatedIds = find(everTreated);

        if ~isempty(opts.cohortSize)
            sizes = round(opts.cohortSize(:)');
            if numel(sizes) ~= K
                error('did:genDIDdata:cohortSizeLen','cohortSize length must equal number of cohorts (%d).', K);
            end
            if sum(sizes) ~= nTreated
                error('did:genDIDdata:cohortSizeSum','sum(cohortSize) must equal number of treated units (%d).', nTreated);
            end
            bag = repelem(cand, sizes);
            bag = bag(randperm(numel(bag)));
            g_unit(treatedIds) = bag;

        else
            if isempty(opts.cohortProbs)
                pr = ones(1,K) / K;
            else
                pr = opts.cohortProbs(:)'; pr = pr / sum(pr);
                if numel(pr) ~= K
                    error('did:genDIDdata:cohortProbsLen','cohortProbs length must equal number of cohorts (%d).', K);
                end
            end
            u = rand(nTreated,1);
            edges = [0, cumsum(pr)]; edges(end) = 1;
            bins = discretize(u, edges); % 1..K
            g_unit(treatedIds) = cand(bins);
        end

        onOffMode = false;  % not used further in these modes
end

% Expand to panel
id   = repelem(Id, numPeriods);
time = repmat(Tm,  numIds, 1);
g    = repelem(g_unit, numPeriods);
unit_has_treat = repelem(everTreated, numPeriods);

% Treatment indicator (on/off)
D = (time >= g) & unit_has_treat;   % on from g
if exist('onOffMode','var') && onOffMode
    % Per-unit end time h_unit: off at/after h_unit (use numPeriods+1 for "never off")
    if ~exist('h_unit','var')
        % legacy scalar endPeriod fallback
        h_unit = repmat(min(max(round(opts.endPeriod), opts.startPeriod+1), numPeriods+1), numIds, 1);
    end
    h = repelem(h_unit, numPeriods);
    D = D & (time < h);
end
D = double(D);

% True treatment effect per observation (ATT)
ATT = zeros(N,1);
switch opts.treatType
    case {"constant"}
        % Homogeneous constant effect once treated, regardless of cohort or event time
        ATT(D==1) = opts.ATT;

    case "constantTime"
        % Constant over post-treatment time within cohort, but level differs by cohort.
        gObs   = g; gObs(~D) = NaN;
        gVals  = unique(gObs(~isnan(gObs)));
        [~, cohIdx] = ismember(gObs, gVals);  % 1..K for treated rows

        if ~isempty(opts.CohortLevels)
            if numel(opts.CohortLevels) ~= numel(gVals)
                error('did:genDIDdata:CohortLevelsLen','CohortLevels must have one entry per cohort (%d).', numel(gVals));
            end
            levelByObs = zeros(N,1);
            good = ~isnan(gObs);
            levelByObs(good) = opts.CohortLevels(cohIdx(good));
        else
            levelByObs = zeros(N,1);
            good = ~isnan(gObs);
            levelByObs(good) = opts.ATT + opts.CohortIncrease * (cohIdx(good) - 1);
        end
        ATT(D==1) = levelByObs(D==1);

    case "timeIncrease"
        % Same baseline + linear growth in event time (k = t - g)
        gObs   = g; gObs(~D) = NaN;
        good = ~isnan(gObs);
        gVals  = unique(gObs(~isnan(gObs)));
        grow = (time - g); grow(g==0) = 0;
        [~, cohIdx] = ismember(gObs, gVals);  % 1..K for treated rows
        ATT(D==1) = opts.ATT + opts.CohortIncrease * (cohIdx(good) - 1) + opts.dynEffect * grow(D==1);

    case "onOff"
        gObs   = g; gObs(~D) = NaN;
        good = ~isnan(gObs);
        gVals  = unique(gObs(~isnan(gObs)));
        grow = (time - g); grow(g==0) = 0;
        [~, cohIdx] = ismember(gObs, gVals);  % 1..K for treated rows
        ATT(D==1) = opts.ATT +opts.CohortIncrease * (cohIdx(good) - 1);


end

% Pre-trend drift (unit-specific slope)
preTerm = zeros(N,1);
slope_i = zeros(numIds,1); % Initialize latent slope

if opts.preTrendType ~= "none"
    % Base random slopes
    slope_i = randn(numIds,1) * opts.preTrendSd;

    % Add systematic difference
    if opts.preTrendMeanTreated ~= 0 || opts.preTrendMeanControl ~= 0
        if opts.preTrendType == "divergentControls"
            % Split controls into Good (Parallel) and Bad (Divergent)
            % Good Controls match Treated Trend
            % Bad Controls match Control Trend (Divergent)
            ctrlIdx = find(~everTreated);
            nC = numel(ctrlIdx);

            % Divergent Share determines Bad Controls
            nBad  = round(nC * opts.preTrendDivergentShare);
            nGood = nC - nBad;

            permC = randperm(nC);
            badC  = ctrlIdx(permC(1:nBad));
            goodC = ctrlIdx(permC(nBad+1:end));

            slope_i(everTreated) = slope_i(everTreated) + opts.preTrendMeanTreated;
            slope_i(goodC)       = slope_i(goodC)       + opts.preTrendMeanTreated; % Parallel
            slope_i(badC)        = slope_i(badC)        + opts.preTrendMeanControl; % Divergent
        else
            % Standard: All controls get MeanControl
            slope_i(everTreated)  = slope_i(everTreated)  + opts.preTrendMeanTreated;
            slope_i(~everTreated) = slope_i(~everTreated) + opts.preTrendMeanControl;
        end
    end

    slope   = repelem(slope_i, numPeriods);
    switch opts.preTrendType
        case "unitLinear"
            if isempty(opts.preCenterTime)
                t0 = mean(1:numPeriods);
            else
                t0 = opts.preCenterTime;
            end
            preTerm = slope .* (time - t0);
        case "preOnly"
            preTerm = slope .* min(0, time - g);  % only before adoption
        case "divergentControls"
            if isempty(opts.preCenterTime)
                t0 = mean(1:numPeriods);
            else
                t0 = opts.preCenterTime;
            end
            preTerm = slope .* (time - t0); % Same linear trend shape as unitLinear
    end
end

% (Covariates X generation moved up)

% Outcome components
y_base  = repelem(i_FE, numPeriods) + repmat(t_FE, numIds, 1) + preTerm;

% Trend Effect from X (confounding trend)
y_trendX = zeros(N,1);
if opts.TrendEffectX ~= 0 && ~isempty(X)
    % Effect of X1 on linear time trend
    % X is N x K. Use X(:,1).
    % Time is 1..T
    % Effect = gamma * X_i * t
    y_trendX = opts.TrendEffectX * X(:,1) .* time;
end

% Always make yX N×1
yX = zeros(N,1);
if xK > 0
    yX = X * betaX.';   % N×xK * xK×1 -> N×1
end

epsilon = opts.meanError + randn(N,1) * opts.errorStd;
y = y_base + yX + y_trendX + ATT + epsilon;

% Event time (never-treated centered at startPeriod)
eventTime = time - g;
evNever   = (g==0);
eventTime(evNever) = time(evNever) - opts.startPeriod;

% Build output table
varNames = ["id","time","y","D","g","eventTime","i_FE","t_FE","everTreated","ATT","latent_slope"];
T = table(id, time, y, D, g, eventTime, ...
    repelem(i_FE,numPeriods), repmat(t_FE,numIds,1), repelem(double(everTreated),numPeriods), ATT, ...
    repelem(slope_i,numPeriods), ...
    'VariableNames', varNames);

% Append covariates if any
if xK > 0
    for k = 1:xK
        T.("x"+k) = X(:,k);
    end
end
end
