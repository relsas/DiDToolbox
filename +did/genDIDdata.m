function T = genDIDdata(numPeriods,numIds,percTreated, opts)
% DID.GENDIDDATA  Simulate panel data for DiD with staggered / on-off treatment,
%                 optional covariates, pre-trend drift, and cohort controls.
%
% Usage (legacy signature + name–values):
%   T = did.genDIDdata(numPeriods, numIds, percTreated, ...
%       ATT=2, startPeriod=5, endPeriod=8, ...
%       treatType="constant"|"constantTime"|"timeIncrease"|"onOff", dynEffect=0.2, ...
%       meanError=0.5, errorStd=1.5, Seed=NaN, ...
%       % Covariates
%       xNum=0, xUnitStd=1, xTimeStd=1, xShockStd=1, betaX=0, ...
%       % Pre-trend
%       preTrendType="none"|"unitLinear"|"preOnly", preTrendSd=0, preCenterTime=[], ...
%       % Cohort control (staggered modes only)
%       cohortTimes=[], numCohorts=NaN, cohortProbs=[], cohortSize=[], treatedNum=NaN, ...
%       % Cohort-level ATT options (constantTime only)
%       CohortIncrease=0, CohortLevels=[])
%
% Output columns:
%   id, time, y, D, g (first-treat time; 0=never), eventTime,
%   i_FE, t_FE, everTreated, ATT, x1..xK (if xNum>0)
% 
% 
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 09/30/2025
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
    opts.preTrendType (1,1) string {mustBeMember(opts.preTrendType,["none","unitLinear","preOnly"])} = "none"
    opts.preTrendSd (1,1) double {mustBeNonnegative} = 0
    opts.preCenterTime double = []

    % Cohort controls (staggered only)
    opts.cohortTimes double = []   % explicit adoption times to use
    opts.numCohorts double = NaN   % if set, auto-pick this many evenly spaced times
    opts.cohortProbs double = []   % length == numel(cohortTimes) probabilities
    opts.cohortSize  double = []   % length == numel(cohortTimes) exact counts; sum == #treated
    opts.treatedNum double = NaN   % optional: force exact number treated

    % Cohort-level ATT controls (for constantTime)
    opts.CohortIncrease (1,1) double = 0.5
    opts.CohortLevels double = []
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

% Treatment assignment at unit level (everTreated)
if ~isnan(opts.treatedNum)
    % exact count requested by user
    nTreated = round(opts.treatedNum);
    if nTreated < 0 || nTreated > numIds
        error('did:genDIDdata:treatedNum','treatedNum must be in [0, %d].', numIds);
    end
    everTreated = false(numIds,1);
    everTreated(randperm(numIds, nTreated)) = true;
else
    everTreated = rand(numIds,1) < percTreated;   % Binomial draw
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
        cand = [];
        K = 1;
        userWantsStaggered = (~isempty(opts.cohortTimes)) || (~isnan(opts.numCohorts) && opts.numCohorts>=1) ...
            || ~isempty(opts.cohortSize) || ~isempty(opts.cohortProbs);

        if userWantsStaggered
            if ~isempty(opts.cohortTimes)
                cand = unique(round(opts.cohortTimes(:)'));
                cand = cand(cand >= opts.startPeriod & cand <= numPeriods);
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
            cand = cand(cand >= opts.startPeriod & cand <= numPeriods);
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
if opts.preTrendType ~= "none" && opts.preTrendSd > 0
    slope_i = randn(numIds,1) * opts.preTrendSd;
    slope   = repelem(slope_i, numPeriods);
    switch opts.preTrendType
        case "unitLinear"
            t0 = isempty(opts.preCenterTime) * mean(1:numPeriods) + ...
                ~isempty(opts.preCenterTime) * opts.preCenterTime;
            preTerm = slope .* (time - t0);
        case "preOnly"
            preTerm = slope .* min(0, time - g);  % only before adoption
    end
end

% Covariates X (xNum of them)
xK = opts.xNum;
X  = [];
if xK > 0
    if isscalar(opts.betaX)
        betaX = repmat(opts.betaX, 1, xK);
    else
        betaX = opts.betaX(:)'; if numel(betaX) ~= xK, error('did:genDIDdata:betaXlen','betaX must be scalar or length xNum.'); end
    end
    X = zeros(N, xK);
    for k = 1:xK
        u_i  = randn(numIds,1)    * opts.xUnitStd;   % unit FE component
        v_t  = randn(numPeriods,1)* opts.xTimeStd;   % time FE component
        e_it = randn(N,1)         * opts.xShockStd;  % idio shock
        X(:,k) = repelem(u_i,numPeriods) + repmat(v_t,numIds,1) + e_it;
    end
else
    betaX = [];
end

% Outcome components
y_base  = repelem(i_FE, numPeriods) + repmat(t_FE, numIds, 1) + preTerm;

% Always make yX N×1
yX = zeros(N,1);
if xK > 0
    yX = X * betaX.';   % N×xK * xK×1 -> N×1
end

epsilon = opts.meanError + randn(N,1) * opts.errorStd;
y = y_base + yX + ATT + epsilon;

% Event time (never-treated centered at startPeriod)
eventTime = time - g;
evNever   = (g==0);
eventTime(evNever) = time(evNever) - opts.startPeriod;

% Build output table
varNames = ["id","time","y","D","g","eventTime","i_FE","t_FE","everTreated","ATT"];
T = table(id, time, y, D, g, eventTime, ...
    repelem(i_FE,numPeriods), repmat(t_FE,numIds,1), repelem(double(everTreated),numPeriods), ATT, ...
    'VariableNames', varNames);

% Append covariates if any
if xK > 0
    for k = 1:xK
        T.("x"+k) = X(:,k);
    end
end
end
