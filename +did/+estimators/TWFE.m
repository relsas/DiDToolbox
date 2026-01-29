classdef TWFE < did.estimators.Estimator
    % Two-way FE (unit + time): point estimate only; VCOV is delegated.
    %
    % Efficient within-estimator (absorption) algebraically equivalent to
    % LSDV with unit and time dummies:
    %
    %   y_it = alpha_i + lambda_t + beta * D_it + gamma' X_it + eps_it
    %
    % We use an APM-style iterative demeaning procedure (similar in spirit
    % to Stata's absorb()) to partial out unit and time fixed effects from
    % y and the regressors (Frisch–Waugh–Lovell).
    % --------------------------------------------------------------------
    % Dr. Ralf Elsas-Nicolle, LMU Germany / Last Change: 11/18/2025
    % --------------------------------------------------------------------

    properties
        Covariates    string = string.empty(1,0)
        Display       (1,1) logical = true
        details       (1,1) logical = true
        absorbMethod  (1,1) string  = "halperin"
        absorbUseGPU  (1,1) logical = false
        absorbTol     (1,1) double  = 1e-8
        absorbMaxIter (1,1) double  = 1000

        % pretrend diagnostic mode ---
        preTrend      (1,1) logical = false      % if true: run pre-trend diag instead of standard ATT
        preTrendMode  (1,1) string  = "pooled_preMinG"   % "pooled_preMinG", "pooled_notYetTreated", "cohortwise_notYetTreated", "event", "event"
        preStart      double = NaN              % first calendar time in window (optional)
        preEnd        double = NaN              % last calendar time in window (optional)

        baseline      double = NaN              % baseline time for interactions (optional)

        TimeEffects   (1,1) logical = true      % If true (default), absorb Time FE. If false, Unit FE only.
        Weights       (:,1) double = []         % Optional weights for WLS estimation
    end


    methods
        function obj = TWFE(varargin)
            for k = 1:2:numel(varargin)
                key = string(varargin{k});
                val = varargin{k+1};
                switch lower(key)
                    case 'covariates',    obj.Covariates   = string(val);
                    case 'display',       obj.Display      = logical(val);
                    case 'details',       obj.details      = logical(val);
                    case 'timeeffects',   obj.TimeEffects  = logical(val);

                    case 'pretrend',      obj.preTrend     = logical(val);
                    case 'pretrendmode',  obj.preTrendMode = string(val);
                    case 'prestart',      obj.preStart     = val;
                    case 'preend',        obj.preEnd       = val;
                    case 'baseline',      obj.baseline     = val;

                    case 'absorbmethod',  obj.absorbMethod  = string(val);
                    case 'absorbusegpu',  obj.absorbUseGPU  = logical(val);
                    case 'absorbtol',     obj.absorbTol     = val;
                    case 'absorbmaxiter', obj.absorbMaxIter = val;
                    case 'weights',       obj.Weights       = double(val);
                end
            end
        end


        function res = fit(obj, ds)
            % Trust ds: already validated by Dataset.fromTable
            T      = ds.T;
            idVar  = ds.idVar;
            tVar   = ds.timeVar;
            yVar   = ds.yVar;
            dVar   = ds.dVar;

            % Raw variables
            y   = T.(yVar);
            id  = T.(idVar);
            tt  = T.(tVar);
            D   = double(T.(dVar));

            N = numel(y);

            % Optional covariates (in levels)
            Xc = [];
            covNames = strings(0,1);
            for vn = obj.Covariates(:).'
                if ismember(vn, T.Properties.VariableNames)
                    xv = T.(vn);
                    if islogical(xv), xv = double(xv); end
                    if ~isfloat(xv),  xv = double(xv); end
                    Xc = [Xc, xv(:)]; %#ok<AGROW>
                    covNames(end+1,1) = vn; %#ok<AGROW>
                end
            end



            % ---------------------------------------------------------------------
            % BRANCH 1: standard TWFE (ATT of D with unit+time FE)
            % ---------------------------------------------------------------------
            if ~obj.preTrend
                % Two-way FE: unit and time
                id_c = categorical(id);
                tt_c = categorical(tt);
                cats_i = categories(id_c);
                cats_t = categories(tt_c);

                % Base regressors: treatment first, then covariates
                Z = [D, Xc];                   %#ok<NASGU>
                K = size(Z,2);                 %#ok<NASGU>

                % ----- Absorb unit + time FE via iterative within transform -----
                Zall = [y, D, Xc];             % N x (1+K)
                [g_i, ~] = findgroups(categorical(T.(idVar)));
                [g_t, ~] = findgroups(categorical(T.(tVar)));

                % Determine FE groups to absorb
                % Determine FE groups to absorb
                if obj.TimeEffects
                    groups = {g_i, g_t};
                    tfe_dropped = string(missing);
                else
                    groups = {g_i};
                    tfe_dropped = "Time FE Dropped (User Request)";
                end

                % Handle Weights
                w = obj.Weights;
                useWLS = ~isempty(w);
                if useWLS
                    if numel(w) ~= N
                        % If weights passed but size mismatch, try matching by ID/Time?
                        % For now, assume strict match or crash.
                        % But commonly, Weights might be passed via .Weights property correctly
                        % if user set it up.
                        if numel(w)==0
                            useWLS=false;
                        else
                            error('did:TWFE:WeightSize','Weights dimensions (%d) do not match data (%d).', numel(w), N);
                        end
                    end
                else
                    w = [];
                end

                [Zall_tilde, infoAbs] = absorbAPM(Zall, groups, obj.absorbMaxIter, obj.absorbTol, ...
                    'method',  obj.absorbMethod, ...
                    'details', obj.details && obj.Display, ...
                    'freq',    20, ...
                    'useGPU',  obj.absorbUseGPU, ...
                    'weights', w); %#ok<NASGU>

                y_tilde = Zall_tilde(:,1);
                X_tilde = Zall_tilde(:,2:end);

                if any(~isfinite(y_tilde)) || any(~isfinite(X_tilde),'all')
                    error('did:TWFE:AbsorbNaN', ...
                        'Non-finite values after within-transformation.');
                end

                % ----- OLS / WLS on within-transformed data -----
                if useWLS
                    % Weighted OLS on transformed data
                    % Note: For FWL with weights, one runs WLS on demeaned data
                    % where demeaning was also weighted.
                    b = lscov(X_tilde, y_tilde, w);

                    resid = y_tilde - X_tilde * b;
                    % SSR and SST should be weighted
                    ssr = sum(w .* (resid.^2), 'omitnan');
                    y_mean_w = sum(w .* y_tilde, 'omitnan') / sum(w,'omitnan'); % should be approx 0 if FE absorbed intercept
                    sst = sum(w .* ((y_tilde - y_mean_w).^2), 'omitnan');
                    r2  = 1 - ssr/sst;

                else
                    b     = X_tilde \ y_tilde;
                    resid = y_tilde - X_tilde * b;
                    ssr   = sum(resid.^2, 'omitnan');
                    sst   = sum((y_tilde - mean(y_tilde, 'omitnan')).^2, 'omitnan');
                    r2    = 1 - ssr/sst;
                end

                idxD = 1;                 % ATT is the first regressor (D)

                % Names: treatment, then covariates
                colnames = {};
                colnames{end+1,1} = char(dVar);
                for kk = 1:numel(covNames)
                    colnames{end+1,1} = char(covNames(kk)); %#ok<AGROW>
                end
                names = string(colnames(:));

                % Pack result (SE/t/p left NaN; VCOV engine will fill)
                res = struct();
                res.Method = "TWFE";
                res.Name   = "TWFE (unit+time FE; within / absorption, no FD)";
                res.ATT    = b(idxD);
                res.beta   = b(idxD);
                res.SE     = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

                res.coef = table(names(:), b, ...
                    'VariableNames', {'Name','Estimate'});

                res.Method = "TWFE";
                if obj.TimeEffects
                    res.Name   = "TWFE (unit+time FE; within / absorption, no FD)";
                else
                    res.Name   = "FixedEffects (Unit FE Only; Time FE Dropped)";
                end
                res.ATT    = b(idxD);
                res.beta   = b(idxD);
                res.SE     = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

                res.coef = table(names(:), b, ...
                    'VariableNames', {'Name','Estimate'});

                res.summaryTable = table("ATT(D)"', "ATT"', b(idxD), NaN, NaN, NaN, ...
                    'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue'});

                % ----- Diagnostics for VCOV engines -----
                ncov = size(Xc,2);

                res.Diagnostics = struct('design', struct( ...
                    'X',          X_tilde, ...
                    'y',          y_tilde, ...
                    'idxD',       1, ...
                    'idxCovars',  2:(1+ncov), ...
                    'idxFE_i',    [], ...
                    'idxFE_t',    [], ...
                    'names',      names, ...
                    'covarNames', covNames, ...
                    'droppedTimeFE', tfe_dropped));

                res.Details = struct( ...
                    'Nobs',         N, ...
                    'Nunits',       numel(cats_i), ...
                    'Ntime',        numel(cats_t), ...
                    'R2',           r2, ...
                    'DroppedTimeFE', tfe_dropped, ...
                    'NoIntercept',  true);

                if obj.details && obj.Display
                    fprintf('[TWFE] ATT(%s) = %.6f\n', dVar, res.ATT);
                end

                return;
            end

            % ---------------------------------------------------------------------
            % BRANCH 2: pre-trend diagnostic TWFE (various modes)
            % ---------------------------------------------------------------------

            % 1. First treatment time per unit and ever-treated indicator
            [uid, ~, g_i_all] = unique(id);
            nUnits = numel(uid);

            G_i = inf(nUnits,1);  % first treatment time per unit
            for ii = 1:nUnits
                idx_i = (g_i_all == ii);
                Di    = D(idx_i);
                ti    = tt(idx_i);
                if any(Di == 1)
                    G_i(ii) = min(ti(Di == 1));
                end
            end

            everTreat_i = isfinite(G_i);
            if ~any(everTreat_i)
                error('did:TWFE_pretrend:NoTreated', ...
                    'No treated units found (no unit ever has D=1).');
            end

            everTreat_full = everTreat_i(g_i_all);
            G_full         = G_i(g_i_all);

            minG = min(G_i(everTreat_i));   % earliest treatment time in sample

            % Normalize mode string: lowercase, no spaces
            mode = lower(strrep(obj.preTrendMode, " ", ""));
            treatedCohorts = unique(G_i(everTreat_i));
            treatedCohorts = treatedCohorts(:)';      % make it a row vector

            switch mode

                % -------------------------------------------------------------
                % MODE A: pooled_preMinG
                % Ever-treated vs never-treated, using periods strictly before
                % the earliest treatment date in the sample.
                % -------------------------------------------------------------
                case 'pooled_preming'

                    % Define pre-period window and baseline
                    if isnan(obj.preEnd)
                        if isnumeric(minG)
                            preEnd = minG - 1;
                        else
                            error('did:TWFE_pretrend:NonNumericTime', ...
                                ['Non-numeric time variable; please set preEnd explicitly ', ...
                                '(e.g., last pre-treatment period).']);
                        end
                    else
                        preEnd = obj.preEnd;
                    end

                    if isnan(obj.preStart)
                        preStart = min(tt);
                    else
                        preStart = obj.preStart;
                    end

                    if isnan(obj.baseline)
                        baseline = preStart;
                    else
                        baseline = obj.baseline;
                    end

                    % Restrict to pre-period sample (never-treated or treated-but-before-G_i)
                    inPreWindow = (tt >= preStart) & (tt <= preEnd);
                    inPreForTr  = (~everTreat_full) | (everTreat_full & (tt < G_full));
                    isPre       = inPreWindow & inPreForTr;

                    if ~any(isPre)
                        error('did:TWFE_pretrend:NoPreSample', ...
                            'No observations in pre-period window [%g,%g].', preStart, preEnd);
                    end

                    y_pre  = y(isPre);
                    id_pre = id(isPre);
                    tt_pre = tt(isPre);
                    ev_pre = everTreat_full(isPre);

                    if isempty(Xc)
                        Xc_pre = zeros(nnz(isPre), 0);
                    else
                        Xc_pre = Xc(isPre,:);
                    end

                    N_pre  = numel(y_pre);

                    tVals = unique(tt_pre);
                    if isempty(tVals)
                        error('did:TWFE_pretrend:EmptyPreTimes', ...
                            'Pre-period sample has no distinct time values.');
                    end

                    if ~any(tVals == baseline)
                        warning('did:TWFE_pretrend:BaselineNotInPre', ...
                            'Baseline time not in pre sample. Using max pre time as baseline.');
                        baseline = max(tVals);
                    end

                    tNonBase = tVals(tVals ~= baseline);
                    K_gamma  = numel(tNonBase);

                    if K_gamma == 0
                        error('did:TWFE_pretrend:NoNonBaseline', ...
                            'No non-baseline pre times found for pretrend interactions.');
                    end

                    % Build everTreat × time dummies
                    Z_gamma = zeros(N_pre, K_gamma);
                    for k = 1:K_gamma
                        Z_gamma(:,k) = ev_pre & (tt_pre == tNonBase(k));
                    end

                    ncov  = size(Xc_pre,2);
                    Z_pre = [Z_gamma, Xc_pre];   % N_pre x (K_gamma+ncov)

                    % Names (full set before dropping)
                    names = strings(K_gamma + ncov, 1);
                    for k = 1:K_gamma
                        names(k) = sprintf('isTreated x pre_t=%g', tNonBase(k));
                    end
                    for j = 1:ncov
                        names(K_gamma + j) = char(covNames(j));
                    end

                    % Absorb unit + time FE
                    [g_i_pre, ~] = findgroups(categorical(id_pre));
                    [g_t_pre, ~] = findgroups(categorical(tt_pre));

                    Zall = [y_pre, Z_pre];

                    [Zall_tilde, infoAbs] = absorbAPM(Zall, {g_i_pre, g_t_pre}, ...
                        obj.absorbMaxIter, obj.absorbTol, ...
                        'method',  obj.absorbMethod, ...
                        'details', obj.details, ...
                        'freq',    20, ...
                        'useGPU',  obj.absorbUseGPU); %#ok<NASGU>

                    y_tilde = Zall_tilde(:,1);
                    X_tilde = Zall_tilde(:,2:end);

                    if any(~isfinite(y_tilde)) || any(~isfinite(X_tilde),'all')
                        error('did:TWFE_pretrend:AbsorbNaN', ...
                            'Non-finite values after within-transformation (pretrend).');
                    end

                    % Rank-safe OLS with bookkeeping
                    [b_full, keepCols, droppedCols, rankInfo] = safeOLSPretrend( ...
                        X_tilde, y_tilde, names, obj.details, 'pooled_preMinG');

                    % Reduced design / names for diagnostics
                    X_tilde_red = X_tilde(:, keepCols);
                    names_red   = names(keepCols);

                    % Which of the kept columns are gamma vs covariates?
                    idxGamma_red  = find(keepCols <= K_gamma);
                    idxCovars_red = find(keepCols >  K_gamma);

                    gammaIdx_orig   = keepCols(idxGamma_red);
                    gammaNames_keep = names(gammaIdx_orig);
                    gamma_hat_keep  = b_full(gammaIdx_orig);

                    % Dropped gamma and covariates (for reporting)
                    droppedGammaIdx   = intersect(droppedCols, 1:K_gamma);
                    droppedGammaNames = names(droppedGammaIdx);
                    droppedCovIdx     = intersect(droppedCols, (K_gamma+1):(K_gamma+ncov));
                    droppedCovNames   = names(droppedCovIdx);

                    % Pack result
                    res = struct();
                    res.Method = "TWFE_pretrend";
                    res.Name   = "TWFE pretrend pooled (ever-treated × time, pre-min(G))";

                    res.ATT  = NaN;
                    res.beta = NaN;
                    res.SE   = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

                    % Full coefficient table (only for kept columns)
                    res.coef = table(names_red(:), b_full(keepCols), ...
                        'VariableNames', {'Name','Estimate'});

                    % SummaryTable: only gamma coefficients that survived
                    res.summaryTable = table( ...
                        gammaNames_keep, ...
                        repelem("pretrend_gamma_pooled_preMinG", numel(gammaNames_keep))', ...
                        gamma_hat_keep, ...
                        NaN(numel(gammaNames_keep),1), ...
                        NaN(numel(gammaNames_keep),1), ...
                        NaN(numel(gammaNames_keep),1), ...
                        'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue'});

                    % Diagnostics design embedded in full sample (reduced columns only)
                    N_total = numel(y);
                    K_total = size(X_tilde_red,2);

                    X_full = nan(N_total, K_total);
                    y_full = nan(N_total, 1);
                    X_full(isPre,:) = X_tilde_red;
                    y_full(isPre)   = y_tilde;

                    res.Diagnostics = struct('design', struct( ...
                        'X',          X_full, ...
                        'y',          y_full, ...
                        'idxD',       [], ...                 % no ATT here
                        'idxGamma',   idxGamma_red, ...       % positions in reduced design
                        'idxCovars',  idxCovars_red, ...
                        'idxFE_i',    [], ...
                        'idxFE_t',    [], ...
                        'names',      names_red, ...
                        'covarNames', covNames, ...
                        'droppedTimeFE', string(missing), ...
                        'sampleMask', isPre));

                    res.Details = struct( ...
                        'Mode',             "pooled_preMinG", ...
                        'Nobs',             N_pre, ...
                        'Nunits',           numel(unique(id_pre)), ...
                        'Ntime',            numel(tVals), ...
                        'PreStart',         preStart, ...
                        'PreEnd',           preEnd, ...
                        'Baseline',         baseline, ...
                        'DroppedGamma',     droppedGammaNames, ...
                        'DroppedCovariates',droppedCovNames, ...
                        'RankInfo',         rankInfo, ...
                        'DroppedTimeFE',    string(missing), ...
                        'NoIntercept',      true);

                    if obj.details
                        fprintf('[TWFE pretrend pooled_preMinG] %d pre-times (excl. baseline), N_pre=%d, Nunits_pre=%d\n', ...
                            numel(tNonBase), N_pre, res.Details.Nunits);
                    end

                    % -------------------------------------------------------------
                    % MODE B: pooled_notYetTreated
                    % Ever-treated vs never-treated, but using only not-yet-treated
                    % observations (D=0), potentially over the full horizon.
                    % -------------------------------------------------------------
                case 'pooled_notyettreated'

                    % Window and baseline
                    if isnan(obj.preEnd)
                        preEnd = max(tt);
                    else
                        preEnd = obj.preEnd;
                    end

                    if isnan(obj.preStart)
                        preStart = min(tt);
                    else
                        preStart = obj.preStart;
                    end

                    if isnan(obj.baseline)
                        baseline = preStart;
                    else
                        baseline = obj.baseline;
                    end

                    % Not-yet-treated observations in window
                    isPre = (tt >= preStart) & (tt <= preEnd) & (D == 0);

                    if ~any(isPre)
                        error('did:TWFE_pretrend:NoPreSample_notYet', ...
                            'No not-yet-treated observations in window [%g,%g].', preStart, preEnd);
                    end

                    y_pre  = y(isPre);
                    id_pre = id(isPre);
                    tt_pre = tt(isPre);
                    ev_pre = everTreat_full(isPre);

                    if isempty(Xc)
                        Xc_pre = zeros(nnz(isPre), 0);
                    else
                        Xc_pre = Xc(isPre,:);
                    end

                    N_pre  = numel(y_pre);
                    tVals  = unique(tt_pre);

                    if isempty(tVals)
                        error('did:TWFE_pretrend:EmptyPreTimes_notYet', ...
                            'Not-yet-treated sample has no distinct time values.');
                    end

                    if ~any(tVals == baseline)
                        warning('did:TWFE_pretrend:BaselineNotInPre_notYet', ...
                            'Baseline time not in sample. Using min time as baseline.');
                        baseline = min(tVals);
                    end

                    tNonBase = tVals(tVals ~= baseline);
                    K_gamma  = numel(tNonBase);

                    if K_gamma == 0
                        error('did:TWFE_pretrend:NoNonBaseline_notYet', ...
                            'No non-baseline times found for pretrend interactions.');
                    end

                    % everTreat × time (while untreated)
                    Z_gamma = zeros(N_pre, K_gamma);
                    for k = 1:K_gamma
                        Z_gamma(:,k) = ev_pre & (tt_pre == tNonBase(k));
                    end

                    ncov  = size(Xc_pre,2);
                    Z_pre = [Z_gamma, Xc_pre];

                    names = strings(K_gamma + ncov, 1);
                    for k = 1:K_gamma
                        names(k) = sprintf('isTreated x t=%g (D=0)', tNonBase(k));
                    end
                    for j = 1:ncov
                        names(K_gamma + j) = char(covNames(j));
                    end

                    [g_i_pre, ~] = findgroups(categorical(id_pre));
                    [g_t_pre, ~] = findgroups(categorical(tt_pre));

                    Zall = [y_pre, Z_pre];

                    [Zall_tilde, infoAbs] = absorbAPM(Zall, {g_i_pre, g_t_pre}, ...
                        obj.absorbMaxIter, obj.absorbTol, ...
                        'method',  obj.absorbMethod, ...
                        'details', obj.details, ...
                        'freq',    20, ...
                        'useGPU',  obj.absorbUseGPU); %#ok<NASGU>

                    y_tilde = Zall_tilde(:,1);
                    X_tilde = Zall_tilde(:,2:end);

                    if any(~isfinite(y_tilde)) || any(~isfinite(X_tilde),'all')
                        error('did:TWFE_pretrend:AbsorbNaN_notYet', ...
                            'Non-finite values after within-transformation (pretrend, not-yet-treated).');
                    end

                    % Rank-safe OLS
                    [b_full, keepCols, droppedCols, rankInfo] = safeOLSPretrend( ...
                        X_tilde, y_tilde, names, obj.details, 'pooled_notYetTreated');

                    X_tilde_red = X_tilde(:, keepCols);
                    names_red   = names(keepCols);

                    idxGamma_red  = find(keepCols <= K_gamma);
                    idxCovars_red = find(keepCols >  K_gamma);

                    gammaIdx_orig   = keepCols(idxGamma_red);
                    gammaNames_keep = names(gammaIdx_orig);
                    gamma_hat_keep  = b_full(gammaIdx_orig);

                    droppedGammaIdx   = intersect(droppedCols, 1:K_gamma);
                    droppedGammaNames = names(droppedGammaIdx);
                    droppedCovIdx     = intersect(droppedCols, (K_gamma+1):(K_gamma+ncov));
                    droppedCovNames   = names(droppedCovIdx);

                    res = struct();
                    res.Method = "TWFE_pretrend";
                    res.Name   = "TWFE pretrend pooled (ever-treated × time, not-yet-treated sample)";

                    res.ATT  = NaN;
                    res.beta = NaN;
                    res.SE   = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

                    res.coef = table(names_red(:), b_full(keepCols), ...
                        'VariableNames', {'Name','Estimate'});

                    res.summaryTable = table( ...
                        gammaNames_keep, ...
                        repelem("pretrend_gamma_pooled_notYetTreated", numel(gammaNames_keep))', ...
                        gamma_hat_keep, ...
                        NaN(numel(gammaNames_keep),1), ...
                        NaN(numel(gammaNames_keep),1), ...
                        NaN(numel(gammaNames_keep),1), ...
                        'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue'});

                    N_total = numel(y);
                    K_total = size(X_tilde_red,2);

                    X_full = nan(N_total, K_total);
                    y_full = nan(N_total, 1);
                    X_full(isPre,:) = X_tilde_red;
                    y_full(isPre)   = y_tilde;

                    res.Diagnostics = struct('design', struct( ...
                        'X',          X_full, ...
                        'y',          y_full, ...
                        'idxD',       [], ...
                        'idxGamma',   idxGamma_red, ...
                        'idxCovars',  idxCovars_red, ...
                        'idxFE_i',    [], ...
                        'idxFE_t',    [], ...
                        'names',      names_red, ...
                        'covarNames', covNames, ...
                        'droppedTimeFE', string(missing), ...
                        'sampleMask', isPre));

                    res.Details = struct( ...
                        'Mode',             "pooled_notYetTreated", ...
                        'Nobs',             N_pre, ...
                        'Nunits',           numel(unique(id_pre)), ...
                        'Ntime',            numel(tVals), ...
                        'PreStart',         preStart, ...
                        'PreEnd',           preEnd, ...
                        'Baseline',         baseline, ...
                        'DroppedGamma',     droppedGammaNames, ...
                        'DroppedCovariates',droppedCovNames, ...
                        'RankInfo',         rankInfo, ...
                        'DroppedTimeFE',    string(missing), ...
                        'NoIntercept',      true);

                    if obj.details
                        fprintf('[TWFE pretrend pooled_notYetTreated] %d times (excl. baseline), N_pre=%d, Nunits_pre=%d\n', ...
                            numel(tNonBase), N_pre, res.Details.Nunits);
                    end

                    % -------------------------------------------------------------
                    % MODE C: cohortwise_notYetTreated
                    % Cohort-specific pretrends using not-yet-treated observations.
                    % Cluster-robust SE computed inside TWFE (per cohort).
                    % -------------------------------------------------------------
                case 'cohortwise_notyettreated'

                    % Identify treated cohorts (first treatment time among ever-treated units)
                    treatedCohorts = unique(G_i(everTreat_i));
                    treatedCohorts = treatedCohorts(:)';   % row vector

                    % Extract cluster spec from attached VCOV engine (if any)
                    clusters    = string.empty;
                    smallSample = true;
                    if isprop(obj, 'VcovEngine') && ~isempty(obj.VcovEngine) ...
                            && isa(obj.VcovEngine, 'did.vcov.Clustered')
                        clusters    = obj.VcovEngine.clusters;
                        smallSample = obj.VcovEngine.smallSample;
                    end

                    Tfull = ds.T;

                    % Resolve cluster variables over full dataset
                    if isempty(clusters)
                        % Default: cluster by id
                        c1_all = Tfull.(ds.idVar);
                        c2_all = [];
                    else
                        c1_all = Tfull.(clusters(1));
                        if numel(clusters) >= 2
                            c2_all = Tfull.(clusters(2));
                        else
                            c2_all = [];
                        end
                    end

                    % Storage for all cohort×time gammas
                    allNames     = strings(0,1);
                    allEstimates = [];
                    allSE        = [];
                    allT         = [];
                    allP         = [];
                    allNobs      = [];

                    totalN      = 0;
                    cohortsUsed = [];
                    droppedList = {};   % each entry: table(Cohort, Name)

                    % Storage for per-cohort joint tests
                    cohVec    = [];
                    FVec      = [];
                    pVec      = [];
                    df1Vec    = [];
                    df2Vec    = [];
                    nGammaVec = [];
                    NobsVec   = [];

                    % -------- loop over cohorts ----------
                    for cc = 1:numel(treatedCohorts)
                        gCoh = treatedCohorts(cc);

                        % Candidate sample: times t<gCoh and D==0
                        isUntreatedPre = (tt < gCoh) & (D == 0);

                        % Cohort indicator at unit level
                        isCohortUnit = (G_i == gCoh);
                        isCohortObs  = isCohortUnit(g_i_all);

                        % Control units: never-treated or treated later than gCoh
                        isControlUnit = (~everTreat_i) | (G_i > gCoh);
                        isControlObs  = isControlUnit(g_i_all);

                        isSample_g = isUntreatedPre & (isCohortObs | isControlObs);

                        if ~any(isSample_g)
                            if obj.details
                                fprintf('[TWFE pretrend cohortwise] Cohort g=%g: no not-yet-treated pre sample; skipped.\n', gCoh);
                            end
                            continue;
                        end

                        y_g   = y(isSample_g);
                        id_g  = id(isSample_g);
                        tt_g  = tt(isSample_g);

                        if isempty(Xc)
                            Xc_g = zeros(nnz(isSample_g), 0);
                        else
                            Xc_g = Xc(isSample_g,:);
                        end

                        totalN      = totalN + numel(y_g);
                        cohortsUsed = [cohortsUsed; gCoh]; %#ok<AGROW>

                        tVals_g = unique(tt_g);
                        if isempty(tVals_g)
                            continue;
                        end

                        % Baseline for this cohort
                        if isnan(obj.baseline)
                            baseline_g = min(tVals_g);
                        else
                            baseline_g = obj.baseline;
                            if ~any(tVals_g == baseline_g)
                                if obj.details
                                    fprintf('[TWFE pretrend cohortwise] Baseline %g not in cohort g=%g pre-times. Using min pre time instead.\n', ...
                                        baseline_g, gCoh);
                                end
                                baseline_g = min(tVals_g);
                            end
                        end

                        tNonBase_g = tVals_g(tVals_g ~= baseline_g);
                        K_gamma_g  = numel(tNonBase_g);

                        if K_gamma_g == 0
                            if obj.details
                                fprintf('[TWFE pretrend cohortwise] Cohort g=%g: no non-baseline pre times; skipped.\n', gCoh);
                            end
                            continue;
                        end

                        N_g = numel(y_g);

                        % Cohort indicator at observation level (in sample)
                        isCohObs_g = isCohortObs(isSample_g);

                        % Build cohort×time dummies
                        Z_gamma_g = zeros(N_g, K_gamma_g);
                        for k = 1:K_gamma_g
                            Z_gamma_g(:,k) = isCohObs_g & (tt_g == tNonBase_g(k));
                        end

                        ncov_g  = size(Xc_g,2);
                        Z_pre_g = [Z_gamma_g, Xc_g];

                        names_g = strings(K_gamma_g + ncov_g, 1);
                        for k = 1:K_gamma_g
                            names_g(k) = sprintf('cohort=%g x pre_t=%g', gCoh, tNonBase_g(k));
                        end
                        for j = 1:ncov_g
                            names_g(K_gamma_g + j) = char(covNames(j));
                        end

                        [g_i_pre, ~] = findgroups(categorical(id_g));
                        [g_t_pre, ~] = findgroups(categorical(tt_g));

                        Zall_g = [y_g, Z_pre_g];

                        [Zall_tilde_g, infoAbs] = absorbAPM(Zall_g, {g_i_pre, g_t_pre}, ...
                            obj.absorbMaxIter, obj.absorbTol, ...
                            'method',  obj.absorbMethod, ...
                            'details', obj.details, ...
                            'freq',    20, ...
                            'useGPU',  obj.absorbUseGPU); %#ok<NASGU>

                        y_tilde_g = Zall_tilde_g(:,1);
                        X_tilde_g = Zall_tilde_g(:,2:end);

                        if any(~isfinite(y_tilde_g)) || any(~isfinite(X_tilde_g),'all')
                            error('did:TWFE_pretrend:AbsorbNaN_cohortwise', ...
                                'Non-finite values after within-transformation (cohortwise pretrend).');
                        end

                        % Handle collinearity in this cohort's design
                        [b_full_g, keepCols_g, droppedCols_g, rankInfo_g] = safeOLSPretrend( ...
                            X_tilde_g, y_tilde_g, names_g, obj.details, sprintf('cohort g=%g', gCoh)); %#ok<NASGU>

                        % Reduced design matching the identified coefficients
                        X_red_g = X_tilde_g(:, keepCols_g);
                        b_red_g = b_full_g(keepCols_g);

                        % Identify which kept columns are gamma vs covariates
                        idxGammaRed_g   = find(keepCols_g <= K_gamma_g);   % positions in reduced design
                        gammaKeepCols_g = keepCols_g(keepCols_g <= K_gamma_g);  % original positions
                        gammaNames_keep_g = names_g(gammaKeepCols_g);
                        gamma_hat_keep_g  = b_full_g(gammaKeepCols_g);

                        % -------- Cluster-robust VCOV for this cohort --------
                        % cluster vectors restricted to this cohort sample
                        c1_g = c1_all(isSample_g);
                        c2_g = [];
                        if exist('c2_all','var') && ~isempty(c2_all)
                            c2_g = c2_all(isSample_g);
                        end

                        XtX_g    = X_red_g' * X_red_g;
                        e_g      = y_tilde_g - X_red_g * b_red_g;
                        XtXinv_g = XtX_g \ eye(size(X_red_g,2));

                        [g1_g, G1_g] = localGroupIds(c1_g);
                        if isempty(g1_g)
                            % no usable clusters, skip SEs
                            se_gamma_keep_g = NaN(size(gamma_hat_keep_g));
                            t_gamma_keep_g  = NaN(size(gamma_hat_keep_g));
                            p_gamma_keep_g  = NaN(size(gamma_hat_keep_g));
                            df_g            = NaN;
                        else
                            S1_g = localMeat(X_red_g, e_g, g1_g);

                            if isempty(c2_g)
                                S_g  = S1_g;
                                df_g = max(1, G1_g - 1);
                            else
                                [g2_g, G2_g] = localGroupIds(c2_g);
                                S2_g  = localMeat(X_red_g, e_g, g2_g);
                                g12_g = findgroups(g1_g, g2_g);
                                S12_g = localMeat(X_red_g, e_g, g12_g);
                                S_g   = S1_g + S2_g - S12_g;
                                df_g  = max(1, min(G1_g, G2_g) - 1);
                            end

                            V_g = XtXinv_g * S_g * XtXinv_g;

                            % Small-sample correction as in Clustered
                            if smallSample
                                N_g_eff = size(X_red_g,1);
                                p_g     = size(X_red_g,2);
                                if isempty(c2_g)
                                    adj = (G1_g/(G1_g-1)) * ((N_g_eff-1)/(N_g_eff-p_g));
                                else
                                    Gmin_g = min(G1_g, G2_g);
                                    adj    = (Gmin_g/(Gmin_g-1)) * ((N_g_eff-1)/(N_g_eff-p_g));
                                end
                                if isfinite(adj) && adj > 0
                                    V_g = adj * V_g;
                                end
                            end

                            % Cluster-robust SEs and t-stats for gamma
                            se_full_g = sqrt(max(diag(V_g), 0));
                            se_gamma_keep_g = se_full_g(idxGammaRed_g);

                            t_gamma_keep_g = gamma_hat_keep_g ./ max(se_gamma_keep_g, eps);
                            p_gamma_keep_g = 2*tcdf(-abs(t_gamma_keep_g), df_g);
                        end

                        % --- Per-cohort joint test H0: all gamma_g(.) = 0 ---
                        % Use sum of squared t-stats (with cluster-robust SEs)
                        valid = isfinite(t_gamma_keep_g);
                        if any(valid) && df_g > 0
                            t_use = t_gamma_keep_g(valid);
                            df1_g = numel(t_use);
                            df2_g = df_g;
                            stat_g = sum(t_use.^2);
                            F_g    = stat_g / max(df1_g, eps);
                            pF_g   = 1 - fcdf(F_g, df1_g, df2_g);
                        else
                            df1_g = sum(valid);
                            df2_g = df_g;
                            F_g   = NaN;
                            pF_g  = NaN;
                        end

                        % --------------------------------------------------
                        % Aggregate across cohorts
                        allNames     = [allNames;     gammaNames_keep_g];
                        allEstimates = [allEstimates; gamma_hat_keep_g];
                        allSE        = [allSE;        se_gamma_keep_g];
                        allT         = [allT;         t_gamma_keep_g];
                        allP         = [allP;         p_gamma_keep_g];
                        allNobs      = [allNobs;      repmat(N_g, numel(gamma_hat_keep_g), 1)];

                        % Store cohortwise joint test
                        cohVec    = [cohVec;    gCoh];
                        FVec      = [FVec;      F_g];
                        pVec      = [pVec;      pF_g];
                        df1Vec    = [df1Vec;    df1_g];
                        df2Vec    = [df2Vec;    df2_g];
                        nGammaVec = [nGammaVec; numel(gamma_hat_keep_g)];
                        NobsVec   = [NobsVec;   N_g];

                        % Record dropped cohort×time interactions
                        droppedGammaIdx_g   = intersect(droppedCols_g, 1:K_gamma_g);
                        droppedGammaNames_g = names_g(droppedGammaIdx_g);

                        if ~isempty(droppedGammaNames_g)
                            droppedList{end+1,1} = table( ...
                                repmat(gCoh, numel(droppedGammaNames_g), 1), ...
                                droppedGammaNames_g, ...
                                'VariableNames', {'Cohort','Name'}); %#ok<AGROW>
                        end

                        if obj.details
                            fprintf('[TWFE pretrend cohortwise_notYetTreated] Cohort g=%g: %d pre-times, N_pre=%d\n', ...
                                gCoh, K_gamma_g, N_g);
                        end
                    end

                    if isempty(allNames)
                        error('did:TWFE_pretrend:NoCohortPretrend', ...
                            'No cohort had usable not-yet-treated pre-period observations.');
                    end

                    % Combine dropped interactions, if any
                    if ~isempty(droppedList)
                        droppedGammaTable = vertcat(droppedList{:});
                    else
                        droppedGammaTable = table([], strings(0,1), 'VariableNames', {'Cohort','Name'});
                    end

                    res = struct();
                    res.Method = "TWFE_pretrend";
                    res.Name   = "TWFE pretrend cohortwise (cohort × time, not-yet-treated sample)";

                    res.ATT  = NaN;
                    res.beta = NaN;
                    res.SE   = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

                    % For cohortwise mode, coef = name + estimate only
                    res.coef = table(allNames, allEstimates, ...
                        'VariableNames', {'Name','Estimate'});

                    % SummaryTable: with clustered SE, t, p, Nobs
                    res.summaryTable = table( ...
                        allNames, ...
                        repelem("cohort_pretrend_gamma_notYetTreated", numel(allNames))', ...
                        allEstimates, ...
                        allSE, ...
                        allT, ...
                        allP, ...
                        allNobs, ...
                        'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue','Nobs'});

                    % No unified design matrix here (cohortwise regressions),
                    % so let the external VCOV engine skip this spec.
                    res.Diagnostics = struct('design', struct( ...
                        'X',          [], ...
                        'y',          [], ...
                        'idxD',       [], ...
                        'idxGamma',   [], ...
                        'idxCovars',  [], ...
                        'idxFE_i',    [], ...
                        'idxFE_t',    [], ...
                        'names',      strings(0,1), ...
                        'covarNames', covNames, ...
                        'droppedTimeFE', string(missing), ...
                        'sampleMask', []));

                    res.Details = struct( ...
                        'Mode',                 "cohortwise_notYetTreated", ...
                        'Nobs',                 totalN, ...
                        'Ncohorts',             numel(unique(cohortsUsed)), ...
                        'CohortsUsed',          unique(cohortsUsed), ...
                        'DroppedGammaCohortwise', droppedGammaTable, ...
                        'PreStart',             NaN, ...
                        'PreEnd',               NaN, ...
                        'Baseline',             obj.baseline, ...
                        'DroppedTimeFE',        string(missing), ...
                        'NoIntercept',          true);

                    % ---- Cohortwise pretrend diagnostic summary ----
                    res.PretrendDiag = struct( ...
                        'type',     "twfe_pretrend_cohortwise", ...
                        'cohort',   cohVec, ...
                        'F',        FVec, ...
                        'pValue',   pVec, ...
                        'df1',      df1Vec, ...
                        'df2',      df2Vec, ...
                        'nGamma',   nGammaVec, ...
                        'Nobs',     NobsVec);


                case 'event'

                    % Window and baseline
                    if isnan(obj.preEnd)
                        preEnd = max(tt);
                    else
                        preEnd = obj.preEnd;
                    end

                    if isnan(obj.preStart)
                        preStart = min(tt);
                    else
                        preStart = obj.preStart;
                    end

                    % Baseline: default to k = -1 (standard reference)
                    if isnan(obj.baseline)
                        baseK = -1;
                    else
                        baseK = obj.baseline;
                    end

                    % Calculate Event Time k = t - G_i
                    k = nan(size(tt));
                    treatedIdx = isfinite(G_full);
                    k(treatedIdx) = tt(treatedIdx) - G_full(treatedIdx);

                    % Sample Restriction: D=0 (not yet treated) inside window
                    isSample = (tt >= preStart) & (tt <= preEnd) & (D == 0);

                    if ~any(isSample)
                        error('did:TWFE_pretrend:NoSample_event', ...
                            'No observations in window [%g,%g] with D=0.', preStart, preEnd);
                    end

                    y_sub  = y(isSample);
                    id_sub = id(isSample);
                    tt_sub = tt(isSample);
                    k_sub  = k(isSample);

                    if isempty(Xc)
                        Xc_sub = zeros(nnz(isSample), 0);
                    else
                        Xc_sub = Xc(isSample,:);
                    end

                    % Identify unique event times present in pre-sample for treated units
                    k_present = unique(k_sub(~isnan(k_sub)));

                    % Baseline Logic
                    % If obj.baseline == -Inf, use min(k) (user request)
                    % If obj.baseline == NaN, use -1 (standard)
                    if isnan(obj.baseline)
                        baseK = -1;
                    elseif isinf(obj.baseline) && obj.baseline < 0
                        baseK = min(k_present);
                    else
                        baseK = obj.baseline;
                    end

                    % Exclude baseline
                    k_reg = k_present(k_present ~= baseK);

                    if isempty(k_reg)
                        error('did:TWFE_pretrend:NoEventTimes', 'No non-baseline event times found in pre-sample.');
                    end

                    K_gamma = numel(k_reg);
                    N_sub   = numel(y_sub);

                    % Build Dummies: I(k == lead)
                    Z_gamma = zeros(N_sub, K_gamma);
                    for j = 1:K_gamma
                        targetK = k_reg(j);
                        hits = (k_sub == targetK);
                        Z_gamma(:, j) = hits;
                    end

                    ncov = size(Xc_sub, 2);
                    Z_pre = [Z_gamma, Xc_sub];

                    names = strings(K_gamma + ncov, 1);
                    for j = 1:K_gamma
                        names(j) = sprintf('Event k=%g', k_reg(j));
                    end
                    for j = 1:ncov
                        names(K_gamma + j) = char(covNames(j));
                    end

                    % FE Absorb
                    [g_i_sub, ~] = findgroups(categorical(id_sub));
                    [g_t_sub, ~] = findgroups(categorical(tt_sub));

                    Zall_sub = [y_sub, Z_pre];

                    [Zall_tilde, ~] = absorbAPM(Zall_sub, {g_i_sub, g_t_sub}, ...
                        obj.absorbMaxIter, obj.absorbTol, ...
                        'method',  obj.absorbMethod, ...
                        'details', obj.details, ...
                        'freq',    20, ...
                        'useGPU',  obj.absorbUseGPU); %#ok<NASGU>

                    y_tilde = Zall_tilde(:,1);
                    X_tilde = Zall_tilde(:,2:end);

                    if any(~isfinite(y_tilde)) || any(~isfinite(X_tilde),'all')
                        error('did:TWFE_pretrend:AbsorbNaN_event', ...
                            'Non-finite values after within-transformation (event mode).');
                    end

                    % Safe OLS
                    [b_full, keepCols, droppedCols, rankInfo] = safeOLSPretrend( ...
                        X_tilde, y_tilde, names, obj.details, 'event');

                    X_tilde_red = X_tilde(:, keepCols);
                    names_red   = names(keepCols);

                    idxGamma_red = find(keepCols <= K_gamma);
                    idxCov_red   = find(keepCols > K_gamma);

                    gammaIdx_orig   = keepCols(idxGamma_red);
                    gammaNames_keep = names(gammaIdx_orig);
                    gamma_hat_keep  = b_full(gammaIdx_orig);

                    droppedGammaNames = names(intersect(droppedCols, 1:K_gamma));
                    droppedCovNames   = names(intersect(droppedCols, (K_gamma+1):end));

                    % Output Construction
                    res = struct();
                    res.Method = "TWFE_pretrend";
                    res.Name   = "TWFE pretrend (Event Time Leads)";
                    res.ATT    = NaN; res.beta = NaN; res.SE = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

                    res.coef = table(names_red(:), b_full(keepCols), 'VariableNames', {'Name','Estimate'});

                    res.summaryTable = table( ...
                        gammaNames_keep, ...
                        repelem("pretrend_event", numel(gammaNames_keep))', ...
                        gamma_hat_keep, ...
                        NaN(numel(gammaNames_keep),1), NaN(numel(gammaNames_keep),1), NaN(numel(gammaNames_keep),1), ...
                        'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue'});

                    N_total = numel(y);
                    K_total = size(X_tilde_red,2);

                    X_full = nan(N_total, K_total);
                    y_full = nan(N_total, 1);
                    X_full(isSample,:) = X_tilde_red;
                    y_full(isSample)   = y_tilde;

                    res.Diagnostics = struct('design', struct( ...
                        'X',          X_full, ...
                        'y',          y_full, ...
                        'idxD',       [], ...
                        'idxGamma',   idxGamma_red, ...
                        'idxCovars',  idxCov_red, ...
                        'idxFE_i',    [], ...
                        'idxFE_t',    [], ...
                        'names',      names_red, ...
                        'covarNames', covNames, ...
                        'droppedTimeFE', string(missing), ...
                        'sampleMask', isSample));

                    res.Details = struct( ...
                        'Mode',             "event", ...
                        'Nobs',             N_sub, ...
                        'Nunits',           numel(unique(id_sub)), ...
                        'PreStart',         preStart, ...
                        'PreEnd',           preEnd, ...
                        'BaselineK',        baseK, ...
                        'DroppedGamma',     droppedGammaNames, ...
                        'DroppedCovariates',droppedCovNames, ...
                        'RankInfo',         rankInfo, ...
                        'NoIntercept',      true);

                    if obj.details
                        fprintf('[TWFE pretrend event] %d leads estimated (baseline k=%g excluded).\n', numel(gammaNames_keep), baseK);
                    end

                otherwise
                    error('did:TWFE_pretrend:BadMode', ...
                        'Unknown preTrendMode "%s". Valid options: "pooled_preMinG", "pooled_notYetTreated", "cohortwise_notYetTreated", "event".', ...
                        obj.preTrendMode);
            end
        end
    end
end


% ================== file-local helper ==================
function [Ztilde, info] = absorbAPM(Z, gList, maxIter, tol, opts)
% absorbAPM  Multi-way absorption (APM) for fixed effects.
%
%   [Ztilde, info] = absorbAPM(Z, gList, maxIter, tol, opts)
%
% REQUIRED INPUTS (positional)
%   Z       : N×K matrix of stacked variables (y, D, X…)
%   gList   : cell array {g1, g2, ..., gR}, each g_r is N×1 grouping ids
%   maxIter : maximum # iterations
%   tol     : stopping tolerance on relative change
%
% NAME–VALUE OPTIONS (opts struct)
%   opts.method   : "halperin" (default) or "cimmino"
%   opts.details  : logical (default false)
%   opts.freq     : how often to print diagnostics (default 20)
%   opts.useGPU   : logical (default false)
%
% OUTPUTS
%   Ztilde        : N×K within-transformed matrix
%   info          : convergence diagnostics

arguments
    Z double
    gList cell
    maxIter (1,1) double {mustBePositive,mustBeInteger}
    tol (1,1) double {mustBePositive}

    opts.method (1,1) string = "halperin"
    opts.details (1,1) logical = false
    opts.freq (1,1) double {mustBePositive,mustBeInteger} = 20
    opts.useGPU (1,1) logical = false
    opts.weights double = []
end

method  = lower(opts.method);
details = opts.details;
freq    = opts.freq;
useGPU  = opts.useGPU;
weights = opts.weights;

if ~(method=="halperin" || method=="cimmino")
    error("absorbAPM:BadMethod", ...
        'opts.method must be "halperin" or "cimmino".');
end

% ------------------------------------------------------------
% Setup
% ------------------------------------------------------------
Ztilde = Z;
[N,K] = size(Ztilde);
R     = numel(gList);

if N==0 || K==0 || R==0
    info = struct('converged',true,'iter',0,'relChange',0);
    return;
end

% Normalize groups to 1..G_r
gCell = cell(R,1);
G     = zeros(R,1);

for r = 1:R
    [~,~,gi] = unique(gList{r});
    gCell{r} = gi;
    G(r)     = max(gi);
end

% ------------------------------------------------------------
% GPU handling
% ------------------------------------------------------------
if useGPU
    Ztilde = gpuArray(Ztilde);
end

% ------------------------------------------------------------
% Orthogonalization: subtract global column means
% ------------------------------------------------------------
Ztilde = Ztilde - mean(Ztilde,1);

% ------------------------------------------------------------
% Precompute column indices for accumarray
% ------------------------------------------------------------
% Z(:) = [Z(:,1); Z(:,2); ...; Z(:,K)]
colIdx = repelem((1:K)', N);   % N*K × 1

% ------------------------------------------------------------
% Main APM loop
% ------------------------------------------------------------
relChange = NaN;

for it = 1:maxIter
    Zold = Ztilde;

    switch method

        % ====================================================
        % HALPERIN (sequential projections)
        % ====================================================
        case "halperin"
            for r = 1:R
                g  = gCell{r};    % N×1
                Gr = G(r);

                % expand to N*K × 1, matching Z(:)
                g_expanded = repmat(g, K, 1);

                % work on CPU for group sums
                Zcpu = gather(Ztilde);

                if isempty(weights)
                    % Unweighted (original)
                    sum_r = accumarray([g_expanded, colIdx], ...
                        Zcpu(:), [Gr, K], @sum, 0, true);
                    cnt_r  = accumarray(g, 1, [Gr,1]);
                    mean_r = sum_r ./ cnt_r;
                else
                    % Weighted
                    w_cpu = gather(weights);
                    % Expand weights to match Z(:) ? No, we need sum(Z*w).
                    % To use accumarray efficiently:
                    % Value is Zcpu(:) .* repmat(w_cpu, K, 1)

                    w_expanded = repmat(w_cpu, K, 1);

                    sum_r = accumarray([g_expanded, colIdx], ...
                        Zcpu(:) .* w_expanded, [Gr, K], @sum, 0, true);

                    sum_w_r = accumarray(g, w_cpu, [Gr,1]);
                    mean_r  = sum_r ./ sum_w_r;
                end

                % subtract means
                if useGPU
                    mean_mat = gpuArray(mean_r(g,:));
                else
                    mean_mat = mean_r(g,:);
                end
                Ztilde = Ztilde - mean_mat;
            end

            % ====================================================
            % CIMMINO (mean of projections)
            % ====================================================
        case "cimmino"
            if useGPU
                Zsum = gpuArray.zeros(N,K,"like",Ztilde);
            else
                Zsum = zeros(N,K,"like",Ztilde);
            end

            for r = 1:R
                g  = gCell{r};
                Gr = G(r);

                g_expanded = repmat(g, K, 1);
                Zcpu       = gather(Ztilde);

                if isempty(weights)
                    sum_r = accumarray([g_expanded, colIdx], ...
                        Zcpu(:), [Gr, K], @sum, 0, true);

                    cnt_r  = accumarray(g,1,[Gr,1]);
                    mean_r = sum_r ./ cnt_r;
                else
                    w_cpu = gather(weights);
                    w_expanded = repmat(w_cpu, K, 1);

                    sum_r = accumarray([g_expanded, colIdx], ...
                        Zcpu(:).*w_expanded, [Gr, K], @sum, 0, true);

                    sum_w_r = accumarray(g, w_cpu, [Gr,1]);
                    mean_r  = sum_r ./ sum_w_r;
                end

                if useGPU
                    mean_mat = gpuArray(mean_r(g,:));
                else
                    mean_mat = mean_r(g,:);
                end

                Zproj = Ztilde - mean_mat;
                Zsum  = Zsum + Zproj;
            end

            Ztilde = Zsum ./ R;
    end

    % --------------------------------------------------------
    % Convergence diagnostics
    % --------------------------------------------------------
    diffNorm = norm(Ztilde(:) - Zold(:));
    baseNorm = max(1, norm(Ztilde(:)));
    relChange = diffNorm / baseNorm;

    if details && (it==1 || mod(it,freq)==0)
        fprintf("  [absorbAPM] iter %4d: relChange = %.3e\n", it, relChange);
    end

    if relChange < tol
        if details
            fprintf("  [absorbAPM] CONVERGED at iter %d (relChange %.3e < tol %.1e)\n", ...
                it, relChange, tol);
        end
        break
    end
end

converged = relChange < tol;

if ~converged && details
    fprintf("  [absorbAPM] WARNING: maxIter=%d reached, relChange=%.3e (tol=%.1e)\n", ...
        maxIter, relChange, tol);
end

if useGPU
    Ztilde = gather(Ztilde);
end

info = struct();
info.converged = converged;
info.iter      = min(it,maxIter);
info.relChange = relChange;
info.method    = method;
info.FEs       = R;
info.maxIter   = maxIter;
info.tol       = tol;
end


function [b_full, keepCols, droppedCols, rankInfo] = safeOLSPretrend(X_tilde, y_tilde, names, details, context)
% safeOLSPretrend  Run OLS after explicitly handling rank deficiency.
%
% INPUTS
%   X_tilde  : N x K design matrix (after FE absorption)
%   y_tilde  : N x 1 outcome
%   names    : K x 1 string array with regressor names
%   details  : logical, if true: print dropped columns
%   context  : char/string label for messages (e.g. 'pooled_notYetTreated')
%
% OUTPUTS
%   b_full      : K x 1 vector of coefficients (NaN for dropped cols)
%   keepCols    : indices of columns retained in the regression
%   droppedCols : indices of columns dropped due to collinearity
%   rankInfo    : struct with rank diagnostics

[~,K] = size(X_tilde);
if K == 0
    error('did:TWFE_pretrend:NoRegressors', ...
        'No regressors in pretrend specification (%s).', context);
end

% QR with column pivoting to identify independent columns.
% Use 'vector' so p is a permutation of 1:K with X(:,p) = Q*R.
[~,~,p] = qr(X_tilde, 'vector');
r = rank(X_tilde);

% Independent columns = first r columns of X(:,p), mapped back to original indices
keepCols    = sort(p(1:r));        % original indices of independent columns
droppedCols = setdiff(1:K, keepCols);

if ~isempty(droppedCols) && details
    fprintf('[TWFE pretrend %s] Dropped %d col(s) due to collinearity:\n', ...
        context, numel(droppedCols));
    for j = droppedCols
        fprintf('   - %s\n', names(j));
    end
end

% Re-run OLS on independent columns only
X_reduced = X_tilde(:, keepCols);
b_reduced = X_reduced \ y_tilde;   % now full rank

% Map back to full K x 1, NaN for dropped columns
b_full = NaN(K,1);
b_full(keepCols) = b_reduced;

rankInfo = struct( ...
    'K',           K, ...
    'rank',        r, ...
    'keepCols',    keepCols, ...
    'droppedCols', droppedCols);
end

function [g, G] = localGroupIds(cl)
% Convert an arbitrary cluster vector to consecutive integer group ids.
if isempty(cl), g = []; G = []; return; end
if iscategorical(cl)
    [~,~,g] = unique(cl, 'stable');
else
    [~,~,g] = unique(string(cl), 'stable');
end
G = max(g);
end

function S = localMeat(X, e, g)
% Sum of (Xg' * eg) * (Xg' * eg)' over clusters
p = size(X,2);
S = zeros(p);
G = max(g);
for k = 1:G
    sel = (g==k);
    Xg  = X(sel,:);
    eg  = e(sel);
    v   = Xg' * eg;
    S   = S + (v * v.');
end
end
