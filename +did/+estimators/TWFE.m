classdef TWFE < did.estimators.Estimator
    % Two-way FE (unit + time): point estimate only; VCOV is delegated.

    properties
        Covariates string = string.empty(1,0)
        Verbose    (1,1) logical = true
    end

    methods
        function obj = TWFE(varargin)
            for k = 1:2:numel(varargin)
                key = string(varargin{k}); val = varargin{k+1};
                switch lower(key)
                    case 'covariates', obj.Covariates = string(val);
                    case 'verbose',    obj.Verbose    = logical(val);
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

            y  = T.(yVar);
            id = categorical(T.(idVar));
            tt = categorical(T.(tVar));
            D  = double(T.(dVar));

            % Optional covariates
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

            % FE design
            cats_i = categories(id);
            cats_t = categories(tt);
            Ii = dummyvar(id);            % all unit FEs
            Ft = dummyvar(tt);            % time FEs
            droppedTime = string(missing);
            if size(Ft,2) >= 1
                droppedTime = string(cats_t{1});
                Ft(:,1) = [];             % drop one time FE
            end

            % Build X and names in lockstep
            X = []; colnames = {};

            X = [X, D(:)];                      colnames{end+1,1} = char(dVar);     % treatment first
            if ~isempty(Xc)
                X = [X, Xc];
                for k = 1:numel(covNames), colnames{end+1,1} = char(covNames(k)); end
            end
            if ~isempty(Ii)
                X = [X, Ii];
                for k = 1:numel(cats_i), colnames{end+1,1} = sprintf('FE_i_%s', char(cats_i{k})); end
            end
            if ~isempty(Ft)
                X = [X, Ft];
                for k = 2:numel(cats_t), colnames{end+1,1} = sprintf('FE_t_%s', char(cats_t{k})); end
            end

            names = string(colnames(:));
            if size(X,2) ~= numel(names)
                error('did:TWFE:NameMismatch','Internal: names (%d) != X columns (%d).', numel(names), size(X,2));
            end

            % Drop rows with any missing
            ok = isfinite(y); if ~isempty(X), ok = ok & all(isfinite(X),2); end
            y = y(ok); X = X(ok,:);

            % OLS point estimates only
            b    = X \ y;
            idxD = 1;                 % ATT is the first regressor (D)

            % Pack result (SE/t/p left NaN; VCOV engine will fill)
            res = struct();
            res.Method = "TWFE";
            res.Name   = "TWFE (unit+time FE; no intercept; drop 1 time FE)";
            res.ATT    = b(idxD);
            res.beta   = b(idxD);
            res.SE     = NaN; res.t = NaN; res.p = NaN; res.df = NaN;

            res.coef         = table(names(:), b, 'VariableNames', {'Name','Estimate'});
            res.summaryTable = table("ATT(D)"', "ATT"', b(idxD), NaN, NaN, NaN, ...
                'VariableNames', {'Name','Effect','Estimate','SE','tStat','pValue'});

            % Expose design for VCOV engines
            ncov = size(Xc,2); nI = size(Ii,2); nT = size(Ft,2);
            res.Diagnostics = struct('design', struct( ...
                'X',X,'y',y, ...
                'idxD',1, ...
                'idxCovars', 2:(1+ncov), ...
                'idxFE_i',    (2+ncov):(1+ncov+nI), ...
                'idxFE_t',    (2+ncov+nI):(1+ncov+nI+nT), ...
                'names',names, ...
                'covarNames',covNames, ...
                'droppedTimeFE', droppedTime));

            res.Details = struct('Nobs',numel(y), ...
                                 'Nunits',numel(cats_i), ...
                                 'Ntime', numel(cats_t), ...
                                 'DroppedTimeFE', droppedTime, ...
                                 'NoIntercept', true);

            if obj.Verbose
                fprintf('[TWFE] ATT(%s) = %.6f\n', dVar, res.ATT);
            end
        end
    end
end
