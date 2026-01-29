classdef Clustered < did.vcov.VcovEngine
    % Cluster-robust VCOV (one- or two-way). Decorates estimator results.
    %
    % Usage (via factory):
    %   eng = did.vcov.Clustered(clusters=["id","time"], smallSample=true);
    %   res = eng.decorate(res, ds);
    % -----------------------------------------------------------------------------------
    %  Dr. Ralf Elsas-Nicolle, LMU, Germany / Last change: 11/17/2025
    % -----------------------------------------------------------------------------------
    properties
        clusters (1,:) string = string.empty   % e.g., ["id"] or ["id","time"]
        smallSample (1,1) logical = true
        Display (1,1) logical=true
    end

    methods
        function obj = Clustered(args)
            % NV-pairs only
            arguments
                args.clusters (1,:) string = string.empty
                args.smallSample (1,1) logical = true
                args.Display (1,1) logical=true
            end
            obj.clusters    = args.clusters;
            obj.smallSample = args.smallSample;
            obj.Display = args.Display;
        end

        function res = decorate(obj, res, ds)
    % Compute clustered VCOV and attach SE/t/p.
    arguments
        obj
        res
        ds
    end

    % ----- Pull design from estimator -----
    X = []; y = []; idxD = [];
    D = [];
    if isfield(res,'Diagnostics') && isfield(res.Diagnostics,'design')
        D = res.Diagnostics.design;
        if isstruct(D) && all(isfield(D, ["X","y","idxD"]))
            X     = D.X;
            y     = D.y;
            idxD  = D.idxD;
        end
    end

    % Detect TWFE pretrend diagnostic
    methodStr   = "";
    isPretrend  = false;
    if isfield(res,'Method')
        methodStr  = string(res.Method);
        isPretrend = startsWith(methodStr, "TWFE_pretrend");
    end

    % Need at least a design; for non-pretrend we also need idxD
    if isempty(X) || isempty(y) || (~isPretrend && (isempty(idxD) || all(~isfinite(idxD))))
        % No usable design => cannot compute VCOV
        return
    end

    % ----- Drop rows with missing values in design -----
    ok = isfinite(y);
    if ~isempty(X)
        ok = ok & all(isfinite(X),2);
    end
    if ~any(ok)
        % Nothing usable
        return
    end
    ok_idx = ok;              % keep mask to apply to clusters
    X = X(ok_idx,:);
    y = y(ok_idx);

    % ----- Resolve cluster variables from dataset -----
    if ~isprop(ds,'T'); return; end
    T = ds.T;

    % Sanity check: design rows must match dataset rows
    if height(T) ~= numel(ok_idx)
        error('did:Clustered:DimensionMismatch', ...
            'Design (X,y) and dataset T must have the same number of rows.');
    end

    if isempty(obj.clusters)
        % Default to id cluster if available
        if isprop(ds,'idVar')
            c1 = T.(ds.idVar);
        else
            return
        end
        c2 = [];
        clusters_used = string(ds.idVar);
    else
        c1 = T.(obj.clusters(1));
        clusters_used = obj.clusters(1);
        c2 = [];
        if numel(obj.clusters) >= 2
            c2 = T.(obj.clusters(2));
            clusters_used = obj.clusters(1:2);
        end
    end

    % Apply same sample mask to clusters
    c1 = c1(ok_idx);
    if ~isempty(c2)
        c2 = c2(ok_idx);
    end

    % ----- OLS fit on provided design -----
    XtX    = X' * X;
    b      = XtX \ (X' * y);
    e      = y - X*b;
    XtXinv = XtX \ eye(size(X,2));

    % ----- Cluster-robust VCOV (CGM) -----
    [g1, G1] = groupIds(c1);
    if isempty(g1)
        return
    end
    S1 = meat(X, e, g1);

    if isempty(c2)
        S  = S1;
        df = max(1, G1 - 1);
    else
        [g2, G2] = groupIds(c2);
        S2  = meat(X, e, g2);
        g12 = findgroups(g1, g2);
        S12 = meat(X, e, g12);
        S   = S1 + S2 - S12;               % Cameron–Gelbach–Miller
        df  = max(1, min(G1, G2) - 1);
    end

    V = XtXinv * S * XtXinv;

    % Small-sample correction
    if obj.smallSample
        N = size(X,1); p = size(X,2);
        if isempty(c2)
            adj = (G1/(G1-1)) * ((N-1)/(N-p));
        else
            Gmin = min(G1, G2);
            adj  = (Gmin/(Gmin-1)) * ((N-1)/(N-p));
        end
        if isfinite(adj) && adj > 0
            V = adj * V;
        end
    end

    % ====== Inference for ALL coefficients ======
    se_all = sqrt(max(diag(V),0));
    t_all  = b ./ max(se_all, eps);
    p_all  = 2*tcdf(-abs(t_all), df);

    % Full coef table (ALL regressors)
    names = [];
    if isfield(D,'names'), names = D.names(:); end
    if isempty(names)
        names = "x" + (1:size(X,2));
    end
    coefTab = table(names, b, se_all, t_all, p_all, ...
        'VariableNames', {'Name','Estimate','SE','tStat','pValue'});
    res.coef = coefTab;

    % ====== ATT-focused scalars (for non-pretrend cases) ======
    if ~isPretrend
        se_att = se_all(idxD);
        beta   = b(idxD);
        if isfield(res,'ATT') && ~isempty(res.ATT)
            beta = res.ATT;    % prefer estimator's ATT if present
        end
        t_att = beta / max(se_att, eps);
        p_att = 2*tcdf(-abs(t_att), df);

        res.SE = se_att; res.t = t_att; res.p = p_att; res.df = df;
    else
        % In pretrend mode there is no ATT scalar
        res.SE = NaN; res.t = NaN; res.p = NaN; res.df = df;
    end

    % ====== summaryTable selection ======
    if isPretrend
        % For pretrend diagnostic: keep ALL pretrend coefficients
        if isfield(D,'idxGamma') && ~isempty(D.idxGamma)
            idxGamma = D.idxGamma(:)';              % positions in design X
            idxGamma = idxGamma(idxGamma >= 1 & idxGamma <= numel(b));
            idxKeep  = idxGamma;
        else
            % Fallback: keep all regressors
            idxKeep = 1:size(X,2);
        end
    else
        % Default behaviour: ATT + declared covariates
        idxKeep = idxD;
        if isfield(D,'idxCovars') && ~isempty(D.idxCovars)
            idxKeep = [idxKeep, D.idxCovars(:)'];  % ATT first, then covariates
        end
        idxKeep = idxKeep(isfinite(idxKeep) & idxKeep>0);
        idxKeep = unique(idxKeep, 'stable');
    end

    % ====== Joint Wald/F test for pretrend coefficients ======
    if isPretrend && ~isempty(D) && isfield(D,'idxGamma') && ~isempty(D.idxGamma)
        idxGamma = D.idxGamma(:);
        idxGamma = idxGamma(idxGamma >= 1 & idxGamma <= numel(b));
        idxGamma = unique(idxGamma, 'stable');

        if ~isempty(idxGamma)
            t_gamma = t_all(idxGamma);

            % Wald statistic: sum of squared t's
            stat  = sum(t_gamma.^2);
            df1   = numel(idxGamma);    % number of restrictions
            df2   = df;                 % cluster-based df from VCOV
            Fstat = stat / max(df1, eps);
            pF    = 1 - fcdf(Fstat, df1, df2);

            res.PretrendDiag = struct( ...
                'type',     "twfe_pretrend", ...
                'idxGamma', idxGamma, ...
                'tStats',   t_gamma, ...
                'Wald',     stat, ...
                'F',        Fstat, ...
                'df1',      df1, ...
                'df2',      df2, ...
                'pValue',   pF);
        end
    end

    % Trim summaryTable accordingly
    res.summaryTable = coefTab(idxKeep, :);
    
    % Add Nobs to pretrend summaryTable (number of observations used in design)
    if isPretrend
        Nobs = size(X, 1);
        res.summaryTable.Nobs = repmat(Nobs, height(res.summaryTable), 1);
    end

    % ----- Attach VCOV info -----
    info = struct();
    info.type        = "clustered";
    info.clusters    = clusters_used;
    info.smallSample = obj.smallSample;
    info.N           = size(X,1);
    info.p           = size(X,2);
    info.G1          = G1;
    if exist('G2','var'); info.G2 = G2; end
    info.df          = df;
    info.matrix      = V;

    res.Vcov = info;
    res.vcov = V;  % convenience mirror
end
end
end

% ================= file-local helpers =================

function [g, G] = groupIds(cl)
% Convert an arbitrary cluster vector to consecutive integer group ids.
if isempty(cl), g = []; G = []; return; end
if iscategorical(cl)
    [~,~,g] = unique(cl, 'stable');
else
    % robust for numeric, string, char
    [~,~,g] = unique(string(cl), 'stable');
end
G = max(g);
end

function S = meat(X, e, g)
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
