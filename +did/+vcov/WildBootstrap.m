classdef WildBootstrap < did.vcov.VcovEngine
    % One-way wild cluster bootstrap for the main effect (DiD/TWFE).
    % Uses estimator's design (X, y, idxD). No full covariance matrix is produced.

    properties
        cluster (1,:) string = string.empty
        B (1,1) double {mustBeInteger, mustBePositive} = 999
        multiplier (1,1) string = "mammen"
        studentize (1,1) logical = true
        rngSeed double = NaN
        smallSample (1,1) logical = true
    end

    methods
        function obj = WildBootstrap(args)
            arguments
                args.cluster (1,:) string = string.empty
                args.B (1,1) double {mustBeInteger, mustBePositive} = 999
                args.multiplier (1,1) string = "mammen"
                args.studentize (1,1) logical = true
                args.rngSeed double = NaN
                args.smallSample (1,1) logical = true
            end
            obj.cluster     = args.cluster;
            obj.B           = args.B;
            obj.multiplier  = lower(args.multiplier);
            obj.studentize  = args.studentize;
            obj.rngSeed     = args.rngSeed;
            obj.smallSample = args.smallSample;
        end

        function res = decorate(obj, res, ds)
            arguments
                obj
                res
                ds
            end

            % ---- pull design ----
            if ~isfield(res,'Diagnostics') || ~isfield(res.Diagnostics,'design'), return; end
            D = res.Diagnostics.design;
            if ~isstruct(D) || ~all(isfield(D,["X","y","idxD"])), return; end

            X = D.X; y = D.y; idxD = D.idxD;
            if isempty(X) || isempty(y) || isempty(idxD), return; end
            [N, p] = size(X);

            % ---- cluster vector ----
            if ~isprop(ds,'T'), return; end
            Ttbl = ds.T;
            if isempty(obj.cluster)
                if ~isprop(ds,'idVar'), return; end
                c = Ttbl.(ds.idVar);
            else
                if ~ismember(obj.cluster, string(Ttbl.Properties.VariableNames)), return; end
                c = Ttbl.(obj.cluster);
            end
            [~,~,g] = unique(c,'stable');
            G = max(g);

            % ---- original OLS & clustered SE for t0 ----
            XtX    = X' * X;
            XtXinv = XtX \ eye(p);
            bhat   = XtX \ (X' * y);
            ehat   = y - X*bhat;

            Vclust = did.vcov.WildBootstrap.clusterVar_(X, ehat, XtXinv, g, obj.smallSample);
            se_hat = sqrt(max(Vclust(idxD,idxD), 0));
            if ~(isfinite(se_hat) && se_hat > eps)
                res.SE = NaN; res.t = NaN; res.p = NaN; res.df = NaN;
                info = struct('type',"wild",'note',"SE too small or non-finite to studentize");
                res.Vcov = info;
                return
            end
            t0 = bhat(idxD) / se_hat;

            % ---- RNG ----
            if ~isnan(obj.rngSeed), rng(obj.rngSeed); end

            % ---- precompute cluster index lists for speed ----
            I = accumarray(g, (1:N)', [G 1], @(v){v});

            % ---- bootstrap loop ----
            tstar = zeros(obj.B,1);
            bstar = zeros(obj.B,1);
            bstar_all = zeros(obj.B, p);   % optional per-coeff SEs

            for b = 1:obj.B
                w = did.vcov.WildBootstrap.drawMultipliers_(obj.multiplier, G);

                % e* = w_g * ehat within each cluster
                e_star = zeros(N,1);
                for gg = 1:G
                    idx = I{gg};
                    e_star(idx) = ehat(idx) * w(gg);
                end
                y_star = X*bhat + e_star;

                % Refit on the same X
                b_b   = XtX \ (X' * y_star);
                bstar_all(b,:) = b_b.';      % track all coefficients
                e_b   = y_star - X*b_b;

                % Clustered SE in bootstrap sample (same g)
                V_b  = did.vcov.WildBootstrap.clusterVar_(X, e_b, XtXinv, g, obj.smallSample);
                se_b = sqrt(max(V_b(idxD,idxD), 0));

                bstar(b) = b_b(idxD);
                if obj.studentize && se_b > eps
                    tstar(b) = (b_b(idxD) - bhat(idxD)) / se_b;
                else
                    tstar(b) = (b_b(idxD) - bhat(idxD));
                end
            end

            % ---- inference ----
            if obj.studentize
                pBoot = mean(abs(tstar) >= abs(t0));
                seBoot = std(bstar - bhat(idxD), 0, 1);
            else
                pBoot = mean(abs(tstar)/se_hat >= abs(t0));
                seBoot = std(bstar - bhat(idxD), 0, 1);
            end

            % Optional: per-coefficient SEs (handy for covariates)
            seBoot_all = std(bstar_all - bhat.', 0, 1);          % 1×p
            res.vcov = diag(seBoot_all.^2);    % <- diagonal vcov for ALL coefficients
            vcovExtras = struct('se_all', seBoot_all(:));
            if isfield(D,'names') && numel(D.names)==p
                vcovExtras.se_by_name = table(D.names(:), seBoot_all(:), ...
                    'VariableNames', {'Name','SE_boot'});
            end

            % ---- Output ----
            res.SE = seBoot;
            res.t  = t0;
            res.p  = pBoot;
            res.df = NaN;

            info = struct();
            info.type        = "wild";
            if isempty(obj.cluster)
                info.cluster = "<ds.idVar>";
            else
                info.cluster = obj.cluster;
            end
            info.B           = obj.B;
            info.multiplier  = obj.multiplier;
            info.studentize  = obj.studentize;
            info.smallSample = obj.smallSample;
            info.t0          = t0;
            info.pBoot       = pBoot;
            info.seBoot      = seBoot;
            info.G           = G;

            % Merge extras
            f = fieldnames(vcovExtras);
            for k = 1:numel(f), info.(f{k}) = vcovExtras.(f{k}); end
            res.Vcov = info;
        end
    end

    methods (Static, Access=private)
        function w = drawMultipliers_(method, G)
            switch method
                case "mammen"
                    % Two-point Mammen distribution
                    p = (sqrt(5)+1)/(2*sqrt(5));
                    u = rand(G,1);
                    w = zeros(G,1);
                    w(u <= 1-p) = (1 - sqrt(5))/2;
                    w(u  > 1-p) = (1 + sqrt(5))/2;
                otherwise % "rademacher"
                    w = 2*(rand(G,1)>0.5) - 1;
            end
        end

        function V = clusterVar_(X, e, XtXinv, g, smallSample)
            % Classic CRV1 with optional small-sample adjustment
            p = size(X,2);
            S = zeros(p);
            G = max(g);
            for kk = 1:G
                idx = (g==kk);
                Xg = X(idx,:); eg = e(idx);
                v  = Xg' * eg;
                S  = S + (v * v.');
            end
            V = XtXinv * S * XtXinv;
            if smallSample
                N = size(X,1);
                adj = (G/(G-1)) * ((N-1)/(N-p));
                if isfinite(adj) && adj > 0
                    V = adj * V;
                end
            end
        end
    end
end
