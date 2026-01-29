function est = getEstimator(method, opts)
% Build an estimator object from method + opts (NV-pairs).
% --------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Germany / Last Change: 11/16/2025
% --------------------------------------------------------------------
arguments
    method (1,1) string
    opts struct
end

m = lower(method);
switch m

    case {"twfe"}
        allowed = ["Covariates","details","absorbMethod","absorbUseGPU","absorbTol","absorbMaxIter","preTrend","preTrendMode","Display"];
        s  = did.utils.keepfields(opts, allowed);
        nv = did.utils.struct2nv(s);
        est = did.estimators.TWFE(nv{:});

    case {"bjs","borusyak"}
        allowed = ["Covariates","Horizons","SEMethod","BootReps","Seed","Display","useParallel"];
        s  = did.utils.keepfields(opts, allowed);
        if opts.SEMethod~="LOO"
            opts.SEMethod="LOO";
        end
        nv = did.utils.struct2nv(s);
        est = did.estimators.BJS(nv{:});

    case {"bacon","goodman-bacon"}
        allowed = ["Display"];
        s = did.utils.keepfields(opts, allowed);
        nv = did.utils.struct2nv(s);
        est = did.estimators.Bacon(nv{:});

    case {"ch","didm","did_m"}
        if isfield(opts,'rngSeed') && ~isfield(opts,'Seed'),      opts.Seed = opts.rngSeed; end
        if isfield(opts,'Weight')  && ~isfield(opts,'WeightVar'), opts.WeightVar = opts.Weight; end
        if ~isfield(opts,'Seed') || ~isfinite(opts.Seed), opts.Seed = []; end

        allowed = ["B","ComputePlacebo","Seed","WeightVar","Display","Details", ...
            "Covariates","CovarSample","Horizons","Placebos"];
        s  = did.utils.keepfields(opts, allowed);
        nv = did.utils.struct2nv(s);
        est = did.estimators.CH(nv{:});

    case {"cs","cs2021"}
        if isfield(opts,'studentize') && ~isfield(opts,'Studentize'), opts.Studentize = opts.studentize; end
        if isfield(opts,'multiplier') && ~isfield(opts,'Multiplier'), opts.Multiplier = opts.multiplier; end


        if isfield(opts,'B'), opts.B = max(0, floor(double(opts.B))); end
        if ~isfield(opts,'Seed') || ~isfinite(opts.Seed) || opts.Seed<0 || opts.Seed>=2^32||isnan(opts.Seed)
            opts.Seed = randi([1,1e7],1,1); else, opts.Seed = double(opts.Seed); end

        allowed = ["Approach","Comparison","Delta","Covariates","WeightVar","Display","Details", ...
            "SEMethod","B","Seed","Multiplier","Studentize","ClusterVar","ClusterVar2", ...
            "SmallSample","UseParallel","CrossFit","Kfolds","StratifyFoldsBy","Weighting", "PreTrendBase"];
        s  = did.utils.keepfields(opts, allowed);

        % Robust property assignment
        est = did.estimators.CS();
        fn  = fieldnames(s);
        for i = 1:numel(fn)
            pname = fn{i};
            if isprop(est, pname), est.(pname) = s.(pname); end
        end

    case {"iw","sa"}
        if isfield(opts,'multiplier')  && ~isfield(opts,'Multiplier'),opts.Multiplier = opts.multiplier; end

        if isfield(opts,'B'), opts.B = max(0, floor(double(opts.B))); end
        if ~isfield(opts,'Seed') || ~isfinite(opts.Seed) || opts.Seed<0 || opts.Seed>=2^32 || isnan(opts.Seed)
            opts.Seed = randi([1,1e7],1,1); else, opts.Seed = double(opts.Seed); end
        if ~isfield(opts,'SEMethod'),   opts.SEMethod = "multiplier"; end
        if ~isfield(opts,'Comparison'), opts.Comparison = "notyet";   end
        if ~isfield(opts,'Delta'),      opts.Delta = 0;               end

        allowed = ["idVar","timeVar","yVar","dVar","WeightVar", ...
            "Delta","Comparison","SEMethod","B","Seed","Multiplier","Studentize", ...
            "ClusterVar","ClusterVar2","Display","Weighting"];
        s  = did.utils.keepfields(opts, allowed);
        nv = did.utils.struct2nv(s);
        est = did.estimators.IW(nv{:});

    case {"sdid","sc","syntheticcontrol"}
        allowed = ["Weights","Regularization","Solver","SEMethod","B","Display","useParallel"];
        s = did.utils.keepfields(opts, allowed);
        nv = did.utils.struct2nv(s);
        est = did.estimators.SDID(nv{:});
        if strcmpi(method,"sc") || strcmpi(method,"syntheticcontrol")
            est.Weights = "SC";
        end

    case "cem"
        allowed = ["Covariates","nBins","Robust","Display","details","Estimator"];
        s = did.utils.keepfields(opts, allowed);
        % map opts.details -> Details property if needed (constructor is case insensitive usually?)
        % But CEM.m has property "Details".
        if isfield(s, 'details'), s.Details = s.details; s = rmfield(s,'details'); end

        nv = did.utils.struct2nv(s);
        est = did.estimators.CEM(nv{:});

    otherwise
        error('did:getEstimator:UnknownMethod','Unknown method "%s".', method);
end
end
