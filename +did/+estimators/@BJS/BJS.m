classdef BJS < did.estimators.Estimator
    % did.estimators.BJS  Wrapper for Borusyak–Jaravel–Spiess (2024) imputation.
    % 
    % 
    % ------------------------------------------------------------------------
    % Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
    % Last change: 09/30/2025
    % ------------------------------------------------------------------------
    properties
        % Forwarded 1:1 to bjs_imputation
        Covariates string = string.empty(1,0)
        Horizons   double = []
        SEMethod   string  = "LOO"    % "LOO" | "BootstrapUnit" | "None"
        BootReps   double  = 0
        Seed       double  = randi([1,1e7],1,1)
        Display    logical = true
        useParallel (1,1) double = 1
    end

    methods
        function obj = BJS(varargin)
            % Simple name–value parser
            if ~isempty(varargin)
                if mod(numel(varargin),2)~=0
                    error('did:BJS:NameValue','Constructor expects name–value pairs.');
                end
                
                for k = 1:2:numel(varargin)
                    name  = lower(string(varargin{k}));
                    value = varargin{k+1};
                    switch name
                        case 'covariates', obj.Covariates = string(value);
                        case 'horizons',   obj.Horizons   = double(value);
                        case 'semethod',   obj.SEMethod   = "LOO";
                        case 'bootreps',   obj.BootReps   = double(value);
                        case 'seed',       obj.Seed       = double(value);
                        case 'display',    obj.Display    = logical(value);
                        case 'useparallel', obj.useParallel = double(value);
                        otherwise
                            warning('did:BJS:UnknownOption','Ignoring unknown option "%s".', name);
                    end
                end
            end
            
        end

        function res = fit(obj, ds)
    [T, idVar, timeVar, yVar, dVar] = obj.unpackDs(ds);

    % Run your function (unchanged behavior)
    res = did.estimators.bjs_imputation(T, ...
        idVar=idVar, timeVar=timeVar, yVar=yVar, dVar=dVar, ...
        Covariates=obj.Covariates, Horizons=obj.Horizons, ...
        SEMethod=obj.SEMethod, BootReps=obj.BootReps, Seed=obj.Seed, ...
        Display=obj.Display, useParallel=obj.useParallel);

    % --- Normalize for did.utils.makeSummaryTable ---
    [coef, vcov, mainIdx] = did.estimators.BJS.packCoefForSummary(res);

    res.Method = "BJS";
    res.coef   = coef;          % must have variables: Name, Estimate
    res.vcov   = vcov;          % diagonal from BJS SEs so SEs show up
    res.df     = Inf;           % normal approx for p-values

    % Minimal design info so keepMainOnly etc. works
    res.Diagnostics.design = struct('names', coef.Name, 'idxD', mainIdx);
end


    end

    methods (Access = private)
        function [T, idVar, timeVar, yVar, dVar] = unpackDs(~, ds)
            % ---- 1) Find the table inside 'ds' robustly ----
            T = [];

            % (a) Fast path: typical field names
            candidates = {'Table','Tbl','T','Data','Dataset','Raw','RawTable','table','tbl','data'};
            for k = 1:numel(candidates)
                nm = candidates{k};
                if isprop(ds, nm)
                    try
                        val = ds.(nm);
                        if istable(val), T = val; break; end
                    catch
                    end
                end
            end

            % (b) Reflect over public properties to find any table
            if isempty(T)
                try
                    mc = metaclass(ds);
                    for p = reshape(mc.PropertyList,1,[])
                        if ~p.Hidden && strcmpi(p.GetAccess,'public')
                            try
                                val = ds.(p.Name);
                                if istable(val), T = val; break; end
                            catch
                            end
                        end
                    end
                catch
                end
            end

            % (c) Try common getter methods
            if isempty(T)
                getters = {'getTable','toTable','asTable','table','export','exportTable'};
                for k = 1:numel(getters)
                    g = getters{k};
                    if ismethod(ds, g)
                        try
                            val = ds.(g)();
                            if istable(val), T = val; break; end
                        catch
                        end
                    end
                end
            end

            if isempty(T)
                error('did:BJS:DatasetNoTable', ['Could not locate a table inside ds. ',
                    'Expose a table as a public property (e.g., .Table) or implement getTable()/toTable().']);
            end

            % ---- 2) Infer variable names ----
            V = string(T.Properties.VariableNames);

            % (a) Preferred: standardized names produced by did.Dataset.fromTable
            if all(ismember(["i","t","y","D"], V))
                idVar="i"; timeVar="t"; yVar="y"; dVar="D"; return
            end

            % (b) Common lowercase
            if all(ismember(["id","time","y","d"], V))
                idVar="id"; timeVar="time"; yVar="y"; dVar="d"; return
            end

            % (c) If ds exposes explicit names, try a few likely properties
            nameProps = { ...
                "idVarName","timeVarName","yVarName","dVarName", ...
                "idVar","timeVar","yVar","dVar", ...
                "IdVarName","TimeVarName","YVarName","DVarName" ...
                };
            vals = strings(1,4); found = false(1,4);
            for k = 1:numel(nameProps)
                nm = nameProps{k};
                if isprop(ds, nm)
                    try
                        val = string(ds.(nm));
                        switch lower(nm)
                            case {'idvarname','idvar'},   vals(1)=val; found(1)=true;
                            case {'timevarname','timevar'}, vals(2)=val; found(2)=true;
                            case {'yvarname','yvar'},     vals(3)=val; found(3)=true;
                            case {'dvarname','dvar'},     vals(4)=val; found(4)=true;
                        end
                    catch
                    end
                end
            end
            if all(found) && all(ismember(vals, V))
                idVar=vals(1); timeVar=vals(2); yVar=vals(3); dVar=vals(4); return
            end

            % (d) Give a precise, helpful error
            error('did:BJS:DatasetVars', ['Could not infer variable names from dataset table. ',
                'Expected standardized columns (i,t,y,D) or (id,time,y,d). ',
                'Available columns: %s'], strjoin(V, ', '));
        end
    end

   methods (Static)
    function [coef, vcov, mainIdx] = packCoefForSummary(res)
        % Build rows for makeSummaryTable: Name, Estimate
        names = strings(0,1);
        ests  = zeros(0,1);
        ses   = nan(0,1);

        % ---- Overall ATT ----
        mainIdx = [];  % we will set to the overall row if present
        if isfield(res,'ATT_overall') && ~isempty(res.ATT_overall)
            r = res.ATT_overall;
            names(end+1,1) = "ATT_overall";
            ests (end+1,1) = double(r.ATT);
            ses  (end+1,1) = double(r.SE);
            mainIdx = 1;
        end

        % ---- ATT by horizon ----
        if isfield(res,'ATT_by_horizon') && ~isempty(res.ATT_by_horizon)
            A = res.ATT_by_horizon;
            for i = 1:height(A)
                k    = double(A.k(i));
                names(end+1,1) = "ATT_k=" + string(k);
                ests (end+1,1) = double(A.ATT_k(i));
                ses  (end+1,1) = double(A.SE(i));
            end
            % If no explicit main yet, prefer k==0 as main effect
            if isempty(mainIdx)
                idx0 = find(double(A.k)==0, 1);
                if ~isempty(idx0)
                    mainIdx = 1 + idx0;  % 1 for overall row offset
                end
            end
        end

        % ---- Optional: cohort-balanced ATT(k) ----
        if isfield(res,'ATT_balanced') && ~isempty(res.ATT_balanced)
            B = res.ATT_balanced;
            for i = 1:height(B)
                k    = double(B.k(i));
                names(end+1,1) = "ATT_bal_k=" + string(k);
                ests (end+1,1) = double(B.ATT_k_bal(i));
                ses  (end+1,1) = double(B.SE(i));
            end
        end

        % Fallback main index
        if isempty(mainIdx)
            mainIdx = 1;
        end

        % Output table and diagonal vcov (so makeSummaryTable can compute SE/t/p)
        coef = table(names, ests, 'VariableNames', ["Name","Estimate"]);
        vcov = diag(ses.^2);
    end
end


end
