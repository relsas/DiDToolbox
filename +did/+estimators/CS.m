classdef CS < did.estimators.Estimator
    % DID.ESTIMATORS.CS  Callaway & Sant'Anna (2021) DiD with multiple periods/timing.

    properties
        % Controls
        Approach (1,1) string = "unconditional"
        Comparison (1,1) string = "never"
        Delta (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        Covariates string = string.empty(1,0)
        WeightVar (1,1) string = ""
        Display (1,1) logical = true
        details (1,1) logical = false

        % Inference controls
        SEMethod (1,1) string = "multiplier"
        B (1,1) double {mustBeInteger,mustBeNonnegative} = 0
        Seed (1,1) double = randi([1,1e7],1,1)
        Multiplier (1,1) string = "rademacher"
        Studentize (1,1) logical = true
        ClusterVar string = string.empty(1,0)
        ClusterVar2 string = string.empty(1,0)
        SmallSample (1,1) logical = true
        CrossFit (1,1) logical = false
        Kfolds (1,1) double = 5
        StratifyFoldsBy (1,1) string = "none";
        % Aggregation
        Weighting string = "treatedObs"
        PreTrendBase (1,1) string = "universal" % "universal" (g-1) or "varying" (t-1)
    end

    methods
        function obj = CS(varargin)
            %CS Callaway & Sant'Anna (2021) Estimator
            %
            % obj = CS('Name',Value,...)

            % Parse Name-Value pairs
            for k = 1:2:numel(varargin)
                key = string(varargin{k});
                if k+1 > numel(varargin)
                    error('did:CS:MissingValue', 'Missing value for option "%s".', key);
                end
                val = varargin{k+1};

                % Case-insensitive property setting
                try
                    obj.(key) = val;
                catch
                    % Try matching property names case-insensitively
                    props = string(properties(obj));
                    idx = find(strcmpi(key, props), 1);
                    if ~isempty(idx)
                        obj.(props(idx)) = val;
                    else
                        error('did:CS:UnknownOption', 'Unknown option "%s".', key);
                    end
                end
            end
        end

        function out = fit(obj, ds_or_T, varargin)
            % Accept did.Dataset OR (table + name–values idVar=..., timeVar=..., ...)
            nv = struct(varargin{:});

            % --- Resolve inputs once
            [T,idVar,timeVar,yVar,dVar,wVar,~] = did.adapters.resolveInputs(ds_or_T, nv);

            % Prefer a Dataset's weightVar, else fall back to obj.WeightVar if present
            % Is a weight variable already resolved from ds/adapter?
            hasWVar   = ~isempty(wVar) && all(strlength(wVar) > 0);

            % If not, try the estimator's WeightVar property (if it names an existing column)
            hasObjW   = ~isempty(obj.WeightVar) && all(strlength(obj.WeightVar) > 0);
            if ~hasWVar && hasObjW && any(strcmp(obj.WeightVar, T.Properties.VariableNames))
                wVar = obj.WeightVar;
                hasWVar = true;
            end


            S = struct( ...
                'idVar',       idVar, ...
                'timeVar',     timeVar, ...
                'yVar',        yVar, ...
                'dVar',        dVar, ...
                'WeightVar',   wVar, ...
                'Covariates',  obj.Covariates, ...
                'Approach',    obj.Approach, ...
                'Comparison',  obj.Comparison, ...
                'Delta',       obj.Delta, ...
                'SEMethod',    obj.SEMethod, ...
                'B',           obj.B, ...
                'Seed',        obj.Seed, ...
                'Multiplier',  obj.Multiplier, ...
                'Studentize',  obj.Studentize, ...
                'ClusterVar',  obj.ClusterVar, ...
                'ClusterVar2', obj.ClusterVar2, ...
                'SmallSample', obj.SmallSample, ...
                'Display',       obj.Display, ...
                'details',     obj.details, ...
                'CrossFit',    obj.CrossFit, ...
                'Kfolds',      obj.Kfolds, ...
                'StratifyFoldsBy', obj.StratifyFoldsBy, ...
                'Weighting', obj.Weighting, ...
                'PreTrendBase', obj.PreTrendBase ...
                );

            nv = did.utils.struct2nv(S);
            out = did.cs_estimator(T, nv{:});

            % coef fallback as before
            if isfield(out,'summaryTable')
                co = out.summaryTable;
                want = {'Name','Estimate','SE','tStat','pValue'};
                have = want(ismember(want,co.Properties.VariableNames));
                out.coef = co(:,have);
            end
        end
    end
end
