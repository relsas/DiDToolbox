classdef Wooldridge < did.estimators.Estimator
    % WOOLDRIDGE  Two-way Mundlak (pooled OLS) estimator for staggered DiD.
    % Wrapper for did.wooldridge_TB.

    properties
        % Options forwarded to wooldridge_TB
        Cluster (1,1) logical = true
        Imputation (1,1) string = "mundlak"  % specific to 2-way mundlak approach
        DropTimeBase (1,1) double = 1
        EventStudy (1,1) logical = false     % If true, estimates full dynamic leads/lags

        % Standard options
        Display (1,1) logical = true
    end

    methods
        function obj = Wooldridge(varargin)
            % WOOLDRIDGE Constructor handles Name-Value pairs
            for k = 1:2:numel(varargin)
                name = string(varargin{k});
                if k+1 > numel(varargin)
                    error('did:Wooldridge:MissingValue','Missing value for option "%s".',name);
                end
                val = varargin{k+1};

                if isprop(obj, name)
                    obj.(name) = val;
                end
                % Also store in general Options struct for flexibility
                obj.Options.(name) = val;
            end
            obj.Name = "Wooldridge (Mundlak)";
        end

        function res = fit(obj, ds)
            % Fits the estimator using the dataset ds
            arguments
                obj
                ds
            end

            % Check requirements
            if ~ismember('g', ds.T.Properties.VariableNames)
                error('did:Wooldridge:MissingG','Wooldridge estimator requires a "g" column (group/cohort time).');
            end

            % Build options for lower-level call
            opts = struct();
            opts.yVar = ds.yVar;
            opts.idVar = ds.idVar;
            opts.timeVar = ds.timeVar;
            opts.gVar = "g"; % enforced by check above

            % MAPPING: wooldridge_TB args are {clusters, ...} not {Cluster}
            opts.Display = obj.Display;
            opts.DropTimeBase = obj.DropTimeBase;

            % Check for specific cluster setting in Options, otherwise default
            if isfield(obj.Options, 'clusters')
                opts.clusters = obj.Options.clusters;
            elseif isfield(obj.Options, 'ClusterVar')
                opts.clusters = obj.Options.ClusterVar;
            else
                opts.clusters = []; % defaults to id in wooldridge_TB
            end

            % Note: wooldridge_TB currently handles clustering internally.
            % If obj.Cluster == false, strictly speaking we should enforce homoskedastic,
            % but wooldridge_TB doesn't support it yet.
            % For now, we rely on wooldridge_TB defaults.
            if isprop(obj, 'Covariates')
                % If class had a Covariates property, we would pass it.
                % For now, pass if in Options
                if isfield(obj.Options, 'Covariates')
                    opts.Covariates = obj.Options.Covariates;
                end
            end

            % Call implementation
            % Note: wooldridge_TB accepts (tbl, Name, Val...)
            % We convert struct opts to name-value pairs
            nvPairs = namedargs2cell(opts);
            out = did.wooldridge_TB(ds.T, nvPairs{:});

            % Convert output to uniform result struct if needed?
            % wooldridge_TB returns a struct with .summaryTable, etc.
            % Estimator.fit is expected to return something compatible with did.Model.
            res = out;
        end
    end
end
