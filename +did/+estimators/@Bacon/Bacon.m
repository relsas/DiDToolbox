classdef Bacon < did.estimators.Estimator
    % did.estimators.Bacon  Wrapper for Goodman–Bacon decomposition.
    % Calls your canonical function: did.bacon_decomp(T, args).

    properties
        Display logical = true
    end

    methods
        function obj = Bacon(varargin)
            % Name–value parser (allows Display=true/false)
            if ~isempty(varargin)
                if mod(numel(varargin),2)~=0
                    error('did:Bacon:NameValue','Constructor expects name–value pairs.');
                end
                for k = 1:2:numel(varargin)
                    name = lower(string(varargin{k}));
                    val  = varargin{k+1};
                    switch name
                        case 'display', obj.Display = logical(val);
                        otherwise, warning('did:Bacon:UnknownOption','Ignoring option "%s".', name);
                    end
                end
            end
        end

        function res = fit(obj, ds)

            % Ensure a cohort variable g (first-treat time index; 0 = never)
            [T2, gVar] = did.utils.ensureGvar(ds.T, ds.idVar, ds.timeVar, ds.dVar);
            nCohort = length(find(unique(T2.(gVar))~=0));
            if nCohort<=1
                error("There are too few cohorts to run the Bacon decomposition: nCohort= "+nCohort);
            end
            % Do Bacon
            out = did.bacon_decomp(T2, ...
                idvar=ds.idVar, tvar=ds.timeVar, gvar=gVar, yvar=ds.yVar, ...
                isPrint=obj.Display);

            % Normalize for did.Model/makeSummaryTable
            res = out;                 % keep legacy fields: Pairs, ByType, etc.
            res.Method = "BACON";

            est = NaN;
            if isfield(out,'BaconTWFE'), est = out.BaconTWFE; end
            if isnan(est) && isfield(out,'TWFE'),  est = out.TWFE;  end
            if isnan(est) && isfield(out,'bacon'), est = out.bacon; end

            res.coef = table("BaconTWFE"', double(est)', 'VariableNames', ["Name","Estimate"]);
            res.vcov = [];      % no SEs from decomposition
            res.df   = NaN;
            res.Diagnostics.design = struct('names', res.coef.Name, 'idxD', 1);
        end
    end

    
end
