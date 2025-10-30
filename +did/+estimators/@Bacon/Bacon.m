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
            [T2, gVar] = obj.ensureGvar(ds.T, ds.idVar, ds.timeVar, ds.dVar);

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

    methods (Access = private)

        function [T2, gVar] = ensureGvar(~, T, idVar, timeVar, dVar)
            % Build integer time index and first-treat cohort g (0 = never)
            T2 = T;
            [~,~,t_idx] = unique(T2.(timeVar), 'stable');
            T2.t_int = double(t_idx);   % <-- single underscore (exists)

            uniqueTime = unique(T2.(timeVar));
            D = T2.(dVar)==1;
            if any(D)
                % summarize MIN over the variable that actually exists: "t_int"
                adopt = groupsummary(T2(D,:), idVar, "min", "t_int");

                % Rename min_t_int -> g_int (robust to column order)
                nv = adopt.Properties.VariableNames;
                ix = find(strcmpi(nv,'min_t_int'), 1);
                if ~isempty(ix)
                    adopt.Properties.VariableNames{ix} = 'g_int';
                else
                    adopt.Properties.VariableNames{end} = 'g_int';
                end

                % Join back and build g_bacon (0 for never-treated)
                T2 = outerjoin(T2, adopt(:,[idVar,"g_int"]), ...
                    "Keys", idVar, "MergeKeys", true, "Type","left");
                T2.g_bacon(~isnan(T2.g_int)) = uniqueTime(T2.g_int(~isnan(T2.g_int)));

                T2.g_int = [];
            else
                T2.g_bacon = zeros(height(T2),1);
            end

            % Cleanup temp column and return gVar name
            T2.t_int = [];
            gVar = "g_bacon";
        end
    end
end
