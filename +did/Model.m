classdef Model
    properties
        ds
        est  % must be did.estimators.Estimator
    end

    methods
        function obj = Model(ds, est)
            arguments
                ds
                est did.estimators.Estimator
            end
            obj.ds  = ds;
            obj.est = est;
        end

        function res = fit(obj)
            % Run estimator (point estimates)
            res = obj.est.fit(obj.ds);

            % ---- Apply VCOV engine if present ----
            if isprop(obj.est, 'VcovEngine') && ~isempty(obj.est.VcovEngine)
                try
                    eng = obj.est.VcovEngine;

                    % Auto-upgrade common mis-set values (strings/structs) into objects
                    if isstring(eng) || ischar(eng)
                        eng = did.factories.getVcov(struct('vcov',string(eng)));
                        obj.est.VcovEngine = eng;
                    elseif isstruct(eng) && isfield(eng,'vcov')
                        eng = did.factories.getVcov(eng);
                        obj.est.VcovEngine = eng;
                    end

                    % Only call if this is truly an object with 'decorate'
                    if isobject(eng) && any(strcmp(methods(eng), 'decorate'))
                        res = eng.decorate(res, obj.ds);
                    else
                        warning('[did.Model] VCOV engine is not an object with a decorate method (type: %s). Skipping.', class(eng));
                    end
                catch ME
                    warning('[did.Model] VCOV decoration failed: %s', '%s',ME.message);
                end
            end


            % ---- Normalize VCOV to a standard field ----
            % Prefer res.Vcov.matrix (from VCOV engines), mirror to res.vcov for convenience
            if isfield(res, 'Vcov') && isstruct(res.Vcov) && isfield(res.Vcov, 'matrix')
                res.vcov = res.Vcov.matrix;
            elseif ~isfield(res,'vcov')
                % leave as-is; user may be using a no-op VCOV
                res.vcov = [];
            end

            try
                if isfield(res, 'summaryTable') && istable(res.summaryTable)
                    res.summaryTable = res.summaryTable;                 % use estimator-supplied table
                else
                    res.summaryTable = did.utils.makeSummaryTable(res);  % fallback builder
                    
                end
            catch ME
                warning('[did.Model] summaryTable not created: %s','%s', ME.message);
                summaryTable = table(); %#ok<NASGU>
            end

        end
    end
end


