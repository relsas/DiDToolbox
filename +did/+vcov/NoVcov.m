classdef NoVcov < did.vcov.VcovEngine
    methods
        function res = decorate(~, res, ~)
            if ~isfield(res,'Vcov'), res.Vcov = struct('type','none'); end
            if ~isfield(res,'SE'), res.SE = NaN; end
            if ~isfield(res,'t'),  res.t  = NaN; end
            if ~isfield(res,'p'),  res.p  = NaN; end
            if ~isfield(res,'df'), res.df = NaN; end
        end
    end
end
