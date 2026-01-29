classdef IW < did.estimators.Estimator
% did.estimators.IW
% Sun & Abraham (2021) Interaction-Weighted (IW) wrapper.
% Calls did.iw_estimator(...), accepting either a did.Dataset or a table.
%
% Typical use (Dataset workflow):
%   ds  = did.Dataset.fromTable(T, idVar="id", timeVar="time", yVar="y", dVar="D", weightVar="w");
%   est = did.estimators.IW(Comparison="notyet", SEMethod="multiplier", B=999);
%   mdl = did.Model(ds, est);
%   out = mdl.fit();

    properties (Access=private)
        opts struct
    end

    methods
        function obj = IW(varargin)
            % Allow: IW(struct)  OR  IW('Name',Value,...)
            if nargin==1 && isstruct(varargin{1})
                S = varargin{1};
            else
                if mod(nargin,2)~=0
                    error('did.estimators.IW:BadArgs','Constructor expects a struct or name–value pairs.');
                end
                S = struct();
                for i = 1:2:nargin
                    key = string(varargin{i});
                    S.(char(key)) = varargin{i+1};
                end
            end

            % ---- Aliases (compat) ----
            if isfield(S,'multiplier')  && ~isfield(S,'Multiplier'), S.Multiplier = S.multiplier;  end
            if isfield(S,'Detail')      && ~isfield(S,'Display'),      S.Display      = S.Detail;      end

            % ---- Defaults / sanitization ----
            if ~isfield(S,'SEMethod'),   S.SEMethod   = "multiplier"; end
            if ~isfield(S,'Comparison'), S.Comparison = "notyet";     end
            if ~isfield(S,'Delta'),      S.Delta      = 0;            end
            if ~isfield(S,'Weighting'),  S.Weighting  = "cohortShare"; end  % classic SA IW
            if isfield(S,'B'), S.B = max(0, floor(double(S.B))); end
            if ~isfield(S,'Seed') || ~isfinite(S.Seed) || S.Seed<0 || S.Seed>=2^32 || isnan(S.Seed)
                S.Seed = randi([1,1e7],1,1); else, S.Seed = double(S.Seed); end
            if ~isfield(S,'Studentize'), S.Studentize = true; end
            if ~isfield(S,'Display'),      S.Display      = true; end

            obj.opts = S;
        end

        function out = fit(obj, T)
            % Accept either a did.Dataset or a table; unwrap Dataset to table + var names
            S = obj.opts;

            if isa(T, 'did.Dataset')
                % Prefer validated, canonical names from Dataset
                S.idVar   = string(T.idVar);
                S.timeVar = string(T.timeVar);
                S.yVar    = string(T.yVar);
                S.dVar    = string(T.dVar);

                % Optional: weight & clusters if present in the Dataset
                if ~isempty(T.weightVar),  S.WeightVar  = string(T.weightVar);  end
                if isprop(T,'clusterVar')  && ~isempty(T.clusterVar),  S.ClusterVar  = string(T.clusterVar);  end
                if isprop(T,'clusterVar2') && ~isempty(T.clusterVar2), S.ClusterVar2 = string(T.clusterVar2); end

                TT = T.T;   % the validated table
            else
                % Fallback: assume user will provide names or defaults
                TT = T;
            end

            % Call functional implementation (it can also accept a Dataset, but we pass the table + names)
            nv  = did.utils.struct2nv(S);
            out = did.iw_estimator(TT, nv{:});
        end
    end
end
