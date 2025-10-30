classdef Dataset
    % DID.DATASET  Validated panel + var-name metadata (no guessing).
    %
    % Create once, reuse everywhere:
    %   ds = did.Dataset.fromTable(T, idVar="id", timeVar="time", yVar="y", dVar="D", weightVar="w");
    %
    % Helpful methods:
    %   k   = ds.get("eventTime");   % k = t_int - g (NaN for never-treated)
    %   t   = ds.get("t_int");       % stable integer time index (1..nT)
    %   g   = ds.get("g");           % first-treat time in t_int units; 0 = never
    %   S   = ds.describe(true);     % reuse did.dataDesc to print a summary
    %   ds2 = ds.subset(mask);       % filtered new Dataset (no mutation)
    %   ds3 = ds.withWeights("w2");  % switch weight variable (validated)
    %   TT  = ds.materialize("t_int","g","eventTime","weights"); % table view incl. helpers

    properties (SetAccess=immutable)
        T table
        idVar   (1,1) string
        timeVar (1,1) string
        yVar    (1,1) string
        dVar    (1,1) string
        weightVar string = string.empty

        % Canonical names for materialized helper columns
        cohortVarName    (1,1) string = "gVar"       % <-- NEW: cohort column name
        eventTimeVarName (1,1) string = "eventTime"  % already conceptually there
    end


    properties (Access=private)
        Cache struct = struct();   % lazy computed helpers: t_int, g, eventTime
        Info  struct = struct();   % stash: dataDesc, provenance flags, etc.
    end

    methods (Static)
        function ds = fromTable(T, varargin)
            % Supports BOTH:
            %   did.Dataset.fromTable(T, idVar="id", timeVar="time", yVar="y", dVar="D", ...)
            %   did.Dataset.fromTable(T, struct('idVar',"id", 'timeVar',"time", ...))

            arguments
                T table
            end
            arguments (Repeating)
                varargin
            end

            % ---- Normalize to an 'opts' struct ----
            if numel(varargin) == 1 && isstruct(varargin{1})
                opts = varargin{1};                     % struct path
            else
                if mod(numel(varargin),2) ~= 0
                    error('did:Dataset:fromTable:NVpairs', ...
                        'Name–value pairs must come in key/value pairs.');
                end
                opts = struct();                        % NV-pair path
                for k = 1:2:numel(varargin)
                    key = string(varargin{k});
                    val = varargin{k+1};
                    opts.(key) = val;
                end
            end

            % ---- Defaults ----
            if ~isfield(opts,'idVar'),            opts.idVar = "";            end
            if ~isfield(opts,'timeVar'),          opts.timeVar = "";          end
            if ~isfield(opts,'yVar'),             opts.yVar = "";             end
            if ~isfield(opts,'dVar'),             opts.dVar = "";             end
            if ~isfield(opts,'weightVar'),        opts.weightVar = string.empty; end
            if ~isfield(opts,'computeEventTime'), opts.computeEventTime = true;   end
            if ~isfield(opts,'describe'),         opts.describe = true;          end
            if ~isfield(opts,'Print'),         opts.Print = true;          end

            % ---- Require the core four ----
            if strlength(string(opts.idVar))==0 || strlength(string(opts.timeVar))==0 || ...
                    strlength(string(opts.yVar))==0  || strlength(string(opts.dVar))==0
                error('did:Dataset:fromTable:MissingNV', ...
                    'idVar, timeVar, yVar, and dVar must be provided.');
            end
            % 1) Lightweight structural validation
            did.utils.validatePanel(T, opts.idVar, opts.timeVar, opts.yVar, opts.dVar, opts.weightVar);

            % 2) Optional rich diagnostics (reuses your existing function)
            dsDesc = [];
            if opts.describe && opts.Print
                dsDesc = did.dataDesc(T, ...
                    idVar=opts.idVar, timeVar=opts.timeVar, dVar=opts.dVar, ...
                    yVar=opts.yVar, Display=true);
            end

            % 3) Construct (no mutation of T)
            ds = did.Dataset(T, opts.idVar, opts.timeVar, opts.yVar, opts.dVar, opts.weightVar);
            ds.Info.Description       = dsDesc;
            ds.Info.CreatedOn         = datetime('now');
            ds.Info.ComputeEventTime  = opts.computeEventTime;  % informational
        end
    end

    methods
        function v = get(ds, what)
            % GET  Access lazy canonical helpers: "t_int", "g", "eventTime"
            key = string(what);
            switch key
                case "t_int"
                    if ~isfield(ds.Cache,'t_int')
                        ds.Cache.t_int = did.utils.timeInt(ds.T, ds.timeVar);
                    end
                    v = ds.Cache.t_int;

                case "g"
                    if ~isfield(ds.Cache,'g')
                        % Build via the exact same helper used elsewhere so results match
                        [T2, gname] = did.utils.ensureGvar(ds.T, ds.idVar, ds.timeVar, ds.dVar);
                        ds.Cache.g = T2.(gname);  % numeric: time label (or 0)
                    end
                    v = ds.Cache.g;

                case "eventTime"
                    if ~isfield(ds.Cache,'eventTime')
                        t = ds.get("t_int");
                        g = ds.get("g");
                        ds.Cache.eventTime = did.utils.eventTimeFrom(ds.T, ds.idVar, t, g);
                    end
                    v = ds.Cache.eventTime;

                otherwise
                    error('did:Dataset:get:UnknownKey','Unknown key "%s".', what);
            end
        end

        function S = describe(ds, Display)
            % DESCRIBE  Re-run (or print) the data overview using did.dataDesc.
            if nargin<2, Display=false; end
            S = did.dataDesc(ds.T, idVar=ds.idVar, timeVar=ds.timeVar, dVar=ds.dVar, ...
                yVar=ds.yVar, Display=Display);
        end

        function ds2 = subset(ds, rowMask)
            % SUBSET  Return a new Dataset filtered by a logical mask (no mutation).
            validateattributes(rowMask, {'logical'},{'vector','numel',height(ds.T)});
            T2 = ds.T(rowMask,:);
            ds2 = did.Dataset.fromTable(T2, ...
                idVar=ds.idVar, timeVar=ds.timeVar, yVar=ds.yVar, dVar=ds.dVar, ...
                weightVar=ds.weightVar, describe=false);
        end

        function ds2 = withWeights(ds, wname)
            % WITHWEIGHTS  Return a new Dataset using another weight column.
            wname = string(wname);
            if strlength(wname)>0 && ~any(strcmp(wname, ds.T.Properties.VariableNames))
                error('did:Dataset:withWeights:NoVar','Weight variable "%s" not found.', wname);
            end
            ds2 = did.Dataset.fromTable(ds.T, ...
                idVar=ds.idVar, timeVar=ds.timeVar, yVar=ds.yVar, dVar=ds.dVar, ...
                weightVar=wname, describe=false);
        end

        function TT = materialize(ds, varargin)
            % MATERIALIZE  Return a table view with requested helpers appended.
            % Example: TT = ds.materialize("t_int","gVar","eventTime","weights");

            TT = ds.T;
            for j = 1:numel(varargin)
                key = string(varargin{j});
                switch key
                    case "weights"
                        % weightVar (if any) is already in TT; nothing to do.

                    case {"t_int","tInt"}   % expose stable time index (helper)
                        TT.t_int = did.utils.timeInt(ds.T, ds.timeVar);

                    case {"g","gVar","cohort"}
                        % Compute once via the cache and append with the *canonical* name
                        g = ds.get("g");            % numeric: time label or 0
                        TT.(ds.cohortVarName) = g;  % e.g., 'gVar'

                    case {"eventTime","k","K"}
                        k = ds.get("eventTime");
                        TT.(ds.eventTimeVarName) = k;

                    otherwise
                        error('did:Dataset:materialize:UnknownKey','Unknown key "%s".', key);
                end
            end
        end

    end

    methods (Access=private)
        function obj = Dataset(T,idVar,timeVar,yVar,dVar,wVar)
            obj.T        = T;
            obj.idVar    = idVar;
            obj.timeVar  = timeVar;
            obj.yVar     = yVar;
            obj.dVar     = dVar;
            obj.weightVar = wVar;
        end
    end
end
