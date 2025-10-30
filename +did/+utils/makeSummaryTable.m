function S = makeSummaryTable(res, args)
% MAKESUMMARYTABLE  Build a compact, filtered summary table.
% Usage examples:
%   S = did.utils.makeSummaryTable(res);                         % default: hide FE & const
%   S = did.utils.makeSummaryTable(res, exclude=["^const$"]);    % hide only intercept
%   S = did.utils.makeSummaryTable(res, keepMainOnly=true);      % only D (and covars if requested)

arguments
    res struct
    args.exclude (1,:) string = ["^FE_i_","^FE_t_","^const$"]  % regex patterns to drop
    args.include (1,:) string = string.empty                   % regex to *keep* (applied after exclude)
    args.keepMainOnly (1,1) logical = false                    % keep only D (and covars if include says so)
    args.addCovars (1,:) string = string.empty                 % if keepMainOnly=true, also keep these by name
end

% --- Special-case: BACON should summarize by component type ---
if isfield(res,'method') && strcmpi(string(res.method), "BACON") ...
        && isfield(res,'ByType') && istable(res.ByType)
    S = res.ByType;
    return
end



% --- get base coef names & estimates ---
assert(isfield(res,'coef') && istable(res.coef), ...
    'makeSummaryTable: res.coef (table) required');

names = string(res.coef.Name);
est   = res.coef.Estimate;

% --- standard errors: prefer vcov, else fall back to res.SE for main effect only ---
if isfield(res,'vcov') && ~isempty(res.vcov)
    seAll = sqrt(max(0, diag(res.vcov)));
else
    seAll = nan(size(est));
    if isfield(res,'Diagnostics') && isfield(res.Diagnostics,'design')
        idxD = res.Diagnostics.design.idxD;
        if isfield(res,'SE') && ~isempty(res.SE)
            seAll(idxD) = res.SE;
        end
    end
end

% --- degrees of freedom ---
df = NaN;
if isfield(res,'df') && ~isempty(res.df), df = res.df;
elseif isfield(res,'Vcov') && isstruct(res.Vcov) && isfield(res.Vcov,'df'), df = res.Vcov.df;
end

% --- t / p ---
tAll = est ./ seAll;
if isfinite(df)
    pAll = 2*tcdf(-abs(tAll), df);
else
    pAll = 2*normcdf(-abs(tAll));  % large-sample fallback
end

T = table(names, est, seAll, tAll, pAll, 'VariableNames', ...
    {'Name','Estimate','SE','tStat','pValue'});

% --- figure out the main effect (D) name once ---
mainName = "";
try
    mainName = string(res.Diagnostics.design.names(res.Diagnostics.design.idxD));
catch
end

% --- filtering rules ---
mask = true(height(T),1);

% Convert to cellstr so regexp returns a cell array even for 1 row
nameCell = cellstr(T.Name);

% 1) default exclude (FE & const)
for pat = args.exclude
    mask = mask & cellfun(@isempty, regexp(nameCell, pat));
end

% 2) include patterns (re-add anything user insists on)
if ~isempty(args.include)
    kee = false(height(T),1);
    for pat = args.include
        kee = kee | ~cellfun(@isempty, regexp(nameCell, pat));
    end
    mask = mask | kee;
end

% 3) keepMainOnly option
if args.keepMainOnly
    keep = (T.Name == mainName);
    if ~isempty(args.addCovars)
        keep = keep | ismember(T.Name, args.addCovars);
    end
    mask = keep;
end

S = T(mask, :);

% (Optional) sort: main first, then covars, then time FE if any slipped through
if any(S.Name == mainName)
    ord = 1:height(S);
    mainIdx = find(S.Name == mainName, 1);
    ord = [mainIdx, setdiff(ord, mainIdx, 'stable')];
    S = S(ord, :);
end
end

