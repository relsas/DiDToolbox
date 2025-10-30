function [T2, gVar] = ensureGvar(T, idVar, timeVar, dVar)
%ENSUREGVAR  Add first-treatment cohort variable 'gVar' (time label or 0).
%   [T2,gVar] = did.utils.ensureGvar(T, idVar, timeVar, dVar)
%
% Output:
%   - T2  : T with a new column 'gVar' giving the calendar time of first
%           treatment for each unit in every row; 0 for never-treated.
%   - gVar: the string "gVar" (the name of the new column).
%
% Notes:
%   * Works for numeric or datetime time variables. For datetime, we store
%     the YEAR number in gVar (e.g., 2006). If you prefer full datetimes,
%     you can adjust the casting below.

T2 = T;
gVar = "gVar";

% Stable integer time index
[~,~,t_int] = unique(T2.(timeVar), 'stable');
t_int = double(t_int);
t_vals = unique(T2.(timeVar), 'stable');

% Per-unit: first index with D==1
[gid, ~] = findgroups(T2.(idVar));
D = T2.(dVar);
if islogical(D), D = double(D); end
if ~isfloat(D),  D = double(D); end

first_idx = splitapply(@(tt,dd) local_first_on_idx(tt,dd), t_int, D, gid); % NaN if never
first_idx_per_row = first_idx(gid);  % align to rows

% Map first-on index -> time label (numeric)
g_label = zeros(height(T2),1);                 % default 0 = never
mask = ~isnan(first_idx_per_row);
if any(mask)
    idx = first_idx_per_row(mask);
    if isnumeric(t_vals)
        g_label(mask) = t_vals(idx);
    elseif isdatetime(t_vals)
        % store year number to keep 'gVar' numeric as required by estimators
        g_label(mask) = year(t_vals(idx));
    else
        % Strings/categorical -> fallback to position index mapped onto 1..T,
        % but keep never-treated as 0 to match estimator expectations.
        g_label(mask) = double(idx);
    end
end

T2.(gVar) = g_label;

% (No temp columns kept; function is non-destructive to original T.)
end

% ---- local helper ----
function a = local_first_on_idx(ti, d)
% Return the first time index where d==1; NaN if never
[ts, ix] = sort(ti);
ds = d(ix);
j = find(ds==1, 1, 'first');
if isempty(j), a = NaN; else, a = ts(j); end
end
