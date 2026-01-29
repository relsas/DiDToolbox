function S = did_standardize(res, varargin)
% Map various estimator outputs -> unified schema:
%   S.es        : table(e, Estimate, [SE|LB|UB])
%   S.byCohort  : table(g, Estimate, [SE|LB|UB])
%   S.calendar  : table(t, Estimate, [SE|LB|UB])
%   S.support   : table(k, nObs|nCohorts)
%   S.chCohortLines : (CH only) table(t, DID_plus, DID_minus)

ip = inputParser;
addParameter(ip,'Alpha',0.95,@(x) isnumeric(x) && x>0 && x<1);
parse(ip,varargin{:});
Alpha = ip.Results.Alpha;
crit = norminv(0.5+Alpha/2);

S = struct('es',[],'byCohort',[],'calendar',[],'overall',[], ...
    'support',[],'chCohortLines',[],'V_es',[]);

% --------- IW (HGardner IW) ----------
% Accept Estimator tags like "SA", "SA_IW", and prefer res.es / res.Aggregates.es
if isfield(res,'Estimator') && any(contains(upper(string(res.Estimator)), ["IW","Gardner"]))
    % Event-time
    if isfield(res,'es') && istable(res.es) && ~isempty(res.es)
        S.es = ensure_bands_(rename_if_needed_(res.es,'e','Estimate'), crit);
    elseif isfield(res,'Aggregates') && isfield(res.Aggregates,'es') && istable(res.Aggregates.es)
        S.es = ensure_bands_(rename_if_needed_(res.Aggregates.es,'e','Estimate'), crit);
    end
    % By-cohort
    if isfield(res,'byCohort') && istable(res.byCohort) && ~isempty(res.byCohort)
        S.byCohort = ensure_bands_(rename_if_needed_(res.byCohort,'g','Estimate'), crit);
    elseif isfield(res,'Aggregates') && isfield(res.Aggregates,'byCohort') && istable(res.Aggregates.byCohort)
        S.byCohort = ensure_bands_(rename_if_needed_(res.Aggregates.byCohort,'g','Estimate'), crit);
    end
    % Calendar-time
    if isfield(res,'calendar') && istable(res.calendar) && ~isempty(res.calendar)
        S.calendar = ensure_bands_(rename_if_needed_(res.calendar,'t','Estimate'), crit);
    elseif isfield(res,'Aggregates') && isfield(res.Aggregates,'calendar') && istable(res.Aggregates.calendar)
        S.calendar = ensure_bands_(rename_if_needed_(res.Aggregates.calendar,'t','Estimate'), crit);
    end
    % Overall (if present)
    if isfield(res,'overall') && isstruct(res.overall) && isfield(res.overall,'Estimate')
        S.overall = res.overall;
    end
    return
end

% --------- CS ----------
if isfield(res,'Estimator') && any(strcmpi(string(res.Estimator),["CS2021","CS"]))
    if isfield(res,'Aggregates')
        A = res.Aggregates;
        if istable(A.es),       S.es       = ensure_bands_(rename_if_needed_(A.es,'e','Estimate'), crit); end
        if istable(A.byCohort), S.byCohort = ensure_bands_(rename_if_needed_(A.byCohort,'g','Estimate'), crit); end
        if istable(A.calendar), S.calendar = ensure_bands_(rename_if_needed_(A.calendar,'t','Estimate'), crit); end
    end
    if isfield(res,'details') && isfield(res.details,'supportByEventTime')
        S.support = res.details.supportByEventTime;
        if ~ismember('k', S.support.Properties.VariableNames)
            S.support.Properties.VariableNames{1} = 'k';
        end
    end
    return
end

% --------- Wooldridge_TB ----------
if isfield(res,'details')
    D = res.details;
    if isfield(D,'attByEventTime') && istable(D.attByEventTime)
        T = D.attByEventTime;
        % Map EventTime -> e
        if ismember('EventTime', T.Properties.VariableNames)
            T.Properties.VariableNames{'EventTime'} = 'e';
        else
            T.Properties.VariableNames{1} = 'e';
        end

        % Map Estimate
        if ismember('Estimate', T.Properties.VariableNames)
            % Good
        elseif ismember('ATT_hat_mean', T.Properties.VariableNames)
            T.Estimate = T.ATT_hat_mean;
        end

        S.es = table(T.e, T.Estimate, 'VariableNames', {'e','Estimate'});
        if ismember('SE', T.Properties.VariableNames)
            S.es.SE = T.SE;
            S.es = ensure_bands_(S.es, crit);
        end
    end
    if isfield(D,'ATTbyCohort') && istable(D.ATTbyCohort)
        T = D.ATTbyCohort; % expects Cohort, ATT(k), [SE]
        if ismember('Cohort', T.Properties.VariableNames) && ismember('ATT(k)', T.Properties.VariableNames)
            S.byCohort = table(T.Cohort, T.("ATT(k)"), 'VariableNames', {'g','Estimate'});
            if ismember('SE', T.Properties.VariableNames)
                S.byCohort.SE = T.SE; S.byCohort = ensure_bands_(S.byCohort, crit);
            end
        end
    end
    if isfield(D,'supportByEventTime') && istable(D.supportByEventTime)
        S.support = D.supportByEventTime;
        if ~ismember('k', S.support.Properties.VariableNames)
            S.support.Properties.VariableNames{1} = 'k';
        end
    end
    if isfield(D,'V_es') && ~isempty(D.V_es)
        S.V_es = D.V_es;
    end
    return
end

% --------- CH (Chaisemartin–Haultfoeuille) ----------
if isfield(res,'Estimator') && any(strcmpi(string(res.Estimator),["CH2020","CH"]))
    % Plot only ATT_by_Cohort (time lines for DID_plus/-minus).
    if isfield(res,'ATT_by_Cohort') && istable(res.ATT_by_Cohort)
        U = res.ATT_by_Cohort;
        % Expect columns: t, DID_plus, DID_minus (some may be missing)
        cols = U.Properties.VariableNames;
        if ~ismember('t', cols)
            % If the first column is time, rename it
            U.Properties.VariableNames{1} = 't';
        end
        S.chCohortLines = U(:, intersect({'t','DID_plus','DID_minus'}, U.Properties.VariableNames));
    end
    return
end

% --------- BJS ----------
if isfield(res,'ATT_by_horizon') || isfield(res,'ATT_overall') || isfield(res,'ATT_cohort_obs')
    if isfield(res,'ATT_by_horizon') && istable(res.ATT_by_horizon)
        H = res.ATT_by_horizon;
        S.es = table(H.k, H.ATT_k, 'VariableNames', {'e','Estimate'});
        if ismember('SE', H.Properties.VariableNames), S.es.SE = H.SE; end
        if all(ismember({'CI_lo','CI_hi'}, H.Properties.VariableNames))
            S.es.LB = H.CI_lo; S.es.UB = H.CI_hi;
        else
            S.es = ensure_bands_(S.es, crit);
        end
        if ismember('N_k', H.Properties.VariableNames)
            S.support = table(H.k, H.N_k, 'VariableNames', {'k','nObs'});
        end
    end
    if isfield(res,'ATT_cohort_obs') && istable(res.ATT_cohort_obs) && ~isempty(res.ATT_cohort_obs)
        C = res.ATT_cohort_obs;
        if all(ismember({'cohort','ATT_obs'}, C.Properties.VariableNames))
            S.byCohort = table(C.cohort, C.ATT_obs, 'VariableNames', {'g','Estimate'});
            if ismember('SE_obs', C.Properties.VariableNames), S.byCohort.SE = C.SE_obs; end
            if all(ismember({'CI_lo_obs','CI_hi_obs'}, C.Properties.VariableNames))
                S.byCohort.LB = C.CI_lo_obs; S.byCohort.UB = C.CI_hi_obs;
            else
                S.byCohort = ensure_bands_(S.byCohort, crit);
            end
        end
    end
    return
end

% --------- Generic fallback (if any) ----------
if isfield(res, 'es')       && istable(res.es),        S.es       = ensure_bands_(res.es, crit); end
if isfield(res, 'byCohort') && istable(res.byCohort),  S.byCohort = ensure_bands_(res.byCohort, crit); end
if isfield(res, 'calendar') && istable(res.calendar),  S.calendar = ensure_bands_(res.calendar, crit); end
if isfield(res, 'support')  && istable(res.support),   S.support  = res.support; end
end

%% ---- Utilities ----
function T = rename_if_needed_(T, keyCol, estCol)
% Ensure T has a specific index column name (keyCol) and an Estimate column
if ~istable(T) || isempty(T), return; end
vn = T.Properties.VariableNames;

% Key column: if missing, rename first column to keyCol
if ~ismember(keyCol, vn) && ~isempty(vn)
    T.Properties.VariableNames{1} = keyCol;
    vn = T.Properties.VariableNames;
end

% Estimate column: if 'Estimate' missing, map a common candidate
if nargin<2 || isempty(estCol), estCol = 'Estimate'; end
vn = T.Properties.VariableNames;
if ~ismember('Estimate', vn)
    cand = intersect(vn, {'theta','ATT','ATT_hat_mean','estimate','ATT_k','Estimate'});
    if ~isempty(cand)
        T.Estimate = T.(cand{1});
    end
end
end

function T = ensure_bands_(T, crit)
if isempty(T) || ~istable(T), return; end
vn = T.Properties.VariableNames;
hasLB = ismember('LB',vn); hasUB = ismember('UB',vn);
if (~hasLB || ~hasUB) && ismember('SE',vn) && ismember('Estimate',vn)
    T.LB = T.Estimate - crit*T.SE;
    T.UB = T.Estimate + crit*T.SE;
end
end
