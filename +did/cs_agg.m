function [summaryTable, Agg] = cs_agg(summaryTable_in, ATT, g_list, t_list, G, att_b, Sinfo, varargin)
% CS_AGG Aggregation logic for Callaway & Sant'Anna estimator
%
% [summaryTable, Agg] = did.cs_agg(summaryTable_in, ATT, g_list, t_list, G, att_b, Sinfo, 'Name', Val...)

ip = inputParser;
addParameter(ip,'TimeValues',[]);
addParameter(ip,'Delta',0);
% NEW params:
addParameter(ip,'CellWeights',[]);         % length R: Wtreat per (g,t)
addParameter(ip,'Weighting',"cohortShare");
parse(ip,varargin{:});
TimeValues = ip.Results.TimeValues;
Delta      = ip.Results.Delta;
CellWts    = ip.Results.CellWeights(:);
Weighting  = lower(string(ip.Results.Weighting));

labelFrom = @(idx) map_labels_(idx, TimeValues);

summaryTable = summaryTable_in;
summaryTable.gYear   = labelFrom(summaryTable.g);
summaryTable.tYear   = labelFrom(summaryTable.t);
summaryTable.refYear = labelFrom(summaryTable.g - Delta - 1);

% Filter out reference period rows (Estimate=0, SE=NaN usually)
% (DISABLED: We need the reference row for Honest DiD anchoring)
% if any(isRef)
%     summaryTable(isRef, :) = [];
% end

% cohort-share weights (existing behavior)
Gfinite = G(~isinf(G));
[uG,~,ic] = unique(Gfinite,'stable');
countG = accumarray(ic, 1, [numel(uG),1], @sum, 0);
wG = countG / max(sum(countG),1);
wGmap = containers.Map(num2cell(uG), num2cell(wG));
wg = @(gv) arrayfun(@(x) wGmap(x), gv);

% --- NEW: select a weight accessor w_fun(I) that returns a vector of weights per subset of cells I
% * cohortShare: w_fun(I) = wg(g_list(I))
% * treatedObs : w_fun(I) = CellWts(I)
w_fun = @(I) wg(g_list(I));   % default
if Weighting == "treatedobs"
    if isempty(CellWts) || numel(CellWts) ~= numel(ATT)
        warning('did:cs:treatedObsWeightsMissing','CellWeights missing or wrong length; falling back to cohortShare.');
    else
        w_fun = @(I) CellWts(I);
    end
end

% event-time  (NOW: pre *and* post)
e_list   = t_list - g_list;

%%%% CHANGED: include all event times, so e<0 show up as genuine pre-trends.
Evals    = unique(e_list);

theta_es = NaN(numel(Evals),1);
for k=1:numel(Evals)
    ek = Evals(k);
    I  = find(e_list==ek);
    if ~isempty(I)
        w = w_fun(I);
        theta_es(k) = nansum( w .* ATT(I) ) / max(nansum(w), eps);
    end
end

% calendar-time
Tuniq   = unique(t_list);
theta_c = NaN(numel(Tuniq),1);
for k=1:numel(Tuniq)
    tk = Tuniq(k);
    I  = find(t_list==tk & g_list<=tk);
    if ~isempty(I)
        w = w_fun(I);
        theta_c(k) = nansum( w .* ATT(I) ) / max(nansum(w), eps);
    end
end

% by-cohort (keep your original simple mean; can be swapped to weighted if desired)
Guniq    = unique(g_list);
theta_g  = NaN(numel(Guniq),1);
for k=1:numel(Guniq)
    gk = Guniq(k);
    I  = find(g_list==gk & t_list>=gk);
    theta_g(k) = mean(ATT(I),'omitnan');
end

% overall
pick_post   = (t_list>=g_list);
I_all       = find(pick_post);
w_all_cells = w_fun(I_all);
theta_OW = nansum( w_all_cells .* ATT(I_all) ) / max(nansum(w_all_cells),eps);

% inference for aggregates
if ~isempty(att_b)
    [es_tab, es_crit]   = agg_infer_(theta_es, build_draws_(Evals,  e_list, g_list, t_list, att_b, w_fun, 'event'));
    [cal_tab, cal_crit] = agg_infer_(theta_c,  build_draws_(Tuniq,  t_list, g_list, t_list, att_b, w_fun, 'calendar'));
    [sel_tab, sel_crit] = agg_infer_(theta_g,  build_draws_(Guniq,  g_list, t_list, t_list, att_b, w_fun, 'cohort'));
    [ow_tab,  ow_crit]  = agg_infer_(theta_OW, build_draws_(NaN,    [],     g_list, t_list, att_b, w_fun, 'overall', pick_post));
    Agg.es       = addvars(table(Evals,'VariableNames',{'e'}), es_tab.Estimate, es_tab.SE, es_tab.tStat, es_tab.pValue, es_tab.LB, es_tab.UB, ...
        'NewVariableNames',{'Estimate','SE','tStat','pValue','LB','UB'});
    Agg.calendar = addvars(table(Tuniq,'VariableNames',{'t'}), cal_tab.Estimate, cal_tab.SE, cal_tab.tStat, cal_tab.pValue, cal_tab.LB, cal_tab.UB, ...
        'NewVariableNames',{'Estimate','SE','tStat','pValue','LB','UB'});
    Agg.byCohort = addvars(table(Guniq,'VariableNames',{'g'}), sel_tab.Estimate, sel_tab.SE, sel_tab.tStat, sel_tab.pValue, sel_tab.LB, sel_tab.UB, ...
        'NewVariableNames',{'Estimate','SE','tStat','pValue','LB','UB'});
    Agg.overall  = ow_tab;
    Agg.Crit     = struct('es',es_crit,'calendar',cal_crit,'byCohort',sel_crit,'overall',ow_crit);
else
    if ~isempty(Sinfo)
        W_es   = weight_matrix_(Evals,  e_list, g_list, t_list, w_fun, 'event');
        W_cal  = weight_matrix_(Tuniq,  t_list, g_list, t_list, w_fun, 'calendar');
        W_sel  = weight_matrix_(Guniq,  g_list, t_list, t_list, w_fun, 'cohort');
        W_over = weight_matrix_(NaN,    [],     g_list, t_list, w_fun, 'overall', pick_post);

        [SE_es, t_es, p_es, LB_es, UB_es] = var_from_S_lincombo_(theta_es, W_es,  Sinfo);
        [SE_c,  t_c,  p_c,  LB_c,  UB_c ] = var_from_S_lincombo_(theta_c,  W_cal, Sinfo);
        [SE_g,  t_g,  p_g,  LB_g,  UB_g ] = var_from_S_lincombo_(theta_g,  W_sel, Sinfo);
        [SE_ow, t_ow, p_ow, LB_ow, UB_ow] = var_from_S_lincombo_(theta_OW, W_over, Sinfo);

        Agg.es       = table(Evals, theta_es, SE_es, t_es, p_es, LB_es, UB_es, ...
            'VariableNames',{'e','Estimate','SE','tStat','pValue','LB','UB'});
        Agg.calendar = table(Tuniq, theta_c,  SE_c,  t_c,  p_c,  LB_c,  UB_c, ...
            'VariableNames',{'t','Estimate','SE','tStat','pValue','LB','UB'});
        Agg.byCohort = table(Guniq, theta_g,  SE_g,  t_g,  p_g,  LB_g,  UB_g, ...
            'VariableNames',{'g','Estimate','SE','tStat','pValue','LB','UB'});
        Agg.overall  = struct('Estimate',theta_OW,'SE',SE_ow,'tStat',t_ow,'pValue',p_ow,'LB',LB_ow,'UB',UB_ow,'crit',1.96);
        Agg.Crit     = struct('es',1.96,'calendar',1.96,'byCohort',1.96,'overall',1.96);
    else
        Agg.es       = table(Evals, theta_es, 'VariableNames',{'e','Estimate'});
        Agg.calendar = table(Tuniq, theta_c,  'VariableNames',{'t','Estimate'});
        Agg.byCohort = table(Guniq, theta_g,  'VariableNames',{'g','Estimate'});
        Agg.overall  = struct('Estimate',theta_OW,'SE',NaN,'tStat',NaN,'pValue',NaN,'LB',NaN,'UB',NaN,'crit',NaN);
        Agg.Crit     = struct('es',NaN,'calendar',NaN,'byCohort',NaN,'overall',NaN);
    end
end

% labels on aggregates
if istable(Agg.byCohort) && ismember('g', Agg.byCohort.Properties.VariableNames)
    Agg.byCohort.gYear = labelFrom(Agg.byCohort.g);
end
if istable(Agg.calendar) && ismember('t', Agg.calendar.Properties.VariableNames)
    Agg.calendar.tYear = labelFrom(Agg.calendar.t);
end
end

% ====== nested helpers (scoped) ======
function outL = map_labels_(idx, L)
if isempty(L), outL = idx; return; end
L = L(:);
outL = repmat(missingLike_(L), size(idx));
valid = idx >= 1 & idx <= numel(L) & isfinite(idx);
if any(valid)
    if isdatetime(L)
        tmp = NaT(size(idx)); tmp(valid) = L(idx(valid)); outL = tmp;
    elseif isstring(L) || iscellstr(L)
        tmp = strings(size(idx)); tmp(valid) = string(L(idx(valid))); outL = tmp;
    elseif iscategorical(L)
        tmp = categorical(missing(size(idx))); %#ok<CTPCT>
        tmp(valid) = L(idx(valid)); outL = tmp;
    else
        tmp = NaN(size(idx)); tmp(valid) = double(L(idx(valid))); outL = tmp;
    end
end
end

function m = missingLike_(L)
if isdatetime(L),      m = NaT;
elseif isstring(L),    m = string(missing);
elseif iscategorical(L), m = categorical(missing);
else,                  m = NaN;
end
end


function A = weight_matrix_(axisVals, a_list, b_list, t_list, w_fun, kind, pick_post)
% Builds linear-combo matrix A so that theta = A * ATT (for variance via S)
% a_list / b_list correspond to the first / second index that define cells.
R = numel(b_list);  %#ok<NASGU> % (only used implicitly through find on conditions)
if nargin==7 && strcmpi(string(kind),'overall')
    % overall: single row
    A = zeros(1, numel(t_list));
    I = find(pick_post);
    if ~isempty(I)
        w = w_fun(I);
        A(1, I) = w / max(sum(w), eps);
    end
    return
end

knd = "event"; if nargin>=6 && ~isempty(kind), knd = lower(string(kind)); end
Rcells = numel(t_list);
A = zeros(0, Rcells); % will set below

switch knd
    case "calendar"
        Tuniq = axisVals;
        A = zeros(numel(Tuniq), Rcells);
        for k=1:numel(Tuniq)
            tk = Tuniq(k);
            I  = find(t_list==tk & b_list<=tk);
            if isempty(I), continue; end
            w = w_fun(I);
            A(k,I) = w / max(sum(w), eps);
        end
    case "cohort"
        Guniq = axisVals;
        A = zeros(numel(Guniq), Rcells);
        for k=1:numel(Guniq)
            gk = Guniq(k);
            I  = find(b_list==gk & t_list>=gk);
            if isempty(I), continue; end
            % keep simple average for "byCohort" as in original
            A(k,I) = 1/numel(I);
        end
    otherwise % "event"
        Evals  = axisVals;
        e_list = t_list - b_list;
        A = zeros(numel(Evals), Rcells);
        for k=1:numel(Evals)
            ek = Evals(k);
            I  = find(e_list==ek);
            if isempty(I), continue; end
            w = w_fun(I);
            A(k,I) = w / max(sum(w), eps);
        end
end
end

function draws = build_draws_(axisVals, a_list, b_list, t_list, att_b, w_fun, kind, pick_post)
B = size(att_b,2);

% overall (single-number) case
if nargin>=8 && strcmpi(string(kind),'overall')
    draws = NaN(1,B);
    I = find(pick_post);
    if isempty(I), return; end
    w = w_fun(I);
    for b=1:B
        draws(1,b) = nansum(w .* att_b(I,b)) / max(sum(w), eps);
    end
    return
end

knd = "event";
if nargin>=7 && ~isempty(kind), knd = lower(string(kind)); end

switch knd
    case "calendar"
        Tuniq = axisVals;
        draws = NaN(numel(Tuniq), B);
        for k=1:numel(Tuniq)
            tk = Tuniq(k);
            I  = find(t_list==tk & b_list<=tk);
            if isempty(I), continue; end
            w = w_fun(I);
            for b=1:B
                draws(k,b) = nansum(w .* att_b(I,b)) / max(sum(w), eps);
            end
        end

    case "cohort"
        Guniq = axisVals;
        draws = NaN(numel(Guniq), B);
        for k=1:numel(Guniq)
            gk = Guniq(k);
            I  = find(b_list==gk & t_list>=gk);
            if isempty(I), continue; end
            for b=1:B
                draws(k,b) = mean(att_b(I,b), 'omitnan');  % keep simple average (as original)
            end
        end

    otherwise  % "event"
        Evals  = axisVals;
        e_list = t_list - b_list;   % event-time offsets
        draws  = NaN(numel(Evals), B);
        for k=1:numel(Evals)
            ek = Evals(k);
            I  = find(e_list==ek);
            if isempty(I), continue; end
            w = w_fun(I);
            for b=1:B
                draws(k,b) = nansum(w .* att_b(I,b)) / max(sum(w), eps);
            end
        end
end
end

function [SEv, tSv, pV, LB, UB] = var_from_S_lincombo_(theta, A, Sinfo)
K = size(A,1); SEv = NaN(K,1); tSv = NaN(K,1); pV = NaN(K,1); LB = NaN(K,1); UB = NaN(K,1);
if isempty(Sinfo), return; end
crit = 1.96;
if Sinfo.kind=="clustered"
    S1 = Sinfo.S1; C1 = Sinfo.C1;
    S1c = S1 - mean(S1,2);
    for j=1:K
        a = A(j,:).';
        s = (a.' * S1c);
        Var = (C1/(C1-1)) * sum(s.^2);
        se = sqrt(Var); SEv(j)=se;
        if isfinite(se) && se>0
            tSv(j) = theta(j)/se;
            pV(j)  = 2*(1-normcdf(abs(tSv(j))));
            LB(j)  = theta(j) - crit*se;
            UB(j)  = theta(j) + crit*se;
        end
    end
else
    S1 = Sinfo.S1; C1 = Sinfo.C1;
    S2 = Sinfo.S2; C2 = Sinfo.C2;
    S12= Sinfo.S12; C12= Sinfo.C12;
    S1c  = S1  - mean(S1, 2);
    S2c  = S2  - mean(S2, 2);
    S12c = S12 - mean(S12,2);
    for j=1:K
        a = A(j,:).';
        s1  = (a.' * S1c);   V1  = (C1 /(C1 -1)) * sum(s1.^2);
        s2  = (a.' * S2c);   V2  = (C2 /(C2 -1)) * sum(s2.^2);
        s12 = (a.' * S12c);  V12 = (C12/(C12-1)) * sum(s12.^2);
        Var = V1 + V2 - V12;
        se = sqrt(Var); SEv(j)=se;
        if isfinite(se) && se>0
            tSv(j) = theta(j)/se;
            pV(j)  = 2*(1-normcdf(abs(tSv(j))));
            LB(j)  = theta(j) - crit*se;
            UB(j)  = theta(j) + crit*se;
        end
    end
end
end

function [tab, crit] = agg_infer_(theta, draws)
K = numel(theta);
if isempty(draws)
    tab = table(theta, NaN(K,1), NaN(K,1), NaN(K,1), NaN(K,1), NaN(K,1), ...
        'VariableNames', {'Estimate','SE','tStat','pValue','LB','UB'});
    crit = NaN; return;
end
dc = draws - mean(draws,2);
SE = std(dc,0,2);
tStat = theta ./ max(SE, eps);
Z = abs((draws - theta) ./ max(SE, eps));
thr = abs(tStat);
pValue = mean(Z >= thr, 2);
maxZ = max(Z, [], 1);
crit = did.quantile(maxZ, 0.95);
LB = theta - crit .* SE;
UB = theta + crit .* SE;
tab = table(theta, SE, tStat, pValue, LB, UB, ...
    'VariableNames', {'Estimate','SE','tStat','pValue','LB','UB'});
end
