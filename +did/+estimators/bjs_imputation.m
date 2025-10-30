function out = bjs_imputation(T, opts)
% BJS_IMPUTATION  Borusyak–Jaravel–Spiess (2024, RES) imputation estimator
%   out = did.bjs_imputation(T, idVar="id", timeVar="time", yVar="y", dVar="D", ...)
%
% Required (name–value inside 'opts' via arguments block):
%   idVar, timeVar, yVar, dVar
%
% Optional:
%   Covariates : string array of covariate names (default = [])
%   Horizons   : double vector of event times k (default = all observed among treated)
%   SEMethod   : "LOO" (default), "BootstrapUnit", or "None"
%   BootReps   : integer >=0 (used if SEMethod="BootstrapUnit")
%   Seed       : RNG seed for bootstrap (default = NaN)
%   Display    : logical, print brief summary (default = true)
%
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/03/2025
% ------------------------------------------------------------------------
arguments
    T table
    opts.idVar (1,1) string
    opts.timeVar (1,1) string
    opts.yVar (1,1) string
    opts.dVar (1,1) string
    opts.Covariates string = string.empty(1,0)
    opts.Horizons double = []
    opts.SEMethod (1,1) string {mustBeMember(opts.SEMethod,["LOO","BootstrapUnit","None"])} = "LOO"
    opts.BootReps (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.Seed (1,1) double = randi([1,1e7],1,1)
    opts.Display (1,1) logical = true
    opts.useParallel double = 1
end

idVar   = opts.idVar;
timeVar = opts.timeVar;
yVar    = opts.yVar;
dVar    = opts.dVar;
covars  = opts.Covariates;

% ----- Basic checks
varsNeeded = [idVar, timeVar, yVar, dVar, covars];
missingVars = varsNeeded(~ismember(varsNeeded, string(T.Properties.VariableNames)));
if ~isempty(missingVars)
    error('Missing variables in T: %s', strjoin(missingVars, ', '));
end
if ~isnumeric(T.(dVar)) && ~islogical(T.(dVar))
    error('Treatment variable %s must be numeric/logical 0/1.', dVar);
end

% ----- Working copy and event-time scaffolding
T2 = T;
[~,~,t_idx] = unique(T2.(timeVar), 'stable');
T2.t_int = double(t_idx);

% Categorical id & time for FE regression
T2.(idVar)   = categorical(T2.(idVar));
T2.(timeVar) = categorical(T2.(timeVar));

% Untreated indicator
D  = T2.(dVar)==1;
U0 = ~D;  % ALL untreated: never-treated + not-yet-treated
if ~any(U0)
    error('No untreated observations (D==0). Cannot fit untreated outcome model.');
end
T2.D = D;                 % used later in cohort summary

% ----- Build regression formula: y ~ 1 + id + time + covariates
useCov = ~isempty(covars);
if useCov
    rhs = "1 + " + idVar + " + " + timeVar + " + " + strjoin(covars, " + ");
else
    rhs = "1 + " + idVar + " + " + timeVar;
end
form = yVar + " ~ " + rhs;

% Which covariates are categorical?
isCatCov = false(size(covars));
for c = 1:numel(covars)
    isCatCov(c) = iscategorical(T2.(covars(c)));
    if isCatCov(c)
        T2.(covars(c)) = categorical(T2.(covars(c)));
    end
end
catVars = [idVar, timeVar, covars(isCatCov)];

% ======================
% Step A: fit on ALL D==0 with cleaned categories
% ======================
T_unt = T2(U0,:);
T_unt.(idVar)   = removecats(T_unt.(idVar));
T_unt.(timeVar) = removecats(T_unt.(timeVar));
for c = find(isCatCov)
    T_unt.(covars(c)) = removecats(T_unt.(covars(c)));
end
mdl = fitlm(T_unt, char(form), 'CategoricalVars', cellstr(catVars));

% ======================
% Step B: predict only on identified rows; compute tau only when treated & identified
% ======================
idOK   = ismember(T2.(idVar),   categories(T_unt.(idVar)));
timeOK = ismember(T2.(timeVar), categories(T_unt.(timeVar)));
predOK = idOK & timeOK;
if any(isCatCov)
    for c = find(isCatCov)
        predOK = predOK & ismember(T2.(covars(c)), categories(T_unt.(covars(c))));
    end
end
T2.predOK = predOK;

Y = T2.(yVar);
Y0hat = NaN(height(T2),1);
Y0hat(predOK) = predict(mdl, T2(predOK,:));

tauHat = NaN(height(T2),1);
tauHat(D & predOK) = Y(D & predOK) - Y0hat(D & predOK);

% ======================
% Event time & adoption time
% ======================
if any(D)
    adoptTbl = groupsummary(T2(D,:), idVar, "min", "t_int");
    idxVar = strcmpi(adoptTbl.Properties.VariableNames, 'min_t_int');
    if any(idxVar)
        adoptTbl.Properties.VariableNames(idxVar) = {'adoptTime'};
    else
        adoptTbl.Properties.VariableNames(width(adoptTbl)) = {'adoptTime'};
    end
    T2 = outerjoin(T2, adoptTbl(:,[idVar,"adoptTime"]), "Keys", idVar, "MergeKeys", true, "Type", "left");
else
    T2.adoptTime = NaN(height(T2),1);
end

k = NaN(height(T2),1);
hasG = ~isnan(T2.adoptTime);
k(hasG) = T2.t_int(hasG) - T2.adoptTime(hasG);
T2.k = k;
T2.cohort = T2.adoptTime;   % <<< minimal fix 
% ======================
% Horizons to report
% ======================
if isempty(opts.Horizons)
    kObs = unique(T2.k(D & ~isnan(T2.k)));
    Horizons = sort(kObs(:))';
else
    Horizons = unique(sort(opts.Horizons(:)'));
end

% ======================
% Aggregation: overall ATT and ATT(k)
% ======================
N_treated = nnz(D);
N_treated_ident = nnz(D & predOK & ~isnan(tauHat));

ATT_overall_point = NaN;
if N_treated>0
    ATT_overall_point = mean(tauHat(D), 'omitnan');
end

ATTk = nan(size(Horizons));
Nk   = zeros(size(Horizons));
for j = 1:numel(Horizons)
    Kj = Horizons(j);
    idx = D & (T2.k==Kj);
    Nk(j) = nnz(idx);
    if Nk(j)>0
        ATTk(j) = mean(tauHat(idx), 'omitnan');
    end
end

T2.tauHat = tauHat;


% ======================
% Inference
% ======================
SE_overall = NaN; CI_overall = [NaN NaN];
SE_k = NaN(size(Horizons)); CI_k = nan(numel(Horizons),2);
SE_k_bal = []; CI_k_bal = [];

uid = categories(T2.(idVar));
G   = numel(uid);

% ---- (A) LOO jackknife (default) ----
if opts.SEMethod=="LOO"
    if G<2
        warning('SEMethod="LOO" requires at least 2 units. SEs will be NaN.');
    else
        att_overall_lo = nan(G,1);
        attk_lo        = nan(G, numel(Horizons));
        attkbal_lo     = nan(G, numel(Horizons));
        hasBalInfo = false;

        % parfor for LOO
        attk_temp=cell(G,1);
        att_temp=cell(G,1);
        pool = gcp('nocreate');
        if ~isempty(opts.useParallel) && isempty(pool)
            parpool(opts.useParallel);
            parfor g = 1:G
                [attk_temp{g},att_temp{g}] = calcATTk(T2,idVar,timeVar,yVar,D, isCatCov, form, catVars, Horizons,uid(g));
            end
        elseif ~isempty(pool)
            parfor g = 1:G
                [attk_temp{g},att_temp{g}] = calcATTk(T2,idVar,timeVar,yVar,D, isCatCov, form, catVars, Horizons,uid(g));
            end
        else
            for g = 1:G
                [attk_temp{g},att_temp{g}] = calcATTk(T2,idVar,timeVar,yVar,D, isCatCov, form, catVars, Horizons,uid(g));
            end
        end

      
        attk_lo = vertcat(attk_temp{:});
        att_overall_lo = vertcat(att_temp{:});
        jack_var = @(v) jackknife_var_nan(v);

        v = att_overall_lo;
        JV = jack_var(v);
        if ~isnan(JV)
            SE_overall = sqrt(JV);
            CI_overall = [ATT_overall_point - 1.96*SE_overall, ATT_overall_point + 1.96*SE_overall];
        end

        for j = 1:numel(Horizons)
            v = attk_lo(:,j);
            JV = jack_var(v);
            if ~isnan(JV)
                SE_k(j) = sqrt(JV);
                CI_k(j,:) = [ATTk(j) - 1.96*SE_k(j), ATTk(j) + 1.96*SE_k(j)];
            end
        end

      end
end

% ---- (B) Unit bootstrap (optional) ----
if opts.SEMethod=="BootstrapUnit" && opts.BootReps>0
    if ~isnan(opts.Seed), rng(opts.Seed); end
    att_overall_b = nan(opts.BootReps,1);
    attk_b        = nan(opts.BootReps, numel(Horizons));
    
    for b = 1:opts.BootReps
        draw = randsample(G, G, true);
        selIds = uid(draw);
        idxB = ismember(T2.(idVar), selIds);
        TB = T2(idxB,:);
        DB = D(idxB);

        isUntB = ~DB;
        if ~any(isUntB) || ~any(DB), continue; end

        TB_unt = TB(isUntB,:);
        TB_unt.(idVar)   = removecats(TB_unt.(idVar));
        TB_unt.(timeVar) = removecats(TB_unt.(timeVar));
        for c = find(isCatCov), TB_unt.(covars(c)) = removecats(TB_unt.(covars(c))); end
        mdlB = fitlm(TB_unt, char(form), 'CategoricalVars', cellstr(catVars));

        idOK_b   = ismember(TB.(idVar),   categories(TB_unt.(idVar)));
        timeOK_b = ismember(TB.(timeVar), categories(TB_unt.(timeVar)));
        predOK_b = idOK_b & timeOK_b;
        if any(isCatCov)
            for c = find(isCatCov)
                predOK_b = predOK_b & ismember(TB.(covars(c)), categories(TB_unt.(covars(c))));
            end
        end

        Y0b = NaN(height(TB),1);
        Y0b(predOK_b) = predict(mdlB, TB(predOK_b,:));
        tauB = NaN(height(TB),1);
        tauB(DB & predOK_b) = TB.(yVar)(DB & predOK_b) - Y0b(DB & predOK_b);

        att_overall_b(b) = mean(tauB(DB), 'omitnan');

        for j = 1:numel(Horizons)
            Kj = Horizons(j);
            idx = DB & (TB.k==Kj);
            if any(idx)
                attk_b(b,j) = mean(tauB(idx), 'omitnan');
            end
        end

       
    end

    SE_overall = std(att_overall_b, 'omitnan');
    if ~all(isnan(att_overall_b)), CI_overall = prctile(att_overall_b, [2.5 97.5]); end
    for j = 1:numel(Horizons)
        SE_k(j) = std(attk_b(:,j), 'omitnan');
        if ~all(isnan(attk_b(:,j))), CI_k(j,:) = prctile(attk_b(:,j), [2.5 97.5]); end
    end
    
end

% ======================
% Pack outputs
% ======================
ATT_overall = table(N_treated, N_treated_ident, ATT_overall_point, SE_overall, CI_overall(:,1), CI_overall(:,2), ...
    'VariableNames', ["N_treated","N_treated_ident","ATT","SE","CI_lo","CI_hi"]);

ATT_by_horizon = table(Horizons(:), Nk(:), ATTk(:), SE_k(:), CI_k(:,1), CI_k(:,2), ...
    'VariableNames', ["k","N_k","ATT_k","SE","CI_lo","CI_hi"]);



%% ATT by cohort (observation–weighted over treated cells)
W = T2;
W = W(W.D==1 & ~isnan(W.cohort) & ~isnan(W.k), :);

att_mean = groupsummary(W, "cohort", "mean", "tauHat");
att_mean.Properties.VariableNames(end) = "ATT_obs";
att_sd   = groupsummary(W, "cohort", "std",  "tauHat");
att_sd.Properties.VariableNames(end)   = "SD_obs";
att_n    = groupcounts(W, "cohort");
att_n.Properties.VariableNames(end)    = "N_obs";

ATT_cohort_obs = outerjoin(att_mean(:,["cohort","ATT_obs"]), att_sd(:,["cohort","SD_obs"]), ...
    "Keys","cohort","MergeKeys",true);
ATT_cohort_obs = outerjoin(ATT_cohort_obs, att_n(:,["cohort","GroupCount"]), ...
    "Keys","cohort","MergeKeys",true);
ATT_cohort_obs.Properties.VariableNames(end) = "N_obs";
ATT_cohort_obs.SE_obs = ATT_cohort_obs.SD_obs ./ sqrt(max(ATT_cohort_obs.N_obs,1));
ATT_cohort_obs.CI_lo_obs = ATT_cohort_obs.ATT_obs + tinv(0.025, max(ATT_cohort_obs.N_obs-1,1)) .* ATT_cohort_obs.SE_obs;
ATT_cohort_obs.CI_hi_obs = ATT_cohort_obs.ATT_obs + tinv(0.975, max(ATT_cohort_obs.N_obs-1,1)) .* ATT_cohort_obs.SE_obs;
ATT_cohort_obs = sortrows(ATT_cohort_obs,"cohort");

% Output
out = struct();
out.Y0hat          = Y0hat;
out.tauHat         = tauHat;
out.adoptTime      = T2.adoptTime;
out.eventTime      = T2.k;
out.ATT_overall    = ATT_overall;
out.ATT_by_horizon = ATT_by_horizon;
out.ATT_cohort_obs = ATT_cohort_obs;
out.ModelInfo      = struct('NumUntreated', nnz(U0), ...
    'NumUnits', numel(categories(T2.(idVar))), ...
    'NumTimes', numel(categories(T2.(timeVar))), ...
    'Formula', char(form), 'Covariates', covars, ...
    'NumTreated', N_treated, ...
    'NumTreatedIdentified', N_treated_ident, ...
    'CoverageShare', N_treated_ident / max(1,N_treated));
out.Options        = opts;

% ----- Display
if opts.Display
    fprintf('\n');
    fprintf('[BJS] Untreated-only fit on ALL D==0: N_untreated = %d | Units = %d | Periods = %d\n', ...
        out.ModelInfo.NumUntreated, out.ModelInfo.NumUnits, out.ModelInfo.NumTimes);
    fprintf('[BJS] Treated identified: %d/%d (%.1f%%)\n', ...
        out.ModelInfo.NumTreatedIdentified, out.ModelInfo.NumTreated, 100*out.ModelInfo.CoverageShare);
    fprintf('[BJS] Overall ATT = %.4f', out.ATT_overall.ATT);
    if ~isnan(out.ATT_overall.SE)
        fprintf(' (SE=%.4f; 95%% CI [%.4f, %.4f])', out.ATT_overall.SE, out.ATT_overall.CI_lo, out.ATT_overall.CI_hi);
    end
    fprintf('\n');
    disp("ATT by Cohort");
    disp(ATT_cohort_obs);
end
end

% ======================
% Local helper: jackknife variance with NaN-safe handling
% ======================
function JV = jackknife_var_nan(v)
ok = ~isnan(v);
m = sum(ok);
if m < 2
    JV = NaN; return;
end
vm = mean(v(ok));
JV = (m-1)/m * sum( (v(ok) - vm).^2 );
end


function  [attk_lo,att_overall_lo] = calcATTk(T2,idVar,timeVar,yVar,D, isCatCov,form, catVars, Horizons, uid)
keep = T2.(idVar) ~= uid;
Tg = T2(keep,:);
Dg = D(keep);

if ~any(~Dg) || ~any(Dg), return; end

Tg_unt = Tg(~Dg,:);
Tg_unt.(idVar)   = removecats(Tg_unt.(idVar));
Tg_unt.(timeVar) = removecats(Tg_unt.(timeVar));
for c = find(isCatCov)
    Tg_unt.(covars(c)) = removecats(Tg_unt.(covars(c)));
end
mdl_g = fitlm(Tg_unt, char(form), 'CategoricalVars', cellstr(catVars));

idOK_g   = ismember(Tg.(idVar),   categories(Tg_unt.(idVar)));
timeOK_g = ismember(Tg.(timeVar), categories(Tg_unt.(timeVar)));
predOK_g = idOK_g & timeOK_g;
if any(isCatCov)
    for c = find(isCatCov)
        predOK_g = predOK_g & ismember(Tg.(covars(c)), categories(Tg_unt.(covars(c))));
    end
end

Y0g = NaN(height(Tg),1);
Y0g(predOK_g) = predict(mdl_g, Tg(predOK_g,:));
tau_g = NaN(height(Tg),1);
tau_g(Dg & predOK_g) = Tg.(yVar)(Dg & predOK_g) - Y0g(Dg & predOK_g);

att_overall_lo = mean(tau_g(Dg), 'omitnan');

kg = Tg.k;
for j = 1:numel(Horizons)
    Kj = Horizons(j);
    idx = Dg & (kg==Kj);
    if any(idx)
        attk_lo(1,j) = mean(tau_g(idx), 'omitnan');
    else
        attk_lo(1,j)=0;
    end
end

end
