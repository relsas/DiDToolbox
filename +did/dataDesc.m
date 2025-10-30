function out = dataDesc(Data_or_ds, opts)
%DID.DATADESC  Quick data overview for DiD panels (periods, units, cohorts, leavers, Y descriptives).
%
% Accepts either:
%   out = did.dataDesc(ds, Display=true)              % ds is did.Dataset
%   out = did.dataDesc(T, idVar="id", timeVar="time", dVar="D", yVar="y")
%
% Outputs (struct):
%   .Overview, .Cohorts, .PeriodStats, .IdStats, .Y_overall, .Y_by_cohort, .Options
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/09/2025
% ------------------------------------------------------------------------

arguments
    Data_or_ds
    % When a TABLE is passed, these are required (except yVar)
    opts.idVar   (1,1) string = ""
    opts.timeVar (1,1) string = ""
    opts.dVar    (1,1) string = ""
    opts.yVar    string = ""
    opts.Display (1,1) logical = true
end

% -------- Resolve inputs from Dataset or Table --------
if isa(Data_or_ds, 'did.Dataset')
    ds      = Data_or_ds;
    T       = ds.T;
    idVar   = ds.idVar;
    timeVar = ds.timeVar;
    dVar    = ds.dVar;
    % yVar: prefer explicit opts.yVar if provided; else pick ds.yVar if present
    if strlength(opts.yVar) > 0
        yVar = string(opts.yVar);
    elseif isprop(ds,'yVar') && strlength(ds.yVar)>0
        yVar = ds.yVar;
    else
        yVar = "";
    end
    Display = opts.Display;
else
    % Legacy table path
    T       = Data_or_ds;
    idVar   = opts.idVar;
    timeVar = opts.timeVar;
    dVar    = opts.dVar;
    yVar    = string(opts.yVar);
    Display = opts.Display;

    if strlength(idVar)==0 || strlength(timeVar)==0 || strlength(dVar)==0
        error('did:dataDesc:MissingVarNames', ...
            'When passing a table, idVar/timeVar/dVar must be provided.');
    end
end

% Ensure sorted (no mutation of external object; we work on a copy)
T = sortrows(T,[idVar, timeVar]);

% ----- Basic checks (lightweight — Dataset.fromTable already validates) -----
varsNeeded = [idVar, timeVar, dVar];
missingVars = varsNeeded(~ismember(varsNeeded, string(T.Properties.VariableNames)));
if ~isempty(missingVars)
    error('did:dataDesc:MissingVars','Missing variables in T: %s', strjoin(missingVars, ', '));
end

% Coerce D to numeric 0/1 if logical; otherwise accept numeric 0/1
D = T.(dVar);
if islogical(D), D = double(D); end
if ~isnumeric(D)
    error('did:dataDesc:Dtype','%s must be numeric/logical 0/1.', dVar);
end

% ----- Working copy & stable integer time index -----
T2 = T;
[~,~,t_idx] = unique(T2.(timeVar), 'stable');
T2.t_int = double(t_idx);

% Canonical time values in order
t_values = unique(T2.(timeVar), 'stable');
nT = numel(t_values);

% Units
uid = unique(T2.(idVar));
nU  = numel(uid);

% Panel balance by unit
[gId, ~] = findgroups(T2.(idVar));
nObsByUnit  = splitapply(@(x) numel(x),          T2.t_int, gId);
nTimeByUnit = splitapply(@(x) numel(unique(x)),  T2.t_int, gId);
isBalancedUnit = (nTimeByUnit == nT);
nBalanced = sum(isBalancedUnit);

% Ever treated, adoption (first-on), first-off (on→off)
adopt_int    = splitapply(@localFirstOn,  T2.t_int, D, gId);
firstOff_int = splitapply(@localFirstOff, T2.t_int, D, gId);
everTreated  = ~isnan(adopt_int);
isLeaver     = ~isnan(firstOff_int);

% Map ints → original time values (preserve type)
adopt_time    = localMapTime(adopt_int,    t_values);
firstOff_time = localMapTime(firstOff_int, t_values);

% ---- Id-level table
IdStats = table;
IdStats.(idVar)       = splitapply(@(x) x(1), T2.(idVar), gId);
IdStats.nObs          = nObsByUnit;
IdStats.adopt_int     = adopt_int;
IdStats.adopt_time    = adopt_time;
IdStats.firstOff_int  = firstOff_int;
IdStats.firstOff_time = firstOff_time;
IdStats.everTreated   = everTreated;
IdStats.isLeaver      = isLeaver;

% ---- Cohorts (first treatment) + "never"
cohort_key = adopt_int;   % NaN = never
[Gcoh, Cvals] = findgroups(cohort_key);
N_by_cohort  = splitapply(@numel, cohort_key, Gcoh);

Cohorts = table;
Cohorts.cohort_int  = Cvals;
Cohorts.cohort_time = localMapTime(Cvals, t_values);

% Labels for display
lab = strings(height(Cohorts),1);
for i=1:height(Cohorts)
    if isnan(Cohorts.cohort_int(i))
        lab(i) = "never";
    else
        lab(i) = string(Cohorts.cohort_time(i));
    end
end
Cohorts.cohort_label = lab;

Cohorts.N_units = N_by_cohort;
Cohorts.Share   = Cohorts.N_units / nU;

% Leavers by cohort (ignore "never" for share calc)
coh_leavers = splitapply(@(a,l) sum(l(~isnan(a))), adopt_int, isLeaver, Gcoh);
Cohorts.N_leavers = coh_leavers;
Cohorts.ShareLeaversAmongCohort = NaN(height(Cohorts),1);
ixTreatedCoh = ~isnan(Cohorts.cohort_int);
Cohorts.ShareLeaversAmongCohort(ixTreatedCoh) = Cohorts.N_leavers(ixTreatedCoh) ./ max(Cohorts.N_units(ixTreatedCoh),1);

% Convenience aliases (explicit first treat time columns)
Cohorts.FirstTreat_t_int = Cohorts.cohort_int;
Cohorts.FirstTreat_time  = Cohorts.cohort_time;

% ---- Period-level stats
[gT, tIntKeys] = findgroups(T2.t_int);
UnitsObserved  = splitapply(@numel, T2.(idVar), gT);
TreatedUnits   = splitapply(@(x) sum(x==1), D, gT);

% Bring adoption per id onto rows for period computations
AdoptById = table;
AdoptById.(idVar) = IdStats.(idVar);
AdoptById.adopt_int = IdStats.adopt_int;

T2 = outerjoin(T2, AdoptById, "Keys", idVar, "MergeKeys", true, "Type","left");
NotYetUnits   = splitapply(@(a,t) sum(~isnan(a) & a>t),           T2.adopt_int, T2.t_int, gT);
NeverUnits    = splitapply(@(a)   sum(isnan(a)),                  T2.adopt_int, gT);
OffAfterOn    = splitapply(@(a,d,t) sum(~isnan(a) & a<=t & d==0), T2.adopt_int, D, T2.t_int, gT);

PeriodStats = table;
PeriodStats.time         = localMapTime(tIntKeys, t_values);
PeriodStats.t_int        = tIntKeys;
PeriodStats.UnitsObs     = UnitsObserved;
PeriodStats.Treated      = TreatedUnits;
PeriodStats.NotYet       = NotYetUnits;
PeriodStats.Never        = NeverUnits;
PeriodStats.OffAfterOn   = OffAfterOn;
PeriodStats.TreatedShare = PeriodStats.Treated ./ max(PeriodStats.UnitsObs,1);

% ---- Overview
Overview = table;
Overview.TimeStart     = t_values(1);
Overview.TimeEnd       = t_values(end);
Overview.NumPeriods    = nT;
Overview.NumUnits      = nU;
Overview.Nobs          = height(T2);
Overview.BalancedUnits = nBalanced;
Overview.BalancedShare = nBalanced / max(nU,1);
Overview.EverTreated   = sum(everTreated);
Overview.NeverTreated  = sum(~everTreated);
Overview.Leavers       = sum(isLeaver);
Overview.LeaversShareAmongEver = Overview.Leavers / max(Overview.EverTreated,1);

% ---- Y descriptives (optional)
Y_overall   = table();
Y_by_cohort = table();
if strlength(yVar)>0 && any(strcmp(string(T.Properties.VariableNames), yVar))
    y = T.(yVar);

    % Overall -> struct -> table
    s_overall = localSummStats(y);
    Y_overall = struct2table(s_overall);

    % By cohort (including "never" as 0 for grouping)
    T3 = outerjoin(AdoptById(:, [idVar "adopt_int"]), T(:, [idVar yVar]), ...
        "Keys", idVar, "MergeKeys", true, "Type","left");

    adopt_key = T3.adopt_int;
    adopt_key(isnan(adopt_key)) = 0;          % 0 == never-treated (groupable)

    [G, K] = findgroups(adopt_key);
    s = splitapply(@localSummStats, T3.(yVar), G);   % struct array

    % Map cohort key -> time (0 => typed missing)
    cohort_time = localMapTimeZero(K, t_values);

    Y_by_cohort = table( ...
        K, cohort_time, ...
        [s.N]', [s.Mean]', [s.SD]', [s.Min]', [s.P25]', [s.Median]', [s.P75]', [s.Max]', ...
        'VariableNames', ["cohort_int","cohort_time","N","Mean","SD","Min","P25","Median","P75","Max"]);

    % Human-friendly label
    lab2 = strings(height(Y_by_cohort),1);
    for i=1:height(Y_by_cohort)
        if Y_by_cohort.cohort_int(i)==0
            lab2(i) = "never";
        else
            lab2(i) = string(Y_by_cohort.cohort_time(i));
        end
    end
    Y_by_cohort.cohort_label = lab2;
    Y_by_cohort = movevars(Y_by_cohort, 'cohort_label', 'Before', 'N');
end

% ---- Pack outputs
out = struct();
out.Overview    = Overview;
out.Cohorts     = sortrows(Cohorts, {'FirstTreat_t_int'}, {'ascend'});
out.PeriodStats = sortrows(PeriodStats, 't_int');
out.IdStats     = IdStats;
out.Y_overall   = Y_overall;

% Guarded sort for optional table
if ~isempty(Y_by_cohort) && any(strcmp('cohort_int', Y_by_cohort.Properties.VariableNames))
    out.Y_by_cohort = sortrows(Y_by_cohort, {'cohort_int'});
else
    out.Y_by_cohort = Y_by_cohort;  % empty or already fine
end

out.Options = struct('idVar',idVar, 'timeVar',timeVar, 'dVar',dVar, 'yVar',yVar, 'Display',Display);

% ---- Display
if Display
    fprintf('\n[DataDesc] Units = %d | Periods = %d | Nobs = %d | Balanced Units = %d (%.1f%%)\n', ...
        Overview.NumUnits, Overview.NumPeriods, Overview.Nobs, ...
        Overview.BalancedUnits, 100*Overview.BalancedShare);
    fprintf('[DataDesc] Ever-treated = %d (%.1f%%) | Never-treated = %d (%.1f%%) | Leavers = %d (%.1f%% of ever)\n', ...
        Overview.EverTreated, 100*Overview.EverTreated/Overview.NumUnits, ...
        Overview.NeverTreated, 100*Overview.NeverTreated/Overview.NumUnits, ...
        Overview.Leavers, 100*Overview.LeaversShareAmongEver);

    % Top cohorts (first few)
    tmp = out.Cohorts;
    if height(tmp)>0
        k = min(6,height(tmp));
        disp("Cohorts (first 6 rows):");
        disp(tmp(1:k, ["cohort_label","N_units","Share","N_leavers","ShareLeaversAmongCohort","FirstTreat_time"]));
    end

    % Y summary
    if ~isempty(Y_overall)
        disp("Outcome Y (overall) – N/Mean/SD/Min/P25/Median/P75/Max:");
        disp(Y_overall);
        if ~isempty(out.Y_by_cohort)
            disp("Outcome Y by cohort (including 'never'):");
            disp(out.Y_by_cohort);
        end
    end
end

% ======== Local helpers ========
    function a = localFirstOn(t, d)
        % first t where d==1
        [ts, ix] = sort(t);
        ds = d(ix);
        j = find(ds==1, 1, 'first');
        if isempty(j), a = NaN; else, a = ts(j); end
    end

    function f = localFirstOff(t, d)
        % first t where unit turns off after being on at least once
        [ts, ix] = sort(t);
        ds = d(ix);
        j = find(ds==1, 1, 'first');
        if isempty(j)
            f = NaN;
        else
            k = find(ds((j+1):end)==0, 1, 'first');
            if isempty(k), f = NaN; else, f = ts(j+k); end
        end
    end

    function tv = localMapTime(ti, tvals)
        % Map integer times (ti) back to original type in tvals
        tv = repmat(tvals(1), size(ti));  % preallocate with correct type
        for ii = 1:numel(ti)
            if isnan(ti(ii))
                % typed missing
                if isdatetime(tvals)
                    tv(ii) = NaT;
                elseif isnumeric(tvals)
                    tv(ii) = NaN;
                elseif isstring(tvals) || ischar(tvals) || iscellstr(tvals)
                    tv(ii) = string(missing);
                else
                    tv(ii) = tvals(1); % fallback
                end
            else
                tv(ii) = tvals(ti(ii));
            end
        end
    end

    function S = localSummStats(v)
        % Return a struct so splitapply can build a struct array cleanly
        v = v(:);
        v = v(~isnan(v));
        if isempty(v)
            S = struct('N',0,'Mean',NaN,'SD',NaN,'Min',NaN,'P25',NaN,'Median',NaN,'P75',NaN,'Max',NaN);
        else
            q = prctile(v, [25 50 75]);
            S = struct( ...
                'N',      numel(v), ...
                'Mean',   mean(v), ...
                'SD',     std(v), ...
                'Min',    min(v), ...
                'P25',    q(1), ...
                'Median', q(2), ...
                'P75',    q(3), ...
                'Max',    max(v));
        end
    end

    function tv = localMapTimeZero(k, tvals)
        % Map cohort_int (0 = never) to typed time (missing for 0)
        tv = repmat(tvals(1), size(k));  % preallocate with correct type
        for ii = 1:numel(k)
            if k(ii) <= 0 || isnan(k(ii))
                if isdatetime(tvals)
                    tv(ii) = NaT;
                elseif isnumeric(tvals)
                    tv(ii) = NaN;
                elseif isstring(tvals) || ischar(tvals) || iscellstr(tvals)
                    tv(ii) = string(missing);
                else
                    tv(ii) = tvals(1); % fallback
                end
            else
                tv(ii) = tvals(k(ii));
            end
        end
    end
end
