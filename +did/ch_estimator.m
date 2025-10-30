function out = ch_estimator(Data, opts)
%DID.CH_ESTIMATOR  CH.m is DID_M wrapper (joiners/leavers, placebo, cluster bootstrap).
%  
% The DID_M estimator targets , the Average Treatment Effect (ATE) across all "switching cells"—groups whose treatment status changes between $t-1$ and $t$. This is calculated as a weighted average of two specific DiD terms:
%     ◦ $DID_{+,t}$: Compares joiners (treatment $0 \to 1$) to groups remaining untreated.
%     ◦ $DID_{-,t}$: Compares leavers (treatment $1 \to 0$) to groups remaining treated.
% $DID_M$ is valid under heterogeneous treatment effects and inherently avoids the source of negative weighting bias by strictly 
% contrasting switching groups against stable groups. In staggered adoption designs, where treatment is typically irreversible, 
% $DID_M$ is a weighted average only of the $DID_{+,t}$ estimators.
% 
% out = did.ch_estimator(Data=tbl, GroupVar="id", TimeVar="time", YVar="y", DVar="D", ...
%                         WeightVar="", B=299, ComputePlacebo=true, Seed=42, ...
%                         Print=true, Details=false, Covariates=[], CovarSample="D0")
%
%  Returns:
%    out.Summary  table with rows: DID_M, Joiners-only, Leavers-only, (Placebo)
%    plus: DIDM, SE, DIDM_joiners/leavers (+ SE), ByT, Placebo, Diagnostics.
% 
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/03/2025  
% ------------------------------------------------------------------------
arguments
    Data {mustBeA(Data,"table")}
    opts.idVar (1,1) string = "id"
    opts.timeVar  (1,1) string = "time"
    opts.yVar     (1,1) string = "y"
    opts.dVar     (1,1) string = "D"
    opts.WeightVar (1,1) string = ""
    opts.B (1,1) double {mustBeInteger, mustBeNonnegative} = 100
    opts.ComputePlacebo (1,1) logical = true
    opts.Seed (1,1) double = NaN
    opts.Print (1,1) logical = true
    opts.Details (1,1) logical = false
    % NEW: FWL covariates
    opts.Covariates string = string.empty(1,0)
    opts.CovarSample (1,1) string {mustBeMember(opts.CovarSample,["D0","never","all"])} = "D0"
end

T = Data;

% --- weights
hasW = ~isempty(opts.WeightVar) && ismember(opts.WeightVar, string(T.Properties.VariableNames));
if hasW
    w = double(T.(opts.WeightVar));
    w(~isfinite(w)) = NaN;
else
    w = ones(height(T),1);
end

% --- covariates
Xnames = string.empty(1,0);
if ~isempty(opts.Covariates)
    Xnames = string(opts.Covariates(:)');
    missX = Xnames(~ismember(Xnames, string(T.Properties.VariableNames)));
    assert(isempty(missX), 'did.ch_estimator:MissingCovariates', ...
        'Covariate(s) not found: %s', strjoin(missX, ', '));
end

% --- micro table for the core (include X as X__name columns)
GT = table(T.(opts.idVar), T.(opts.timeVar), double(T.(opts.yVar)), ...
           double(T.(opts.dVar)), double(w), ...
    'VariableNames', {'G','T','Y','D','N'});
for k = 1:numel(Xnames)
    GT.('X__' + Xnames(k)) = double(T.(Xnames(k)));
end
GT = sortrows(GT, {'G','T'});


% crude df for p-values (t approx)
nClusters = numel(unique(GT.G));
df = max(10, nClusters-1);

% --- run core
ob = ch_core_(GT=GT, B=opts.B, ComputePlacebo=opts.ComputePlacebo, Seed=opts.Seed, ...
              Xnames=Xnames, CovarSample=opts.CovarSample);

% --- assemble
out = ob;
out.Estimator = "CH2020";
out.Call = struct('idVar',opts.idVar, 'timeVar',opts.timeVar, 'yVar',opts.yVar, ...
                  'dVar',opts.dVar, 'WeightVar',opts.WeightVar, 'B',opts.B, ...
                  'ComputePlacebo',opts.ComputePlacebo, 'Seed',opts.Seed, ...
                  'Covariates',{Xnames}, 'CovarSample',opts.CovarSample);

% summary table
Names = ["DID_M"; "Joiners-only"; "Leavers-only"];
Est   = [ob.DIDM; ob.DIDM_joiners; ob.DIDM_leavers];
SEs   = [ob.SE;   ob.SE_joiners;  ob.SE_leavers];
tval  = Est ./ SEs;
pval  = 2*tcdf(-abs(tval), df);

if isfield(ob,'OverallCW')
    Names = [Names; "Overall (cohort-weighted)"];
    Est   = [Est;   ob.OverallCW];
    SEs   = [SEs;   ob.SECW];
    tval  = [tval;  Est(end)/SEs(end)];
    pval  = [pval;  2*tcdf(-abs(tval(end)), df)];
end

if isfield(ob,'Placebo') && ~isempty(ob.Placebo) && isfield(ob.Placebo,'DIDM_pl')
    Names = [Names; "Placebo DID_M"];
    Est   = [Est;   ob.Placebo.DIDM_pl];
    SEs   = [SEs;   ob.Placebo.SE_pl];
    tval  = [tval;  ob.Placebo.DIDM_pl / ob.Placebo.SE_pl];
    pval  = [pval;  2*tcdf(-abs(tval(end)), df)];
end

out.summaryTable = table(Names, Est, SEs, tval, pval, ...
    'VariableNames', {'Effect','Estimate','SE','t','p'});

% pretty print
if opts.Print
    fprintf('\n=== DID_M summary (table) ===\n');
    disp(out.summaryTable);

    fprintf('\n=== DID_M by Cohort (table) ===\n');
    disp(out.ATT_by_Cohort);

end
end

% ---------- core (now FWL-aware + bootstrap re-estimates γ) ----------
function out = ch_core_(args)
arguments
    args.GT table
    args.B (1,1) double {mustBeInteger, mustBeNonnegative} = 399
    args.ComputePlacebo (1,1) logical = true
    args.Seed (1,1) double = NaN
    args.Xnames string = string.empty(1,0)
    args.CovarSample (1,1) string {mustBeMember(args.CovarSample,["D0","never","all"])} = "D0"
end

GT0 = args.GT;

% --- FWL residualization if covariates present ---
if ~isempty(args.Xnames)
    % Build X from X__* columns
    X = zeros(height(GT0), numel(args.Xnames));
    for k = 1:numel(args.Xnames)
        X(:,k) = double(GT0.("X__" + args.Xnames(k)));
    end
    y = double(GT0.Y);
    w = double(GT0.N);

    % Choose sample to estimate gamma
    switch args.CovarSample
        case "D0"
            mask = (GT0.D==0);
        case "never"
            [ggrp,~] = findgroups(GT0.G);
            ever1 = splitapply(@(d) any(d==1), GT0.D, ggrp);
            ever1 = ever1(ggrp);
            mask = ~ever1;
        case "all"
            mask = true(height(GT0),1);
    end
    mask = mask & isfinite(y) & all(isfinite(X),2) & isfinite(w) & (w>0);

    if nnz(mask) >= size(X,2)
        sw = sqrt(w(mask));
        Xw = X(mask,:).*sw;  yw = y(mask).*sw;
        beta = Xw\yw;  if any(~isfinite(beta)), beta = pinv(Xw'*Xw)*(Xw'*yw); end
        GT0.Y = y - X*beta;  % residualized outcome
    end
end

% --- collapse to (G,T)
[Ggt, gk, tk] = findgroups(GT0.G, GT0.T);
wmean2 = @(x,w) sum(x.*w,'omitnan')./sum(w,'omitnan');
Ygt  = splitapply(wmean2, GT0.Y, GT0.N, Ggt);
Dbar = splitapply(wmean2, GT0.D, GT0.N, Ggt);
Ngt  = splitapply(@(w) sum(w,'omitnan'), GT0.N, Ggt);

GT = table(gk, tk, Ygt, Dbar, Ngt, 'VariableNames',{'G','T','Y','Dbar','N'});

% (new) enforce order by G then T before computing lags
GT = sortrows(GT, {'G','T'});

if any(GT.Dbar~=0 & GT.Dbar~=1)
    warning('did.ch_estimator:NotSharp','Non-sharp cells detected; rounding D.');
end
GT.D = round(GT.Dbar);

% --- lags within G
% Preallocate
GT.Ylag1 = NaN(height(GT),1);
GT.Ylag2 = NaN(height(GT),1);
GT.Dlag1 = NaN(height(GT),1);
GT.Dlag2 = NaN(height(GT),1);

uG = unique(GT.G);
for gi = 1:numel(uG)
    idx = find(GT.G==uG(gi));              % rows for this group, already sorted by T
    y = GT.Y(idx); d = GT.D(idx);
    if numel(idx) >= 1
        GT.Ylag1(idx) = [NaN; y(1:end-1)];
        GT.Dlag1(idx) = [NaN; d(1:end-1)];
    end
    if numel(idx) >= 2
        GT.Ylag2(idx) = [NaN; NaN; y(1:end-2)];
        GT.Dlag2(idx) = [NaN; NaN; d(1:end-2)];
    end
end


GT.dY     = GT.Y - GT.Ylag1;
GT.dYprev = GT.Ylag1 - GT.Ylag2;

isJoiner_t = (GT.D==1 & GT.Dlag1==0);
isLeaver_t = (GT.D==0 & GT.Dlag1==1);
isStab0_t  = (GT.D==0 & GT.Dlag1==0);
isStab1_t  = (GT.D==1 & GT.Dlag1==1);

isJoiner_pl = (GT.D==1 & GT.Dlag1==0 & GT.Dlag2==0);
isLeaver_pl = (GT.D==0 & GT.Dlag1==1 & GT.Dlag2==1);
isStab0_pl  = (GT.D==0 & GT.Dlag1==0 & GT.Dlag2==0);
isStab1_pl  = (GT.D==1 & GT.Dlag1==1 & GT.Dlag2==1);

[Gt, ut] = findgroups(GT.T);
sumW  = @(x) sum(x,'omitnan');
wmean = @(x,w) sum(x.*w,'omitnan')./sum(w,'omitnan');

N10t = splitapply(@(w) sumW(w), GT.N.*isJoiner_t, Gt);
N01t = splitapply(@(w) sumW(w), GT.N.*isLeaver_t, Gt);
N00t = splitapply(@(w) sumW(w), GT.N.*isStab0_t,  Gt);
N11t = splitapply(@(w) sumW(w), GT.N.*isStab1_t,  Gt);

dY_join   = splitapply(@(x,w) wmean(x,w), GT.dY, GT.N.*isJoiner_t, Gt);
dY_stab0  = splitapply(@(x,w) wmean(x,w), GT.dY, GT.N.*isStab0_t,  Gt);
dY_stab1  = splitapply(@(x,w) wmean(x,w), GT.dY, GT.N.*isStab1_t,  Gt);
dY_leave  = splitapply(@(x,w) wmean(x,w), GT.dY, GT.N.*isLeaver_t, Gt);

DID_plus_t  = dY_join  - dY_stab0;
DID_minus_t = dY_stab1 - dY_leave;

ByT = table(ut, N10t, N01t, N00t, N11t, DID_plus_t, DID_minus_t, ...
            'VariableNames', {'t','N10','N01','N00','N11','DID_plus','DID_minus'});

% clean missing
ByT.N10(isnan(ByT.DID_plus))  = 0;   ByT.DID_plus(isnan(ByT.DID_plus))   = 0;
ByT.N01(isnan(ByT.DID_minus)) = 0;   ByT.DID_minus(isnan(ByT.DID_minus)) = 0;

% === Helper: cohort-weighted overall from micro GT0 and ByT ===
cohortWeighted = compute_cohort_weighted_overall_(GT0, ByT);

% (store now; SE comes from bootstrap below)
OverallCW = cohortWeighted;
SECW = NaN;


NS = sum(ByT.N10 + ByT.N01);
wP = ByT.N10 / max(NS, eps);
wM = ByT.N01 / max(NS, eps);

DIDM         = sum(wP .* ByT.DID_plus + wM .* ByT.DID_minus);
DIDM_joiners = sum((ByT.N10/max(sum(ByT.N10),eps)) .* ByT.DID_plus);
DIDM_leavers = sum((ByT.N01/max(sum(ByT.N01),eps)) .* ByT.DID_minus);

% placebo
DIDM_pl = NaN; SEpl = NaN; Placebo = struct([]);
if args.ComputePlacebo
    N100t = splitapply(@(w) sumW(w), GT.N.*isJoiner_pl, Gt);
    N011t = splitapply(@(w) sumW(w), GT.N.*isLeaver_pl, Gt);
    dYprev_join  = splitapply(@(x,w) wmean(x,w), GT.dYprev, GT.N.*isJoiner_pl, Gt);
    dYprev_stab0 = splitapply(@(x,w) wmean(x,w), GT.dYprev, GT.N.*isStab0_pl,  Gt);
    dYprev_leave = splitapply(@(x,w) wmean(x,w), GT.dYprev, GT.N.*isLeaver_pl, Gt);
    dYprev_stab1 = splitapply(@(x,w) wmean(x,w), GT.dYprev, GT.N.*isStab1_pl,  Gt);

    DID_plus_pl_t  = dYprev_join - dYprev_stab0;
    DID_minus_pl_t = dYprev_stab1 - dYprev_leave;

    mP = isfinite(DID_plus_pl_t)  & (N100t>0);
    mM = isfinite(DID_minus_pl_t) & (N011t>0);
    N100t(~mP) = 0; DID_plus_pl_t(~mP) = 0;
    N011t(~mM) = 0; DID_minus_pl_t(~mM) = 0;

    NplS = sum(N100t + N011t);
    if NplS>0
        wPpl = N100t / NplS;
        wMpl = N011t / NplS;
        DIDM_pl = sum(wPpl .* DID_plus_pl_t + wMpl .* DID_minus_pl_t);
    end

    Placebo = struct('ByT',table(ut, N100t, N011t, DID_plus_pl_t, DID_minus_pl_t, ...
                        'VariableNames',{'t','N100','N011','DID_plus_pl','DID_minus_pl'}), ...
                     'DIDM_pl',DIDM_pl);
end

Diagnostics = struct('Times',ByT.t, 'HasStable0',ByT.N00>0, 'HasStable1',ByT.N11>0, ...
                     'NSwitchers',NS, 'SharpCellFraction',mean(GT.Dbar==0 | GT.Dbar==1));

% --- cluster bootstrap (re-estimate FWL each draw)
SE=NaN; SEJ=NaN; SEL=NaN; SEpl=NaN;
B = args.B;
if B>0
    rng(args.Seed,'twister');
    [Gid,uG] = findgroups(GT0.G);
    groupRows = accumarray(Gid, (1:height(GT0))', [numel(uG) 1], @(ix){ix});
    nG = numel(uG);

   b_est  = nan(B,1); bJ = nan(B,1); bL = nan(B,1); bPl = nan(B,1); bCW = nan(B,1);
for b = 1:B
    draw = randi(nG, nG, 1);
    rows = vertcat(groupRows{draw});
    Tb   = GT0(rows,:);
    ob   = ch_core_(GT=Tb, B=0, ComputePlacebo=args.ComputePlacebo, Seed=args.Seed, ...
                    Xnames=args.Xnames, CovarSample=args.CovarSample);
    b_est(b) = ob.DIDM; 
    bJ(b)    = ob.DIDM_joiners; 
    bL(b)    = ob.DIDM_leavers;
    if args.ComputePlacebo && isfield(ob,'Placebo'), bPl(b) = ob.Placebo.DIDM_pl; end
    if isfield(ob,'OverallCW'), bCW(b) = ob.OverallCW; end
end
SE   = std(b_est,0,1,'omitnan');
SEJ  = std(bJ,0,1,'omitnan');
SEL  = std(bL,0,1,'omitnan');
SEpl = std(bPl,0,1,'omitnan');
SECW = std(bCW,0,1,'omitnan');

end


%% ATT by Cohort
% cohort efffect entry / exit
wer =find(ByT.DID_plus~=0);
ATT_Plus_Cohort = ByT(wer,["t","DID_plus"]);
wer =find(ByT.DID_minus~=0);
ATT_Minus_Cohort = ByT(wer,["t","DID_minus"]);

ATT_by_Cohort = outerjoin(ATT_Plus_Cohort,ATT_Minus_Cohort,'MergeKeys',true);
%% Pack outputs
out = struct();
out.DIDM = DIDM; out.SE = SE;
out.DIDM_joiners = DIDM_joiners; out.SE_joiners = SEJ;
out.DIDM_leavers = DIDM_leavers; out.SE_leavers = SEL;
out.ByT = ByT; out.Diagnostics = Diagnostics;
out.Placebo = struct('Enabled',args.ComputePlacebo,'DIDM_pl',DIDM_pl,'SE_pl',SEpl,'Details',Placebo);
out.ATT_by_Cohort = ATT_by_Cohort;
out.OverallCW = OverallCW;      % cohort-weighted overall
out.SECW      = SECW;           % bootstrap SE for OverallCW

    function val = compute_cohort_weighted_overall_(GTmicro, ByTtab)
    % Build maps t -> DID_plus / DID_minus
    tvals = ByTtab.t;
    mapP = containers.Map(num2cell(tvals), num2cell(ByTtab.DID_plus));
    mapM = containers.Map(num2cell(tvals), num2cell(ByTtab.DID_minus));

    % For each unit, find first 0->1 (join) time and first 1->0 (leave) time
    [gU,~,gi] = unique(GTmicro.G,'stable');
    Tu = GTmicro.T; Du = GTmicro.D;
    nU = numel(gU);

    joinTimes  = NaN(nU,1);
    leaveTimes = NaN(nU,1);

    % indices per unit for fast loop
    idxPer = accumarray(gi, (1:height(GTmicro))', [nU 1], @(ix){ix});
    for u = 1:nU
        ix = idxPer{u};
        [~,ord] = sort(Tu(ix),'ascend');
        ix = ix(ord);
        Dseq = Du(ix);
        Tseq = Tu(ix);

        % first 0->1
        ch = diff(Dseq);
        jpos = find(ch==1,1,'first');
        if ~isempty(jpos), joinTimes(u) = Tseq(jpos+1); end

        % first 1->0
        lpos = find(ch==-1,1,'first');
        if ~isempty(lpos), leaveTimes(u) = Tseq(lpos+1); end
    end

    % Counts per cohort (ignore NaNs)
    if all(isnan(joinTimes)) && all(isnan(leaveTimes))
        val = NaN; return;
    end
    % counts per time
    [uj,~,cj] = unique(joinTimes(~isnan(joinTimes)),'stable');
    cntJ = accumarray(cj,1);
    [ul,~,cl] = unique(leaveTimes(~isnan(leaveTimes)),'stable');
    cntL = accumarray(cl,1);

    SJ = sum(cntJ); SL = sum(cntL); S = SJ + SL;
    if S==0, val = NaN; return; end

    % cohort-weighted average over joiner and leaver cohorts
    agg = 0.0;
    if ~isempty(cntJ)
        for k=1:numel(uj)
            tk = uj(k);
            if isKey(mapP, tk)
                agg = agg + (cntJ(k)/S) * mapP(tk);
            end
        end
    end
    if ~isempty(cntL)
        for k=1:numel(ul)
            tk = ul(k);
            if isKey(mapM, tk)
                agg = agg + (cntL(k)/S) * mapM(tk);
            end
        end
    end
    val = agg;
end

end


