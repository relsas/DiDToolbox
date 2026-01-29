function did_plot(res_in, what, varargin)
% DID_PLOT  Unified, light-themed plots for staggered DiD results.
%   did_plot(res, "event")
%   did_plot({res1,res2}, "cohort")
%   did_plot(resCH, "cohort")   % CH: plots DID_plus & DID_minus lines
%
% Supported 'what': "event" | "cohort" | "calendar" | "support"
% (Note: "overall" intentionally does nothing per spec.)
%
% ------------------------------------------------------------------------
% Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
% Last change: 10/03/2025
% ------------------------------------------------------------------------
ip = inputParser;
addRequired(ip,'res_in');
addRequired(ip,'what',@(s) any(strcmpi(string(s),["event","cohort","calendar","support","trends"])));
addParameter(ip,'Alpha',0.95,@(x) isnumeric(x) && x>0 && x<1);
addParameter(ip,'Legend',true,@islogical);
addParameter(ip,'ZeroLine',true,@islogical);
addParameter(ip,'Labels',{},@(x) iscellstr(x) || isstring(x));
addParameter(ip,'XAxis',"calendar",@(x) any(strcmpi(string(x),["calendar","event"])));
parse(ip,res_in,what,varargin{:});
P = ip.Results;

% Normalize input to a cell array
if ~iscell(res_in), res_in = {res_in}; end
K = numel(res_in);

% Standardize
S = cell(1,K);
L = strings(1,K);
userLabels = string(P.Labels);
for k = 1:K
    if isa(res_in{k}, 'did.Dataset')
        S{k} = res_in{k}; % Pass thorugh datasets for "trends" mode
    else
        S{k} = did.did_standardize(res_in{k}, 'Alpha', P.Alpha);
    end
    if ~isempty(userLabels) && numel(userLabels) >= k && strlength(userLabels(k))>0
        L(k) = userLabels(k);
    elseif isfield(res_in{k},'Method') && ~isempty(res_in{k}.Method)
        L(k) = string(res_in{k}.Method);
    elseif isfield(res_in{k},'Estimator') && ~isempty(res_in{k}.Estimator)
        L(k) = string(res_in{k}.Estimator);
    else
        L(k) = "Series " + k;
    end
end

% Figure title suffix with estimator tags (unique, comma-separated)
tags = unique(L);
titleSuffix = " — " + strjoin(tags, ", ");

switch lower(string(P.what))
    case "event"
        hasAny = any(cellfun(@(s) isfield(s,'es') && istable(s.es) && ~isempty(s.es), S));
        if ~hasAny, fprintf('No data for event plot \n'); return; end
        f = figure('Color','w','Name', "Event-time ATT" + titleSuffix); %#ok<NASGU>
        ax = gca; hold(ax,'on'); set_light_axes_(ax);

        for k=1:K
            if ~isfield(S{k},'es') || isempty(S{k}.es), continue; end
            T = S{k}.es;
            x = T.e; y = T.Estimate;

            if ismember('LB', T.Properties.VariableNames) && ismember('UB', T.Properties.VariableNames)
                lo = T.LB; hi = T.UB;
                errorbar(ax, x, y, y-lo, hi-y, 'o-','LineWidth',1.4,'CapSize',0,'DisplayName',L(k));
            elseif ismember('SE', T.Properties.VariableNames)
                crit = norminv(0.5+P.Alpha/2);
                se = T.SE; errorbar(ax, x, y, crit*se, crit*se, 'o-','LineWidth',1.4,'CapSize',0,'DisplayName',L(k));
            else
                plot(ax, x, y, 'o-','LineWidth',1.6,'DisplayName',L(k));
            end
        end
        grid(ax,'on'); xlabel(ax,'Event time e'); ylabel(ax,'ATT(e)');
        if P.ZeroLine, yline(ax,0,'k:'); xline(ax,0,'k:'); end
        if P.Legend, legend(ax,'Location','best'); end

    case "cohort"
        % Special CH handling: plot lines for DID_plus / DID_minus from ATT_by_Cohort
        hasCH    = any(cellfun(@(s) isfield(s,'chCohortLines') && istable(s.chCohortLines) && ~isempty(s.chCohortLines), S));
        hasOther = any(cellfun(@(s) isfield(s,'byCohort')      && istable(s.byCohort)      && ~isempty(s.byCohort), S));
        if ~(hasCH || hasOther), fprintf('No data for cohort plot \n'); return; end
        f = figure('Color','w','Name', "By-cohort / cohort-lines" + titleSuffix); %#ok<NASGU>
        ax = gca; hold(ax,'on'); set_light_axes_(ax);

        % 1) CH cohort-lines (DID_plus / DID_minus)
        if hasCH
            for k=1:K
                if ~isfield(S{k},'chCohortLines') || isempty(S{k}.chCohortLines), continue; end
                U = S{k}.chCohortLines;
                if ismember('DID_plus', U.Properties.VariableNames) && any(isfinite(U.DID_plus))
                    plot(ax, U.t, U.DID_plus, '-o','LineWidth',1.6,'DisplayName', L(k)+" (DID^+)");
                end
                if ismember('DID_minus', U.Properties.VariableNames) && any(isfinite(U.DID_minus))
                    plot(ax, U.t, U.DID_minus, '--s','LineWidth',1.6,'DisplayName', L(k)+" (DID^-)");
                end
            end
        end

        % 2) Non-CH by-cohort points with error bars
        if hasOther
            for k=1:K
                if ~isfield(S{k},'byCohort') || isempty(S{k}.byCohort), continue; end
                T = S{k}.byCohort;
                x = T.g; y = T.Estimate;
                if ismember('LB', T.Properties.VariableNames) && ismember('UB', T.Properties.VariableNames)
                    lo=T.LB; hi=T.UB; errorbar(ax,x,y,y-lo,hi-y,'o','LineWidth',1.4,'CapSize',0,'DisplayName',L(k));
                elseif ismember('SE', T.Properties.VariableNames)
                    crit = norminv(0.5+P.Alpha/2);
                    se=T.SE; errorbar(ax,x,y,crit*se,crit*se,'o','LineWidth',1.4,'CapSize',0,'DisplayName',L(k));
                else
                    plot(ax,x,y,'o','LineWidth',1.6,'DisplayName',L(k));
                end
            end
        end

        grid(ax,'on');
        if hasCH && ~hasOther
            xlabel(ax,'Time t'); ylabel(ax,'DID^+(t) / DID^-(t)');
        else
            xlabel(ax,'Cohort g (first treatment time)'); ylabel(ax,'ATT by cohort');
        end
        if P.ZeroLine, yline(ax,0,'k:'); end
        if P.Legend, legend(ax,'Location','best'); end

    case "calendar"
        hasAny = any(cellfun(@(s) isfield(s,'calendar') && istable(s.calendar) && ~isempty(s.calendar), S));
        if ~hasAny, fprintf('No data for calendar plot \n'); return; end
        f = figure('Color','w','Name', "Calendar-time ATT" + titleSuffix); %#ok<NASGU>
        ax = gca; hold(ax,'on'); set_light_axes_(ax);

        for k=1:K
            if ~isfield(S{k},'calendar') || isempty(S{k}.calendar), continue; end
            T = S{k}.calendar; x = T.t; y = T.Estimate;
            if ismember('LB', T.Properties.VariableNames) && ismember('UB', T.Properties.VariableNames)
                lo=T.LB; hi=T.UB; errorbar(ax,x,y,y-lo,hi-y,'-o','LineWidth',1.4,'CapSize',0,'DisplayName',L(k));
            elseif ismember('SE', T.Properties.VariableNames)
                crit = norminv(0.5+P.Alpha/2);
                se=T.SE; errorbar(ax,x,y,crit*se,crit*se,'-o','LineWidth',1.4,'CapSize',0,'DisplayName',L(k));
            else
                plot(ax,x,y,'-o','LineWidth',1.6,'DisplayName',L(k));
            end
        end
        grid(ax,'on'); xlabel(ax,'Calendar time t'); ylabel(ax,'ATT_c(t)');
        if P.ZeroLine, yline(ax,0,'k:'); end
        if P.Legend, legend(ax,'Location','best'); end

    case "support"
        hasAny = any(cellfun(@(s) isfield(s,'support') && istable(s.support) && ~isempty(s.support), S));
        if ~hasAny, fprintf('No data for support plot \n'); return; end
        f = figure('Color','w','Name', "Support by event time" + titleSuffix); %#ok<NASGU>
        ax = gca; hold(ax,'on'); set_light_axes_(ax);

        for k=1:K
            if ~isfield(S{k},'support') || isempty(S{k}.support), continue; end
            U = S{k}.support;
            if ismember('nObs', U.Properties.VariableNames)
                bar(ax, U.k, U.nObs, 'FaceAlpha',0.25, 'EdgeColor','none','DisplayName',L(k));
            elseif ismember('nCohorts', U.Properties.VariableNames)
                bar(ax, U.k, U.nCohorts, 'FaceAlpha',0.25, 'EdgeColor','none','DisplayName',L(k));
            else
                plot(ax, U.k, U{:,2},'-o','LineWidth',1.6,'DisplayName',L(k));
            end
        end
        grid(ax,'on'); xlabel(ax,'Event time e'); ylabel(ax,'Support');
        if P.ZeroLine, yline(ax,0,'k:'); xline(ax,0,'k:'); end
        if P.Legend, legend(ax,'Location','best'); end

    case "trends"
        % New case: Plot raw average trends by cohort (including never-treated)
        % Expects res_in to contain did.Dataset objects (not standardized structs)
        hasAny = any(cellfun(@(s) isa(s,'did.Dataset'), res_in));

        if ~hasAny
            % Check if it's the standardized struct but maybe we can't plot trends from summary alone?
            % Actually, trends require raw data.
            fprintf('Error: "trends" plot requires passing a did.Dataset object.\n');
            return;
        end

        f = figure('Color','w','Name', "Raw Trends by Cohort"); %#ok<NASGU>
        ax = gca; hold(ax,'on'); set_light_axes_(ax);

        % Color palette (MATLAB default approx)
        colors = lines(7);

        for k=1:K
            ds = res_in{k};
            if ~isa(ds,'did.Dataset'), continue; end

            % Materialize cohort info if not present
            g = ds.get("g");
            if strcmpi(P.XAxis, "event")
                t = ds.get("eventTime");
            else
                t = ds.get("t_int");
            end
            y = ds.T.(ds.yVar);

            % Aggregate mean(y) by (t, g)
            % Create temporary table for aggregation
            tmp = table(t, g, y, 'VariableNames', ["t", "g", "y"]);
            stats = groupsummary(tmp, ["t", "g"], "mean", "y");

            gVals = unique(stats.g);

            for i = 1:numel(gVals)
                g_val = gVals(i);
                sub = stats(stats.g == g_val, :);

                % Style
                if g_val == 0
                    % Never treated: Grey, dashed
                    lStyle = '--';
                    mStyle = 'o';
                    col    = [0.5 0.5 0.5];
                    name   = "Never Treated";
                    lw     = 2;
                else
                    % Treated: Solid circles
                    lStyle = '-';
                    mStyle = 'o';
                    col    = colors(mod(i-1,7)+1, :);
                    name   = sprintf("Cohort g=%g", g_val);
                    lw     = 1.5;
                end

                plot(ax, sub.t, sub.mean_y, ...
                    'LineStyle', lStyle, 'Marker', mStyle, ...
                    'Color', col, 'LineWidth', lw, ...
                    'DisplayName', name);

                % Add vertical line at treatment start (if treated)
                if strcmpi(P.XAxis, "event")
                    if i==1, xline(ax, -0.5, ':', 'Color','k', 'LineWidth', 1, 'HandleVisibility','off'); end
                else
                    if g_val > 0
                        xline(ax, g_val - 0.5, ':', 'Color', col, 'LineWidth', 1, 'HandleVisibility','off');
                    end
                end
            end
        end

        grid(ax,'on');
        if strcmpi(P.XAxis, "event")
            xlabel(ax,'Event Time e');
        else
            xlabel(ax,'Calendar Time t');
        end
        ylabel(ax,'Average Outcome Y');
        title(ax, ["Parallel Trends Inspection (Raw Means)", "X-Axis: " + P.XAxis]);
        if P.Legend, legend(ax,'Location','best'); end
end
end

%% ---- Light theme helper ----
function set_light_axes_(ax)
set(ax,'Color','w'); grid(ax,'on'); box(ax,'on');
ax.GridAlpha = 0.15; ax.LineWidth = 1.0;
end

