classdef SDID < did.estimators.Estimator
    % did.estimators.SDID  Synthetic Difference-in-Differences (and Synthetic Control)
    %
    % Implements the SDID estimator from Arkhangelsky et al. (2021).
    % Also supports pure Synthetic Control (SC) by setting Weights='unit'.
    %
    % Properties:
    %   Weights        : "SDID" (default), "SC" (Unit only), "Time" (Time only)
    %   Regularization : Ridge penalty parameter (zeta). Default 1e-6.
    %   Solver         : "auto", "quadprog", "frank-wolfe"
    %   SEMethod       : "Placebo" (default) or "Jackknife"
    %   B              : Number of bootstrap/placebo replications (default 50)
    %
    % ------------------------------------------------------------------------

    properties
        Weights    string = "SDID"
        Regularization double = NaN % NaN = Auto (Data-driven)
        Solver     string = "auto"
        SEMethod   string = "Placebo"
        B          double = 50
        Display    logical = true
        UseParallel logical = true
    end

    methods
        function obj = SDID(varargin)
            % Parser
            if ~isempty(varargin)
                for k = 1:2:numel(varargin)
                    name  = lower(string(varargin{k}));
                    value = varargin{k+1};
                    switch name
                        case 'weights',    obj.Weights = string(value);
                        case 'regularization', obj.Regularization = double(value);
                        case 'solver',     obj.Solver = string(value);
                        case 'semethod',   obj.SEMethod = string(value);
                        case 'b',          obj.B = double(value);
                        case 'display',    obj.Display = logical(value);
                        case 'useparallel', obj.UseParallel = logical(value);
                    end
                end
            end
        end

        function res = fit(obj, ds)
            % 1. Data Prep (Pivot to Wide)
            [Y, N_units, T_periods, id_map, time_map, D_mat] = obj.pivotData(ds);

            % 2. Main Estimation (Aggregated over Cohorts)
            [tau, se, aggreg_details] = obj.estimate_aggregate(Y, D_mat, id_map, time_map);

            % 3. Result Construction
            res.Method = "SDID";
            res.tau = tau;
            res.se  = se;
            res.WeightsType = obj.Weights;

            % Store aggregated or representative weights?
            % For staggered, we have multiple omegas.
            % We will store them in Details.
            res.Details = aggreg_details;

            % Construct Cohort Table
            if isfield(aggreg_details, 'Cohorts') && ~isempty(aggreg_details.Cohorts)
                nC = numel(aggreg_details.Cohorts);
                Cohort   = strings(nC,1);
                Estimate = zeros(nC,1);
                SE       = zeros(nC,1);
                Weight   = zeros(nC,1);
                N        = zeros(nC,1);

                for k=1:nC
                    c = aggreg_details.Cohorts(k);
                    Cohort(k)   = string(c.G);
                    Estimate(k) = c.Tau;
                    if isfield(c, 'SE')
                        SE(k) = c.SE;
                    else
                        SE(k) = NaN;
                    end
                    Weight(k) = c.Weight;
                    N(k)      = numel(c.Units);
                end

                tStat  = Estimate ./ SE;
                pValue = 2 * (1 - normcdf(abs(tStat)));

                res.CohortTable = table(Cohort, N, Estimate, SE, tStat, pValue, Weight);
            end

            % Backward compatibility for single block (store first cohort's weights if only 1?)
            if numel(aggreg_details.Cohorts) == 1
                res.omega  = aggreg_details.Cohorts(1).Omega;
                res.lambda = aggreg_details.Cohorts(1).Lambda;
            else
                res.omega  = []; % Ambiguous for staggered
                res.lambda = [];
            end

            % Visualization Data
            % Reconstruct Aggregate Curves (weighted by cohort size)
            if numel(aggreg_details.Cohorts) > 0
                T = size(Y, 2);
                y_tr_agg = zeros(T, 1);
                y_sc_agg = zeros(T, 1);
                total_N = 0;

                for k = 1:numel(aggreg_details.Cohorts)
                    c = aggreg_details.Cohorts(k);
                    ng = length(c.Units);

                    % Cohort Treated Average
                    y_tr_g = mean(Y(c.Units, :), 1)'; % (T x 1)

                    % Cohort Synthetic Control
                    % Y_co * omega_g
                    % Note: Y_co is Y(c.Controls, :)
                    Y_co_g = Y(c.Controls, :);
                    y_sc_g = Y_co_g' * c.Omega;

                    % Accumulate
                    y_tr_agg = y_tr_agg + y_tr_g * ng;
                    y_sc_agg = y_sc_agg + y_sc_g * ng;
                    total_N  = total_N + ng;
                end

                if total_N > 0
                    y_tr_agg = y_tr_agg / total_N;
                    y_sc_agg = y_sc_agg / total_N;
                end

                res.Curves.Treated   = y_tr_agg;
                res.Curves.Synthetic = y_sc_agg;
                res.Diagnostics.ControlIds = []; % Ambiguous for staggered
                res.Diagnostics.Omega = [];
                res.Diagnostics.is_pre = []; % Varies by cohort
            else
                res.Curves.Treated = [];
                res.Curves.Synthetic = [];
            end

            res.Curves.Time = time_map;

            % Summary Table
            nm = "SDID";
            if strcmpi(obj.Weights, "SC"), nm="SynthControl"; end

            res.coef = table(nm, tau, se, tau/se, 2*normcdf(-abs(tau/se)), 'VariableNames',{'Name','Estimate','SE','tStat','pValue'});
            res.summaryTable = res.coef;

            if obj.Display
                disp(res.summaryTable);
                if numel(aggreg_details.Cohorts) > 1
                    fprintf('Aggregated over %d cohorts.\n', numel(aggreg_details.Cohorts));
                end
            end
        end

        function [tau_att, se_att, details] = estimate_aggregate(obj, Y, D_mat, id_map, time_map)
            % Core estimation logic handling Staggered Adoption

            [N, T] = size(Y);

            % Identify Never-Treated (Control Pool)
            % Units that are NEVER treated in the window
            is_ever_treated = any(D_mat == 1, 2);
            control_indices = find(~is_ever_treated);

            if isempty(control_indices)
                error('did:SDID:NoControls', 'No never-treated units found. SDID requires a clean control pool.');
            end

            % Identify Cohorts among Treated
            treated_indices = find(is_ever_treated);

            % For each treated unit, find first treatment time
            start_times = zeros(length(treated_indices), 1);
            for i = 1:length(treated_indices)
                idx = treated_indices(i);
                t_event = find(D_mat(idx, :) == 1, 1, 'first');
                if isempty(t_event), t_event = Inf; end
                start_times(i) = t_event;
            end

            % Filter out weird cases?
            valid_trt = isfinite(start_times);
            treated_indices = treated_indices(valid_trt);
            start_times = start_times(valid_trt);

            unique_starts = unique(start_times);

            % --- Point Estimate Loop ---
            cohort_res = struct();
            weights_att = [];
            taus = [];

            % Static params
            wType = obj.Weights;
            solv  = obj.Solver;

            % Auto-Regularization Pre-calc
            sigma_hat = 0;
            if isnan(obj.Regularization)
                % diff(Y) over time for controls
                % Y_co is Y(control_indices, :)
                Y_co_pool = Y(control_indices, :);
                dY = diff(Y_co_pool, 1, 2);
                sigma_hat = std(dY(:), 'omitnan');
                if sigma_hat == 0, sigma_hat = 1e-6; end
            end

            for g_idx = 1:length(unique_starts)
                g_time = unique_starts(g_idx);

                % Units in this cohort
                units_in_g = treated_indices(start_times == g_time);

                % Resolve Zeta
                if isnan(obj.Regularization)
                    ng = length(units_in_g);
                    % T_post relative to g_time
                    % We need to know T_post construction below
                    t_post_count = (T - g_time + 1);
                    % Or use is_post definition below which is derived from g_time
                    zeta = sigma_hat * (ng * t_post_count)^0.25;
                else
                    zeta = obj.Regularization;
                end

                % Define Pre/Post relative to g_time
                is_post = false(T, 1);
                is_post(g_time:end) = true;
                is_pre  = ~is_post;

                if sum(is_pre) < 1 || sum(is_post) < 1
                    continue; % Skip if no pre or post periods
                end

                % Data for Estimation
                Y_trt = Y(units_in_g, :);
                Y_co  = Y(control_indices, :);

                % Target: Average of Cohort
                Y_pre_trt_avg = mean(Y_trt(:, is_pre), 1)';
                Y_pre_co      = Y_co(:, is_pre)';

                % Estimate
                [tau_g, omega_g, lambda_g] = did.estimators.SDID.estimate_sdid_static(...
                    Y_pre_trt_avg, Y_pre_co, Y_co, Y_trt, is_pre, is_post, wType, zeta, solv);

                % Weight for Aggregation (Arkhangelsky Appendix A)
                % Weight = N_g * T_post_g ?
                % Actually usually weighted by N_g (if we care about unit-level ATT)
                % or N_g * T_post (if we care about person-period ATT).
                % Standard DiD usually targets person-period.
                ng = length(units_in_g);
                tp = sum(is_post);
                wg = ng * tp;

                % Store
                cR.G = g_time;
                cR.Tau = tau_g;
                cR.Weight = wg;
                cR.Omega = omega_g;
                cR.Lambda = lambda_g;
                cR.Units = units_in_g;
                cR.Controls = control_indices;
                cR.is_pre = is_pre;

                cohort_res(g_idx).Result = cR; % Assign to struct array

                taus(end+1) = tau_g;
                weights_att(end+1) = wg;
            end

            % Aggregate ATT
            if isempty(taus)
                tau_att = NaN;
            else
                tau_att = sum(taus .* weights_att) / sum(weights_att);
            end

            details.Cohorts = [cohort_res.Result]; % flatten
            details.WeightsATT = weights_att;

            % --- Inference (Bootstrap/Placebo) ---
            % We must replicate the *aggregation* procedure.
            se_att = NaN;

            if obj.B > 0 && ismember(lower(obj.SEMethod), ["placebo","clustered"])
                % Placebo on Controls:
                % We assign the *actual* treatment histories (start times) to
                % random sets of control units.

                N_co = length(control_indices);
                N_tr = length(treated_indices);

                % If we don't have enough controls to mimic the treated size,
                % we might have issues. But usually N_co >> N_tr in SC.
                % If N_co is small, we sample with replacement?
                % Arkhangelsky typically assumes N_co is large.

                % Auto-Regularization Pre-calc for Inference
                % We need to use the SAME regularization schedule as in the point estimate.
                % Sigma_hat is roughly constant (based on full control pool deviation).
                sigma_hat_inf = 0;
                if isnan(obj.Regularization)
                    diff_Y_inf = diff(Y_co_pool, 1, 2);
                    sigma_hat_inf = std(diff_Y_inf(:), 'omitnan');
                    if sigma_hat_inf == 0, sigma_hat_inf = 1e-6; end
                end

                % Pre-calculate cohort structure (with Zeta)
                cohort_struct = struct('g_time', {}, 'n_units', {}, 'is_pre', {}, 'is_post', {}, 'weight', {}, 'zeta', {});
                for k = 1:length(details.Cohorts)
                    c = details.Cohorts(k);
                    cs.g_time = c.G;
                    cs.n_units = length(c.Units);
                    cs.is_pre = c.is_pre;
                    cs.is_post = ~c.is_pre;
                    cs.weight  = c.Weight;

                    % Zeta Calculation matching point estimate
                    if isnan(obj.Regularization)
                        % t_post_count is sum(cs.is_post)
                        t_post_inf = sum(cs.is_post);
                        cs.zeta = sigma_hat_inf * (cs.n_units * t_post_inf)^0.25;
                    else
                        cs.zeta = obj.Regularization;
                    end

                    cohort_struct(k) = cs;
                end


                n_coh = length(cohort_struct);
                % Y_co_pool = Y(control_indices, :); % Already defined above? No, check scope.
                % control_indices defined in outer scope.

                taus_boot_agg = NaN(obj.B, 1);
                taus_boot_coh = NaN(obj.B, n_coh);

                wType = obj.Weights;
                % zeta  = obj.Regularization; % REMOVED: Use cs.zeta
                solv  = obj.Solver;

                use_p = obj.UseParallel;
                num_B = obj.B;

                if use_p
                    parfor b = 1:num_B
                        try
                            % Simplified Placebo Logic ...
                            t_agg_b = 0;
                            w_agg_b = 0;
                            t_coh_b = zeros(1, n_coh);

                            for k = 1:n_coh
                                cs = cohort_struct(k);
                                perm = randperm(N_co);
                                idx_trt = perm(1:cs.n_units);
                                idx_co  = perm(cs.n_units+1:end);

                                Y_p_trt = Y_co_pool(idx_trt, :);
                                Y_p_co  = Y_co_pool(idx_co, :);

                                y_p_target = mean(Y_p_trt(:, cs.is_pre), 1)';
                                y_p_co_pre = Y_p_co(:, cs.is_pre)';

                                % Use cs.zeta
                                [t_p, ~, ~] = did.estimators.SDID.estimate_sdid_static(...
                                    y_p_target, y_p_co_pre, Y_p_co, Y_p_trt, ...
                                    cs.is_pre, cs.is_post, wType, cs.zeta, solv);

                                t_coh_b(k) = t_p;
                                t_agg_b = t_agg_b + t_p * cs.weight;
                                w_agg_b = w_agg_b + cs.weight;
                            end

                            taus_boot_agg(b) = t_agg_b / w_agg_b;
                            taus_boot_coh(b, :) = t_coh_b;
                        catch
                        end
                    end
                else
                    % Serial Loop
                    for b = 1:num_B
                        try
                            t_agg_b = 0;
                            w_agg_b = 0;
                            t_coh_b = zeros(1, n_coh);

                            for k = 1:n_coh
                                cs = cohort_struct(k);
                                perm = randperm(N_co);
                                idx_trt = perm(1:cs.n_units);
                                idx_co  = perm(cs.n_units+1:end);

                                Y_p_trt = Y_co_pool(idx_trt, :);
                                Y_p_co  = Y_co_pool(idx_co, :);

                                y_p_target = mean(Y_p_trt(:, cs.is_pre), 1)';
                                y_p_co_pre = Y_p_co(:, cs.is_pre)';

                                [t_p, ~, ~] = did.estimators.SDID.estimate_sdid_static(...
                                    y_p_target, y_p_co_pre, Y_p_co, Y_p_trt, ...
                                    cs.is_pre, cs.is_post, wType, cs.zeta, solv);

                                t_coh_b(k) = t_p;
                                t_agg_b = t_agg_b + t_p * cs.weight;
                                w_agg_b = w_agg_b + cs.weight;
                            end

                            taus_boot_agg(b) = t_agg_b / w_agg_b;
                            taus_boot_coh(b, :) = t_coh_b;
                        catch
                        end
                    end
                end

                if any(~isnan(taus_boot_agg))
                    se_att = std(taus_boot_agg, 'omitnan');

                    % Cohort SEs
                    se_coh = std(taus_boot_coh, 'omitnan');
                    for k = 1:n_coh
                        details.Cohorts(k).SE = se_coh(k);
                    end
                else
                    for k = 1:n_coh
                        details.Cohorts(k).SE = NaN;
                    end
                end
            end
        end
        function [Y, N, T, id_rec, time_rec, D_mat] = pivotData(obj, ds)
            % Convert ds table to Wide Matrix Y (N x T)
            % Assumes balanced panel for now or fills NaNs? SDID needs balanced usually.

            tab = ds.T; % Access underlying table

            % Factorize ID and Time
            [id_idx, id_rec] = findgroups(tab.(ds.idVar));
            [t_idx, time_rec] = findgroups(tab.(ds.timeVar));

            N = max(id_idx);
            T = max(t_idx);

            Y = NaN(N, T);
            D_mat = zeros(N, T);

            % Linear indexing
            idx = sub2ind([N, T], id_idx, t_idx);

            Y(idx) = tab.(ds.yVar);
            D_mat(idx) = tab.(ds.dVar);

            % Filling NaNs?
            % For now, error if NaNs in Y
            if any(isnan(Y), 'all')
                warning('did:SDID:Unbalanced', 'Panel is unbalanced. Filling missing with 0 (dangerous) or interpolating? Keeping NaNs implies solver failure.');
                Y(isnan(Y)) = 0; % TODO: Improving handling
            end
        end

        function [tau, omega, lambda] = estimate_sdid(obj, Y_pre_trt_avg, Y_pre_co, Y_co, Y_trt, is_pre, is_post)
            % Wrapper for static method to allow class property access easier
            [tau, omega, lambda] = did.estimators.SDID.estimate_sdid_static(...
                Y_pre_trt_avg, Y_pre_co, Y_co, Y_trt, is_pre, is_post, ...
                obj.Weights, obj.Regularization, obj.Solver);
        end
    end

    methods(Static)
        function plot(res)
            % SDID.plot(res)
            % Visualizes the results of an SDID/SC estimation.
            % Panels:
            % 1. Trends: Treated vs Synthetic
            % 2. Gap: Difference
            % 3. Unit Weights (Top contributors)

            if ~isfield(res, 'Method') || res.Method ~= "SDID"
                error('did:SDID:plot:InvalidInput', 'Input must be a result struct from SDID or SC.');
            end

            if ~isfield(res, 'Curves')
                error('did:SDID:plot:NoData', 'Result struct does not contain curve data. Re-run fit.');
            end

            % Unpack
            T_axis = res.Curves.Time; % Time vector
            Y_tr   = res.Curves.Treated;
            Y_sc   = res.Curves.Synthetic;
            omega  = res.Diagnostics.Omega;
            ids    = res.Diagnostics.ControlIds;
            is_pre = res.Diagnostics.is_pre;

            % Identify treatment start (first post period)
            % For Staggered, is_pre is empty -> x_line = NaN
            if isempty(is_pre)
                x_line = NaN;
            else
                t_start_idx = find(~is_pre, 1, 'first');
                if ~isempty(t_start_idx)
                    x_line = T_axis(t_start_idx);
                else
                    x_line = NaN;
                end
            end

            figure('Name', sprintf('SDID Results (%s)', res.WeightsType), 'Color', 'w');
            tiledlayout(2,2, 'TileSpacing','compact');

            % --- 1. Trends ---
            nexttile([1 2]);
            hold on;
            plot(T_axis, Y_tr, '-k', 'LineWidth', 1.5, 'DisplayName', 'Treated');
            plot(T_axis, Y_sc, '--r', 'LineWidth', 1.5, 'DisplayName', 'Synthetic');
            if ~isnan(x_line)
                xline(x_line, ':', 'Treatment', 'LabelVerticalAlignment','bottom');
            end
            legend('Location','best');
            title('Trends: Treated vs Synthetic');
            xlabel('Time'); ylabel('Outcome');
            grid on;
            hold off;


            % --- 2. Gap ---
            nexttile;
            gap = Y_tr - Y_sc;
            yline(0, '-k');
            hold on;
            plot(T_axis, gap, '-b', 'LineWidth', 1.2);
            if ~isnan(x_line)
                xline(x_line, ':');
            end
            title('Gap (Treated - Synthetic)');
            xlabel('Time');
            grid on;

            % --- 3. Weights ---
            nexttile;

            if isempty(omega)
                % Staggered case: Weights are per-cohort, not aggregated easily yet.
                axis off;
                text(0.5, 0.5, {'Weights not available', 'for Staggered Aggregation'}, ...
                    'HorizontalAlignment', 'center', 'Units', 'normalized');
                title('Top Unit Weights');
            else
                % Sort weights
                [w_sort, idx] = sort(omega, 'descend');

                % Show top 10 or those > 0.01
                mask = w_sort > 0.01;
                if sum(mask) < 5
                    mask(1:min(5, length(mask))) = true; % Show at least 5
                end

                w_top = w_sort(mask);
                id_top = ids(idx(mask));

                % Convert ids to string for categorical axis
                str_ids = string(id_top);

                bar(w_top, 'FaceColor', [0.2 0.4 0.6]);
                title('Top Unit Weights');
                set(gca, 'XTick', 1:length(w_top), 'XTickLabel', str_ids, 'XTickLabelRotation', 45);
                ylabel('Weight');
            end
        end

        function [tau, omega, lambda] = estimate_sdid_static(Y_pre_trt_avg, Y_pre_co, Y_co, Y_trt, is_pre, is_post, wType, zeta, solv)
            % 1. Unit Weights (SC)
            % min || Y_pre_trt - Y_pre_co * omega ||

            N_co = size(Y_co, 1);
            T_pre = sum(is_pre);
            T_post = sum(is_post);

            if strcmpi(wType, "Time") % Only time weights, Unit weights = 1/N
                omega = ones(N_co, 1) / N_co;
            else
                % Solve Simplex
                % A = Y_pre_co (T_pre x N_co)
                % b = Y_pre_trt_avg (T_pre x 1)
                omega = did.optimization.solveSimplex(Y_pre_co, Y_pre_trt_avg, 'Zeta', zeta, 'Solver', solv);
            end

            % 2. Time Weights (SDID)
            % min || Y_pre_co_avg - Y_pre_co' * lambda || ??
            % Actually: we want to match the Control Trends to be constant?
            % Arkhangelsky: Match Y_post_co (average) using pre periods.
            % Target: constant?

            if strcmpi(wType, "SC") || strcmpi(wType, "Unit")
                % SC: Lambda = 1/T_post for post, 0 otherwise?
                % Pure SC just compares Post Trt to Post Synth.
                % Implicit lambda is uniform on post.
                lambda = zeros(length(is_pre), 1);
                lambda(is_post) = 1/T_post;
            else
                % SDID Time Weights
                % Target: Constant 1 (intercept)? Or
                % "Find lambda_pre such that Y_co(pre)' * lambda_pre approx Y_co(post)' * 1/T_post"

                % Y_co is (N_co x T).
                % Y_co_pre = Y_co(:, is_pre)  (N_co x T_pre)
                % Y_co_post = Y_co(:, is_post) (N_co x T_post)

                Y_post_mean = mean(Y_co(:, is_post), 2); % (N_co x 1)

                % A = Y_co(:, is_pre)'  (T_pre x N_co)? No.
                % We want lambda (T_pre x 1) such that:
                % Y_co(:, is_pre) * lambda approx Y_post_mean

                % A = Y_co(:, is_pre)  (N_co x T_pre)
                % b = Y_post_mean      (N_co x 1)

                lambda_pre = did.optimization.solveSimplex(Y_co(:, is_pre), Y_post_mean, 'Zeta', zeta, 'Solver', solv);

                lambda = zeros(length(is_pre), 1);
                lambda(is_pre)  = lambda_pre;
                lambda(is_post) = 1/T_post; % Fixed for post periods? Wait, usually we compare simple avgs in post.
            end

            % 3. Estimate
            % Tau = (Y_trt(post) - Y_trt(pre)*lambda) - (Y_co(post)*omega - Y_co(pre)*omega*lambda) ??
            % Simplified:
            % Diff-in-Diffs with weights.

            % Weighted Means
            mu_trt_post = mean(mean(Y_trt(:, is_post)));

            % mu_trt_pre (weighted by lambda)
            % Y_trt is (N_trt x T). Average over N_trt first?
            y_trt_vec = mean(Y_trt, 1)'; % (T x 1)
            mu_trt_pre = sum(y_trt_vec(is_pre) .* lambda(is_pre));

            % Controls
            % Y_co is (N_co x T).
            % Weight by omega first (Synth Control Unit)
            y_sc_vec = Y_co' * omega; % (T x 1)

            mu_co_post = mean(y_sc_vec(is_post)); % (Simple mean if lambda_post uniform)
            mu_co_pre  = sum(y_sc_vec(is_pre) .* lambda(is_pre));

            % SC special case: if lambda_pre is not optimizing, just mean?
            % If SC, lambda(is_pre) should be 0? No, SC is Y_trt_post - Y_sc_post.
            % So SC implies "DiD" term is 0?
            % Arkhangelsky: SC is just Y_trt_post - Y_sc_post.
            % This implies mu_trt_pre and mu_co_pre are ignored or cancel out?
            % Actually SC assumes pre-trends match perfectly so difference is 0.

            if strcmpi(wType, "SC")
                tau = mu_trt_post - mu_co_post;
            else
                % SDID (Double Diff)
                tau = (mu_trt_post - mu_trt_pre) - (mu_co_post - mu_co_pre);
            end
        end
    end
end
