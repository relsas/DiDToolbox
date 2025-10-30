classdef baconDecomp
    % Bacon decomposition
    % 
    % ------------------------------------------------------------------------
    % Dr. Ralf Elsas-Nicolle, LMU Munich, Germany
    % Last change: 09/30/2025
    % ------------------------------------------------------------------------

    properties
        data table
        timeVar string
        idVar string
        treatVar string
        outcomeVar string
        treatmentTimes double

        % Computed outputs stored for convenience
        Results table
        ResultsAgg table
    end

    methods
        function obj = baconDecomp(T,args)
            arguments
                T table
                args.timeVar (1,1) string
                args.idVar (1,1) string
                args.treatVar (1,1) string
                args.outcomeVar (1,1) string
                
            end

            obj.data       = T;
            obj.timeVar    = args.timeVar;
            obj.idVar      = args.idVar;
            obj.treatVar   = args.treatVar;
            obj.outcomeVar = args.outcomeVar;

            % Replace missing treatment times with 0 (never-treated)
            werNaN = find(isnan(obj.data.(obj.treatVar)));
            if ~isempty(werNaN)
                obj.data.(obj.treatVar)(werNaN) = 0;
            end

            % Unique treatment times (exclude never-treated)
            obj.treatmentTimes = unique(obj.data.(obj.treatVar));
            obj.treatmentTimes = obj.treatmentTimes(obj.treatmentTimes > 0);

            % Compute and store results immediately
            [res, resAgg] = obj.decompose();
            obj.Results   = res;
            obj.ResultsAgg = resAgg;

        end

        function [results, results_aggregate] = decompose(obj)
            % Initialize containers
            n_comparisons = nchoosek(length(obj.treatmentTimes), 2) * 2 + length(obj.treatmentTimes);
            type     = cell(n_comparisons, 1);
            group1   = zeros(n_comparisons, 1);
            group2   = zeros(n_comparisons, 1);
            weight   = zeros(n_comparisons, 1);
            estimate = zeros(n_comparisons, 1);

            idx = 1;

            % Total variance of D (indicator ever-treated-at-time) per observation
            VD = var(obj.data.(obj.treatVar) > 0);

            % Treated vs Never Treated
            for i = 1:length(obj.treatmentTimes)
                time = obj.treatmentTimes(i);
                [est, ~]   = obj.computeTreatedVsNever(time);
                weight(idx) = obj.computeTreatedVsNeverWeight(time, VD);

                type{idx}   = 'Treated vs Never Treated';
                group1(idx) = time;
                group2(idx) = 0;
                estimate(idx) = est;
                idx = idx + 1;
            end

            % Earlier vs Later Treated
            for i = 1:length(obj.treatmentTimes)
                for j = i+1:length(obj.treatmentTimes)
                    time1 = obj.treatmentTimes(i);
                    time2 = obj.treatmentTimes(j);
                    [est, ~]   = obj.computeTwoByTwo(time1, time2);
                    weight(idx) = obj.computeEarlierVsLaterWeight(time1, time2, VD);

                    type{idx}     = 'Treated earlier vs later';
                    group1(idx)   = time1;
                    group2(idx)   = time2;
                    estimate(idx) = est;
                    idx = idx + 1;
                end
            end

            % Later vs Earlier Treated
            for i = 1:length(obj.treatmentTimes)
                for j = 1:i-1
                    time1 = obj.treatmentTimes(i);
                    time2 = obj.treatmentTimes(j);
                    [est, ~]   = obj.computeTwoByTwo(time1, time2);
                    weight(idx) = obj.computeLaterVsEarlierWeight(time1, time2, VD);

                    type{idx}     = 'Treated later vs earlier';
                    group1(idx)   = time1;
                    group2(idx)   = time2;
                    estimate(idx) = est;
                    idx = idx + 1;
                end
            end

            % Normalize weights
            weight = weight / sum(weight);

            % Results table of all 2x2s
            results = table(type, group1, group2, weight, estimate);

            % Aggregate by type
            g_sum  = groupsummary(results, "type", ["sum"], ["weight"]);
            % Weighted mean of estimates using weights (groupsummary supports function handles with multiple vars in recent MATLAB releases)
            g_sum2 = groupsummary(results, "type", @(w,e) (w' * e) ./ sum(w), {"weight","estimate"});

            results_aggregate = g_sum(:, ["type","sum_weight"]);
            results_aggregate.estimate = g_sum2.fun1_weight_estimate;
            results_aggregate.Properties.VariableNames = ["Type","Weight","Estimate"];
        end

        function [estimate, weight] = computeTreatedVsNever(obj, treatTime)
            % 2x2 DID: treated-at-time vs never-treated
            treated = obj.data.(obj.treatVar) == treatTime;
            control = obj.data.(obj.treatVar) == 0;

            prePeriod  = obj.data.(obj.timeVar) < treatTime;
            postPeriod = obj.data.(obj.timeVar) >= treatTime;

            y = obj.data.(obj.outcomeVar);
            treated_pre   = mean(y(treated & prePeriod), 'omitnan');
            treated_post  = mean(y(treated & postPeriod), 'omitnan');
            control_pre   = mean(y(control & prePeriod), 'omitnan');
            control_post  = mean(y(control & postPeriod), 'omitnan');

            estimate = (treated_post - treated_pre) - (control_post - control_pre);
            weight   = 0; % computed elsewhere
        end

        function [estimate, weight] = computeTwoByTwo(obj, time1, time2)
            % 2x2 DID: earlier vs later treated (order-sensitive logic preserved)
            group1 = obj.data.(obj.treatVar) == time1;
            group2 = obj.data.(obj.treatVar) == time2;

            if time1 < time2
                prePeriod     = obj.data.(obj.timeVar) < time1;
                betweenPeriod = (obj.data.(obj.timeVar) >= time1) & (obj.data.(obj.timeVar) < time2);
            else
                prePeriod     = obj.data.(obj.timeVar) >= time1;
                betweenPeriod = (obj.data.(obj.timeVar) >= time2) & (obj.data.(obj.timeVar) < time1);
            end

            y = obj.data.(obj.outcomeVar);
            g1_pre     = mean(y(group1 & prePeriod), 'omitnan');
            g1_between = mean(y(group1 & betweenPeriod), 'omitnan');
            g2_pre     = mean(y(group2 & prePeriod), 'omitnan');
            g2_between = mean(y(group2 & betweenPeriod), 'omitnan');

            if time1 < time2
                estimate = (g1_between - g1_pre) - (g2_between - g2_pre);
            else
                estimate = (g1_pre - g1_between) - (g2_pre - g2_between);
            end
            weight = 0; % computed elsewhere
        end

        function weight = computeTreatedVsNeverWeight(obj, k, VD)
            % Equation (10e)
            treated       = obj.data.(obj.treatVar) == k;
            never_treated = obj.data.(obj.treatVar) == 0;
            isSample      = treated | never_treated;

            nk = sum(treated);
            nu = sum(never_treated);

            post_k = obj.data.(obj.timeVar)(isSample) >= k;
            Dk = mean(post_k);

            nku = nk ./ (nk + nu);
            weight = ((nk + nu)^2 * nku * (1 - nku) * Dk * (1 - Dk)) / VD;
        end

        function weight = computeEarlierVsLaterWeight(obj, k, l, VD)
            % Equation (10f)
            treated_k = obj.data.(obj.treatVar) == k;
            treated_l = obj.data.(obj.treatVar) == l;

            nk = sum(treated_k);
            nl = sum(treated_l);

            post_k = obj.data.(obj.timeVar) >= k;
            post_l = obj.data.(obj.timeVar) >= l;
            Dk = mean(post_k);
            Dl = mean(post_l);

            nkl = nk ./ (nk + nl);
            weight = (((nk + nl) * (1 - Dl))^2) * nkl * (1 - nkl) * ((Dk - Dl) / (1 - Dl)) * ((1 - Dk) / (1 - Dl)) / VD;
        end

        function weight = computeLaterVsEarlierWeight(obj, k, l, VD)
            % Equation (10g)
            treated_k = obj.data.(obj.treatVar) == k;
            treated_l = obj.data.(obj.treatVar) == l;

            nk = sum(treated_k);
            nl = sum(treated_l);

            post_k = obj.data.(obj.timeVar) >= k;
            post_l = obj.data.(obj.timeVar) >= l;
            Dk = mean(post_k);
            Dl = mean(post_l);

            nkl = nk ./ (nk + nl);
            weight = (((nk + nl) * Dl)^2) * nkl * (1 - nkl) * (Dk / Dl) * ((Dl - Dk) / Dl) / VD;
        end
    end
end


