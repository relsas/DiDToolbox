classdef HonestDiD
    % HONESTDID  Sensitivity analysis for Event Studies (Rambachan & Roth 2023).
    % Wraps did.diagnostics.honest_pt to provide Fixed Length Confidence Intervals (FLCI)
    % and sensitivity plots.
    %
    % Usage:
    %   hd = did.diagnostics.HonestDiD(beta, sigma, numPre, numPost);
    %   out = hd.fit(M=0:0.01:2, Target="e=5");
    %   hd.plot(out);
    %

    properties
        Betas       (:,1) double
        Sigma       (:,:) double
        NumPre      (1,1) double
        NumPost     (1,1) double
        EventTimes  (:,1) double

        % Derived
        PreSE       (:,1) double
        PostSE      (:,1) double
    end

    methods
        function obj = HonestDiD(beta, sigma, numPre, numPost, opts)
            arguments
                beta    (:,1) double
                sigma   (:,:) double
                numPre  (1,1) double
                numPost (1,1) double
                opts.EventTimes (:,1) double = []
            end

            if numel(beta) ~= numPre + numPost
                error('HonestDiD:SizeMismatch', 'Beta length must be NumPre + NumPost.');
            end

            obj.Betas   = beta;
            obj.Sigma   = sigma;
            obj.NumPre  = numPre;
            obj.NumPost = numPost;

            % Extract diagonals for PreSE
            se_all = sqrt(diag(sigma));
            obj.PreSE  = se_all(1:numPre);
            obj.PostSE = se_all(numPre+1:end);

            if isempty(opts.EventTimes)
                obj.EventTimes = [(-numPre:-1) 0:(numPost-1)]';
            else
                obj.EventTimes = opts.EventTimes;
            end
        end

        function res = fit(obj, varargin)
            % FIT  Run Honest PT analysis.
            % Options:
            %   M (double) : vector of M values (default 0:0.01:2)
            %   DeltaType (string) : "SD", "RM", "SDRM" (default "SD")
            %   Alpha (double) : significance level (default 0.05)
            %   Target (string or vector) : "e=5", or weight vector l_vec.

            % Display   = opts.Display; (already parsed)

            p = inputParser;
            addParameter(p, 'M', 0:0.05:2);
            addParameter(p, 'DeltaType', "RM");
            addParameter(p, 'Alpha', 0.05);
            addParameter(p, 'PreBand', 1.96); % Default generic
            addParameter(p, 'Target', "last"); % last, average, or e=X
            addParameter(p, 'Display', true);
            parse(p, varargin{:});

            Mvec  = p.Results.M;
            alpha = p.Results.Alpha;
            dtype = p.Results.DeltaType;
            dispOn= p.Results.Display;

            % Construct l_vec from Target
            l_vec = obj.parseTarget(p.Results.Target);

            % Compute scalar SE for target theta
            % Var(l' beta_post) = l' Sigma_post l
            sig_post = obj.Sigma(obj.NumPre+1:end, obj.NumPre+1:end);
            theta_se = sqrt(l_vec' * sig_post * l_vec);

            % Call engine
            eng = did.diagnostics.honest_pt(obj.Betas, obj.NumPre, obj.NumPost, ...
                'M', Mvec, ...
                'DeltaType', dtype, ...
                'l_vec', l_vec, ...
                'eventTimes', obj.EventTimes, ...
                'preSE', obj.PreSE, ...
                'PreBand', p.Results.PreBand, ...
                'ThetaSE', theta_se, ...
                'Display', dispOn, ...
                'StoreDelta', false);

            % Post-process for Inference (Naive FLCI)
            cv = norminv(1 - alpha/2);

            % Bounds of Identified Set
            lb_id = eng.theta_lb;
            ub_id = eng.theta_ub;

            % FLCI (Naive Imbens-Manski style worst-case)
            % CI = [ LB - cv*SE, UB + cv*SE ]
            % Note: R&R FLCI uses a specific CV derived from M, but for fixed M
            % typically we assume bias is worst-case constant and add max noise.
            ci_lo = lb_id - cv * theta_se;
            ci_hi = ub_id + cv * theta_se;

            % Pack results
            res = struct();

            % Filter out NaNs (infeasible M)
            valid_idx = ~isnan(ci_lo) & ~isnan(ci_hi);

            res.M = reshape(Mvec(valid_idx), [], 1);
            res.LB_ID = reshape(lb_id(valid_idx), [], 1);
            res.UB_ID = reshape(ub_id(valid_idx), [], 1);
            res.CI_Lo = reshape(ci_lo(valid_idx), [], 1);
            res.CI_Hi = reshape(ci_hi(valid_idx), [], 1);

            % Create summary table for easy display
            res.summaryTable = table(res.M, res.LB_ID, res.UB_ID, res.CI_Lo, res.CI_Hi, ...
                'VariableNames', {'M', 'LB_ID', 'UB_ID', 'CI_Lo', 'CI_Hi'});

            res.ThetaHat = eng.theta_hat;
            res.OriginalCI = eng.ci95;
            res.Estimator = eng; % full engine output
            res.Alpha = alpha;
        end

        function plot(obj, res, ax)
            % PLOT  Sensitivity plot (M vs Estimates)
            %   obj.plot(res)       -> Creates new figure
            %   obj.plot(res, ax)   -> Plots into specific axes handle

            if nargin < 3 || isempty(ax)
                f = figure('Color','w');
                ax = axes(f);
            end

            hold(ax,'on');

            % Plot Identified Set (Shaded?) or Lines
            % 1. CI Area (Light Gray)
            x = [res.M; flipud(res.M)];
            y = [res.CI_Lo; flipud(res.CI_Hi)];
            valid = ~isnan(y);
            fill(x(valid), y(valid), [0.9 0.9 0.9], 'EdgeColor','none', 'DisplayName', sprintf('%.0f%% Confidence', (1-res.Alpha)*100));

            % 2. Identification Region (Darker Gray)
            y_id = [res.LB_ID; flipud(res.UB_ID)];
            valid_id = ~isnan(y_id);
            fill(x(valid_id), y_id(valid_id), [0.7 0.7 0.7], 'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName', 'Identified Set');

            % 3. Original Estimate Line
            yline(ax, res.ThetaHat, '--k', 'Original Estimate');
            yline(ax, 0, '-', 'Color', [0.5 0 0]);

            xlabel(ax, 'Sensitivity Parameter M');
            ylabel(ax, 'Theta (Treatment Effect)');
            title(ax, 'Honest DiD Sensitivity Analysis');
            legend(ax, 'Location','Best');
            grid(ax, 'on');
        end
    end

    methods (Access=private)
        function l_vec = parseTarget(obj, target)
            % Returns NumPost x 1 vector
            if isnumeric(target)
                l_vec = target(:);
                if numel(l_vec) ~= obj.NumPost
                    error('Target vector length mismatch');
                end
                return;
            end

            % String parsing
            switch string(target)
                case "last"
                    l_vec = zeros(obj.NumPost, 1);
                    l_vec(end) = 1;
                case "average"
                    l_vec = ones(obj.NumPost, 1) ./ obj.NumPost;
                case "first"
                    l_vec = zeros(obj.NumPost, 1);
                    l_vec(1) = 1;
                otherwise
                    error('Unknown target: %s', target);
            end
        end
    end



    methods (Static)
        function obj = fromFit(res)
            % FROMFIT  Factory: Create from did.fit result.
            %   hd = did.diagnostics.HonestDiD.fromFit(res);
            %
            % Supports:
            %   1. CS result struct (must have Aggregates.es)
            %   2. Table with 'e'/'EventTime', 'Estimate', 'SE' columns.

            beta = []; sigma = []; eventTimes = [];

            % Case 1: CS Result Struct
            if isstruct(res) && isfield(res, 'Aggregates') && isfield(res.Aggregates, 'es')
                T = res.Aggregates.es;
                % CS table usually has 'e' or 'EventTime'
                if ismember('e', T.Properties.VariableNames)
                    ename = 'e';
                elseif ismember('EventTime', T.Properties.VariableNames)
                    ename = 'EventTime';
                else
                    error('HonestDiD:fromFit:CS_NoEventTime','CS Aggregates.es table missing "e" column.');
                end

                T = sortrows(T, ename);
                eventTimes = T.(ename);
                beta = T.Estimate;

                % SE / Sigma
                if ismember('VCov', T.Properties.VariableNames)
                    se = T.SE;
                    sigma = diag(se.^2);
                elseif ismember('SE', T.Properties.VariableNames)
                    se = T.SE;
                    sigma = diag(se.^2);
                else
                    error('HonestDiD:fromFit:NoSE','No SE in CS/ES table.');
                end

                % Case 1b: Result struct with 'EventStudy' (Wooldridge/New DID_M)
            elseif isstruct(res) && isfield(res, 'EventStudy')
                T = res.EventStudy;
                if ismember('EventTime', T.Properties.VariableNames)
                    ename = 'EventTime';
                elseif ismember('e', T.Properties.VariableNames)
                    ename = 'e';
                else
                    error('HonestDiD:fromFit:NoEventTime','EventStudy table missing "EventTime" column.');
                end
                T = sortrows(T, ename);
                eventTimes = T.(ename);
                beta = T.Estimate;
                se   = T.SE;
                sigma = diag(se.^2);

                % Case 1c: Wooldridge details.attByEventTime (Legacy/Specific)
            elseif isstruct(res) && isfield(res,'details') && isfield(res.details,'attByEventTime')
                T = res.details.attByEventTime;
                T = sortrows(T, 'EventTime');
                eventTimes = T.EventTime;
                beta = T.Estimate;
                se   = T.SE;
                sigma = diag(se.^2);

                % Case 2: Generic Table
            elseif istable(res)
                T = res;
                if ismember('EventTime', T.Properties.VariableNames)
                    T = sortrows(T, 'EventTime');
                    eventTimes = T.EventTime;
                    beta = T.Estimate;
                    se = T.SE;
                    sigma = diag(se.^2);
                elseif ismember('e', T.Properties.VariableNames)
                    T = sortrows(T, 'e');
                    eventTimes = T.e;
                    beta = T.Estimate;
                    se = T.SE;
                    sigma = diag(se.^2);
                else
                    error('HonestDiD:fromFit:UnknownTable', ...
                        'Table must have "EventTime" (or "e"), "Estimate", "SE".');
                end
            else
                error('HonestDiD:fromFit:InvalidInput', ...
                    'Input must be CS result struct or EventStudy table.');
            end

            % Remove NaN values (e.g. incomplete event study bins)
            validMask = ~isnan(beta) & ~isnan(se) & ~isnan(eventTimes);
            if any(~validMask)
                warning('HonestDiD:fromFit:NaNsRemoved', ...
                    'Removed %d periods with NaN estimates or SEs.', sum(~validMask));
                beta = beta(validMask);
                sigma = sigma(validMask, validMask);
                eventTimes = eventTimes(validMask);
            end

            % Determine NumPre / NumPost

            % Determine NumPre / NumPost
            numPre  = sum(eventTimes < 0);
            numPost = sum(eventTimes >= 0);

            if numPost == 0
                error('HonestDiD:fromFit:NoPost','No post-period coefficients found (e>=0).');
            end

            obj = did.diagnostics.HonestDiD(beta, sigma, numPre, numPost, EventTimes=eventTimes);
        end
    end
end
