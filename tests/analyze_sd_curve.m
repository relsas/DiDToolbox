function analyze_sd_curve()
% ANALYZE_SD_CURVE  Investigate the shape of the sensitivity bounds.
%
%   Purpose:
%   This script calculates the "Sensitivity Curve" of the Honest DiD method
%   given a theoretical set of coefficients (beta) and standard errors.
%   It demonstrates how the identification region [LB(M), UB(M)] grows
%   as we relax the assumption of parallel trends (M).
%
%   Theoretical Insight:
%   - The bounds LB(M) and UB(M) are the solutions to a Linear Program.
%   - They are piecewise linear and convex/concave with respect to M.
%   - Expansion Logic: As M increases, the constraint |Delta^2 delta| <= M
%     loosens, allowing for larger bias terms (delta). Since we minimize/maximize
%     theta over this larger set, the bounds expand monotonically.
% -------------------------------------------------------------------------
clear; close all; clc;

fprintf('--- Setting up Theoretical Event Study ---\n');
% Setup: 3 pre-periods (-3, -2, -1), 3 post-periods (0, 1, 2)
numPre = 3; numPost = 3;
e = [-3 -2 -1 0 1 2]';

% 1. Define Coefficients (Beta)
% Scenario: True Effect = 1 in post-periods. Parallel trends hold exactly (0 in pre).
beta = [0; 0; 0; 1; 1; 1];

% 2. Define Standard Errors (Sigma) for all coefficients
% We assume a constant standard error for each event-study coefficient.
sigma_val = 0.1;
sigma_vec = repmat(sigma_val, numPre + numPost, 1);

% 3. Define Target Parameter (Theta)
% We target the effect in the LAST period (e=2).
% Corresponding weight vector l_vec for post-periods [0, 1, 2]:
l_vec = [0; 0; 1];

% 4. Calculate Implied Standard Error for Theta
% Theta_hat = l_vec' * beta_post
% Var(Theta_hat) = l_vec' * Cov(beta_post) * l_vec
% Assuming independence (diagonal cov) for this theoretical input:
post_sigma = sigma_vec(numPre+1:end);
ThetaSE = sqrt(sum((l_vec .* post_sigma).^2));

fprintf('Betas:   %s\n', mat2str(beta', 2));
fprintf('ThetaSE: %.4f (derived from sigma=%.2f and target weights)\n', ThetaSE, sigma_val);

% Run honest_pt with SDRM
M_grid = 0:0.1:2.0;

fprintf('\n--- Running optimized Honest DiD (honest_pt) ---\n');
out = did.diagnostics.honest_pt(beta, numPre, numPost, ...
    'M', M_grid, ...
    'DeltaType', 'SDRM', ...
    'l_vec', l_vec, ...
    'preSE', sigma_vec(1:numPre), ...
    'ThetaSE', ThetaSE, ...      % Consistent SE passed to function
    'eventTimes', e, ...
    'Display', true);

% Extract bounds
lb = out.theta_lb;
ub = out.theta_ub;
diff_lb = diff(lb);

fprintf('\n--- Shape Analysis of the Sensitivity Curve ---\n');
fprintf('The "Sensitivity Curve" maps M (allowed violation) -> [LB, UB].\n');
fprintf('Notice the expansion is linear or piecewise linear as constraints relax.\n');
fprintf('\n');
fprintf('M       | LB      | UB      | Marginal Cost (Slope of LB)\n');
fprintf('----------------------------------------------------------\n');
for i=1:min(11, numel(M_grid))
    if i < numel(M_grid)
        if i==1
            slope = (lb(i+1)-lb(i))/(M_grid(i+1)-M_grid(i));
        else
            slope = (lb(i)-lb(i-1))/(M_grid(i)-M_grid(i-1));
        end
    else
        slope = NaN;
    end
    fprintf('%0.2f    | %0.4f | %0.4f | %0.4f\n', ...
        M_grid(i), lb(i), ub(i), slope);
end

end
