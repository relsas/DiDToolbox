function [b, yhat, resid] = fastOLS(y, X)
% FASTOLS  Minimal overhead OLS solver (A\B)
%   [b, yhat, resid] = did.utils.fastOLS(y, X)
%
%   Input:
%       y : (N x 1) Outcome vector
%       X : (N x K) Design matrix (must include intercept if desired)
%
%   Output:
%       b     : (K x 1) Coefficients
%       yhat  : (N x 1) Predicted values (optional)
%       resid : (N x 1) Residuals (optional)
%
%   Note: Does not perform data checks (NaNs, rank deficiency) for speed.
%         Caller must ensure inputs are valid double arrays.

% Solve

% Suppress excessive rank deficiency warnings (common in FE models)
warnState = warning('off', 'MATLAB:rankDeficientMatrix');
cleanup = onCleanup(@() warning(warnState));

b = X \ y;

% Optional outputs
if nargout > 1
    yhat = X * b;
end

if nargout > 2
    resid = y - yhat;
end
end
