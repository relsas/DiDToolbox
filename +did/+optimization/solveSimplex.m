function [w, exitflag, out] = solveSimplex(A, b, varargin)
% did.optimization.solveSimplex  Solve min ||Ax - b||^2 + zeta*||x||^2 s.t. sum(x)=1, x>=0
%
% Usage:
%   w = solveSimplex(A, b);
%   w = solveSimplex(A, b, 'Zeta', 1e-6, 'Solver', 'quadprog');
%
% Inputs:
%   A      : (M x N) Design matrix
%   b      : (M x 1) Target vector
%   Zeta   : (scalar) L2 regularization parameter (default 0)
%   Solver : "auto", "quadprog", "frank-wolfe" (default "auto")
%   Tol    : Tolerance (default 1e-6)
%   MaxIter: Max iterations (default 1000 for FW)
%   Display: "off" or "iter"

p = inputParser;
addParameter(p, 'Zeta', 0, @(x) isscalar(x) && x>=0);
addParameter(p, 'Solver', 'auto', @(x) any(validatestring(x, {'auto','quadprog','frank-wolfe'})));
addParameter(p, 'Tol', 1e-6, @isscalar);
addParameter(p, 'MaxIter', 1000, @isscalar);
addParameter(p, 'Display', 'off');
parse(p, varargin{:});
opts = p.Results;

[M, N] = size(A);

% Decide solver
useQP = false;
if strcmpi(opts.Solver, 'quadprog')
    useQP = true;
elseif strcmpi(opts.Solver, 'auto')
    % Use quadprog if available and N is not huge
    if exist('quadprog', 'file')==2 && N < 2000
        useQP = true;
    end
end

if useQP
    [w, exitflag, out] = solve_quadprog(A, b, opts.Zeta, opts);
else
    [w, exitflag, out] = solve_frankwolfe(A, b, opts.Zeta, opts);
end
end

function [w, exitflag, output] = solve_quadprog(A, b, zeta, opts)
% Problem: min (Ax-b)'(Ax-b) + z*x'x
%        = x'A'Ax - 2b'Ax + b'b + z*x'x
%        = 0.5 * x'(2(A'A + zI))x + (-2A'b)'x

H = 2 * (A'*A);
if zeta > 0
    H = H + 2*zeta*eye(size(H,1));
end
f = -2 * (A'*b);

% Constraints: sum(x)=1 => Aeq x = beq
Aeq = ones(1, size(A,2));
beq = 1;
% x >= 0
lb = zeros(size(A,2), 1);

qpOpts = optimoptions('quadprog', 'Display', opts.Display, ...
    'OptimalityTolerance', opts.Tol);

[w, val, exitflag, output] = quadprog(H, f, [], [], Aeq, beq, lb, [], [], qpOpts);
end

function [w, exitflag, output] = solve_frankwolfe(A, b, zeta, opts)
% Frank-Wolfe for Simplex Constraints
% L(w) = ||Aw - b||^2 + zeta ||w||^2
% Grad L(w) = 2 A'(Aw - b) + 2 zeta w

N = size(A, 2);

% 1. Init (Uniform start)
w = ones(N, 1) / N;

exitflag = 0;

if strcmpi(opts.Display, 'iter')
    fprintf('Frank-Wolfe (SDID Simplex)\n');
    fprintf('Iter     ObjVal        Gap\n');
end

% Precompute? No, A is potentially large, A'A is expensive.
% We compute gradients iteratively.

for k = 1:opts.MaxIter
    % Gradient
    % r = Aw - b
    r = A*w - b;
    grad = 2 * (A'*r + zeta*w);

    % Linear Oracle: s = min_{s in Simplex} <s, grad>
    % This is just the unit vector at index of min gradient
    [~, min_idx] = min(grad);
    s = zeros(N,1);
    s(min_idx) = 1;

    % Direction
    d = s - w;

    % Duality Gap (Stopping Criterion)
    % Gap = <grad, w - s>
    gap = grad' * (-d); % Note: d = s-w => -d = w-s

    if strcmpi(opts.Display, 'iter') && mod(k,10)==0
        obj = sum(r.^2) + zeta * sum(w.^2);
        fprintf('%4d  %10.4e  %10.4e\n', k, obj, gap);
    end

    if gap < opts.Tol
        exitflag = 1;
        break;
    end

    % Line Search: min_{gamma in [0,1]} L(w + gamma d)
    % L(w+gd) = || A(w+gd) - b ||^2 + z || w+gd ||^2
    % Let v = Ad.
    % min || r + gamma v ||^2 + z || w + gamma d ||^2
    % = ||r||^2 + 2g <r,v> + g^2 ||v||^2 + z(||w||^2 + 2g <w,d> + g^2 ||d||^2)
    % deriv wrt gamma:
    % 2 <r,v> + 2 g ||v||^2 + z(2 <w,d> + 2 g ||d||^2) = 0
    % g (||v||^2 + z ||d||^2) = - (<r,v> + z <w,d>)
    % gamma_opt = - (<r,v> + z <w,d>) / (||v||^2 + z ||d||^2)

    v = A*d;
    numerator = -(r'*v + zeta * (w'*d));
    denominator = (v'*v + zeta * (d'*d));

    if denominator < 1e-12
        gamma = 0;
    else
        gamma = numerator / denominator;
    end

    % Clip to [0,1]
    gamma = max(0, min(1, gamma));

    % Update
    w = w + gamma * d;
end

output.iterations = k;
output.gap = gap;
end
