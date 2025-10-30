function [X, Xnames] = extractCovariates(T, covarNames, keyvars)
% EXTRACTCOVARIATES  Build numeric covariate matrix X from T.
% covarNames: string array of column names; may be empty
% keyvars (optional): names to forbid (id/time/y/d)
if nargin < 3, keyvars = string.empty; end
Xnames = string(covarNames(:))';
if isempty(Xnames)
    X = []; return;
end

% existence
names = string(T.Properties.VariableNames);
miss  = Xnames(~ismember(Xnames, names));
assert(isempty(miss), 'did:extractCovariates:Missing', ...
    'Covariate(s) not found: %s', strjoin(miss, ', '));

% forbid keys
if ~isempty(keyvars)
    bad = Xnames(ismember(Xnames, keyvars));
    assert(isempty(bad), 'did:extractCovariates:KeyVars', ...
        'Covariates cannot include key vars: %s', strjoin(bad, ', '));
end

% type check + build
X = zeros(height(T), numel(Xnames));
for k = 1:numel(Xnames)
    v = T.(Xnames(k));
    assert(isnumeric(v) || islogical(v), 'did:extractCovariates:Type', ...
        'Covariate %s must be numeric/logical.', Xnames(k));
    X(:,k) = double(v);
end
end
