classdef (Abstract) Estimator < handle
% Base interface for all estimators.


properties
Name (1,1) string = ""
Options struct = struct()
% VcovEngine (1,1) did.vcov.VcovEngine = did.vcov.NoVcov()

VcovEngine = []   % untyped, set later by fit.m
end


methods (Abstract)
res = fit(obj, ds)
end


methods
function obj = setOptions(obj, opts), obj.Options = opts; end %#ok<INUSD>
function tf = supportsEventStudy(obj), tf = false; end
end
end