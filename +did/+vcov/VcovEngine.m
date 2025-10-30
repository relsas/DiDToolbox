classdef (Abstract) VcovEngine < handle
    % Base interface for VCOV engines that decorate estimator results.
    methods (Abstract)
        res = decorate(obj, res, ds)
    end
end