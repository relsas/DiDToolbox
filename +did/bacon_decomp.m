function out = bacon_decomp(T, args)
% BACON_DECOMP  Goodman-Bacon (2021) TWFE decomposition wrapper for the DiD Toolbox.
% 
% Usage (name–value only):
%   out = did.bacon_decomp(T, idvar="id", tvar="time", gvar="treatmentYear", yvar="y")
%
% Inputs (table T must be one row per unit-time):
%   idvar : unit id variable name in T (string)
%   tvar  : time variable name in T (string, numeric or datetime supported)
%   gvar  : treatment time/cohort (0 for never-treated; >0 for first treat time)
%   yvar  : outcome variable name in T (string)
%
% Output struct `out`:
%   .BaconTWFE : scalar, sum(weight .* estimate)
%   .Pairs     : table of all 2x2 components (type, group1, group2, weight, estimate)
%   .ByType    : table aggregated by component type (Type, Weight, Estimate)
%   .Call      : string, reproducible call signature
%   .Vars      : struct of variable names used
%
% Notes
% - Thin wrapper around the OOP class `baconDecomp` 
% 

arguments
    T table
    args.idvar (1,1) string = "id"
    args.tvar  (1,1) string = "time"
    args.gvar  (1,1) string = "treatVar"
    args.yvar  (1,1) string = "y"
    args.isPrint logical = true
end

% Instantiate the class (constructor also computes and stores results)
mdl = did.estimators.baconDecomp(T, timeVar=args.tvar, idVar=args.idvar, treatVar=args.gvar, outcomeVar=args.yvar);

% Collect outputs in toolbox-friendly struct
Pairs  = mdl.Results;
ByType = mdl.ResultsAgg;

% Overall (TWFE) is the weighted average over all 2x2 blocks
BaconTWFE = sum(Pairs.weight .* Pairs.estimate, 'omitnan');

out = struct();
out.BaconTWFE = BaconTWFE;
out.Pairs     = Pairs;
out.ByType    = ByType;
out.Call      = sprintf('did.bacon_decomp(T, idvar="%s", tvar="%s", gvar="%s", yvar="%s")', ...
                        args.idvar, args.tvar, args.gvar, args.yvar);
out.Vars      = struct('idvar',args.idvar, 'tvar',args.tvar, 'gvar',args.gvar, 'yvar',args.yvar);

% print
if args.isPrint
    fprintf('\n');
    disp("---- Bacon: Overall (TWFE) -----");
    disp(out.BaconTWFE)

    disp("Bacon: ---- Comparisons -----");
    out.Pairs.type = string(out.Pairs.type);
    disp(out.Pairs);

    disp("Bacon: ---- Comparisons agg. by Type ----")
    out.ByType.Type = string(out.ByType.Type); 
    disp(out.ByType)
end

end


