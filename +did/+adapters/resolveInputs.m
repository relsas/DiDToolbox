function [T,idVar,timeVar,yVar,dVar,wVar,canon] = resolveInputs(firstArg, nv)
% Resolve either a did.Dataset or (T + names) into a unified set of inputs.
% nv: struct with possible fields idVar, timeVar, yVar, dVar, weightVar

if isa(firstArg, 'did.Dataset')
    ds     = firstArg;
    T      = ds.T;
    idVar  = ds.idVar;   timeVar = ds.timeVar;
    yVar   = ds.yVar;    dVar    = ds.dVar;
    wVar   = ds.weightVar;

    % canonical helpers (lazy via Dataset if available)
    canon.t_int     = ds.get("t_int");
    canon.g         = ds.get("g");
    canon.eventTime = ds.get("eventTime");

else
    % Legacy path: firstArg must be a table
    T      = firstArg;
    idVar  = nv.idVar;   timeVar = nv.timeVar;
    yVar   = nv.yVar;    dVar    = nv.dVar;
    wVar   = "";
    if isfield(nv,'weightVar'), wVar = nv.weightVar; end

    % validate basics with your existing function
    did.utils.validatePanel(T, idVar, timeVar, yVar, dVar, wVar);

    % compute canonical helpers on the fly
    canon.t_int     = did.utils.timeInt(T, timeVar);
    canon.g         = did.utils.firstTreatCohort(T, idVar, timeVar, dVar);
    canon.eventTime = did.utils.eventTimeFrom(T, idVar, canon.t_int, canon.g);
end
end
