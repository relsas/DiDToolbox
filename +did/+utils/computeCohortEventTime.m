function T = computeCohortEventTime(T, idVar, timeVar, dVar)
% Add `cohort` and `eventTime` once, consistently.
arguments
T table
idVar (1,1) string
timeVar (1,1) string
dVar (1,1) string
end


tvals = T.(timeVar);
if isdatetime(tvals)
tnum = days(tvals - min(tvals));
elseif isduration(tvals)
tnum = seconds(tvals - min(tvals));
else
tnum = double(tvals);
end


[g, ~] = findgroups(T.(idVar));
firstTreat = splitapply(@localFirstTreat, T.(dVar), tnum, g);
cohortByObs = firstTreat(g);


T.("cohort") = cohortByObs;
T.("eventTime") = tnum - cohortByObs; % can be negative pre, 0 at first treat, NaN if never treated
end


function ft = localFirstTreat(d, t)
mask = (d ~= 0) & ~isnan(d);
if any(mask)
ft = min(t(mask));
else
ft = NaN;
end
end