function k = eventTimeFrom(T, idVar, t_int, g_cohort)
% k = t - g(id), with k=NaN for never-treated (g==0)
k = nan(height(T),1);
% group rows by id once
[gId, uids] = findgroups(T.(idVar));
% first row per id gets the cohort value; broadcast it
gPerId = splitapply(@(gi) gi(1), g_cohort, gId);  % one g per id row
% Bring back to rows order
gThis = gPerId(gId);
m = (gThis>0) & isfinite(t_int);
k(m) = t_int(m) - gThis(m);
end
