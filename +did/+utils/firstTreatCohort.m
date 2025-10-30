function g = firstTreatCohort(T, idVar, timeVar, dVar)
% First treatment time per unit mapped to integer time; 0 for never
t_int = did.utils.timeInt(T, timeVar);
D = T.(dVar)==1;

% group ids
[gId, ids] = findgroups(T.(idVar));
firstOn = splitapply(@(t,d) localFirstOn(t,d), t_int, D, gId);  % NaN if never

% Expand back to row level by joining on id
adopt = table(ids, firstOn, 'VariableNames', [idVar, "g_int"]);
T2 = outerjoin(T(:, idVar), adopt, "Keys", idVar, "MergeKeys", true, "Type","left");

g = T2.g_int;                 % NaN for never
g(~isfinite(g)) = 0;          % 0 denotes never-treated (canonical)

    function a = localFirstOn(t, d)
        [ts, ix] = sort(t);
        ds = d(ix);
        j = find(ds==1, 1, 'first');
        if isempty(j), a = NaN; else, a = ts(j); end
    end
end
