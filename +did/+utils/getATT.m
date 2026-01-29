function att = getATT(res)
% GETATT  Extract main "overall" effect from a DiD result struct.
% Returns NaN if not found.

att = NaN;

% By estiamtor - get overall effect
if contains(res.Method,["TWFE","BJS"]) 
    att = res.ATT; return;
end

if contains(res.Method,["Wooldridge","IW"]) 
    att = res.overall.Estimate; return;
end

if contains(res.Method,["CH"]) 
    att = res.OverallCW; return;
end

if contains(res.Method,["CS"]) 
    att = res.overall.Aggregates.Estimate; return;
end



% Fallback via did_standardize: average post-event ATT
try
    S = did.did_standardize(res);
    if isfield(S,'overall') && isstruct(S.overall) && isfield(S.overall,'Estimate')
        att = S.overall.Estimate; if isfinite(att), return; end
    end
    if isfield(S,'es') && istable(S.es) && ...
            ismember('Estimate', S.es.Properties.VariableNames) && ...
            ismember('e', S.es.Properties.VariableNames)
        % Example: simple mean of post-event ATT (e >= 0)
        mask = S.es.e >= 0 & isfinite(S.es.Estimate);
        if any(mask)
            att = mean(S.es.Estimate(mask), 'omitnan');
            if isfinite(att), return; end
        end
    end
catch
    % swallow and leave att = NaN
end

end
