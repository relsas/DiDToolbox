function validatePanel(T, idVar, timeVar, yVar, dVar, wVar)
% Validate required variables & types. No renaming, no mutation.
arguments
    T table
    idVar   (1,1) string
    timeVar (1,1) string
    yVar    (1,1) string
    dVar    (1,1) string
    wVar string = string.empty
end

names = string(T.Properties.VariableNames);

% Presence (exact match)
req  = [idVar, timeVar, yVar, dVar];
miss = req(~ismember(req, names));
if ~isempty(miss)
    error("did:validatePanel:MissingVars", ...
        "Missing variables: %s. Available: %s", ...
        strjoin(miss, ", "), strjoin(names, ", "));
end

% Types
if ~isnumeric(T.(yVar)) && ~islogical(T.(yVar))
    error("did:validatePanel:Type", "%s must be numeric or logical.", yVar);
end
if ~isnumeric(T.(dVar)) && ~islogical(T.(dVar))
    error("did:validatePanel:Type", "%s must be numeric or logical.", dVar);
end

% (Optional) weight presence/type if provided
if strlength(wVar)>0
    if ~ismember(wVar, names)
        error("did:validatePanel:MissingWeight", "Weight variable %s not found.", wVar);
    end
    if ~isnumeric(T.(wVar)) && ~islogical(T.(wVar))
        error("did:validatePanel:Type", "Weight %s must be numeric/logical.", wVar);
    end
end
end
