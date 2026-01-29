function v = getVcov(opts)
%GETVCOV Build a VCOV engine object from opts (NV-pairs).
arguments
    opts struct
end

spec = "none";
if isfield(opts,'vcov') && ~isempty(opts.vcov)
    spec = lower(string(opts.vcov));
end

% Clusters: let the engine default to ds.idVar when empty
clusters = string.empty(1,0);
if isfield(opts,'clusters') && ~isempty(opts.clusters)
    clusters = string(opts.clusters);
end

smallSample = true;
if isfield(opts,'smallSample') && ~isempty(opts.smallSample)
    smallSample = logical(opts.smallSample);
end

display = true;
if isfield(opts,'Display') && ~isempty(opts.Display)
    display = logical(opts.Display);
end

% Wild bootstrap options (only if you actually have a WildBootstrap engine)
cluster = string.empty(1,0);
if isfield(opts,'cluster') && ~isempty(opts.cluster)
    cluster = string(opts.cluster);
end

B = 999; 
if isfield(opts,'B') && ~isempty(opts.B)
    B = double(opts.B);
end

mult = "mammen";
if isfield(opts,'multiplier') && ~isempty(opts.multiplier)
    mult = lower(string(opts.multiplier));
end

stud = true;
if isfield(opts,'studentize') && ~isempty(opts.studentize)
    stud = logical(opts.studentize);
end

seed = randi([1,1e7],1,1);
if isfield(opts,'rngSeed') && ~isempty(opts.rngSeed)
    seed = double(opts.rngSeed);
end

switch spec
    case "clustered"
        v = did.vcov.Clustered(clusters=clusters, smallSample=smallSample, Display=display);

    case "wild"
        % Only keep this branch if did.vcov.WildBootstrap exists in your repo.
        v = did.vcov.WildBootstrap( ...
                cluster=cluster, B=B, multiplier=mult, ...
                studentize=stud, rngSeed=seed, smallSample=smallSample);

    otherwise
        v = [];   % no VCOV decoration
end
end

