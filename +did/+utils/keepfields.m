function s = keepfields(s, names)
%KEEPFIELDS Keep only named fields in a scalar struct
all = string(fieldnames(s));
drop = setdiff(all, names);
if ~isempty(drop)
    s = rmfield(s, cellstr(drop));
end
end
