function nv = struct2nv(s)
%STRUCT2NV Convert scalar struct to name–value cell row
if ~isscalar(s), error('struct2nv:ScalarOnly','opts must be scalar'); end
f = fieldnames(s);
v = struct2cell(s);
nv = reshape([f.'; v.'], 1, []);
end
