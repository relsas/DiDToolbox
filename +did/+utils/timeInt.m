function t_int = timeInt(T, timeVar)
% Stable integer index for time (1..nT), preserving original order
[~,~,t_idx] = unique(T.(timeVar), 'stable');
t_int = double(t_idx);
end
