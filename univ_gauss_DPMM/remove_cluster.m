function [params, old_params] = remove_cluster(params, k)

old_params = params;
params.c(k) = [];
params.sums(k) = [];
params.ss(k) = [];
params.z(find(params.z>k)) = params.z(find(params.z>k))-1;
params.K = params.K-1;

end