function [params, old_params] = rem_observation(x, params, i)

old_params = params;
class = params.z(i);
params.c(class) = params.c(class) - 1;
params.sums(class) = params.sums(class) - x(i);
params.ss(class) = params.ss(class) - x(i)^2;

end