function [params, old_params] = add_observation(x, params, i, k)

old_params = params;
if k <= params.K
    params.c(k) = params.c(k) + 1;
    params.sums(k) = params.sums(k) + x(i);
    params.ss(k) = params.ss(k) + x(i)^2;
else
    params.c(k) = 1;
    params.sums(k) = x(i);
    params.ss(k) = x(i)^2;
end
params.z(i) = k;

end