function [hyperpars, old_hyperpars] = get_NIG_hyperpars(hyperpars, params, k)

% Compute the posterior parameters of the Normal-inverse-Gamma distribution
% for class k, given the prior parameters and sufficient statistics.

old_hyperpars = hyperpars;
hyperpars.lambda = old_hyperpars.lambda + params.c(k);
hyperpars.mu = (old_hyperpars.lambda*old_hyperpars.mu + params.sums(k)) / ...
    (hyperpars.lambda);
hyperpars.a = old_hyperpars.a + params.c(k)/2;
hyperpars.b = old_hyperpars.b + 1/2*(old_hyperpars.mu^2*old_hyperpars.lambda + ...
    params.ss(k) - hyperpars.mu^2*hyperpars.lambda);

end