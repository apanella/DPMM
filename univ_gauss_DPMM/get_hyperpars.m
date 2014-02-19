function hyperpars_new = get_NIG_hyperpars(hyperpars, params, k)

hyperpars_new.lambda = hyperpars.lambda + params.c(k);
hyperpars_new.mu = (hyperpars.lambda*hyperpars.mu + params.sums(k)) / ...
    (hyperpars_new.lambda);
hyperpars_new.a = hyperpars.a + params.c(k)/2;
hyperpars_new.b = hyperpars.b + 1/2*(hyperpars.mu^2*hyperpars.lambda + ...
    params.ss - hyperpars_new.mu^2*hyperpars_new.lambda);

end