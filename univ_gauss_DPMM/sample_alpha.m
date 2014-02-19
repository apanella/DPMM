function [params, old_params] = sample_alpha(params, hyperpars, N)

% Sample the concentration parameter alpha, which is given a gamma prior,
% using the method explained in Mike West's 1995 white paper
% "Hyperparameter estimation in Dirichlet process mixture models"

old_params = params;
% Sample auxiliary variable x
x_alpha = randbeta(params.alpha+1, N);
% Sample alpha
prop = (hyperpars.a_alpha+params.K-1) / ...
    (hyperpars.a_alpha+params.K-1 + N * (hyperpars.b_alpha-log(x_alpha)));
if rand < prop
    params.alpha = randgamma(hyperpars.a_alpha+params.K) / ...
        (hyperpars.b_alpha-log(x_alpha));
else
    params.alpha = randgamma(hyperpars.a_alpha+params.K-1) / ...
        (hyperpars.b_alpha-log(x_alpha));
end

end