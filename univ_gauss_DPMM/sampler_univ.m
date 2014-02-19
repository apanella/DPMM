function [params] = sampler_univ(x, T, hyperpars)

% Gibbs sampler for mixtures of univariate Gaussians.
% This code implements a Rao-Blackwellized (collapsed) Gibbs sampler,
% following the algorithm on page 108 of Erik Sudderth's PhD thesis.
%
% Arguments:
% x is a vector of the data points to be clustered
% T is the number of Gibbs iterations to be performed
% hyperpars is a structure containing the parameters of the Bayesian model

%% Variable declaration and initialization
N = numel(x); % length of observation sequence
draw_plots = true; % draw plots at each iteration
var_alpha = true; % sample the concentration parameter
plot_prior = false; % plot prior distribution over parameters
plot_posterior = false; % plot posterior distribution over parameters after each iteration

% Density function of a generalized student-t distribution
gtpdf = @(x, v, mu, sigma) gamma((v+1)/2) / (gamma(v/2) * sqrt(pi*v) * sigma) * ...
    (1 + 1/v * ((x-mu)/sigma)^2)^(-(v+1)/2);

params(1).z = zeros(1,N); % cluster id, for each observation
params(1).c = 0; % counts (how many observations currently assigned to clusters)
params(1).K = 0; % cumber of clusters
params(1).sums = 0; % sum of data points associated to each clsuter
params(1).ss = 0; % sum of squares of data points
% Concentration parameter
if var_alpha
    params(1).alpha = randgamma(hyperpars.a_alpha) / hyperpars.b_alpha;
else
    params(1).alpha = hyperpars.alpha;
end

%% Plot prior distributions, if needed
if plot_prior
    figure(4);
    subplot(1,2,1);
    plot_gamma(hyperpars.a_alpha, hyperpars.b_alpha, 200);
    title('$p(\alpha)$', 'interpreter', 'Latex', 'fontsize', 14);
    subplot(1,2,2);
    plot_NIG(hyperpars.mu, hyperpars.lambda, hyperpars.a, hyperpars.b);
    title('$p(\mu,\sigma^2)$', 'interpreter', 'Latex', 'fontsize', 14);
end

%% Gibbs sampling
for t=2:T
    % Copy parameters
    params(t)  = params(t-1);
    % Compute random order to scan data points
    perm = randperm(N);
    for i=1:N
        idx = perm(i);
        % Update the counts array by removing x from its cluster
        % If the cluster_id is zero, then we are in the first iteration and
        % this is not needed.
        if params(t-1).z(idx) ~= 0
            params(t);
            params_if = rem_observation(x, params(t), idx);
        else
            params_if = params(t);
        end
        % Compute the predictive likelihood for each of the existing clusters
        likelihood = zeros(1,params(t).K+1);
        for k=1:params(t).K
            %params_k = add_observation(x, params_if, idx, k);
            hyperpars_k = get_NIG_hyperpars(hyperpars, params_if, k);
            v_t = 2*hyperpars_k.a;
            mu_t = hyperpars_k.mu;
            sigma_t = sqrt(hyperpars_k.b*(1+hyperpars_k.lambda^-1) / hyperpars_k.a);
            likelihood(k) = gtpdf(x(idx), v_t, mu_t, sigma_t);
        end
        % Compute the predictive likelihood for a new cluster
        %params_k = add_observation(x, params(1), idx, 1)
        %hyperpars_k = get_NIG_hyperpars(hyperpars, params_k, 1)
        hyperpars_k = hyperpars;
        %pause
        v_t = 2*hyperpars_k.a;
        mu_t = hyperpars_k.mu;
        sigma_t = sqrt(hyperpars_k.b*(1+hyperpars_k.lambda^-1) / hyperpars_k.a);
        likelihood(params(t).K+1) = gtpdf(x(idx), v_t, mu_t, sigma_t);
        % Compute the probability of each cluster
        % if K=0, we're picking a new cluster for sure
        if params(t).K == 0
            prob_cluster = 1;
        else
            prob_cluster = [params_if.c, params(t).alpha] .* likelihood;
        end
%         prob_cluster
%         pause
        cum_prob_cluster = cumsum(prob_cluster);
        % Sample the cluster for observation x(idx)
        new_z = 1 + sum(cum_prob_cluster(end)*rand(1) > cum_prob_cluster);
        params(t) = add_observation(x, params_if, idx, new_z);
        % Increment K if a new cluster is chosen
        if new_z == params(t).K+1
            params(t).K = params(t).K+1;
        end
        % Deal with possible empty clusters
        if ~isempty(find(params(t).c == 0,1))
            % Remove from parameter arrays and decrease cluster ids
            params(t) = remove_cluster(params(t), find(params(t).c == 0));
        end
    end
    % Resample alpha, if required
    if var_alpha
       params(t) = sample_alpha(params(t), hyperpars, N);
    end
    % If needed, plot the posterior parameters after each
    % iteration
    if plot_posterior
        figure(5);
        for k=1:params(t).K
            subplot(1,params(t).K,k);
            hyperpars_k = get_NIG_hyperpars(hyperpars, params(t), k);
            plot_NIG(hyperpars_k.mu, hyperpars_k.lambda, hyperpars_k.a, hyperpars_k.b);
            title('$p(\mu_k,\sigma_k^2)$', 'interpreter', 'Latex', 'fontsize', 14);
        end
    end
    
    if draw_plots
        figure(2)
        plot_mixture(x, params(t).z)
        
        figure(3)
        subplot(2,1,1)
        plot([params(1:t).K]);
        title('K');
        subplot(2,1,2)
        plot([params(1:t).alpha]);
        title('alpha');
        drawnow
    end
end

end

