g1 = gaussian(60, 10);
g2 = gaussian(90, 5);
g3 = gaussian(150, 20);
gmm = gaussian_mixture_model(1/4,g1,1/4,g2,1/2,g3);

N=500;
[training_data, training_data_class] = sample(gmm,N);

% training_data = rand(1,5)
% training_data_class = ones(1,5)

% training_data(1,1) = 40;
% training_data(1,2) = 56;
% training_data_class(1,1) = 1;
% training_data_class(1,2) = 2;

figure(1)
plot_mixture(training_data, training_data_class);

hyperpars.alpha = 1; % Concentration parameter; ignored if alpha is given Gamma prior
% Parameters of the NIG(mu, lambda, a, b) base distribution
hyperpars.mu = 75;
hyperpars.lambda = 1;
hyperpars.a = 1;
hyperpars.b = 10;
% Parameters for alpha's Gamma hyperprior
hyperpars.a_alpha = 1; % Shape
hyperpars.b_alpha = 0.01; % Inverse scale

%[class_id, phi, K_record, lP_record, alpha_record] = my_sampler(training_data, 150);
params_est = sampler_univ(training_data,200,hyperpars);
