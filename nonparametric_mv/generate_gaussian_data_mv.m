function [data, true_mixture] = generate_gaussian_data_mv(m, prop, mu, sigma)
% generate samples.

latent_assignment = rand(1, m) > prop;
tmp = sum(latent_assignment);
true_mixture = [tmp / m, 1 - tmp / m];

data = zeros(m, 3);
data_ = randn(3, m);
data(latent_assignment, :) = data_(:, latent_assignment)' * sigma{1} + repmat(mu(1, :), sum(latent_assignment), 1);
data(~latent_assignment, :) = data_(:, ~latent_assignment)' * sigma{2} + repmat(mu(2, :), sum(~latent_assignment), 1);
