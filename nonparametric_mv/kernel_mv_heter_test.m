
clear; close all; 
clc

% generate synthetic data

m = 1000
k = 2;
nfold = 5;
nview = 3;

[mu, sigma, true_mixture, data] = generate_mix_heter(m, k, nview, 0);


X{1} = data(:, 1);
X{2} = data(:, 2);
X{3} = data(:, 3);

%%
kcoeff = [0.1];

% construct kernels we need

options.KernelType = 'Gaussian';

Kcell = cell(1, nview);
for i = 1: nview	
    D = pdist(X{i});
    median_distance = median(D);
    options.t = median_distance * kcoeff;
	Kcell{i} = constructKernel(X{i}, [], options);
end

[cond_opt, prior] = kernel_mv(Kcell, k);

m_test = 200;
x_test = linspace(-10, 15, m_test)';

center_1 = -5; sigma_1 = 1.5;
center_2 = 5; sigma_2 = 3;


for i = 1 :nview
    K_test{i} = constructKernel(x_test, X{i}, options);
    testlik{i} = K_test{i} * cond_opt{i};

    truelik{i}(:, 1) = normpdf(x_test', mu(1, i), sigma(1,i))';
    truelik{i}(:, 2) = gampdf(x_test', mu(2, i), sigma(2,i))';
    figure; hold on
   
    plot(x_test, testlik{i}(:, 1), 'r');
    plot(x_test, testlik{i}(:, 2));
    
    plot(x_test, truelik{i}(:, 1), 'c');
    plot(x_test, truelik{i}(:, 2), 'g');
    hold off;
end


