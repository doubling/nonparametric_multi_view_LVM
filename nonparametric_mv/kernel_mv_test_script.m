
clear; close all; 
clc;

% generate synthetic dataset

m = 3000
k = 2;
nfold = 5;
nview = 3;

mu(1, :) = [0, -5, -2]; 
mu(2, :) = [1, 0.5, 1] * 10;
sigma{1} = diag([1.5, 1, 0.5]);
sigma{2} = diag([4, 3, 2]);
[data, true_mixture] = generate_gaussian_data_mv(m, 0.1, mu, sigma);

X{1} = data(:, 1);
X{2} = data(:, 2);
X{3} = data(:, 3);

%%

% kernel bandwith
kcoeff = [1];

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
    truelik{i}(:, 1) = normpdf(x_test', mu(1, i), sigma{1}(i,i))';
    truelik{i}(:, 2) = normpdf(x_test', mu(2, i), sigma{2}(i,i))';
    % truelik{i}(:, 1) = normpdf(x_test', center_1, sigma_1)';
    % truelik{i}(:, 2) = normpdf(x_test', center_2, sigma_2)';
    figure; hold on
   
    plot(x_test, testlik{i}(:, 1), 'r');
    plot(x_test, testlik{i}(:, 2));
    
    plot(x_test, truelik{i}(:, 1), 'c');
    plot(x_test, truelik{i}(:, 2), 'g');
    hold off;
end


