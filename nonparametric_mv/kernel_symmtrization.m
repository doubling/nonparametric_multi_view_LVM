
function [beta, S] = kernel_symmtrization(R, L, k)

% param K/L: the kernel for view_1 and view_2
% param k: the number of status of latent variables
% beta: n-by-k matrix

[n,m] = size(R);
% R = chol(K);

M = R * L * R';

% here different: 
% Le use function eig with sort
[beta, S] = eigs(M, k);
[sval, sindx] = sort(diag(S), 'descend');
S = diag(sval);
beta = beta(:, sindx);

beta = R\beta;
% S = S / (4 * n^2);
S = sqrt(S) / (2 * n);
