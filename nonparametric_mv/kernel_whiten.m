
function [beta, S] = kernel_whiten(G, K_nk, L_nk, kinv_LK, k)

% param G: kernel matrix for corresponding view
% param K_nk, L_nk: compute H = K_nk inv(L_nk^t K_nk) L_nk^t. 
% param k: number of status of latent variable


[n,m] = size(G);
R = chol(G);

% M = R * H * R';
tmp = K_nk / kinv_LK;
mid = tmp' * G * tmp;
leftmp = R * L_nk;
M = leftmp * mid * leftmp';

[beta, S] = eigs(M, k);
[sval, sindx] = sort(diag(S), 'descend');
S = diag(sval);
beta = beta(:, sindx);

beta = R\beta;
S = sqrt(S) / n;

