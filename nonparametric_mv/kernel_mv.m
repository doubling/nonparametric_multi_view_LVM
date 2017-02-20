
function [cond_opt, prior] = kernel_mv(Kcell, k)

% param Kcell: it is a cell, each component of Kcell is a kernel matrix for 
%          one view observation. Assume here we deal with three-view problem
% param k: the number of status of latent variable
% 
%
%   Written by Bo Dai (bohr.dai@gmail.com)

% assume we are dealing with three-view problem
epsilon = 1e-8;
n = size(Kcell{1}, 1);
nview = length(Kcell);

cond_opt = cell(1, nview);

for i = 1: nview
    triR{i} = chol(Kcell{i} +  epsilon * eye(size(Kcell{i})));
end


for i = 1: nview
    
    % estimate the third moment embedding
	id = 1;
	for j = 1: nview
		if j ~= i
			tri_other{id} = triR{j};
            other_kernel{id} = Kcell{j};
			id = id + 1;
		end
    end
	
    
    [betaK, S] = kernel_symmtrization(tri_other{1}, other_kernel{2}, k);
    partK = other_kernel{1} * betaK;
    
    [betaL, S] = kernel_symmtrization(tri_other{2}, other_kernel{1}, k);
    partL = other_kernel{2} * betaL;
    k_LK = partL' * partK;
    
    % whitening
	[beta, S] = kernel_whiten(Kcell{i} + epsilon * eye(size(Kcell{i})), partK, partL, k_LK, k);
    
	%construct whitened third order embedding
	
    H = partK / k_LK * partL';
    W = diag(1./sqrt(diag(S))) * beta';
	obs_v3 = W * Kcell{i};
	obs_v1 = obs_v3 * H';
	obs_v2 = obs_v3 * H;
	
	lambda = ones(n, 1);
	third_moment = double(full(ktensor(lambda, obs_v1, obs_v2, obs_v3)));
	third_moment = tensor(third_moment + permute(third_moment, [3, 1, 2]) + permute(third_moment, [2, 3, 1]));
    third_moment = 1/(3*n) * third_moment;
    
    % 
	% tensor decomposition
    iter = 50;
    trial = 10;
 	[U, tS] = tpower(third_moment, k, iter, trial);
	cond_opt{i} = beta * sqrt(S) * U;
    cond_opt{i} = abs(cond_opt{i});
    cond_opt{i} = bsxfun(@rdivide, cond_opt{i}, sum(cond_opt{i}));
    tS = 1 ./ tS.^2;
	prior(:, i) = tS ./ sum(tS);
	
end

