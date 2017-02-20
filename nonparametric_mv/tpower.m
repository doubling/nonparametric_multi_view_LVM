
function [U, S] = tpower(tmoment, k, iter, trial)

% param tmoment: empirical estimation of high order moment, here it
%                    should be a tensor
% param k: the number of status of latent variable, a.k.a, the number of
%          top eigenvectors we compute
% param iter: the iteration maximum number
% param trial: the number of trials
%
% this funciton implements the tensor power method for eigenvector and 
% corresponing eigenvalues 
% 
%
%   Written by Bo Dai (bohr.dai@gmail.com)

l = size(tmoment);

if l(1) ~= l(2)
    printf('the size of tensor is not match');
    return
elseif l(2) ~= l(3)
    printf('the size of tensor is not match');
    return
elseif l(1)<k
    printf('k should be less than the length of each mode');
    return
end

U = zeros(l(1), k);
S = zeros(k, 1);

for i = 1: k
    
    % restore the results for each trial
    restraj = zeros(l(1), trial);
    straj = zeros(1, trial);
    
    % several trials
    for j = 1: trial
        
        % random generate the initial vectors which is on unit sphere
        r = randn(l(1),1); 
        r = r ./ sqrt(r'*r);     
        
        % fixed-point ite rarion 
        for m = 1: iter
            r = double(ttv(tmoment, {r,r}, [2,3]));
            r = r ./ sqrt(r'*r);
        end
        
        restraj(:, j) = r;
        straj(j) = ttv(tmoment, {r,r,r}, [1,2,3]);
    end
    
    % find the maximum in straj
    [ms , indx] = max(straj);
    S(i) = max(ms, 0);
    mr = restraj(:, indx);
    U(:, i) = mr;
    
    
    % modify tmoment
    tmoment = tmoment - full(ktensor(ms, mr, mr, mr));
    
end
