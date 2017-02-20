function [mu, sigma, true_mixture, data] = generate_mix_heter(num, k, nview, symmetric)

% param num: the number samples to be generated
% param k: the number of components, k = 2.^l
% param nview: the number of view
% param symmetric: the flag for symmetric distribution or not

if symmetric
%     % k = 2 setting 
%     mu(1, :) = [-2, -2, -2]; % + repmat(rand() * 0.2, 1, 3);
%     
%     %alpha
%     mu(2, :) = [3, 3, 3] ;% + repmat(rand() * 0.2, 1, 3);
%     
%     sigma(1, :) = [1, 1, 1]; % + repmat(rand() * 0.2, 1, 3);
%     
%     %beta
%     sigma(2, :) = [2, 2, 2]; % + repmat(rand() * 0.2, 1, 3);
    
    
    
    mu(1, :) = [-6, -6, -6]; % + repmat(rand() * 0.2, 1, 3);
    
    %alpha
    mu(2, :) = [0.75, 0.75, 0.75] ;% + repmat(rand() * 0.2, 1, 3);
    
    sigma(1, :) = [2, 2, 2]; % + repmat(rand() * 0.2, 1, 3);
    
    %beta
    sigma(2, :) = [3, 3, 3]; % + repmat(rand() * 0.2, 1, 3);

    
else

    % this parameters for k = 2
%     mu(1, :) = [-3, -5, -2]; 
% 
%     %alpha
%     mu(2, :) = [2, 0.5, 3];
% 
%     sigma(1, :) = [1.5, 1, 0.5];
% 
%     %beta
%     sigma(2, :) = [8, 3, 2];

    mu(1, :) = [-3, -5, -2]; 

    %alpha
    % mu(2, :) = [0.75, 0.5, 1] ;
    mu(2, :) = [2, 3, 1] ;
    
    sigma(1, :) = [1.5, 1, 0.5];

    %beta
    sigma(2, :) = [8, 3, 2];

end

gap = mu(2, :) - mu(1, :);

% if mod(k, 2)==0
%     sigma =  repmat(sigma, k/2, 1);
% else
%     sigma = [sigma; sigma(2, :)/2];
% end
%
% for i = 3 : k
%     mu = [mu; mu(i-1, :) + gap];
% end
% 
% if ~symmetric
%     mu = mu + rand(size(mu))*0.2;
%     sigma = sigma + rand(size(sigma))*0.2;
% end

%the proportion for each component in ascending
prop = (1:k)/sum(1:k);

% true_mixture = mnrnd(num, prop);
true_mixture = ceil(prop * num);
true_mixture(end) = true_mixture(end) + (num - sum(true_mixture));

indx = zeros(1, k+1);
for i = 2: k+1
    indx(i) = indx(i-1) + true_mixture(i-1);
end

data = randn(num, nview);

for i = 2: k+1
    
    if mod(i, 2) == 0
        data((indx(i-1)+1):indx(i) ,:) = data((indx(i-1)+1):indx(i) ,:) * diag(sigma(1, :)) ...
            + repmat(mu(1,:), indx(i)-indx(i-1),1) + repmat(gap * (i-2) / 2, indx(i)-indx(i-1), 1);
        
        if i > 3
            mu = [mu; mu(1, :) + gap * (i-2) / 2];
            sigma = [sigma; sigma(1, :)];
        end
        
    else
        for p = 1:3
            data((indx(i-1)+1):indx(i), p) = gamrnd(mu(2, p), sigma(2, p), true_mixture(i-1), 1) + gap(p) * (i-3) / 2;
        end
        
        if i >3
            mu = [mu; mu(2, :)];
            sigma = [sigma; sigma(2, :)];
        end
    end
end

if k == 3
    for p = 1:3
        data((indx(k)+1):indx(k+1), p) = gamrnd(mu(2, p), sigma(2, p), true_mixture(i-1), 1) + gap(p);
    end
    
    mu(end, :) = mu(2, :);
    sigma(end, :) = sigma(2, :);
end

true_mixture = true_mixture ./ num;
