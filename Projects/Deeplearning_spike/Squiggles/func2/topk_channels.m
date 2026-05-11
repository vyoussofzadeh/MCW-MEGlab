function [Xk, idx] = topk_channels(X, k)
% Keep top-k channels by energy (fit from X).
C = size(X{1},1);
eng = zeros(C,1,'double');
for i=1:numel(X), eng = eng + sum(double(X{i}).^2, 2); end
[~,idx] = maxk(eng, k);
Xk = cellfun(@(z) z(idx,:), X, 'UniformOutput', false);
end
