function Xo = instance_robust_z(Xi)
Xo = cell(size(Xi));
for i=1:numel(Xi)
  x = Xi{i}; med = median(x,2); madv = 1.4826*median(abs(x-med),2)+1e-12;
  Xo{i} = (x - med) ./ madv;
end
end
