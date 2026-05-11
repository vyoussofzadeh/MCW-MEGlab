function [Yhat,p] = tta_predict(net, Xcells, cats, fs)
K=5; jitter = round(0.01*fs);   % ±10 ms
p = zeros(numel(Xcells),1);
for i=1:numel(Xcells)
  px = zeros(K,1);
  for k=1:K
    x = Xcells{i};
    j = randi([-jitter,jitter]);
    if j>0, x = [x(:,1+j:end) x(:,end*ones(1,j))];
    elseif j<0, x = [x(:,1*ones(1,-j)) x(:,1:end+j)];
    end
    [~,s] = classify(net, {x}); px(k) = s(:,strcmp(cats,'Spike'));
  end
  p(i) = mean(px);
end
Yhat = categorical((p>=0.5)+1,[1 2],cats);
end
