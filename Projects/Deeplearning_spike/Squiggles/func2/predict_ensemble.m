function [Yhat, pFuse] = predict_ensemble(net1d, net2d, Xcells, Im4d, alpha, cats)
[~, s1] = classify(net1d, Xcells); p1 = s1(:, strcmp(cats,'Spike'));
[~, s2] = classify(net2d, Im4d);   p2 = s2(:, strcmp(cats,'Spike'));
pFuse   = alpha*p1 + (1-alpha)*p2;
Yhat    = categorical( (pFuse>=0.5)+1, [1 2], cats);
end
