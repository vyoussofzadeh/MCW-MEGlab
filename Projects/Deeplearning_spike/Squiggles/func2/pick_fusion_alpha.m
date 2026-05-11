function [alpha, BA, cm] = pick_fusion_alpha(net1d, net2d, Xva2, ImVa, YVal, cats)
% Grid-search alpha in [0..1] to maximize Balanced Accuracy on VAL
alphas = linspace(0,1,51);

% Spike prob from each model
[~, s1] = classify(net1d, Xva2);  p1 = s1(:, strcmp(cats,'Spike'));
[~, s2] = classify(net2d, ImVa);  p2 = s2(:, strcmp(cats,'Spike'));

bestBA=-Inf; bestA=0.5; bestCM=[];
for a = alphas
    p  = a*p1 + (1-a)*p2;                          % fused prob(spike)
    Yp = categorical( (p >= 0.5) + 1, [1 2], cats);% threshold @0.5
    C  = confusionmat(YVal, Yp, 'Order', cats);
    rec = diag(C) ./ max(1,sum(C,2)); BAcur = mean(rec);
    if BAcur > bestBA, bestBA=BAcur; bestA=a; bestCM=C; end
end
alpha = bestA; BA = bestBA; cm = bestCM;
fprintf('Fusion alpha=%.2f  BA=%.3f\n', alpha, BA);
end
