function [Ztr, model] = pca_fit_apply(Ftr, varargin)
% [Ztr, model] = pca_fit_apply(Ftr, 'var',95)
% model has fields: mu, coeff, D, expl

p = inputParser;
addParameter(p,'var',95,@(x) isscalar(x) && x>0 && x<=100);
parse(p,varargin{:});
varTarget = p.Results.var;

% standardize (fit on TRAIN)
mu = mean(Ftr,1,'omitnan');
sd = std(Ftr,0,1,'omitnan'); sd(sd<1e-8)=1;
Fz = (Ftr - mu) ./ sd;

[coeff,score,~,~,expl] = pca(Fz,'Centered',false);  % already centered/scaled
D = find(cumsum(expl)>=varTarget, 1);
if isempty(D), D = size(score,2); end

Ztr = score(:,1:D);  % reduced TRAIN
model = struct('mu',mu,'sd',sd,'coeff',coeff,'D',D,'expl',expl);
end

function Z = pca_apply(F, model)
% Z = pca_apply(F, model)  ? apply TRAIN-fitted PCA to new data
F(~isfinite(F)) = 0;
Fz = (F - model.mu) ./ model.sd;
Z  = Fz * model.coeff(:,1:model.D);
end
