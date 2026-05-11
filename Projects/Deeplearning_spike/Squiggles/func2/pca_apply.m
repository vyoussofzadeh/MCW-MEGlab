function Z = pca_apply(F, model)
% Apply TRAIN-fitted PCA model to new feature matrix F (N×D_raw)
% model.mu, model.sd, model.coeff, model.D come from pca_fit_apply
F(~isfinite(F)) = 0;
Fz = (F - model.mu) ./ model.sd;
Z  = Fz * model.coeff(:, 1:model.D);
end
