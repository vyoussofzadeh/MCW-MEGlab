clc
X = a';
[R,~] = corrcoef(X','rows','pairwise');
imagesc(R); axis square; colorbar; title('Pearson-r across parcels')

% ---------- p-values with effective df ----------
rho1  = diag(corrcoef(X(:,1:end-1)',X(:,2:end)'));     % lag-1 AC
T     = size(X,2);
Neff  = max(   T.*(1-rho1)./(1+rho1) , 3);             % Bartlett
P     = 0*R;
for i = 1:size(X,1)-1
    for j = i+1:size(X,1)
        df     = floor(min(Neff([i j])) - 2);
        tstat  = R(i,j)*sqrt(df/(1-R(i,j)^2));
        P(i,j) = 2*tcdf(-abs(tstat),df);
        P(j,i) = P(i,j);
    end
end



C  = cov(X);                          % p × p
C  = C + 0.05*trace(C)/size(C,1)*eye(size(C));   % Tikhonov
Prec = pinv(C);
pcorr = -Prec ./ sqrt(diag(Prec)*diag(Prec)');   % partial-r

imagesc(abs(pcorr) > 0.05); axis square
title('Conditional-dependency map (|partial-r| > 0.05)')
colorbar


%%
% --- compute Gram (p × p) ---
G = X*X';                    

% --- eigen / SVD for stability ---
[V,D]   = eig((G+G')/2);     % ensure symmetry
D_inv_sqrt = diag(1./sqrt(max(diag(D), 1e-12)));  % protect tiny eigs

% --- whitening transform ---
W = D_inv_sqrt * V';          % p × p

X_sym  = W * X;               % orthogonalised matrix  (p × T)

% sanity check
nearEye = (X_sym'*X_sym)/T;   % should be ~identity


C  = cov(X_sym);                          % p × p
C  = C + 0.05*trace(C)/size(C,1)*eye(size(C));   % Tikhonov
Prec = pinv(C);
pcorr = -Prec ./ sqrt(diag(Prec)*diag(Prec)');   % partial-r

figure, imagesc(abs(pcorr) > 0.05); axis square
title('Conditional-dependency map (|partial-r| > 0.05)')
colorbar