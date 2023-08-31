function unmixing = do_pca(sel_val, n)

numcomponent = n;
Nchans = size(sel_val,2);
dat = sel_val';
C = (dat*dat')./(size(dat,2)-1);

% eigenvalue decomposition (EVD)
[E,D] = eig(C);

% sort eigenvectors in descending order of eigenvalues
d = cat(2,(1:1:Nchans)',diag(D));
d = sortrows(d, -2);

% return the desired number of principal components
unmixing = E(:,d(1:numcomponent,1))';

end