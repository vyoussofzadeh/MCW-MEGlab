function nAB = vy_norm_diff(A,B)

nA = A/max(A);
nB = B/max(B);
nAB = (nA  - nB) ./ (nA + nB);
nAB(nAB<0)=0;