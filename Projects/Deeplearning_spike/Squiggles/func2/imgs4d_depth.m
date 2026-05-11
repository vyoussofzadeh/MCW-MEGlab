function [Im4d, Y1] = imgs4d_depth(Xcells, Y)
N = numel(Xcells);
[C,T,D] = size(Xcells{1});
Im4d = zeros(C, T, D, N, 'single');
for i=1:N
    Xi = Xcells{i};
    assert(isequal(size(Xi), [C T D]), 'imgs4d_depth: size mismatch at %d', i);
    Im4d(:,:,:,i) = single(Xi);
end
Y1 = Y;
end
