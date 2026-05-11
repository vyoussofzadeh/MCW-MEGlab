function [Im4d, Y1] = imgs4dD(Xcells, Y)
N = numel(Xcells); C = size(Xcells{1},1); T = size(Xcells{1},2); D = size(Xcells{1},3);
Im4d = zeros(C, T, D, N, 'single');
for i=1:N
    Xi = Xcells{i};
    assert(isequal(size(Xi), [C T D]), 'imgs4dD: epoch %d size mismatch', i);
    Im4d(:,:,:,i) = single(Xi);
end
Y1 = Y;
end
