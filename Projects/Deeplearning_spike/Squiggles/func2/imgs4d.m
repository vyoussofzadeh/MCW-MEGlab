function [Im4d, Y1] = imgs4d(Xcells, Y)
N = numel(Xcells); C = size(Xcells{1},1); T = size(Xcells{1},2);
Im4d = zeros(C, T, 1, N, 'single');
for i=1:N, Im4d(:,:,1,i) = single(Xcells{i}); end
Y1 = Y;
end
