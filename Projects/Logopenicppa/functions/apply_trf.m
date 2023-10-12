function [new] = apply_trf(transform, old)

[m, n] = size(old);
if m~=3 && n==3
  % each row is one position
  old(:,4) = 1;
  new = old * transform';
  new = new(:,1:3);
elseif m==3 && n~=3
  % each column is one position
  old(4,:) = 1;
  new = transform * old;
  new = new(1:3,:);
else
  % assume that each row is one position
  old(:,4) = 1;
  new = old * transform';
  new = new(:,1:3);
end