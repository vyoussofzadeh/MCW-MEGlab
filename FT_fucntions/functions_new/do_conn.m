function [outsum] = do_conn(mom)

[~, nrpt] = size(mom);
crsspctrm = (mom*mom')./nrpt;
tmp = crsspctrm; crsspctrm = []; crsspctrm(1,:,:) = tmp;

input = crsspctrm;
pownorm = 1;

siz = [size(input) 1];
% crossterms are described by chan_chan_therest
outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));
% outcnt = zeros(siz(2:end));
for j = 1:siz(1)
    if pownorm
        p1  = zeros([siz(2) 1 siz(4:end)]);
        p2  = zeros([1 siz(3) siz(4:end)]);
        for k = 1:siz(2)
            p1(k,1,:,:,:,:) = input(j,k,k,:,:,:,:);
            p2(1,k,:,:,:,:) = input(j,k,k,:,:,:,:);
        end
        p1    = p1(:,ones(1,siz(3)),:,:,:,:);
        p2    = p2(ones(1,siz(2)),:,:,:,:,:);
        denom = sqrt(p1.*p2); clear p1 p2;
    end
    tmp    = abs(reshape(input(j,:,:,:,:,:,:), siz(2:end))./denom); % added this for nan support marvin
    %tmp(isnan(tmp)) = 0; % added for nan support
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
%     outcnt = outcnt + double(~isnan(tmp));
end

% size(outsum)
% figure,imagesc(outsum), colorbar, title('conn (across voxels)');

end