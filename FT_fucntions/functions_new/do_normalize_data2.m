function dataCell = do_normalize_data2(dataCell)


for i = 1:size(dataCell,1)
    dataCell(i,:) = (dataCell(i,:) - min(dataCell(i,:), [], 'all')) / (max(dataCell(i,:), [], 'all') - min(dataCell(i,:), [], 'all'));
end



% function X = do_normalize_data(X, method)
% % method: 'none' (default) | 'demean' | 'minmax'
% 
% if nargin<2, method = 'none'; end
% 
% switch lower(method)
%     case 'demean'
%         X = X - mean(X,2);
%     case 'minmax'
%         X = (X - min(X,[],2)) ./ (max(X,[],2) - min(X,[],2) + eps);
%     otherwise
%         % 'none' => return unchanged
% end
% end
