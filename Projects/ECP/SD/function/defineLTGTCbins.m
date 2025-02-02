function LTGTCcat = defineLTGTCbins(LTGTCraw)
% DEFINELTGTCBINS_BIGGER  Merges narrower 'c_...' categories into bigger bins.
% 
% For example:
%   c_0                  -> '0'
%   c_1, c_2-3, c_4-5    -> '1-5'
%   c_6-10, c_11-20      -> '6-20'
%   c_21-50, c_51-100, c_>100 -> '21plus'
%   MISSING              -> 'Missing'
%
% Any unexpected category is labeled 'Unknown'.

    % 1) Ensure we have a categorical array
    if iscell(LTGTCraw) || isstring(LTGTCraw)
        LTGTCraw = categorical(LTGTCraw);
    end
    if ~iscategorical(LTGTCraw)
        error('LTGTC must be a categorical array with c_0, c_1, c_2-3, etc.');
    end

    % 2) Prepare a new string array of the same size
    newLabels = repmat("Missing", size(LTGTCraw));  % default to "Missing"

    % 3) Assign merged bins for each row
    for i = 1:numel(LTGTCraw)
        oldCat = string(LTGTCraw(i));  % e.g. "c_1", "c_2-3"

        switch oldCat
            case "c_0"
                newLabels(i) = "0";

            case {"c_1","c_2-3","c_4-5"}
                newLabels(i) = "1-5";

            case {"c_6-10","c_11-20"}
                newLabels(i) = "6-20";

            case {"c_21-50","c_51-100","c_>100"}
                newLabels(i) = "21plus";

            case "MISSING"
                newLabels(i) = "Missing";

            otherwise
                % If there's an unexpected category, label it 'Unknown'
                newLabels(i) = "Unknown";
        end
    end

    % 4) Convert the new string array to a categorical
    % Define a final order of bins
    finalOrder = ["0","1-5","6-20","21plus","Missing","Unknown"];
    LTGTCcat = categorical(newLabels, finalOrder, 'Ordinal', false);

    % 5) Remove categories we didn't actually use
    usedCats  = categories(LTGTCcat);
    toRemove  = setdiff(usedCats, usedCats);  % likely none, but this is optional
    LTGTCcat  = removecats(LTGTCcat, toRemove);
end


% function LTGTCcat = defineLTGTCbins(LTGTCraw)
% % DEFINELTGTCBINSLARGE  Merges narrower 'c_...' categories into fewer, bigger bins.
% %
% % For example:
% %   c_0       -> '0'
% %   c_1,c_2-3 -> '1-3'
% %   c_4-5,c_6-10 -> '4-10'
% %   c_11-20   -> '11-20'
% %   c_21-50   -> '21-50'
% %   c_51-100  -> '51-100'
% %   c_>100    -> '>100'
% %   MISSING   -> 'Missing'
% %
% % If your data lacks some categories, they just won't appear in the final set.
% 
%     % 1) Ensure we have a categorical array
%     if iscell(LTGTCraw) || isstring(LTGTCraw)
%         LTGTCraw = categorical(LTGTCraw);
%     end
%     if ~iscategorical(LTGTCraw)
%         error('LTGTC must be a categorical array with c_0, c_1, c_2-3, c_4-5, etc.');
%     end
%     
%     % 2) Prepare a new string array of the same size
%     newLabels = repmat("Missing", size(LTGTCraw));  % default to "Missing"
%     
%     % 3) Assign merged bins for each row
%     for i = 1:numel(LTGTCraw)
%         oldCat = string(LTGTCraw(i));  % e.g. "c_1", "c_2-3"
%         
%         switch oldCat
%             case "c_0"
%                 newLabels(i) = "0";
%             case {"c_1","c_2-3"}
%                 newLabels(i) = "1-3";
%             case {"c_4-5","c_6-10"}
%                 newLabels(i) = "4-10";
%             case "c_11-20"
%                 newLabels(i) = "11-20";
%             case "c_21-50"
%                 newLabels(i) = "21-50";
%             case "c_51-100"
%                 newLabels(i) = "51-100";
%             case "c_>100"
%                 newLabels(i) = ">100";
%             case "MISSING"
%                 newLabels(i) = "Missing";
%             otherwise
%                 % If there's an unexpected category, either keep it or label it "Unknown"
%                 newLabels(i) = "Unknown";
%         end
%     end
%     
%     % 4) Convert the new string array to a categorical
%     % Define a final order of bins
%     finalOrder = ["0","1-3","4-10","11-20","21-50","51-100",">100","Missing","Unknown"];
%     LTGTCcat = categorical(newLabels, finalOrder, 'Ordinal', false);
%     
%     % 5) You might want to remove categories you didn't use
%     LTGTCcat = removecats(LTGTCcat, setdiff(categories(LTGTCcat), categories(LTGTCcat)));
% end


% function LTGTCcat = defineLTGTCbins(LTGTCraw)
% % DEFINELTGTCBINS  Renames each 'c_...' category to a simpler form,
% %   keeping them distinct (no merging).
% %
% % The function tries to rename:
% %   'c_0'     -> '0'
% %   'c_1'     -> '1'
% %   'c_2-3'   -> '2-3'
% %   'c_4-5'   -> '4-5'
% %   'c_6-10'  -> '6-10'
% %   'c_11-20' -> '11-20'
% %   'c_21-50' -> '21-50'
% %   'c_51-100'-> '51-100'
% %   'c_>100'  -> '>100'
% %   'MISSING' -> 'Missing'
% %
% % If your data doesn't contain some of these old categories, they won't appear.
% % Finally, we reorder them in a logical ascending sequence.
% 
%     % 1) Ensure input is categorical
%     if iscell(LTGTCraw) || isstring(LTGTCraw)
%         LTGTCraw = categorical(LTGTCraw);
%     end
%     if ~iscategorical(LTGTCraw)
%         error('LTGTC must be a categorical array with c_0, c_2-3, etc.');
%     end
% 
%     % 2) Current categories (e.g. c_0, c_1, c_2-3, c_4-5, c_6-10, etc.)
%     oldCats = categories(LTGTCraw);
% 
%     % 3) Define parallel arrays for old -> new mapping (1:1)
%     oldNames = { 'c_0','c_1','c_2-3','c_4-5','c_6-10','c_11-20','c_21-50','c_51-100','c_>100','MISSING' };
%     newNames = { '0',  '1',  '2-3',   '4-5',   '6-10',  '11-20',  '21-50',  '51-100', '>100',  'Missing' };
% 
%     % Intersect to see which categories actually appear in the data
%     [~, idxInMap, ~] = intersect(oldNames, oldCats, 'stable');
%     if isempty(idxInMap)
%         % no recognized categories found -> no changes
%         warning('None of the expected LTGTC categories found in the data. Using raw categories.');
%         LTGTCcat = LTGTCraw;
%         return;
%     end
% 
%     renameOld = oldNames(idxInMap);  % e.g. {'c_0','c_6-10','c_1',...}
%     renameNew = newNames(idxInMap);  % e.g. {'0','6-10','1',...}
% 
%     % 4) Rename only the oldCats that exist
%     LTGTCrenamed = renamecats(LTGTCraw, renameOld, renameNew);
% 
%     % 5) Reorder in ascending sequence
%     finalOrder = {'0','1','2-3','4-5','6-10','11-20','21-50','51-100','>100','Missing'};
% 
%     usedCats   = categories(LTGTCrenamed);              % which new names actually appear
%     catInData  = intersect(finalOrder, usedCats, 'stable');  % keep only what's used
% 
%     LTGTCcat   = reordercats(LTGTCrenamed, catInData);
% end