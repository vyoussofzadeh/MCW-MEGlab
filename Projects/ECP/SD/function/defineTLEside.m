function TLEcat = defineTLEside(TLEsideVec)
% DEFINETLESIDE  Ensures the input TLE side data is a categorical
%   variable with recognized categories 'Left', 'Bilateral', 'Right'.
%
% USAGE:
%   TLEcat = defineTLEside(T1_epil_measures.TLEside);
%
% INPUT:
%   TLEsideVec : Nx1 categorical or string array (or cell array of chars)
%                with values like 'Left', 'Right', 'Bilateral', etc.
%
% OUTPUT:
%   TLEcat     : Nx1 categorical array with categories
%                {'Left','Bilateral','Right'}
%
% EXAMPLE:
%   % Suppose T1_epil_measures.TLEside has strings like
%   % {'Left','Left','Right','Bilateral','Right',...}
%   T1_epil_measures.TLEside = defineTLEside(T1_epil_measures.TLEside);

    % 1) Convert to categorical if not already
    if ~iscategorical(TLEsideVec)
        TLEsideVec = categorical(TLEsideVec);
    end

    % 2) Define the valid set of categories
    validCats = {'Left','Bilateral','Right'};

    % 3) Remove any categories not in validCats
    allCats   = categories(TLEsideVec);
    catsToRemove = setdiff(allCats, validCats);
    if ~isempty(catsToRemove)
        TLEsideVec = removecats(TLEsideVec, catsToRemove);
    end

    % 4) Reorder categories to a consistent order
    TLEcat = reordercats(TLEsideVec, validCats);
end
