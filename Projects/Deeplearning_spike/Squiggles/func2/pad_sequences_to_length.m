% function Xpad = pad_sequences_to_length(X, Tfix)
% % Pad/crop each {C x Ti} to {C x Tfix} (right-pad with zeros)
% C = size(X{1},1);
% Xpad = cell(size(X));
% for i=1:numel(X)
%     xi = X{i}; Ti = size(xi,2);
%     if     Ti > Tfix, xi = xi(:,1:Tfix);
%     elseif Ti < Tfix, xi = [xi, zeros(C, Tfix-Ti, 'like', xi)];
%     end
%     Xpad{i} = xi;
% end
% end

function Xp = pad_sequences_to_length(Xcells, Tfix)
% Right-pad/crop to Tfix for a cell array of [C x T] matrices
assert(iscell(Xcells) && ~isempty(Xcells), 'Xcells must be a non-empty cell array.');
C = size(Xcells{1},1);
Xp = cell(size(Xcells));
for i = 1:numel(Xcells)
    x = Xcells{i};
    assert(size(x,1) == C, 'pad_sequences_to_length: epoch %d has C=%d (expected %d).', i, size(x,1), C);
    T = size(x,2);
    if T >= Tfix
        Xp{i} = x(:,1:Tfix);
    else
        Xp{i} = [x, zeros(C, Tfix-T, 'like', x)];  % <-- C rows, not 1
    end
end
end

