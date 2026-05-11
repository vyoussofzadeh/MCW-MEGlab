function Xi = keep_grads(X, refLbl)
% Keep only gradiometers (MEGIN style ...2/3). Use after alignment.
isGrad = endsWith(refLbl,'2') | endsWith(refLbl,'3');
Xi = cellfun(@(z) z(isGrad,:), X, 'UniformOutput', false);
end
