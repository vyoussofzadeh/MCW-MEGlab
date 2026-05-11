function lgraph = setLearnRateFactor(lgraph, layerNames, fac)
for i = 1:numel(layerNames)
    L = lgraph.Layers(strcmp({lgraph.Layers.Name}, layerNames(i)));
    if isempty(L), continue; end
    if isprop(L,'WeightLearnRateFactor'), L.WeightLearnRateFactor = fac; end
    if isprop(L,'BiasLearnRateFactor'),   L.BiasLearnRateFactor   = fac; end
    lgraph = replaceLayer(lgraph, layerNames(i), L);
end
end
