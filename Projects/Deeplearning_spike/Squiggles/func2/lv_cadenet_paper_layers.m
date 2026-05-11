function lg = lv_cadenet_paper_layers(C, L, D, numClasses, clsLayer)
% LV-CadeNet (paper-style, MATLAB-only)
% Input  : imageInputLayer([C L D])  (C=channels, L=time, D=7 maps)
% Stages : (temporal conv -> spatial conv -> pool) × 2 + bottleneck
% Fusion : {GAP, GMP} pooled from each stage ? concat ? FC ? softmax

assert(D>=1,'Depth D must be >=1 (expect D=7: raw+6 LV maps)');

lg = layerGraph();
lg = addLayers(lg, imageInputLayer([C L D],"Normalization","none","Name","in"));

% ===== Encoder: Stage 1 =====
e1 = [
    convolution2dLayer([1 7], 32, "Padding","same","Name","e1_tconv")
    reluLayer("Name","e1_tr")
    convolution2dLayer([7 1], 32, "Padding","same","Name","e1_sconv")
    reluLayer("Name","e1_sr")
    maxPooling2dLayer([2 2], "Stride",[2 2], "Name","e1_pool")
    globalAveragePooling2dLayer("Name","e1_gap")
    globalMaxPooling2dLayer("Name","e1_gmp")];
lg = addLayers(lg, e1);
lg = connectLayers(lg,"in","e1_tconv");

% ===== Encoder: Stage 2 =====
e2 = [
    convolution2dLayer([1 5], 64, "Padding","same","Name","e2_tconv")
    reluLayer("Name","e2_tr")
    convolution2dLayer([5 1], 64, "Padding","same","Name","e2_sconv")
    reluLayer("Name","e2_sr")
    maxPooling2dLayer([2 2], "Stride",[2 2], "Name","e2_pool")
    globalAveragePooling2dLayer("Name","e2_gap")
    globalMaxPooling2dLayer("Name","e2_gmp")];
lg = addLayers(lg, e2);
lg = connectLayers(lg,"e1_pool","e2_tconv");

% ===== Encoder: Stage 3 (bottleneck, no downsample) =====
e3 = [
    convolution2dLayer([1 3], 128, "Padding","same","Name","e3_tconv")
    reluLayer("Name","e3_tr")
    convolution2dLayer([3 1], 128, "Padding","same","Name","e3_sconv")
    reluLayer("Name","e3_sr")
    globalAveragePooling2dLayer("Name","e3_gap")
    globalMaxPooling2dLayer("Name","e3_gmp")];
lg = addLayers(lg, e3);
lg = connectLayers(lg,"e2_pool","e3_tconv");

% ===== Fusion head =====
% Concat pooled features from all stages (GAP+GMP -> 6 inputs)
% Each pooled tensor is 1x1xC_i; we concat along depth (channels)
lg = addLayers(lg, depthConcatenationLayer(6,"Name","concat_pools"));
lg = connectLayers(lg,"e1_gap","concat_pools/in1");
lg = connectLayers(lg,"e1_gmp","concat_pools/in2");
lg = connectLayers(lg,"e2_gap","concat_pools/in3");
lg = connectLayers(lg,"e2_gmp","concat_pools/in4");
lg = connectLayers(lg,"e3_gap","concat_pools/in5");
lg = connectLayers(lg,"e3_gmp","concat_pools/in6");

head = [
    dropoutLayer(0.3,"Name","head_drop")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")];
lg = addLayers(lg, head);
lg = connectLayers(lg,"concat_pools","head_drop");

% Optional weighted class layer
if nargin>=5 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end
