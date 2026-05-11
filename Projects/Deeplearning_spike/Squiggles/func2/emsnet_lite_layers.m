function lg = emsnet_lite_layers(megH, megW, numClasses, clsLayer)
% No custom layers; fast on older MATLAB.

lg = layerGraph();
lg = addLayers(lg, imageInputLayer([megH megW 1], "Normalization","none","Name","in"));

% ----- Local 1D (temporal) branch -----
local = [
    PermuteHW1_to_1WH("perm_to_1WH")                             % [H,W,1]?[1,W,H]
    convolution2dLayer([1 7], 64, "Padding","same","Name","tconv7_shared")
    reluLayer("Name","l_relu1")
    convolution2dLayer([1 1], 64, "Padding","same","Name","l_pointwise")
    reluLayer("Name","l_relu2")
    maxPooling2dLayer([1 2], "Stride",[1 2], "Name","l_pool")
    globalAveragePooling2dLayer("Name","l_gap")
];
lg = addLayers(lg, local);
lg = connectLayers(lg,"in","perm_to_1WH");

% ----- Global 2D (spatio-temporal) branch -----
global2d = [
    convolution2dLayer([7 7], 32, "Padding","same","Name","g_conv1")
    reluLayer("Name","g_relu1")
    maxPooling2dLayer([2 2], "Stride",[2 2], "Name","g_pool1")

    convolution2dLayer([3 3], 64, "Padding","same","Name","g_conv2")
    reluLayer("Name","g_relu2")
    maxPooling2dLayer([2 2], "Stride",[2 2], "Name","g_pool2")

    globalAveragePooling2dLayer("Name","g_gap")
];
lg = addLayers(lg, global2d);
lg = connectLayers(lg,"in","g_conv1");

% ----- Fusion head -----
head = [
    concatenationLayer(1,2,"Name","concat")
    dropoutLayer(0.3,"Name","head_drop")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")
];
lg = addLayers(lg, head);
lg = connectLayers(lg,"l_gap","concat/in1");
lg = connectLayers(lg,"g_gap","concat/in2");

if nargin>=4 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end
