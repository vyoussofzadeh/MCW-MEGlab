function lg = emsnet_layers(megH, megW, numClasses, clsLayer)
% megH: #sensors (e.g., 306), megW: padded time steps, numClasses: 2
% clsLayer: optional weighted classificationLayer('Name','cls',...)

% ----- input -----
lg = layerGraph();
lg = addLayers(lg, imageInputLayer([megH megW 1], "Normalization","none", "Name","in"));

% ================== Local (per-sensor 1D temporal) branch ==================
% Permute to [1 x T x H] so channels=H; depthwise temporal conv; mix with 1x1.
% local = [
%     PermuteHW1_to_1WH("perm_to_1WH")
%     convolution2dLayer([1 7], 64, "Padding","same", "NumGroups", megH, "Name","dw_tconv7") % depthwise per sensor
%     reluLayer("Name","l_relu1")
%     convolution2dLayer([1 1], 64, "Padding","same", "Name","l_pointwise")                  % channel mixing
%     reluLayer("Name","l_relu2")
%     maxPooling2dLayer([1 2], "Stride",[1 2], "Name","l_pool")
%     globalAveragePooling2dLayer("Name","l_gap")];


local = [
    PermuteHW1_to_1WH("perm_to_1WH")                % [H,W,1] -> [1,W,H] (H becomes channels)
    convolution2dLayer([1 7], 64, ...               % shared temporal filters across sensors
    "Padding","same", "Name","tconv7_shared")
    reluLayer("Name","l_relu1")
    convolution2dLayer([1 1], 64, ...               % channel mixing
    "Padding","same", "Name","l_pointwise")
    reluLayer("Name","l_relu2")
    maxPooling2dLayer([1 2], "Stride",[1 2], "Name","l_pool")
    globalAveragePooling2dLayer("Name","l_gap")
    ];
lg = addLayers(lg, local);
lg = connectLayers(lg,"in","perm_to_1WH");

% ================== Global (spatio-temporal 2D) branch ====================
% Locally-connected 2D (unshared) + shared 1x1 conv
locConn = LocallyConnected2dLayer([7 7], 32, [megH megW 1], "lc2d", true);
gb = [
    locConn
    reluLayer("Name","g_relu1")
    convolution2dLayer([1 1], 64, "Padding","same", "Name","g_1x1")
    reluLayer("Name","g_relu2")
    maxPooling2dLayer([2 2], "Stride",[2 2], "Name","g_pool")
    globalAveragePooling2dLayer("Name","g_gap")];
lg = addLayers(lg, gb);
lg = connectLayers(lg,"in","lc2d");

% ================== Fusion head ==================
head = [
    concatenationLayer(1,2,"Name","concat")
    dropoutLayer(0.3,"Name","head_drop")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")];
lg = addLayers(lg, head);

lg = connectLayers(lg,"l_gap","concat/in1");
lg = connectLayers(lg,"g_gap","concat/in2");

% Replace with weighted classifier if provided
if nargin>=4 && ~isempty(clsLayer)
    lg = replaceLayer(lg, "cls", clsLayer);
end
end
