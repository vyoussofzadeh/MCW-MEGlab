function lgraph = emsnet_lite_backbone(C,W)
% Input: imageInputLayer([C W 1]) ; Output tap: layer named "feat_out"
lgraph = layerGraph();
lgraph = addLayers(lgraph, imageInputLayer([C W 1], "Normalization","none", "Name","in"));

% 2D spatio-temporal path
g = [
    convolution2dLayer([7 7], 32, "Padding","same", "Name","g1")
    reluLayer("Name","r1")
    convolution2dLayer([3 3], 64, "Padding","same", "Name","g2")
    reluLayer("Name","r2")
];
lgraph = addLayers(lgraph, g);
lgraph = connectLayers(lgraph, "in", "g1");

% 1D temporal path (shared across sensors)
t = [
    convolution2dLayer([1 7], 32, "Padding","same", "Name","t1")
    reluLayer("Name","tr1")
    convolution2dLayer([1 5], 64, "Padding","same", "Name","t2")
    reluLayer("Name","tr2")
];
lgraph = addLayers(lgraph, t);
lgraph = connectLayers(lgraph, "in", "t1");

% fuse both paths
lgraph = addLayers(lgraph, depthConcatenationLayer(2, "Name","cat"));
lgraph = connectLayers(lgraph, "r2",  "cat/in1");
lgraph = connectLayers(lgraph, "tr2", "cat/in2");

% projection head that preserves spatial size; output node = "feat_out"
h = [
    convolution2dLayer([1 1], 64, "Padding","same", "Name","proj")
    reluLayer("Name","hr")
    convolution2dLayer([3 3], 32, "Padding","same", "Name","up1")
    reluLayer("Name","hr2")
    convolution2dLayer([1 1], 1,  "Padding","same", "Name","feat_out")
];
lgraph = addLayers(lgraph, h);
lgraph = connectLayers(lgraph, "cat", "proj");
end
