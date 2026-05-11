function lg = lv_cadenet_lite_layers(C, W, D, numClasses, clsLayer)
% Fast spatiotemporal encoder with projection -> concat (no size mismatch).
% Input is imageInputLayer([C W D]).

lg = layerGraph();
lg = addLayers(lg, imageInputLayer([C W D], "Normalization","none", "Name","in"));

% ----- Encoder 1 -----
e1 = [
    convolution2dLayer([1 7], 32, "Padding","same","Name","e1_tconv")
    reluLayer("Name","e1_relu1")
    convolution2dLayer([7 1], 32, "Padding","same","Name","e1_sconv")
    reluLayer("Name","e1_relu2")
    maxPooling2dLayer([2 2],"Stride",[2 2],"Name","e1_pool")
    globalAveragePooling2dLayer("Name","e1_gap")
    fullyConnectedLayer(64,"Name","e1_fc")   % <-- project to same size
    reluLayer("Name","e1_prelu")];
lg = addLayers(lg, e1);
lg = connectLayers(lg,"in","e1_tconv");

% ----- Encoder 2 -----
e2 = [
    convolution2dLayer([1 5], 64, "Padding","same","Name","e2_tconv")
    reluLayer("Name","e2_relu1")
    convolution2dLayer([5 1], 64, "Padding","same","Name","e2_sconv")
    reluLayer("Name","e2_relu2")
    maxPooling2dLayer([2 2],"Stride",[2 2],"Name","e2_pool")
    globalAveragePooling2dLayer("Name","e2_gap")
    fullyConnectedLayer(64,"Name","e2_fc")   % <-- same size
    reluLayer("Name","e2_prelu")];
lg = addLayers(lg, e2);
lg = connectLayers(lg,"e1_pool","e2_tconv");

% ----- Encoder 3 -----
e3 = [
    convolution2dLayer([1 3], 128, "Padding","same","Name","e3_tconv")
    reluLayer("Name","e3_relu1")
    convolution2dLayer([3 1], 128, "Padding","same","Name","e3_sconv")
    reluLayer("Name","e3_relu2")
    globalAveragePooling2dLayer("Name","e3_gap")
    fullyConnectedLayer(64,"Name","e3_fc")   % <-- same size
    reluLayer("Name","e3_prelu")];
lg = addLayers(lg, e3);
lg = connectLayers(lg,"e2_pool","e3_tconv");

% ----- Fusion head -----
% Use depthConcatenationLayer to concat along channels (dim 3).
head = [
    depthConcatenationLayer(3,"Name","concat_gaps")  % 3 inputs (64 each) -> 1x1x192
    dropoutLayer(0.3,"Name","head_drop")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")];
lg = addLayers(lg, head);

% Wire the three projected GAPs to concat (each is 1x1x64)
lg = connectLayers(lg,"e1_prelu","concat_gaps/in1");
lg = connectLayers(lg,"e2_prelu","concat_gaps/in2");
lg = connectLayers(lg,"e3_prelu","concat_gaps/in3");

% Optional weighted class layer
if nargin>=5 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end
