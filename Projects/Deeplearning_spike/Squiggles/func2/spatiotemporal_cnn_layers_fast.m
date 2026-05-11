function lg = spatiotemporal_cnn_layers_fast(C, T, numInCh, numClasses, clsLayer)
% Fast Spatio-Temporal 2D CNN over [channels x time]
% C,T       : image height (channels) and width (time bins)
% numInCh   : 1 (MEG only) or 3 (MEG + X/Y maps)
% numClasses: e.g., 2
% clsLayer  : optional weighted classificationLayer named 'cls'

lg = layerGraph();

% ----- Input -----
in = imageInputLayer([C T numInCh], "Normalization","none", "Name","in");
lg = addLayers(lg, in);

% ==================== Block 1: lightweight, strided ======================
% Stride ALL branches along time so cat sizes match.
tconv1 = convolution2dLayer([1 7], 24, "Padding","same", "Stride",[1 2], "Name","tconv1");
t1relu = reluLayer("Name","t1_relu");

sconv1 = convolution2dLayer([7 1], 24, "Padding","same", "Stride",[1 2], "Name","sconv1");
s1relu = reluLayer("Name","s1_relu");

jconv1 = convolution2dLayer([3 3], 24, "Padding","same", "Stride",[1 2], "Name","jconv1");
j1relu = reluLayer("Name","j1_relu");

cat1  = depthConcatenationLayer(3,"Name","cat1");   % -> 72 ch
bn1   = batchNormalizationLayer("Name","bn1");
drop1 = dropoutLayer(0.2,"Name","drop1");

% add
for L = [tconv1 t1relu sconv1 s1relu jconv1 j1relu cat1 bn1 drop1]
    lg = addLayers(lg, L);
end
% wire
lg = connectLayers(lg,"in","tconv1");   lg = connectLayers(lg,"tconv1","t1_relu");
lg = connectLayers(lg,"in","sconv1");   lg = connectLayers(lg,"sconv1","s1_relu");
lg = connectLayers(lg,"in","jconv1");   lg = connectLayers(lg,"jconv1","j1_relu");
lg = connectLayers(lg,"t1_relu","cat1/in1");
lg = connectLayers(lg,"s1_relu","cat1/in2");
lg = connectLayers(lg,"j1_relu","cat1/in3");
lg = connectLayers(lg,"cat1","bn1");
lg = connectLayers(lg,"bn1","drop1");

% ==================== Block 2: anisotropic + pooling =====================
tconv2 = convolution2dLayer([1 5], 32, "Padding","same","Name","tconv2");
t2relu = reluLayer("Name","t2_relu");

sconv2 = convolution2dLayer([5 1], 32, "Padding","same","Name","sconv2");
s2relu = reluLayer("Name","s2_relu");

jconv2 = convolution2dLayer([3 3], 48, "Padding","same","Name","jconv2");
j2relu = reluLayer("Name","j2_relu");

bn2   = batchNormalizationLayer("Name","bn2");
pool2 = maxPooling2dLayer([1 2], "Stride",[1 2], "Name","pool2"); % time /2 again
drop2 = dropoutLayer(0.25,"Name","drop2");

for L = [tconv2 t2relu sconv2 s2relu jconv2 j2relu bn2 pool2 drop2]
    lg = addLayers(lg, L);
end
lg = connectLayers(lg,"drop1","tconv2");
lg = connectLayers(lg,"tconv2","t2_relu");
lg = connectLayers(lg,"t2_relu","sconv2");
lg = connectLayers(lg,"sconv2","s2_relu");
lg = connectLayers(lg,"s2_relu","jconv2");
lg = connectLayers(lg,"jconv2","j2_relu");
lg = connectLayers(lg,"j2_relu","bn2");
lg = connectLayers(lg,"bn2","pool2");
lg = connectLayers(lg,"pool2","drop2");

% ============================ Head =======================================
gap   = globalAveragePooling2dLayer("Name","gap");
hdrop = dropoutLayer(0.30,"Name","head_drop");
hfc   = fullyConnectedLayer(48,"Name","head_fc");
fc    = fullyConnectedLayer(numClasses,"Name","fc");
sm    = softmaxLayer("Name","sm");
cls   = classificationLayer("Name","cls");

for L = [gap hdrop hfc fc sm cls], lg = addLayers(lg, L); end
lg = connectLayers(lg,"drop2","gap");
lg = connectLayers(lg,"gap","head_drop");
lg = connectLayers(lg,"head_drop","head_fc");
lg = connectLayers(lg,"head_fc","fc");
lg = connectLayers(lg,"fc","sm");
lg = connectLayers(lg,"sm","cls");

% Optional: weighted classifier
if nargin>=5 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end
