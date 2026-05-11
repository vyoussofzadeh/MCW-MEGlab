function lg = spatiotemporal_cnn_layers(C, T, numInCh, numClasses, clsLayer)
% Spatio-Temporal 2D CNN over [channels x time]
% C,T      : image height (channels) and width (time bins)
% numInCh  : 1 (MEG only) or 3 (MEG + X/Y spatial maps)
% numClasses: e.g., 2
% clsLayer : optional weighted classificationLayer named 'cls'

lg = layerGraph();

% ----- Input -----
in = imageInputLayer([C T numInCh], "Normalization","none", "Name","in");
lg = addLayers(lg, in);

% ===== Inception-style block: temporal, spatial, joint =====
tconv = convolution2dLayer([1 7], 32, "Padding","same", "Name","tconv1");
tact  = reluLayer("Name","t_relu");

sconv = convolution2dLayer([7 1], 32, "Padding","same", "Name","sconv1");
sact  = reluLayer("Name","s_relu");

jconv = convolution2dLayer([3 3], 32, "Padding","same", "Name","jconv1");
jact  = reluLayer("Name","j_relu");

cat1  = depthConcatenationLayer(3,"Name","cat1");
bn1   = batchNormalizationLayer("Name","bn1");
do1   = dropoutLayer(0.25,"Name","drop1");

% Add individually
lg = addLayers(lg, tconv); lg = addLayers(lg, tact);
lg = addLayers(lg, sconv); lg = addLayers(lg, sact);
lg = addLayers(lg, jconv); lg = addLayers(lg, jact);
lg = addLayers(lg, cat1);  lg = addLayers(lg, bn1); lg = addLayers(lg, do1);

% Connect with guards
lg = connectIfAbsent(lg,"in","tconv1");
lg = connectIfAbsent(lg,"tconv1","t_relu");

lg = connectIfAbsent(lg,"in","sconv1");
lg = connectIfAbsent(lg,"sconv1","s_relu");

lg = connectIfAbsent(lg,"in","jconv1");
lg = connectIfAbsent(lg,"jconv1","j_relu");

lg = connectIfAbsent(lg,"t_relu","cat1/in1");
lg = connectIfAbsent(lg,"s_relu","cat1/in2");
lg = connectIfAbsent(lg,"j_relu","cat1/in3");
lg = connectIfAbsent(lg,"cat1","bn1");
lg = connectIfAbsent(lg,"bn1","drop1");

% ===== Block 2: anisotropic + pooling =====
tconv2 = convolution2dLayer([1 5], 64, "Padding","same", "Name","tconv2");
t2relu = reluLayer("Name","t2_relu");
sconv2 = convolution2dLayer([5 1], 64, "Padding","same", "Name","sconv2");
s2relu = reluLayer("Name","s2_relu");
jconv2 = convolution2dLayer([3 3], 96, "Padding","same", "Name","jconv2");
j2relu = reluLayer("Name","j2_relu");
bn2    = batchNormalizationLayer("Name","bn2");
pool2  = maxPooling2dLayer([1 2], "Stride",[1 2], "Name","pool2");
drop2  = dropoutLayer(0.25,"Name","drop2");

lg = addLayers(lg, tconv2); lg = addLayers(lg, t2relu);
lg = addLayers(lg, sconv2); lg = addLayers(lg, s2relu);
lg = addLayers(lg, jconv2); lg = addLayers(lg, j2relu);
lg = addLayers(lg, bn2);    lg = addLayers(lg, pool2); lg = addLayers(lg, drop2);

lg = connectIfAbsent(lg,"drop1","tconv2");
lg = connectIfAbsent(lg,"tconv2","t2_relu");
lg = connectIfAbsent(lg,"t2_relu","sconv2");
lg = connectIfAbsent(lg,"sconv2","s2_relu");
lg = connectIfAbsent(lg,"s2_relu","jconv2");
lg = connectIfAbsent(lg,"jconv2","j2_relu");
lg = connectIfAbsent(lg,"j2_relu","bn2");
lg = connectIfAbsent(lg,"bn2","pool2");
lg = connectIfAbsent(lg,"pool2","drop2");

% ===== Block 3 =====
tconv3 = convolution2dLayer([1 3], 128,"Padding","same","Name","tconv3");
t3relu = reluLayer("Name","t3_relu");
sconv3 = convolution2dLayer([3 1], 128,"Padding","same","Name","sconv3");
s3relu = reluLayer("Name","s3_relu");
jconv3 = convolution2dLayer([3 3], 192,"Padding","same","Name","jconv3");
j3relu = reluLayer("Name","j3_relu");
bn3    = batchNormalizationLayer("Name","bn3");
pool3  = maxPooling2dLayer([1 2], "Stride",[1 2], "Name","pool3");
drop3  = dropoutLayer(0.30,"Name","drop3");

lg = addLayers(lg, tconv3); lg = addLayers(lg, t3relu);
lg = addLayers(lg, sconv3); lg = addLayers(lg, s3relu);
lg = addLayers(lg, jconv3); lg = addLayers(lg, j3relu);
lg = addLayers(lg, bn3);    lg = addLayers(lg, pool3); lg = addLayers(lg, drop3);

lg = connectIfAbsent(lg,"drop2","tconv3");
lg = connectIfAbsent(lg,"tconv3","t3_relu");
lg = connectIfAbsent(lg,"t3_relu","sconv3");
lg = connectIfAbsent(lg,"sconv3","s3_relu");
lg = connectIfAbsent(lg,"s3_relu","jconv3");
lg = connectIfAbsent(lg,"jconv3","j3_relu");
lg = connectIfAbsent(lg,"j3_relu","bn3");
lg = connectIfAbsent(lg,"bn3","pool3");
lg = connectIfAbsent(lg,"pool3","drop3");

% ===== Head =====
gap    = globalAveragePooling2dLayer("Name","gap");
hdrop  = dropoutLayer(0.40,"Name","head_drop");
hfc    = fullyConnectedLayer(64,"Name","head_fc");
fc     = fullyConnectedLayer(numClasses,"Name","fc");
sm     = softmaxLayer("Name","sm");
cls    = classificationLayer("Name","cls");

lg = addLayers(lg, gap); lg = addLayers(lg, hdrop);
lg = addLayers(lg, hfc); lg = addLayers(lg, fc);
lg = addLayers(lg, sm);  lg = addLayers(lg, cls);

lg = connectIfAbsent(lg,"drop3","gap");
lg = connectIfAbsent(lg,"gap","head_drop");
lg = connectIfAbsent(lg,"head_drop","head_fc");
lg = connectIfAbsent(lg,"head_fc","fc");
lg = connectIfAbsent(lg,"fc","sm");
lg = connectIfAbsent(lg,"sm","cls");

% Weighted classifier (optional)
if nargin>=5 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end

% ---- helper ----
function lg = connectIfAbsent(lg, src, dst)
conn = lg.Connections;
if ~any(strcmp(conn.Source, src) & strcmp(conn.Destination, dst))
    lg = connectLayers(lg, src, dst);
end
end
