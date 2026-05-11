function lgraph = morgoth1d_layers(numChannels, Tfix, args)
% numChannels = #features per time-step (e.g., K after PCA, or 306)
% Tfix        = fixed sequence length (samples)
% args        = struct with fields: dModel, patch, nHeads, nEnc, ffDim, numClasses

% ---- defaults ----
dModel     = getArg(args,'dModel',128);
patch      = getArg(args,'patch',20);
nHeads     = getArg(args,'nHeads',4);
nEnc       = getArg(args,'nEnc',4);
ffDim      = getArg(args,'ffDim',256);
numClasses = getArg(args,'numClasses',2);

nTok = floor(Tfix/patch);  % tokens after patching

lgraph = layerGraph();

% 1) sequence input (features = numChannels)
lgraph = addLayers(lgraph, sequenceInputLayer(numChannels, ...
    "Name","seq_in","Normalization","none","MinLength",Tfix));

% 2) patch embed = conv1d over time with stride=patch ? dModel channels
patchBlock = [
    convolution1dLayer(patch, dModel, ...
        "Stride", patch, ...
        "Padding","same", ...   % < was "none"
        "Name","patch_conv")
    layerNormalizationLayer("Name","patch_ln")
    reluLayer("Name","patch_relu")
];


% patchBlock = [
%     convolution1dLayer(patch, dModel, "Stride",patch, "Padding","none", "Name","patch_conv")
%     layerNormalizationLayer("Name","patch_ln")
%     reluLayer("Name","patch_relu")
% ];
lgraph = addLayers(lgraph, patchBlock);
lgraph = connectLayers(lgraph,"seq_in","patch_conv");

% 3) learnable positional embeddings (dModel x nTok) + add
pos = positionalEmbedding1dLayer(dModel, nTok, "Name","pos");
add = additionLayer(2,"Name","add_pos");
lgraph = addLayers(lgraph,pos);
lgraph = addLayers(lgraph,add);
lgraph = connectLayers(lgraph,"patch_relu","add_pos/in1");
lgraph = connectLayers(lgraph,"pos","add_pos/in2");

% 4) transformer encoder blocks
prev = "add_pos";
for i = 1:nEnc
    blk = encoderBlock(dModel, nHeads, ffDim, i);
    lgraph = addLayers(lgraph, blk.layers);
    % wire residuals inside the stack
    lgraph = connectLayers(lgraph, prev, sprintf("enc%d_mha_norm",i));
    lgraph = connectLayers(lgraph, prev, sprintf("enc%d_add1/in2",i));
    lgraph = connectLayers(lgraph, sprintf("enc%d_add1",i), sprintf("enc%d_ffn_norm",i));
    lgraph = connectLayers(lgraph, sprintf("enc%d_add1",i), sprintf("enc%d_add2/in2",i));
    prev = sprintf("enc%d_add2",i);
end

% 5) token pooling ? classifier head
head = [
    globalAveragePooling1dLayer("Name","tok_pool")   % average over tokens
    layerNormalizationLayer("Name","prehead_ln")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")
];
lgraph = addLayers(lgraph, head);
lgraph = connectLayers(lgraph, prev, "tok_pool");
end

function v = getArg(s, f, def)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = def; end
end

function blk = encoderBlock(dModel, nHeads, ffDim, idx)
% LN ? MHA ? Add, LN ? FFN ? Add
pre1 = layerNormalizationLayer("Name",sprintf("enc%d_mha_norm",idx));
mha  = multiheadAttentionLayer(nHeads, dModel, "Name",sprintf("enc%d_mha",idx), "Dropout",0.1);
add1 = additionLayer(2,"Name",sprintf("enc%d_add1",idx));

pre2 = layerNormalizationLayer("Name",sprintf("enc%d_ffn_norm",idx));
ffn  = [
    fullyConnectedLayer(ffDim,"Name",sprintf("enc%d_ff1",idx))
    reluLayer("Name",sprintf("enc%d_ffrelu",idx))
    dropoutLayer(0.1,"Name",sprintf("enc%d_ffdrop",idx))
    fullyConnectedLayer(dModel,"Name",sprintf("enc%d_ff2",idx))
];
add2 = additionLayer(2,"Name",sprintf("enc%d_add2",idx));

L = layerGraph(pre1);
L = addLayers(L,mha);  L = addLayers(L,add1);
L = addLayers(L,pre2); L = addLayers(L,ffn); L = addLayers(L,add2);

% q=k=v = normalized tokens
L = connectLayers(L, sprintf("enc%d_mha_norm",idx), sprintf("enc%d_mha/query",idx));
L = connectLayers(L, sprintf("enc%d_mha_norm",idx), sprintf("enc%d_mha/key",idx));
L = connectLayers(L, sprintf("enc%d_mha_norm",idx), sprintf("enc%d_mha/value",idx));
L = connectLayers(L, sprintf("enc%d_mha",idx),      sprintf("enc%d_add1/in1",idx));

% FFN + residual
L = connectLayers(L, sprintf("enc%d_add1",idx),     sprintf("enc%d_ffn_norm",idx));
L = connectLayers(L, sprintf("enc%d_ffn_norm",idx), sprintf("enc%d_ff1",idx));
L = connectLayers(L, sprintf("enc%d_ff2",idx),      sprintf("enc%d_add2/in1",idx));

blk.layers = L;
end
