% function lg = spatiotemporal_dualbranch_layers(inputSize, numClasses, clsLayer)
% % inputSize = #channels (e.g., 306)
% % numClasses = 2 (NoSpike/Spike)
% % clsLayer  = optional weighted classificationLayer named 'cls'


function lg = spatiotemporal_dualbranch_layers(inputSize, numClasses, clsLayer)
% Spatio-temporal dual-branch network (temporal BiLSTM + spatial TD FC + GRU)
% inputSize  : #channels (e.g., 306)
% numClasses : #classes (e.g., 2)
% clsLayer   : optional weighted classificationLayer named 'cls'

% ===== Init / Input =====
lg = layerGraph();
in = sequenceInputLayer(inputSize,"Normalization","none","Name","in");
lg = addLayers(lg, in);

% ===== Temporal branch (smaller; stronger dropout) =====
t_bilstm1 = bilstmLayer(64,"OutputMode","last","Name","t_bilstm1");
t_drop    = dropoutLayer(0.30,"Name","t_drop");
t_fc      = fullyConnectedLayer(64,"Name","t_fc");

lg = addLayers(lg, t_bilstm1);
lg = addLayers(lg, t_drop);
lg = addLayers(lg, t_fc);

lg = connectIfAbsent(lg,"in","t_bilstm1");
lg = connectIfAbsent(lg,"t_bilstm1","t_drop");
lg = connectIfAbsent(lg,"t_drop","t_fc");

% ===== Spatial branch (time-distributed mixing + GRU; BN + dropout) =====
fold  = sequenceFoldingLayer("Name","fold");
td_fc = fullyConnectedLayer(96,"Name","td_fc");   % mixes across channels per time-step
td_bn = batchNormalizationLayer("Name","td_bn");
td_rl = reluLayer("Name","td_relu");
td_do = dropoutLayer(0.30,"Name","td_drop");
unf   = sequenceUnfoldingLayer("Name","unfold");
s_gru = gruLayer(48,"OutputMode","last","Name","s_gru");
s_do  = dropoutLayer(0.30,"Name","s_drop");
s_fc  = fullyConnectedLayer(64,"Name","s_fc");

lg = addLayers(lg, fold);
lg = addLayers(lg, td_fc);
lg = addLayers(lg, td_bn);
lg = addLayers(lg, td_rl);
lg = addLayers(lg, td_do);
lg = addLayers(lg, unf);
lg = addLayers(lg, s_gru);
lg = addLayers(lg, s_do);
lg = addLayers(lg, s_fc);

lg = connectIfAbsent(lg,"in","fold/in");
lg = connectIfAbsent(lg,"fold/out","td_fc");
lg = connectIfAbsent(lg,"td_fc","td_bn");
lg = connectIfAbsent(lg,"td_bn","td_relu");
lg = connectIfAbsent(lg,"td_relu","td_drop");
lg = connectIfAbsent(lg,"td_drop","unfold/in");
lg = connectIfAbsent(lg,"fold/miniBatchSize","unfold/miniBatchSize");
lg = connectIfAbsent(lg,"unfold/out","s_gru");
lg = connectIfAbsent(lg,"s_gru","s_drop");
lg = connectIfAbsent(lg,"s_drop","s_fc");

% ===== Fusion + head (BN + bottleneck + dropout) =====
concat    = concatenationLayer(1,2,"Name","concat");
head_bn   = batchNormalizationLayer("Name","head_bn");
head_drop = dropoutLayer(0.40,"Name","head_drop");
head_fc   = fullyConnectedLayer(64,"Name","head_fc");
fc        = fullyConnectedLayer(numClasses,"Name","fc");
sm        = softmaxLayer("Name","sm");
cls       = classificationLayer("Name","cls");

lg = addLayers(lg, concat);
lg = addLayers(lg, head_bn);
lg = addLayers(lg, head_drop);
lg = addLayers(lg, head_fc);
lg = addLayers(lg, fc);
lg = addLayers(lg, sm);
lg = addLayers(lg, cls);

lg = connectIfAbsent(lg,"t_fc","concat/in1");
lg = connectIfAbsent(lg,"s_fc","concat/in2");
lg = connectIfAbsent(lg,"concat","head_bn");
lg = connectIfAbsent(lg,"head_bn","head_drop");
lg = connectIfAbsent(lg,"head_drop","head_fc");
lg = connectIfAbsent(lg,"head_fc","fc");
lg = connectIfAbsent(lg,"fc","sm");
lg = connectIfAbsent(lg,"sm","cls");

% ===== Optional weighted classification layer =====
if nargin>=3 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end

% -------------------------------------------------------------------------
function lg = connectIfAbsent(lg, src, dst)
% Only add the connection if it doesn't already exist (avoids duplicate-edge error)
conn = lg.Connections;
if ~any(strcmp(conn.Source, src) & strcmp(conn.Destination, dst))
    lg = connectLayers(lg, src, dst);
end
end

% % Init
% lg = layerGraph();
% 
% % ===== Input =====
% in = sequenceInputLayer(inputSize, "Normalization","none", "Name","in");
% lg = addLayers(lg, in);
% 
% % ===== Temporal branch: BiLSTM over time =====
% t_bilstm1 = bilstmLayer(128, "OutputMode","last", "Name","t_bilstm1");
% t_drop    = dropoutLayer(0.2, "Name","t_drop");
% t_fc      = fullyConnectedLayer(128, "Name","t_fc");
% 
% lg = addLayers(lg, t_bilstm1);
% lg = addLayers(lg, t_drop);
% lg = addLayers(lg, t_fc);
% 
% lg = connectIfAbsent(lg, "in",        "t_bilstm1");
% lg = connectIfAbsent(lg, "t_bilstm1", "t_drop");
% lg = connectIfAbsent(lg, "t_drop",    "t_fc");
% 
% % ===== Spatial branch: time-distributed channel mixing + GRU =====
% fold  = sequenceFoldingLayer("Name","fold");
% td_fc = fullyConnectedLayer(128,"Name","td_fc");   % mixes across channels per time-step
% td_rl = reluLayer("Name","td_relu");
% unf   = sequenceUnfoldingLayer("Name","unfold");
% s_gru = gruLayer(64,"OutputMode","last","Name","s_gru");
% s_fc  = fullyConnectedLayer(128,"Name","s_fc");
% 
% lg = addLayers(lg, fold);
% lg = addLayers(lg, td_fc);
% lg = addLayers(lg, td_rl);
% lg = addLayers(lg, unf);
% lg = addLayers(lg, s_gru);
% lg = addLayers(lg, s_fc);
% 
% lg = connectIfAbsent(lg, "in",                 "fold/in");
% lg = connectIfAbsent(lg, "fold/out",           "td_fc");
% lg = connectIfAbsent(lg, "td_fc",              "td_relu");
% lg = connectIfAbsent(lg, "td_relu",            "unfold/in");
% lg = connectIfAbsent(lg, "fold/miniBatchSize", "unfold/miniBatchSize");
% lg = connectIfAbsent(lg, "unfold/out",         "s_gru");
% lg = connectIfAbsent(lg, "s_gru",              "s_fc");
% 
% % ===== Fuse branches + head =====
% concat   = concatenationLayer(1,2,"Name","concat");
% head_drop= dropoutLayer(0.3,"Name","head_drop");
% fc       = fullyConnectedLayer(numClasses,"Name","fc");
% sm       = softmaxLayer("Name","sm");
% cls      = classificationLayer("Name","cls");
% 
% lg = addLayers(lg, concat);
% lg = addLayers(lg, head_drop);
% lg = addLayers(lg, fc);
% lg = addLayers(lg, sm);
% lg = addLayers(lg, cls);
% 
% lg = connectIfAbsent(lg, "t_fc",  "concat/in1");
% lg = connectIfAbsent(lg, "s_fc",  "concat/in2");
% lg = connectIfAbsent(lg, "concat","head_drop");
% lg = connectIfAbsent(lg, "head_drop","fc");
% lg = connectIfAbsent(lg, "fc","sm");
% lg = connectIfAbsent(lg, "sm","cls");
% 
% % ===== Optional: weighted classification layer =====
% if nargin >= 3 && ~isempty(clsLayer)
%     lg = replaceLayer(lg, "cls", clsLayer);
% end
% end
% 
% % -------------------------------------------------------------------------
% function lg = connectIfAbsent(lg, src, dst)
% % Only add the connection if it doesn't already exist
% conn = lg.Connections;
% exists = any(strcmp(conn.Source, src) & strcmp(conn.Destination, dst));
% if ~exists
%     lg = connectLayers(lg, src, dst);
% end
% end
