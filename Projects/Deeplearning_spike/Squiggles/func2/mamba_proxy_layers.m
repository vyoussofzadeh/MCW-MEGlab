function lg = mamba_proxy_layers(inputSize, numClasses, clsLayer)
% MAMBA-PROXY: dilated 1D convs + GLU gating + residuals + GRU

lg = layerGraph();

% input + pre-norm
lg = addLayers(lg, sequenceInputLayer(inputSize,'Normalization','none','Name','in'));
lg = addLayers(lg, layerNormalizationLayer('Name','ln0'));
lg = connectLayers(lg,'in','ln0');

% ---- Block A1 (dilation 1) ----
% main (conv + gate -> GLU)
lg = addLayers(lg, convolution1dLayer(7,128,'Padding','same','DilationFactor',1,'Name','convA1'));
lg = addLayers(lg, convolution1dLayer(7,128,'Padding','same','DilationFactor',1,'Name','gateA1'));
lg = addLayers(lg, sigmoidLayer('Name','sigA1'));
lg = addLayers(lg, multiplicationLayer(2,'Name','gluA1'));

% skip projection to match channels (306 -> 128)
lg = addLayers(lg, convolution1dLayer(1,128,'Padding','same','Name','proj0'));

% residual + norm
lg = addLayers(lg, additionLayer(2,'Name','resA1'));
lg = addLayers(lg, layerNormalizationLayer('Name','lnA1'));

% wires A1
lg = connectLayers(lg,'ln0','convA1');
lg = connectLayers(lg,'ln0','gateA1');
lg = connectLayers(lg,'gateA1','sigA1');
lg = connectLayers(lg,'convA1','gluA1/in1');
lg = connectLayers(lg,'sigA1','gluA1/in2');

lg = connectLayers(lg,'ln0','proj0');
lg = connectLayers(lg,'proj0','resA1/in1');   % projected skip
lg = connectLayers(lg,'gluA1','resA1/in2');   % GLU output

lg = connectLayers(lg,'resA1','lnA1');

% ---- Block A2 (dilation 2) ----
lg = addLayers(lg, convolution1dLayer(7,128,'Padding','same','DilationFactor',2,'Name','convA2'));
lg = addLayers(lg, convolution1dLayer(7,128,'Padding','same','DilationFactor',2,'Name','gateA2'));
lg = addLayers(lg, sigmoidLayer('Name','sigA2'));
lg = addLayers(lg, multiplicationLayer(2,'Name','gluA2'));
lg = addLayers(lg, additionLayer(2,'Name','resA2'));
lg = addLayers(lg, layerNormalizationLayer('Name','lnA2'));

lg = connectLayers(lg,'lnA1','convA2');
lg = connectLayers(lg,'lnA1','gateA2');
lg = connectLayers(lg,'gateA2','sigA2');
lg = connectLayers(lg,'convA2','gluA2/in1');
lg = connectLayers(lg,'sigA2','gluA2/in2');
lg = connectLayers(lg,'lnA1','resA2/in1');
lg = connectLayers(lg,'gluA2','resA2/in2');
lg = connectLayers(lg,'resA2','lnA2');

% ---- Block A3 (dilation 4) ----
lg = addLayers(lg, convolution1dLayer(7,128,'Padding','same','DilationFactor',4,'Name','convA3'));
lg = addLayers(lg, convolution1dLayer(7,128,'Padding','same','DilationFactor',4,'Name','gateA3'));
lg = addLayers(lg, sigmoidLayer('Name','sigA3'));
lg = addLayers(lg, multiplicationLayer(2,'Name','gluA3'));
lg = addLayers(lg, additionLayer(2,'Name','resA3'));
lg = addLayers(lg, layerNormalizationLayer('Name','lnA3'));

lg = connectLayers(lg,'lnA2','convA3');
lg = connectLayers(lg,'lnA2','gateA3');
lg = connectLayers(lg,'gateA3','sigA3');
lg = connectLayers(lg,'convA3','gluA3/in1');
lg = connectLayers(lg,'sigA3','gluA3/in2');
lg = connectLayers(lg,'lnA2','resA3/in1');
lg = connectLayers(lg,'gluA3','resA3/in2');
lg = connectLayers(lg,'resA3','lnA3');

% ---- Temporal integrator + head ----
lg = addLayers(lg, gruLayer(64,'OutputMode','last','Name','gru'));
lg = connectLayers(lg,'lnA3','gru');

head = [
    dropoutLayer(0.2,'Name','drop')
    fullyConnectedLayer(numClasses,'Name','fc')
    softmaxLayer('Name','sm')
    classificationLayer('Name','cls')];
lg = addLayers(lg, head);
lg = connectLayers(lg,'gru','drop');

if nargin>=3 && ~isempty(clsLayer)
    lg = replaceLayer(lg,'cls',clsLayer);
end
end
