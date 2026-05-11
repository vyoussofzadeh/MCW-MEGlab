function net = train_temporal_cnn(Wtr, Ytr, Wva, Yva, Tfix, cats, w, useGPU)
layers = [
  sequenceInputLayer(1,'Normalization','none','Name','in')
  convolution1dLayer(7, 32, 'Padding','same'), reluLayer
  convolution1dLayer(5, 32, 'Padding','same'), reluLayer
  maxPooling1dLayer(2,'Stride',2)
  convolution1dLayer(3, 64, 'Padding','same'), reluLayer
  globalAveragePooling1dLayer
  dropoutLayer(0.2)
  fullyConnectedLayer(numel(cats))
  softmaxLayer
  classificationLayer('Classes',cats,'ClassWeights',w)];
opts = trainingOptions('adam', ...
  'InitialLearnRate',3e-4,'MaxEpochs',40,'MiniBatchSize',256, ...
  'SequenceLength',Tfix,'SequencePaddingDirection','right', ...
  'Shuffle','every-epoch', 'ValidationData',{Wva,Yva}, ...
  'ValidationFrequency',max(10,ceil(numel(Wtr)/256)), ...
  'ValidationPatience',8, 'ExecutionEnvironment', ternaryExec(useGPU), ...
  'Plots','training-progress','Verbose',false);
net = trainNetwork(Wtr, Ytr, layers, opts);
end

function net = train_spatial_cnn(ImTr, YTr, ImVa, YVa, cats, w, useGPU)
layers = [
  imageInputLayer([size(ImTr,1) size(ImTr,2) 1], 'Normalization','none')
  convolution2dLayer(5, 32, 'Padding','same'), reluLayer
  maxPooling2dLayer(2,'Stride',2)
  convolution2dLayer(3, 64, 'Padding','same'), reluLayer
  maxPooling2dLayer(2,'Stride',2)
  convolution2dLayer(3, 128,'Padding','same'), reluLayer
  globalAveragePooling2dLayer
  dropoutLayer(0.3)
  fullyConnectedLayer(numel(cats))
  softmaxLayer
  classificationLayer('Classes',cats,'ClassWeights',w)];
opts = trainingOptions('adam', ...
  'InitialLearnRate',3e-4,'MaxEpochs',40,'MiniBatchSize',128, ...
  'Shuffle','every-epoch', 'ValidationData',{ImVa,YVa}, ...
  'ValidationFrequency',max(10,ceil(size(ImTr,4)/128)), ...
  'ValidationPatience',8, 'ExecutionEnvironment', ternaryExec(useGPU), ...
  'Plots','training-progress','Verbose',false);
net = trainNetwork(ImTr, YTr, layers, opts);
end

function e = ternaryExec(flag), e = 'auto'; if flag && gpuDeviceCount>0, e='gpu'; end, end
