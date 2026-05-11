function dlnet = train_with_consistency(trainL_img, trainL_lab, trainU_img, varargin)
% trainL_img: [C W 1 NL], trainU_img: [C W 1 NU], labels categorical
p=inputParser; addParameter(p,'Epochs',30); addParameter(p,'MB',128);
addParameter(p,'LR',3e-4); addParameter(p,'LambdaU',1.0); parse(p,varargin{:});
E=p.Results.Epochs; MB=p.Results.MB; LR=p.Results.LR; lamU=p.Results.LambdaU;

[C,W,~,~]=size(trainL_img);
lgraph = emsnet_lite_backbone(C,W);
% classifier head directly on feats
cls = [
    globalAveragePooling2dLayer("Name","gap")
    fullyConnectedLayer(numel(categories(trainL_lab)),"Name","fc")
    softmaxLayer("Name","sm")];
lgraph = addLayers(lgraph, cls); lgraph = connectLayers(lgraph,"feat_out","gap");
dlnet = dlnetwork(lgraph);

NL = size(trainL_img,4); NU=size(trainU_img,4);
iter = ceil(max(NL,NU)/MB); vel=[]; ag=[]; asq=[];

for epoch=1:E
  for it=1:iter
    % labeled minibatch
    idxL = randi(NL,[1 min(MB,NL)]);
    Xl = trainL_img(:,:,:,idxL); Yl = trainL_lab(idxL);
    % unlabeled minibatch
    idxU = randi(NU,[1 min(MB,NU)]);
    Xu = trainU_img(:,:,:,idxU);
    % two simple augs
    Xu1 = weak_aug(Xu); Xu2 = weak_aug(Xu);

    [loss,grads]= dlfeval(@consistencyLoss, dlnet, Xl, Yl, Xu1, Xu2, lamU);
    [dlnet,vel,ag,asq] = adamupdate(dlnet, grads, vel, ag, asq, (epoch-1)*iter+it, LR,0.9,0.999,1e-8);
  end
  fprintf('Epoch %d/%d done\n',epoch,E);
end
end

function [loss,grads]=consistencyLoss(dlnet, Xl, Yl, Xu1, Xu2, lamU)
dlXl = dlarray(single(Xl),'SSCB');
dlXu1= dlarray(single(Xu1),'SSCB');
dlXu2= dlarray(single(Xu2),'SSCB');

% forward
pl = forward(dlnet, dlXl);   % probs
pu1= forward(dlnet, dlXu1);
pu2= forward(dlnet, dlXu2);

% CE on labeled
T = onehotencode(Yl,1); T = dlarray(single(T),'CB'); % [Classes x B]
plCB = permute(pl,[3 4 2 1]); plCB = squeeze(plCB);  % [Classes x B]
ce = -sum(T.*log(plCB+1e-7),1); ce = mean(ce);

% consistency on unlabeled (MSE between probs)
pu1CB = squeeze(permute(pu1,[3 4 2 1]));
pu2CB = squeeze(permute(pu2,[3 4 2 1]));
cons = mean(sum((pu1CB - pu2CB).^2,1));

loss = ce + lamU*cons;
grads = dlgradient(loss, dlnet.Learnables);
end

function X2 = weak_aug(X)
% light jitter/crop/flip
[C,W,~,B]=size(X); X2=X;
% small time crop (keep 90100% width)
keep = randi([round(0.9*W), W]);
start = randi([1, W-keep+1]);
for b=1:B
    x = X(:,:,1,b);
    xb = zeros(C,W,'like',x);
    seg = x(:,start:start+keep-1);
    xb(:,1:size(seg,2)) = seg;
    % random sign flip
    if rand<0.5, xb=-xb; end
    X2(:,:,1,b) = xb;
end
end
