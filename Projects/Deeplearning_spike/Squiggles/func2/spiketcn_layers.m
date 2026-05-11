function lg = spiketcn_layers(inputSize, numClasses, params, clsLayer)
% Residual TCN (dup-safe wiring) for sequences shaped [C x T] (cell arrays)
% inputSize : #channels
% numClasses: e.g., 2
% params (optional): struct with fields:
%   stemWidth (64), blockWidths ([64 64 96 128]), kernels ([3 7 15]),
%   dilations ([1 2 4 8]), dropout (0.30), useGAPGMP (true)
% clsLayer (optional): weighted classificationLayer named 'cls'

if nargin < 3 || isempty(params), params = struct(); end
stemWidth  = getf(params,'stemWidth',64);
blockW     = getf(params,'blockWidths',[64 64 96 128]);
kernels    = getf(params,'kernels',[3 7 15]);
dilsIn     = getf(params,'dilations',[1 2 4 8]);
dropP      = getf(params,'dropout',0.30);
useGAPGMP  = getf(params,'useGAPGMP',true);

assert(numel(blockW)==numel(dilsIn),'blockWidths and dilations must match.');

lg = layerGraph();

% ----- Input -----
in = sequenceInputLayer(inputSize,"Normalization","none","Name","in");
lg = addLayers(lg,in);

% ----- Stem (add individually) -----
stem_conv = convolution1dLayer(1,stemWidth,"Padding","same","Name","stem_conv","WeightsInitializer","he");
stem_bn   = batchNormalizationLayer("Name","stem_bn");
stem_relu = reluLayer("Name","stem_relu");
stem_drop = dropoutLayer(dropP*0.5,"Name","stem_drop");
lg = addLayers(lg,stem_conv); lg = addLayers(lg,stem_bn);
lg = addLayers(lg,stem_relu); lg = addLayers(lg,stem_drop);
lg = connectIfAbsent(lg,"in","stem_conv");
lg = connectIfAbsent(lg,"stem_conv","stem_bn");
lg = connectIfAbsent(lg,"stem_bn","stem_relu");
lg = connectIfAbsent(lg,"stem_relu","stem_drop");

prevName = "stem_drop";
prevCh   = stemWidth;

% ----- Residual TCN blocks -----
for b = 1:numel(blockW)
    tag = sprintf("b%d",b);
    [lg, prevName, prevCh] = tcnResBlock_dupSafe(lg, prevName, prevCh, blockW(b), kernels, dilsIn(b), dropP, tag);
end

% ----- Head -----
if useGAPGMP
    gap = globalAveragePooling1dLayer("Name","gap");
    gmp = globalMaxPooling1dLayer("Name","gmp");
    cat = concatenationLayer(1,2,"Name","gap_gmp_cat");
    head_drop = dropoutLayer(dropP,"Name","head_drop");
    head_fc   = fullyConnectedLayer(64,"Name","head_fc");
    head_relu = reluLayer("Name","head_relu");
    fc  = fullyConnectedLayer(numClasses,"Name","fc");
    sm  = softmaxLayer("Name","sm");
    cls = classificationLayer("Name","cls");

    for L = [gap gmp cat head_drop head_fc head_relu fc sm cls], lg = addLayers(lg,L); end
    lg = connectIfAbsent(lg,prevName,"gap");
    lg = connectIfAbsent(lg,prevName,"gmp");
    lg = connectIfAbsent(lg,"gap","gap_gmp_cat/in1");
    lg = connectIfAbsent(lg,"gmp","gap_gmp_cat/in2");
    lg = connectIfAbsent(lg,"gap_gmp_cat","head_drop");
else
    gap = globalAveragePooling1dLayer("Name","gap");
    head_drop = dropoutLayer(dropP,"Name","head_drop");
    head_fc   = fullyConnectedLayer(64,"Name","head_fc");
    head_relu = reluLayer("Name","head_relu");
    fc  = fullyConnectedLayer(numClasses,"Name","fc");
    sm  = softmaxLayer("Name","sm");
    cls = classificationLayer("Name","cls");

    for L = [gap head_drop head_fc head_relu fc sm cls], lg = addLayers(lg,L); end
    lg = connectIfAbsent(lg,prevName,"gap");
    lg = connectIfAbsent(lg,"gap","head_drop");
end

lg = connectIfAbsent(lg,"head_drop","head_fc");
lg = connectIfAbsent(lg,"head_fc","head_relu");
lg = connectIfAbsent(lg,"head_relu","fc");
lg = connectIfAbsent(lg,"fc","sm");
lg = connectIfAbsent(lg,"sm","cls");

% Weighted classifier (optional)
if nargin>=4 && ~isempty(clsLayer)
    lg = replaceLayer(lg,"cls",clsLayer);
end
end

% ============================ Sub-functions ==============================

function [lg, outName, outCh] = tcnResBlock_dupSafe(lg, inName, inCh, outFilters, kernels, dilation, dropP, tag)
% Multi-kernel dilated conv branches -> concat -> BN/ReLU/Drop -> 1x1 -> BN
% Skip: identity or 1x1 -> BN
% Add -> ReLU

nBr = numel(kernels);
catName = tag + "_cat";
cat = depthConcatenationLayer(nBr,"Name",catName);
lg = addLayers(lg,cat);

branchOuts = strings(1,nBr);
for i = 1:nBr
    k   = kernels(i);
    cNm = sprintf("%s_conv_%d",tag,k);
    bNm = sprintf("%s_bn_%d",tag,k);
    rNm = sprintf("%s_relu_%d",tag,k);
    conv = convolution1dLayer(k, floor(outFilters/nBr), ...
        "Padding","same","DilationFactor",dilation,"Name",cNm,"WeightsInitializer","he");
    bn   = batchNormalizationLayer("Name",bNm);
    relu = reluLayer("Name",rNm);
    lg = addLayers(lg,conv); lg = addLayers(lg,bn); lg = addLayers(lg,relu);
    lg = connectIfAbsent(lg,inName,cNm);
    lg = connectIfAbsent(lg,cNm,bNm);
    lg = connectIfAbsent(lg,bNm,rNm);
    branchOuts(i) = rNm;
end
for i=1:nBr
    lg = connectIfAbsent(lg,branchOuts(i),char(catName + "/in" + i));
end

post_bn   = tag + "_post_bn";
post_relu = tag + "_post_relu";
post_drop = tag + "_drop";
projNm    = tag + "_proj";
proj_bn   = tag + "_proj_bn";
bnL   = batchNormalizationLayer("Name",post_bn);
rL    = reluLayer("Name",post_relu);
dL    = dropoutLayer(dropP,"Name",post_drop);
pL    = convolution1dLayer(1,outFilters,"Padding","same","Name",projNm,"WeightsInitializer","he");
pBN   = batchNormalizationLayer("Name",proj_bn);
for L=[bnL rL dL pL pBN], lg = addLayers(lg,L); end

lg = connectIfAbsent(lg,catName,post_bn);
lg = connectIfAbsent(lg,post_bn,post_relu);
lg = connectIfAbsent(lg,post_relu,post_drop);
lg = connectIfAbsent(lg,post_drop,projNm);
lg = connectIfAbsent(lg,projNm,proj_bn);

% Skip
if inCh ~= outFilters
    skipNm  = tag + "_skip";
    skipBNm = tag + "_skip_bn";
    sL  = convolution1dLayer(1,outFilters,"Padding","same","Name",skipNm,"WeightsInitializer","he");
    sBN = batchNormalizationLayer("Name",skipBNm);
    lg = addLayers(lg,sL); lg = addLayers(lg,sBN);
    lg = connectIfAbsent(lg,inName,skipNm);
    lg = connectIfAbsent(lg,skipNm,skipBNm);
    skipOut = skipBNm;
else
    skipOut = inName;
end

addNm = tag + "_add";
outNm = tag + "_out";
aL = additionLayer(2,"Name",addNm);
oR = reluLayer("Name",outNm);
lg = addLayers(lg,aL); lg = addLayers(lg,oR);

lg = connectIfAbsent(lg,proj_bn,addNm + "/in1");
lg = connectIfAbsent(lg,skipOut,addNm + "/in2");
lg = connectIfAbsent(lg,addNm,outNm);

outName = char(outNm);
outCh   = outFilters;
end

function val = getf(s,field,default)
if isfield(s,field) && ~isempty(s.(field)), val = s.(field); else, val = default; end
end

function lg = connectIfAbsent(lg, src, dst)
% Guarded connection to avoid "connection already exists"
conn = lg.Connections;
if ~any(strcmp(conn.Source, src) & strcmp(conn.Destination, dst))
    lg = connectLayers(lg, src, dst);
end
end
