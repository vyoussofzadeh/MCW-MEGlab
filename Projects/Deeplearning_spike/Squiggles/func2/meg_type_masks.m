function [magIdx, gradIdx] = meg_type_masks(refLbl)
magIdx  = endsWith(refLbl,'1');
gradIdx = endsWith(refLbl,'2') | endsWith(refLbl,'3');
end

function [W1, topo, eigsOut] = pca_by_type(Xcells, magIdx, gradIdx)
% For each epoch: grads -> PC1 waveform (1×T), mags -> PC1 topography (M×1)
n = numel(Xcells);
W1   = cell(n,1); topo = cell(n,1);
eigsOut = cell(n,1);
for i=1:n
    X = single(Xcells{i});              % C×T
    % grads
    G  = X(gradIdx,:);                  % GxT
    G  = G - mean(G,2);
    [U,~,~] = svd(G,'econ');           % U: GxG
    pc1w = U(:,1)' * G;                % 1×T waveform
    % mags
    M  = X(magIdx,:); M = M - mean(M,2);
    [U2,~,~] = svd(M,'econ');
    pc1t = U2(:,1);                    % M×1 topography
    W1{i}   = pc1w;
    topo{i} = pc1t;
    eigsOut{i} = struct('Ugrad',U(:,1),'Umag',U2(:,1));
end
end

function XYmag = get_channel_xy(refLbl, magIdx, sens)
% Return XY positions (Nx2) for magnetometers in refLbl order (magIdx)
% If FieldTrip sens provided, use it; else fall back to circle layout.
idx = find(magIdx);
Nm = numel(idx);
XYmag = zeros(Nm,2);
if ~isempty(sens) && isfield(sens,'label') && isfield(sens,'chanpos')
    % map labels -> rows
    L = cellstr(sens.label);
    for i=1:Nm
        k = find(strcmp(L, refLbl{idx(i)}),1);
        if ~isempty(k), XYmag(i,:) = sens.chanpos(k,1:2); end
    end
    if all(XYmag(:)==0), XYmag = circular_xy(Nm); end
else
    XYmag = circular_xy(Nm);
end
% normalize to [-1,1]
XYmag = XYmag - mean(XYmag,1);
XYmag = XYmag ./ max( max(abs(XYmag)), 1e-6);
end

function XY = circular_xy(N)
theta = linspace(0,2*pi,N+1); theta(end)=[];
XY = [cos(theta(:)) sin(theta(:))];
end

function I = topo_to_image(topo, XY, outSize)
% topo: M×1, XY: M×2, output: outSize×outSize single
F = scatteredInterpolant(XY(:,1),XY(:,2), topo, 'natural','nearest');
[xg,yg] = meshgrid(linspace(-1,1,outSize), linspace(-1,1,outSize));
I = single(F(xg, yg));
% fill NaNs if outside convex hull
m = isnan(I); if any(m(:)), I(m)=0; end
% zscore image (optional)
I = (I - mean(I(:))) / (std(I(:))+1e-6);
end

function [W_aug, Y_aug] = augment_temporal(Wcells, Y, step_ms, Nsteps, fs, doFlip)
% jitter Spikes by ±Nsteps of 'step_ms'; optional sign-flip for all classes
W_aug = {}; Y_aug = categorical(); 
step  = max(1, round(step_ms*fs/1000));
for i=1:numel(Wcells)
    w = Wcells{i};
    cls = Y(i);
    % always include original
    W_aug{end+1,1}=w; Y_aug(end+1,1)=cls; %#ok<AGROW>
    if cls=="Spike"
        for s=1:Nsteps
            k = step*s;
            W_aug{end+1,1} = circshift(w,[0, k]);  Y_aug(end+1,1)=cls;
            W_aug{end+1,1} = circshift(w,[0,-k]);  Y_aug(end+1,1)=cls;
        end
    end
    if doFlip
        W_aug{end+1,1} = -w;          Y_aug(end+1,1)=cls;
    end
end
end

function [I_aug, Y_aug] = augment_spatial(Icells, Y, rot90on, doFlip)
I_aug = {}; Y_aug = categorical();
for i=1:numel(Icells)
    im = Icells{i}; cls=Y(i);
    I_aug{end+1,1}=im; Y_aug(end+1,1)=cls;
    if rot90on
        I_aug{end+1,1}=rot90(im,1); Y_aug(end+1,1)=cls;
        I_aug{end+1,1}=rot90(im,2); Y_aug(end+1,1)=cls;
        I_aug{end+1,1}=rot90(im,3); Y_aug(end+1,1)=cls;
    end
    if doFlip
        I_aug{end+1,1}=-im; Y_aug(end+1,1)=cls;
    end
end
end

function [Xb, Yb] = balance_pairs(Xcells, Y)
% downsample majority to match minority
cats = categories(Y); c = countcats(Y);
[minN, minIdx] = min(c);
majIdx = 3-minIdx;                 % 1<->2 for binary
idxMin = find(Y==cats{minIdx});
idxMaj = find(Y==cats{majIdx});
idxMaj = idxMaj(randperm(numel(idxMaj), minN));
idx = [idxMin; idxMaj];
Xb = Xcells(idx);
Yb = Y(idx);
end

function [Im4d, Y1] = imgs4d(Icells, Y)
N = numel(Icells); H=size(Icells{1},1); W=size(Icells{1},2);
Im4d = zeros(H, W, 1, N, 'single');
for i=1:N, Im4d(:,:,1,i) = single(Icells{i}); end
Y1 = Y;
end

function Xp = pad_sequences_to_length(Xcells, Tfix)
Xp = cell(size(Xcells));
for i=1:numel(Xcells)
    x=Xcells{i}; t=size(x,2);
    if t>=Tfix, Xp{i}=x(:,1:Tfix);
    else, Xp{i}=[x, zeros(1, Tfix-t, 'like', x)];
    end
end
end
