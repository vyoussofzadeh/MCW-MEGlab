function [ImTr, ImVa, meta] = timelock_to_4d(tl_tr, tl_va, varargin)
% TIMELOCK_TO_4D Convert FieldTrip timelock (rpt×chan×time) to 4D tensors (C×T×D×N).
% 
% [ImTr, ImVa, meta] = timelock_to_4d(tl_tr, tl_va)
% [ImTr, ImVa, meta] = timelock_to_4d(...,'Depth4',true)   % add {raw,d1,d2,env}
% [ImTr, ImVa, meta] = timelock_to_4d(...,'Normalize',true) % per-channel z-norm (fit on TRAIN)
%
% Inputs:
%   tl_tr : timelock struct with keeptrials='yes' (fields: trial [Ntr×C×T], time, label)
%   tl_va : timelock struct with keeptrials='yes' (fields: trial [Nva×C×T], time, label)
%
% Outputs:
%   ImTr : C×T×D×Ntr  (single)
%   ImVa : C×T×D×Nva  (single)
%   meta : struct with fields C, T, D, mu, sd

% ---- parse options
p = inputParser;
p.addParameter('Depth4', false, @(x)islogical(x)||isscalar(x));
p.addParameter('Normalize', false, @(x)islogical(x)||isscalar(x));
p.parse(varargin{:});
wantDepth4 = logical(p.Results.Depth4);
doNorm     = logical(p.Results.Normalize);

% ---- sizes
Ntr = size(tl_tr.trial,1);
C   = size(tl_tr.trial,2);
T   = size(tl_tr.trial,3);
Nva = size(tl_va.trial,1);

% ---- base amplitude (D=1)
ImTr = single(permute(tl_tr.trial, [2 3 4 1]));   % C×T×1×Ntr
ImVa = single(permute(tl_va.trial, [2 3 4 1]));   % C×T×1×Nva

% ---- optional depth features {raw, d1, d2, env}
if wantDepth4
    % TRAIN
    Xraw = ImTr;
    Xd1  = cat(2, diff(ImTr,1,2), ImTr(:,end,:,:));
    Xd2  = cat(2, diff(ImTr,2,2), ImTr(:,end,:,:), ImTr(:,end,:,:));
    Xenv = zeros(C,T,1,Ntr,'single');
    for i = 1:Ntr
        Xi = ImTr(:,:,1,i);                         % C×T
        Xenv(:,:,1,i) = single(abs(hilbert(double(Xi).').').^2);
    end
    ImTr = cat(3, Xraw, Xd1, Xd2, Xenv);

    % VAL
    Xraw = ImVa;
    Xd1  = cat(2, diff(ImVa,1,2), ImVa(:,end,:,:));
    Xd2  = cat(2, diff(ImVa,2,2), ImVa(:,end,:,:), ImVa(:,end,:,:));
    Xenv = zeros(C,T,1,Nva,'single');
    for i = 1:Nva
        Xi = ImVa(:,:,1,i);
        Xenv(:,:,1,i) = single(abs(hilbert(double(Xi).').').^2);
    end
    ImVa = cat(3, Xraw, Xd1, Xd2, Xenv);
end

% ---- optional per-channel normalization (fit on TRAIN, apply to VAL)
mu = []; sd = [];
if doNorm
    Xall = reshape(ImTr, C, []);          % collapse (T*D*Ntr)
    mu = mean(Xall,2);
    sd = std(Xall,0,2); sd(sd<1e-12) = 1;
    for i = 1:size(ImTr,4)
        for d = 1:size(ImTr,3)
            ImTr(:,:,d,i) = (ImTr(:,:,d,i) - mu) ./ sd;
        end
    end
    for i = 1:size(ImVa,4)
        for d = 1:size(ImVa,3)
            ImVa(:,:,d,i) = (ImVa(:,:,d,i) - mu) ./ sd;
        end
    end
end

meta = struct('C',C,'T',T,'D',size(ImTr,3),'mu',mu,'sd',sd);
end
