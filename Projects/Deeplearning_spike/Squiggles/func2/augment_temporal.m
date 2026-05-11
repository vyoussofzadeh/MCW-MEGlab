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