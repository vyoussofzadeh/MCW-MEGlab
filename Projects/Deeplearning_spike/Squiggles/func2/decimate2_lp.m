function Xd = decimate2_lp(X, fs, fcut)
% Low-pass (Butterworth) then 2× decimation by pairwise averaging.
% Handles even/odd T safely.
% X : cell array of C×T
% fs: sampling rate (Hz)
% fcut: cutoff (Hz), e.g., 80

[b,a] = butter(4, fcut/(fs/2), 'low');   % 4th-order low-pass

Xd = cell(size(X));
for i = 1:numel(X)
    xi = X{i};
    if isempty(xi)
        Xd{i} = xi; 
        continue;
    end
    % Filter along time (transpose to T×C for filtfilt)
    xf = filtfilt(b, a, double(xi.')).';   % back to C×T

    T = size(xf,2);
    if T < 2
        % nothing to decimate
        Xd{i} = single(xf);
        continue;
    end

    % Make an even number of columns for pairwise averaging
    Teven = 2*floor(T/2);         % largest even = T
    x1 = xf(:, 1:2:Teven);        % C × (Teven/2)
    x2 = xf(:, 2:2:Teven);        % C × (Teven/2)
    y  = 0.5 * (x1 + x2);         % averaged pairs

    if T > Teven
        % If odd, append the last (unpaired) sample without averaging
        y = [y, xf(:, end)];
    end

    Xd{i} = single(y);
end
end
