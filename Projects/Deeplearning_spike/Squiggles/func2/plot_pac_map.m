% function plot_pac_map(PAC, ip, ia, fs, titleStr)
% Z = PAC.z{ip,ia};            % C×nb
% t = PAC.tbins{ip}/fs;        % seconds
% figure('Color','w'); imagesc(t, 1:size(Z,1), Z); axis tight xy
% xlabel('Time (s)'); ylabel('Channel'); colorbar
% title(sprintf('zPAC  phase [%d] / amp [%d]  %s', ip, ia, titleStr));
% end


function plot_pac_map(PAC, ip, ia, titleStr)
Z = PAC.z{ip,ia};  % C×nb
t = PAC.tbin_sec{ip};
imagesc(t, 1:size(Z,1), Z); axis tight xy
xlabel('time (s) wrt spike onset'); ylabel('channel'); colorbar
title(titleStr);
end

