function plot_spectral_scalars(S)
% S from feats_spectral_epoch
figure('Color','w'); 
subplot(2,1,1);
bar([S.bp_delta S.bp_theta S.bp_alpha S.bp_beta S.bp_gamma]); 
legend('\delta','\theta','\alpha','\beta','\gamma'); grid on
ylabel('Relative power'); title('Relative bandpowers per channel')

subplot(2,1,2);
bar(S.sent); ylabel('Spectral entropy'); xlabel('Channel'); grid on
end
