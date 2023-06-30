spike = ft_read_spike('p029_sort_final_01.nex');

cfg = [];
cfg.spikechannel = {'sig002a_wf', 'sig003a_wf'}; % select only the two single-units
spike = ft_spike_select(cfg, spike);


cfg             = [];
cfg.fsample     = 40000;
cfg.interpolate = 1; % keep the density of samples as is
[wave, spikeCleaned] = ft_spike_waveform(cfg,spike);


for k = [1 2]
  figure,
  sl = squeeze(wave.dof(k,:,:))>1000; % only keep samples with enough spikes
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl)),'k') % factor 10^6 to get microseconds
  hold on

  % plot the standard deviation
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl))+sqrt(squeeze(wave.var(k,:,sl))),'k--')
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl))-sqrt(squeeze(wave.var(k,:,sl))),'k--')

  axis tight
  set(gca,'Box', 'off')
  xlabel('time')
  ylabel('normalized voltage')
end


f = squeeze(wave.avg(k,:,sl));

save('/MEG_data/LAB_MEMBERS/Vahab/Github/tools/megclinic_development/tools/func/spk_temp_ft.mat','f')