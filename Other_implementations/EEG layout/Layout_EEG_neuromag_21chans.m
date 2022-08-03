%% EEG neuromag layout

lay_23 = load('easycapM23.mat');

cfg = [];
cfg.layout = lay_23.lay; % EEG
lay = ft_prepare_layout(cfg);
ft_layoutplot(cfg);

lay_21_neuromag = lay_23;
labels = {'Fp1', 'Fp2', 'F7', 'F3' , 'Fz', 'F4', 'F8', 'FT9', 'FT10', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2'};  

lay_21_neuromag.lay.pos = [];
lay_21_neuromag.lay.width = [];
lay_21_neuromag.lay.height = [];
lay_21_neuromag.lay.label = [];

for i=1:length(labels)
  idx = find(strcmp(lay_23.lay.label, labels{i})==1);
%   idx_label(i) = idx; 
  
  lay_21_neuromag.lay.pos(i,:) = lay_23.lay.pos(idx,:);
  lay_21_neuromag.lay.width(i,:) = lay_23.lay.width(idx);
  lay_21_neuromag.lay.height(i,:) = lay_23.lay.height(idx);
  lay_21_neuromag.lay.label{i} = lay_23.lay.label{idx};

end

cfg = [];
cfg.layout = lay_21_neuromag.lay; % EEG
lay = ft_prepare_layout(cfg);
ft_layoutplot(cfg);


save('neuromag_21.mat','lay_21_neuromag')


%%
labels2 =    {'EEG001','EEG003','EEG008','EEG010','EEG012','EEG014','EEG016',...
'EEG017','EEG024','EEG026','EEG028','EEG030','EEG032','EEG034','EEG044',...
'EEG046','EEG048','EEG050','EEG052','EEG057','EEG059'};

lay_21_neuromag2 = lay_21_neuromag;
for i=1:length(lay_21_neuromag.lay.label)
  lay_21_neuromag2.lay.label{i} = labels2{i};
end

cfg = [];
cfg.layout = lay_21_neuromag2.lay; % EEG
lay = ft_prepare_layout(cfg);
ft_layoutplot(cfg);


save('neuromag_21_2.mat','lay_21_neuromag2')
