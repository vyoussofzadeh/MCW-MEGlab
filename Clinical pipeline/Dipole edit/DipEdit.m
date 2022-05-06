
clc, clear,
Path =  '/MEG_data/epilepsy/garcia_brian/brainstorm_db/data';
Name =  'garcia_brian/Run03_spont_supine_raw_sss_ecgClean_raw/dipoles_Run03_event1002.mat';

fname = fullfile(Path,Name);
Dip_org = load(fname);

L  = length(Dip_org.Dipole);
disp(['L=', num2str(L)]);
% dipin = input('Dipole inputs in range, e.g, [1,200] ');
dipin  = [1, round(L/2)]

Dip_updt = Dip_org;
Dip_updt.DipoleNames = Dip_updt.DipoleNames(1:end-1);
% Dip_updt.DipoleNames = Dip_updt.DipoleNames(dipin(1):dipin(2));
% Dip_updt.Dipole = Dip_updt.Dipole(dipin(1):dipin(2));
Dip_updt.Dipole = Dip_updt.Dipole;
Dip_updt.Comment = 'dipole update';

idx = strfind(fname,'/');
cd(fname(1:idx(end)))
savetag = [fname(idx(end)+1:end-4), '_updt.mat'];
save(fullfile(savetag),'-struct', 'Dip_updt'),
