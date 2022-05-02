

bsdir = '/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/dSPM/Parcellation/HCPMMP1/Mean/parc/HCPMMP1/event 2';
dd = rdir(fullfile (bsdir,'/**/*.pkl'));

clc
sFiles = []; k=1; sFiles_names = [];
for jj=1:length(dd)
    disp(dd(jj).name);
    pause,
    delete(dd(jj).name);
end