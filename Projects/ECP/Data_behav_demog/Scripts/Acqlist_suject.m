clear; clc, close('all'); warning off

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

indir = '/rcc/stor1/projects/ECP/MEG/MEG_Work';

%% Step1,
% - Finding resting-state MEG data
% d = rdir([indir,['/**/tSSS','/ec*Rest*_raw*fif']]);
d = rdir([indir,['/**/tSSS','/ec*oise*run1*.fif']]);
clear subj_all datenumm
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    date{i} = d(i).date;
    datenumm(i) = d(i).datenum;
    Index = strfind(datafile{i}, '/');
    subj_all{i} = datafile{i}(Index(7)+1:Index(8)-1);
end
[a,b] = sort(datenumm);

datainfo = [datafile(b); date(b)]';

%%
T = [];
for i=1:length(b)
    T.sub{i} = subj_all{b(i)};
    tmp = date(b(i));
    T.date{i} = tmp{1};
end

%%
clear T1
T1 = table(T.sub', T.date');
T1.Properties.VariableNames{'Var1'} = 'Subject_name';
T1.Properties.VariableNames{'Var2'} = 'Session_date';

%%
outdir = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/ECP/Subj_demog';
filename = fullfile(outdir,'ECP_sublog.xlsx');
writetable(T1,filename)
filename = fullfile(outdir,'ECP_sublog.csv');
disp(filename)
writetable(T1,filename)
