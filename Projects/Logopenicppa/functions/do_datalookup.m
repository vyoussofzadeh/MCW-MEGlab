function datafile = do_datalookup(cfg_main)

indir = cfg_main.indir;
tag = cfg_main.tag;


%%
% tag = 'Run'; %- Semantic decision

% clc
% if exist(['datalog_',tag,'.mat'], 'file') == 2
%     load(['datalog_',tag,'.mat'])
% else
clear datafolder datafile subj_all sub
datafile_fif = [];
d = rdir([indir,['/**/tsss/*',tag,'*_raw.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    tkz = tokenize(datafile{i}(Index(end)+1:end),'_');
    task{i} = tkz{2};
    run{i} = tkz{1};
    task_run{i} = [tkz{1},'_', tkz{2}];
end
datafile_fif = vertcat(datafile_fif,datafile);
datafile_fif = datafile_fif';
%     save(['datalog_',tag,'.mat'],'datafile_fif','subj_all', 'sub')
% end
% disp(datafile_fif)

%%
datafile = [];
datafile.datafile_fif = datafile_fif;
datafile.task_run = task_run;
datafile.task = task;
datafile.run = run;

