%% BCI dataset, Ulster University & Medical College of Wisconsin

% Script: Process (wPLI connectivity analysis)
% Project: BCI_nonstationarity
% Writtern by: Vahab Youssof Zadeh
% Update: 11/03/2022

clear; clc, close('all'); warning off

%% Analysis flags
flag = [];
flag.plotconn = 0;
flag.plottriggers = 1;
flag.anatomy_check = 1;
flag.dicsanalysis = 0;
flag.connanalysis = 1;
flag.network_exportfig = 2;

%% Paths
addpath('./data')
addpath('./run')
Run_setpath

%%
cd(outdir)

clear d_all sub session
d = rdir(fullfile(outdir,'/**/vs_*.mat'));
for i=1:length(d)
    tkz = tokenize(d(i).name,'/');
    D.df{i} = tkz{end};
    
    tkz1 = tokenize(D.df{i},'_');
    tkz2 = tokenize(tkz1{end},'.mat');
    
    sub{i} = tkz1{end-1};
    session{i} = tkz2{1};
    
    D.ss{i} = [sub{i}, '_', session{i}];
    D.sub{i} = sub{i};
    D.session{i} = session{i};
    D.datafile{i} = d(i).name;
end
disp(D.ss') % Data subject/Session
disp(D.datafile') % Data subject/Session

%%
D_sel=  [];
for i=1:length(D.ss)
    D_sel{i} = [num2str(i), ':', D.ss{i}];
end
disp(D_sel') % Data subject/Session

%% template Anatomy
load('anat')

%% Process
for df = 1:length(D.ss)
    
    subj = D.ss{df};
    disp(['Analyzing: ', subj])
    
    idx = strfind(D.datafile{df},'/');
    D.datafile{df}(1:idx(end)-1);
    
    load(D.datafile{df})
    
    idx = [];
    idx.hand = find(trialinfo == 1); % 'HAND'
    idx.feet = find(trialinfo == 2); % 'FEET'
    idx.word = find(trialinfo == 3); % 'WORD'
    idx.sub  = find(trialinfo == 4); % 'SUB'
    
    disp({'hand', 'feet', 'word', 'sub'})
    disp([length(idx.hand), length(idx.feet), length(idx.word), length(idx.sub)])
    
    vs_roi_tsk = [];
    vs_roi_tsk.hand = vs_roi; vs_roi_tsk.hand.trial = vs_roi.trial(idx.hand); vs_roi_tsk.hand.time = vs_roi.time(idx.hand);
    vs_roi_tsk.feet = vs_roi; vs_roi_tsk.feet.trial = vs_roi.trial(idx.feet); vs_roi_tsk.feet.time = vs_roi.time(idx.feet);
    vs_roi_tsk.word = vs_roi; vs_roi_tsk.word.trial = vs_roi.trial(idx.word); vs_roi_tsk.word.time = vs_roi.time(idx.word);
    vs_roi_tsk.sub  = vs_roi; vs_roi_tsk.sub.trial  = vs_roi.trial(idx.sub);  vs_roi_tsk.sub.time  = vs_roi.time(idx.sub);
    
    disp({'1) hand', '2) feet', '3) word', '4) sub'})
    task_sel = input('select task:');
    
    switch task_sel
        case 1
            vs_in = vs_roi_tsk.hand;
        case 2
            vs_in = vs_roi_tsk.feet;
        case 3
            vs_in = vs_roi_tsk.word;
        case 4
            vs_in = vs_roi_tsk.sub;
    end
    
    disp([num2str(vs_in.time{1}(1)),' to ', num2str(vs_in.time{1}(end)), 'sec'])
    toi = input('time range e.g, [0.4,4]: ');
    
    close all,
    [mn, idx_mn] = min(abs(toi(1) - vs_in.time{1}));
    [mx, idx_mx] = min(abs(toi(2) - vs_in.time{1}));
    
    vs_in_sel = vs_in;
    for i=1:length(vs_in_sel.trial)
        vs_in_sel.trial{i} = vs_in_sel.trial{i}(:,idx_mn:idx_mx);
        vs_in_sel.time{i} = vs_in_sel.time{i}(idx_mn:idx_mx);
    end
    
    cfg = [];
    cfg.flag = flag;
    source_conn = do_connanalysis(cfg, vs_in_sel); %- Conn analysis (wPLI)
    
    cfg = [];
    cfg.par =  'wpli_debiasedspctrm';
    cfg.label = vs_in_sel.label;
    source_net = do_networkanalysis(cfg, source_conn); %- network analysis (network degrees)
    
    cfg = [];
    cfg.par =  'wpli_debiasedspctrm';
    cfg.coor = vs_in_sel.coor;
    cfg.individual_headmodel = individual_headmodel;
    cfg.flag = flag;
    cfg.atlas = atlas;
    cfg.label = vs_in_sel.label;
    do_vis_netconn_analysis(cfg, source_conn) %- mapping
    
end


