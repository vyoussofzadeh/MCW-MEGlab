speechonset_time = tt(idx_good);
speechonset_time = speechonset_time(ia);
L = length(speechonset_time);

speechonset_sample = ipoints_all(idx_good);
speechonset_sample = speechonset_sample(ia);

figure,
plot(speechonset_time, '*'),
hline(thre,'b',['mean:', num2str(thre),'sec']),
box off;
set(gca,'color','none');
title('Speech onset'),
ylabel('Time (sec)');
xlabel('Trials');
set(gca,'Xtick', 1:L,'XtickLabel',1:L);
xlim([1 L]);
set(gca,'FontSize',10,'XTickLabelRotation',90);


%%
ep_data.all = cln_data;
ep_data.pst = [];

%% PST
tinterval = 700;
trl_res = [];
for i = 1:length(cln_data.trial)   
    trl = cln_data.trial{i};
    trl_res.trl{i} = trl(:, speechonset_sample(i)-tinterval:speechonset_sample(i));
    ttrl = cln_data.time{i};
%     trl_res.time{i} = ttrl(:, speechonset_sample(1)-1000:speechonset_sample(1));
    trl_res.time{i} = linspace(-tinterval/1000,0,tinterval+1);
end

ep_data.pst = cln_data;
ep_data.pst.trial = trl_res.trl;
ep_data.pst.time = trl_res.time;


%% BSL
tinterval = 300;
trl_res = [];
for i = 1:length(cln_data.trial)
    trl = cln_data.trial{i};
    idx = find(cln_data.time{1} == 0);     idx1 = find(cln_data.time{1} == -tinterval/1000);
    trl_res.trl{i} = trl(:, idx1:idx);
    ttrl = cln_data.time{i};
    %     trl_res.time{i} = ttrl(:, speechonset_sample(1)-1000:speechonset_sample(1));
    trl_res.time{i} = linspace(-tinterval/1000,0,tinterval+1);
end

ep_data.bsl = cln_data;
ep_data.bsl.trial = trl_res.trl;
ep_data.bsl.time = trl_res.time;


%- Appending data
cfg = [];
ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);

%% Data-covarinace estimation
t_data = vy_timelock(ep_data);

toi = [-700, 0; -300, 0];