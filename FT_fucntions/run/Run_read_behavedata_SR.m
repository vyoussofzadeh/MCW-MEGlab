
switch subj
    case '2233'
        if Clz == 1
            datafile = '/mnt/file1/binder/binder_data/spendl/Clozedata/2233/cloze1/cloze1_raw_defaultHead_sss.fif';
        end
end

%%
hdr = ft_read_header(datafile);
Fsample = hdr.Fs;

Index = strfind(hdr.label,{'STI101'});
Index = find(not(cellfun('isempty',Index)));

clear detTrig
detTrig=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','neuromag_fif','dataformat','neuromag_fif');
detTrig = (detTrig - min(detTrig));
detTrig=bitand(detTrig,255);
figure,plot(detTrig)

detResp=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','digital trigger','dataformat','neuromag_fif');
detResp = detResp - detTrig;
detResp = (detResp - min(detResp));
figure,plot(detResp)
Nsamps=length(detTrig);

%%
if isempty(find(detTrig == 8, 1))
    detTrig = detTrig+8;
end

%%
figure,plot(detTrig)
for i=1:length(sampletrig)
    hold on
    y = ylim; % current y-axis limits
    plot([sampletrig(i) sampletrig(i)],y)
end
%             saveas(gcf, 'trigger_time.png')
disp('check triggers, then enter to proceed!')
% pause,

%% response/reaction time
resp = (detResp - mean(detResp))./max(detTrig);
difResp=diff(resp);
difResp = abs(difResp);

con_id = [];
figure,plot(detTrig)
hold on
plot(resp, 'r')
tmp =[];
for i=1:length(sampletrig)
    hold on
    y = ylim; % current y-axis limits
    plot([sampletrig(i) sampletrig(i)],y)
    tmp(i)  = detTrig(sampletrig(i)+10);
%     disp(detTrig(sampletrig(i)+10))
    switch detTrig(sampletrig(i)+10) % Identifying conditions
        case 33
            con_id = [con_id;1];
        case 34
            con_id = [con_id;2];
        case 35
            con_id = [con_id;3];
        case 36
            con_id = [con_id;4];
        case 17
            con_id = [con_id;-1];
        case 18
            con_id = [con_id;-2];
        case 19
            con_id = [con_id;-3];
        case 20
            con_id = [con_id;-4];
    end
end

%%
con_id = Condition;

%%
figure,plot(detTrig)
hold on
plot(resp, 'r')

tmp = []; tmp2 = [];
k = 1; rs = []; ws =[];
for i=1:length(sampletrig)-1
    y = ylim; % current y-axis limits
    idx1 = sampletrig(i):sampletrig(i+1);
    tmp = difResp(idx1);
    tmp2 = detTrig(idx1);
    idx = find(tmp > 50);
    if length(idx) >= 2 
        rs(k) = idx1(idx(1)); % resposne sample
        plot([rs(k) rs(k)],y);
        %                     plot([sampletrig(i) sampletrig(i)],y)
        %                     rs(k) - sampletrig(i)
        idx2 = find(tmp2 == 8);
        if length(idx1) > 1
            ws(k) = idx1(idx2(2)); % waiting sample
            plot([ws(k) ws(k)],y);
            %                         rs(k) - sampletrig(i)
        else
            disp('issue with trl,', num2str(i))
        end
        k=k+1;
    end
end
%%
%- Reaction time
RT = [(rs - ws)/1e3];
con_id = con_id(1:length(RT));

%%
% idx = find(RT < 0);
% RT(idx) = [];
% con_id(idx)=[];

%%
figure,
subplot 211
bar(RT, 0.4); % (resposne - waiting) / fs
set(gca,'Xtick', 1:length(ws),'XtickLabel',1:length(ws)+1);
title(['Reaction time, mean:', num2str(mean((rs - ws)/1e3)), ' sec'])
xlabel('Trials'), ylabel('time (sec)'),
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
disp([num2str(mean((rs - ws)/1e3)), ' +- ', num2str(std((rs - ws)/1e3))]);
%
subplot 212
bar(con_id, 0.4),
set(gca,'Xtick', 1:length(RT),'XtickLabel',1:length(RT));
xlabel('Trials'), ylabel('time (sec)'),
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
title('data conditions');
saveas(gcf, fullfile('behav', 'Reaction_time.png'));

%%
% RT = (rs - ws)/1e3;
con_id1 = con_id(1:length(RT));

idx = [];
idx.expected = find(con_id1 == 1);
idx.unexpect = find(con_id1 == 2);
idx.anomalou = find(con_id1 == 3);
idx.pseudo   = find(con_id1 == 4);

% idx.Fexpected   = find(con_id1 == -1); 
% idx.Funexpect   = find(con_id1 == -2);
% idx.Fanomalou   = find(con_id1 == -3);
% idx.Fpseudo   = find(con_id1 == -4);

idx.Fexpected   = find(con_id1 == 11); 
idx.Funexpect   = find(con_id1 == 22);
idx.Fanomalou   = find(con_id1 == 33);
idx.Fpseudo   = find(con_id1 == 44);

%% Acc
Acc = [];
Acc.expected = (length(idx.expected) - length(idx.Fexpected)) / (length(idx.expected) + length(idx.Fexpected))*100;
Acc.unexpect = (length(idx.unexpect) - length(idx.Funexpect)) / (length(idx.unexpect) + length(idx.Funexpect))*100;
Acc.anomalou = (length(idx.anomalou) - length(idx.Fanomalou)) / (length(idx.anomalou) + length(idx.Fanomalou))*100;
Acc.pseudo = (length(idx.pseudo) - length(idx.Fpseudo)) / (length(idx.pseudo) + length(idx.Fpseudo))*100;
Acc.all = (Acc.expected + Acc.unexpect + Acc.anomalou + Acc.pseudo)/4;
disp(Acc)

figure,
bar([Acc.expected, Acc.unexpect, Acc.anomalou, Acc.pseudo, Acc.all], 0.4)
set(gca,'Xtick', 1:5,'XtickLabel',{'Expected','Unexpected','Anomalous','Pseudo', 'Mean'});
ylim([0,100])
grid on
set(gca,'color','none');
saveas(gcf, fullfile('behav', 'Accuracy.png'));

%%
RT1.exp = RT(idx.expected);
RT1.unex = RT(idx.unexpect);
RT1.anam = RT(idx.anomalou);
RT1.pseud = RT(idx.pseudo);

figure,
subplot 221, bar(RT1.exp),set(gca,'color','none'); title(['Ex, m: ',num2str(mean(RT1.exp))])
subplot 222, bar(RT1.unex),set(gca,'color','none'); title(['Unex, m: ',num2str(mean(RT1.unex))])
subplot 223, bar(RT1.anam),set(gca,'color','none'); title(['Anam, m: ',num2str(mean(RT1.anam))])
subplot 224, bar(RT1.pseud),set(gca,'color','none'); title(['Pseudo, m: ',num2str(mean(RT1.pseud))])

save(fullfile('behav', 'behavinfo'), 'RT', 'RT1', 'Acc','rs','ws','con_id1');
saveas(gcf, fullfile('behav', 'Reaction_time_cond.png'));

