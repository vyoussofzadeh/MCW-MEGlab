
%%
d = rdir([subjdir,'/cloze',num2str(k),'/*raw_sss.fif']);
if isempty(d)
    d = rdir([subjdir,'/cloze',num2str(k),'/*sss_raw.fif']);
end
if isempty(d)
    d = rdir([subjdir,'/cloze',num2str(k),'/*xtraClean_raw.fif']);
end
if isempty(d)
    d = rdir([subjdir,'/cloze',num2str(k),'/*_sss.fif']);
end

%%
if ~isempty(d)
    clear subj datafolder datafile datafile1
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        subj = datafile{i}(Index(end-2)+1:Index(end-1)-1);
    end
    datafile1 = datafile';
    disp(datafile1)
    if length(datafile1) > 1
        datasel = input('choose which data to analyze, row number:');
    else
        datasel = 1;
    end
    
    %- looking for alternative data location
else
    d = rdir([indir_second,'/*/*cloze',num2str(k),'*sss.fif']);
    
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        subj_all{i} = datafile{i}(Index(end-1)+1:Index(end)-1);
    end
    datafile1 = datafile';
    disp(datafile1)
    if length(datafile1) > 1
        datasel = input('choose which data to analyze, row number:');
        subj = subj_all{datasel};
    else
        datasel = 1;
        subj = subj_all{1};
    end
end
disp([subj, ' and,'])
disp([datafile1{datasel}, ' was selected for the analysis ...'])
disp('============');
    
%% Reading event files
% d1 = rdir([subjdir,'/cloze',num2str(k),'/*cloze1_eventInfo.mat']);
% 
% if ~isempty(d1)
%     load(d1.name)
% end
% 
% %% Reading log file
% % d2 = rdir([subjdir,'/cloze',num2str(k),'/*xtraClean_ave.log']);
% d2 = rdir([subjdir,'/cloze',num2str(k),'/*probe*_ave.log']);
% 
% clear dd
% if ~isempty(d2)
%     if length(d2)>1
%         for j=1:length(d2), dd{j,:}=d2(j).name; end
%         disp(dd)
%         sel = input('select file:');
% %         sel = 1;
%     else
%         sel = 1;
%     end
%     trig = importdata(d2(sel).name);
% end
% 
% Condition = [];
% sampletrig = [];
% for i=1:length(trig)
%     tmp = trig{i};
%     cac = textscan( tmp, '%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f', 'Whitespace',' ');
%     sampletrig(i) = cac{1}; Condition(i) = cac{4};
% end
% disp(d2(sel).name)