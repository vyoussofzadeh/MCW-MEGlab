clear; clc, close('all'); warning off

%%
%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/data/MEG/epilepsy';

%%
addpath(genpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions'));

%%
disp('1: 2019')
disp('2: 2018');
disp('3: 2020');
disp('4: 2021');
Y = input('Year data were acquired: ');
switch Y
    case 1
        yr = '19';
    case 2
        yr = '18';
    case 3
        yr = '20';
    case 4
        yr = '21';
end

d = rdir([indir,['/**/', yr, '*/sessionLog.txt']]);
save([yr,'_sessionlog.mat'],'d');

%%
switch Y
    case {1,2}
        k = 1;
        clear Sub MRN DOB DOS Age DFNM PN spont Name
        
        for i = 1:length(d)
            [pathstr, name] = fileparts(d(i).name);
            filename = d(i).name;
            
            fid = fopen(filename,'rt');
            tmp = textscan(fid,'%s','Delimiter','\n');
            
            if ~isempty(tmp{1,1})
                
                tmp2 = tmp{1,1};
                
                Name{k} = tmp2{1};
                Sub{k}.name = Name{k};
                
                idx = find(contains(tmp2,'MRN')==1);
                if isempty(idx), idx = find(contains(tmp2,'MNR')==1); end
                Index = strfind(tmp2{idx}, ':');
                if isempty(Index), Index = strfind(tmp2{idx}, '='); end
                MRN{k} = tmp2{idx}(Index+1:end);
                Sub{k}.MRN = MRN{k};
                
                idx = find(contains(tmp2,'DOS')==1);
                if isempty(idx), idx = find(contains(tmp2,'DOs')==1); end
                Index = strfind(tmp2{idx}, ':');
                if isempty(Index), Index = strfind(tmp2{idx}, '='); end
                DOS{k} = tmp2{idx}(Index+1:end);
                Sub{k}.DOS = DOS{k};
                DOS_d = datenum(DOS{k},'mm/DD/YYYY');  % note year,month,day
                
                idx = find(contains(tmp2,'DOB')==1);
                Index = strfind(tmp2{idx}, ':');
                if isempty(Index), Index = strfind(tmp2{idx}, '='); end
                DOB{k} = tmp2{idx}(Index+1:end);
                Sub{k}.DOB = DOB{k};
                DOB_d = datenum(DOB{k},'mm/DD/YYYY');  % note year,month,day
                
                Age{k} = datestr(DOS_d-DOB_d,'YY.mm.dd');
                Sub{k}.Age = Age{k};
                
                Index = strfind(tmp2, 'spont');
                spont{k} = double(~isempty(cell2mat(Index)));
                Sub{k}.spont = spont{k};
                
                Index = strfind(tmp2, 'DFNM');
                DFNM{k} = double(~isempty(cell2mat(Index)));
                Sub{k}.DFNM = DFNM{k};
                
                Index = strfind(tmp2, 'PN');
                PN{k} = double(~isempty(cell2mat(Index)));
                Sub{k}.PN = PN{k};
                
                k = k+1;
            end
        end
        
        %%
        % T = cell2table(Sub);
        clear T
        T = table(Name', MRN', DOS', DOB', Age', spont', DFNM', PN');
        T.Properties.VariableNames{'Var1'} = 'Name';
        T.Properties.VariableNames{'Var2'} = 'MRN';
        T.Properties.VariableNames{'Var3'} = 'DOS';
        T.Properties.VariableNames{'Var4'} = 'DOB';
        T.Properties.VariableNames{'Var5'} = 'Age';
        T.Properties.VariableNames{'Var6'} = 'spont';
        T.Properties.VariableNames{'Var7'} = 'DFNM';
        T.Properties.VariableNames{'Var8'} = 'PN';
        
        %%
        [~,~,c]  = unique(T.Name);
        [~,b] = unique(sort(reshape(c,size(T.Name)),2),'rows');
        T1 = T;
        T = T1(b,:);
        
        %%
        filename = ['patientdata',yr,'.xlsx'];
        writetable(T,filename)
        filename = ['patientdata',yr,'.csv'];
        writetable(T,filename)
        
end

