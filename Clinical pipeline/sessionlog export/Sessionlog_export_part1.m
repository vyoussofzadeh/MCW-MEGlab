clear; clc, close('all'); warning off

%%
%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/sessionlog';

%%
% if exist(outdir, 'file') == 0, mkdir(outdir); end

%%
addpath(genpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions'));

%%
% yr = '19';
disp('1: 2019')
disp('2: 2018');
disp('3: 2017');
disp('4: 2020');
disp('5: 2016'); % not working ... (missing data - inconsistant data)
% disp('4: 2016')
% disp('5: 2015');
Y = input('Year data were acquired: ');
switch Y
    case 1
        yr = '19';
    case 2
        yr = '18';
    case 3
        yr = '17';
    case 4
        yr = '20';
    case 5
        yr = '16';
end

savefile = fullfile(outdir,[yr,'_sessionlog.mat']);
if ~exist(savefile,'file')
    d = rdir([indir,['/**/', yr, '*/sessionLog.txt']]);
    save(savefile,'d');
else
    load(savefile)
end

%%
k = 1;
clear Sub MRN DOB DOS Age DFNM PN spont Name

switch yr
    case '16'
        for i = 14:length(d)
            [pathstr, name] = fileparts(d(i).name);
            filename = d(i).name;
            
            
            
            [a,~] = fileparts(filename);
            cd(a),
            hs = [];
            if exist('help_sessionlog.mat','file')
                load('help_sessionlog.mat'),
                disp(hs)
            end
            
            fid = fopen(filename,'rt');
            tmp = textscan(fid,'%s','Delimiter','\n');
            disp(filename)
            
            if ~isempty(tmp{1,1})
                
                tmp2 = tmp{1,1};
                disp(tmp2)
                
                Name{k} = tmp2{1};
                if isempty(Name{k})
                    N = input('enter name of patn (lastn firstn):');
                    Name{k} = N;
                    hs.name = N;
                end
                Sub{k}.name = Name{k};
                
                bad_id = strfind(Name{k}, 'cta');
                
                if isempty(bad_id)
                    
                    idx = find(contains(tmp2,'MRN')==1);
                    if isempty(idx), idx = find(contains(tmp2,'MNR')==1); end
                    if isempty(idx), idx = find(contains(tmp2,'MR')==1); end
                    
                    if isempty(idx)
                        mrn = input('enter mrn:');
                        Sub{k}.MRN = mrn;
                        hs.mrn = mrn;
                    else
                        Index = strfind(tmp2{idx}, ':');
                        if isempty(Index), Index = strfind(tmp2{idx}, '='); end
                        MRN{k} = tmp2{idx}(Index+1:end);
                        Sub{k}.MRN = MRN{k};
                    end
                    
                    idx = find(contains(tmp2,'DOS')==1);
                    if isempty(idx), idx = find(contains(tmp2,'DOs')==1); end
                    if isempty(idx)
                        dos = input('enter DOS:');
                        Sub{k}.DOS = dos;
                        hs.dos = dos;
                    else
                        Index = strfind(tmp2{idx}, ':');
                        if isempty(Index), Index = strfind(tmp2{idx}, '='); end
                        DOS{k} = tmp2{idx}(Index+1:end);
                        Sub{k}.DOS = DOS{k};
                    end
                    DOS_d = datenum(Sub{k}.DOS,'mm/DD/YYYY');  % note year,month,day
                    
                    idx = find(contains(tmp2,'DOB')==1);
                    if isempty(idx)
                        dob = input('enter DOB:');
                        Sub{k}.DOB = dob;
                        hs.dob = dob;
                    else
                        Index = strfind(tmp2{idx(1)}, ':');
                        if isempty(Index), Index = strfind(tmp2{idx(1)}, '='); end
                        DOB{k} = tmp2{idx(1)}(Index+1:end);
                        Sub{k}.DOB = DOB{k};
                    end
                    DOB_d = datenum(Sub{k}.DOB,'mm/DD/YYYY');  % note year,month,day
                    
                    Age{k} = datestr(DOS_d-DOB_d,'YY.mm.dd');
                    Sub{k}.Age = Age{k};
                    
                    Index = strfind(tmp2, 'spont');
                    spont{k} = double(~isempty(cell2mat(Index)));
                    Sub{k}.spont = spont{k};
                    
                    Index = strfind(tmp2, 'DFNAM');
                    DFNM{k} = double(~isempty(cell2mat(Index)));
                    Sub{k}.DFNM = DFNM{k};
                    
                    Index = strfind(tmp2, 'PN');
                    PN{k} = double(~isempty(cell2mat(Index)));
                    Sub{k}.PN = PN{k};
                    
                    save('help_sessionlog.mat','hs')
                    k = k+1;
                end
            end
        end
        
        
    otherwise
        for i = 1:length(d)
            [pathstr, name] = fileparts(d(i).name);
            filename = d(i).name;
            
            fid = fopen(filename,'rt');
            tmp = textscan(fid,'%s','Delimiter','\n');
            disp(filename)
            
            if ~isempty(tmp{1,1})
                
                tmp2 = tmp{1,1};
                
                Name{k} = tmp2{1};
                Sub{k}.name = Name{k};
                
                bad_id = strfind(Name{k}, 'cta');
                
                if isempty(bad_id)
                    
                    idx = find(contains(tmp2,'MRN')==1);
                    if isempty(idx), idx = find(contains(tmp2,'MNR')==1); end
                    if isempty(idx), idx = find(contains(tmp2,'MR')==1); end
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
                    Index = strfind(tmp2{idx(1)}, ':');
                    if isempty(Index), Index = strfind(tmp2{idx(1)}, '='); end
                    DOB{k} = tmp2{idx(1)}(Index+1:end);
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
filename = fullfile(outdir,['patientdata',yr,'.xlsx']);
writetable(T,filename)
filename = fullfile(outdir,['patientdata',yr,'.csv']);
disp(filename)
writetable(T,filename)



