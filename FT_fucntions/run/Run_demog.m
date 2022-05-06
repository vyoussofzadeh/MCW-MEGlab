

load('/data/MEG/Clinical/group_dics/taskcomp/baddata.mat');

T = readtable('master_session.xlsx','sheet', 'patientdata');

%%
for i=1:length(T.Name)
    name1 = strtrim(T.Name{i});
    Index = strfind(name1, ', ');
    if isempty(Index), Index = strfind(name1, ','); end
    if isempty(Index), Index = strfind(name1, ' '); end
    if isempty(Index), Index = strfind(name1, '_'); end
    tmp = strtrim(name1(Index+1:end));
    Index1 = strfind(tmp, ' ');
    if ~isempty(Index1)
        tmp = tmp(1:Index1-1);
        name2{i} = [lower(strtrim(name1(1:Index-1))),'_', lower(tmp)];
    else
        name2{i} = [lower(strtrim(name1(1:Index-1))),'_', lower(strtrim(name1(Index+1:end)))];
    end
end

[C,ia,ib] = intersect(C1, name2,'stable');

age_val = (cell2mat(T.Age(ib)));

for i=1:length(age_val)
    age(i) = str2num(age_val(i,1:2));
end

AGE1 = [num2str(mean(age)),' +- ', num2str(std(age))];
AGE2 = [num2str(min(age)),' to ', num2str(max(age))];
disp(['Age mean-std: ',AGE1])
disp(['Age range: ',AGE2])
disp(age')

disp(C'),
mrn = T.MRN(ib);
disp(mrn)

Leg_sex = ['F'
    'M'
    'F'
    'F'
    'F'
    'M'
    'M'
    'M'
    'M'
    'F'
    'F'
    'M'
    'F'
    'M'
    'M'
    'F'
    'M'
    'M'
    'M'
    'F'
    'F'
    'M'
    'M'
    'F'
    'M'
    'F'
    'F'
    'F'
    'M']; % Derived from EPIC

fmale = length(find(Leg_sex == 'F'));
male = length(find(Leg_sex == 'M'));



