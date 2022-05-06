function LI = vy_laterality(data, idx)


% for i=1:size(data.value,1)
%     tmp =  abs(data.value(i,:));
%     tmp(tmp < thre.*max(tmp))=0;
%     data_thre(i,:) = tmp;
% end
% data.value = data_thre;

% idx = [4:2:20,24:2:26,80:2:90];
% idx = [2:2:90];
% idx = [12:2:20,80:2:88];

% idx = [12:2:20,80:2:88];
% idx(6) = []; idx(5) = [];
% idx = [idx-1, idx];


% idx1 = [];
% idx1.central = [2,20]; 
% idx1.frontal = 4:2:18; 
% idx1.subcor = [22:2:48,78]; 
% idx1.Occ = 50:2:56; 
% idx1.pari= 58:2:70; 
% idx1.temp = 80:2:90;
% 
% % idx = [idx_central,idx_frontal,idx_subcor,idx_Occ,idx_pari,idx_temp];
% idx2         = [idx1.frontal,idx1.temp,idx1.pari];
% idx2         = [idx1.frontal];
% idx2         = [idx1.frontal, idx1.temp];
% idx = idx2;

%% added, 210621
% idx1 = [];
% idx1.central = [2,20];
% idx1.frontal = 4:2:18;
% idx1.subcor = [22:2:48,78];
% idx1.Occ = 50:2:56;
% idx1.pari= 58:2:70;
% idx1.temp = 82:2:90; %idx1.temp = 80:2:90;
% 
% idx = [idx1.central,idx1.frontal,idx1.pari,idx1.temp];
% % idx = [idx1.temp];
% idx = [idx1.frontal,idx1.pari,idx1.temp];

%%

clear rightFT
for i=1:length(idx)
    rightFT{i,:} = data.label{idx(i)};
end
disp('right hemisphre ROIs:')
disp(rightFT)
% m_right  = mean(data.value(:,idx),2);
m_right  = mean(data.value(idx));


% Left FT lobe
idx = idx-1;
% idx = [1:2:90];

clear leftFT
for i=1:length(idx)
    leftFT{i,:} = data.label{idx(i)};
end
disp('left hemisphre ROIs:')
disp(leftFT)

% m_left  = mean(data.value(:,idx),2);
m_left  = mean(data.value(idx));


LI = (m_left - m_right)./ (m_left + m_right);
disp('Laterality index:')
disp(LI)




