function LI = vy_laterality2(data, idx)


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

k = 1;
clear rightFT
for i=2:2:length(idx)
    rightFT{k,:} = data.label{idx(i)}; 
    r_val(:,k) = data.value(:,idx(i));
    k=1+k;
end
disp('right hemisphre ROIs:')
disp(rightFT)
m_right  = mean(r_val,2);

% Left FT lobe
% idx = idx-1;
% idx = [1:2:90];

k=1;
clear leftFT
for i=1:2:length(idx)
    leftFT{k,:} = data.label{idx(i)};
    l_val(:,k) = data.value(:,idx(i));
    k=1+k;
end
disp('left hemisphre ROIs:')
disp(leftFT)

m_left  = mean(l_val,2);

LI = (m_left - m_right)./ (m_left + m_right);
disp('Laterality index:')
disp(LI)




