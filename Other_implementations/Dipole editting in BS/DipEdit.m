test2.DipoleNames = test2.DipoleNames(1:20);
test2.Dipole = test2.Dipole(1:20);
test2.Comment = 'dipole test';
% A = load('dipoles_Run02_Run03_Run.mat');


savetag = 'dipoles_Run02_Run03_Run_edit';
save(fullfile(savetag),'-struct', 'test2'),
