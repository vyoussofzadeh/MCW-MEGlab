clear
load ('colormap.mat')
figure, colormap(A)
a = axes;
c = colorbar(a);
a.Visible = 'off';