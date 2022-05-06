clc;                     % clear screen
clear;                   % clear all the varibles
close all;               % Close all the open figures

t = 0:1/2000:2-1/2000;
x = chirp(t-2,4,1/2,6,'quadratic',100,'convex').*exp(-4*(t-1).^2);
plot(t,x)


[up,lo] = envelope(x);
hold on
plot(t,up,t,lo,'linewidth',1.5)
legend('x','up','lo')
hold off
%% Fourier Transform:
X = fftshift(fft(x));
N = size(t,2);
dt = t(2) - t(1);
fs = 1/dt;                     % FS
% Frequency specifications:
df = fs/N;                      % hertz
f = -fs/2:df:fs/2-df;           % hertz
figure;
plot(f,abs(X)/N);
xlabel('Frequency (in hertz)');
title('Magnitude Response');
hold on
plot(f,unwrap(angle(X))/N);

% H = ifft(X);
% y = sqrt(H.^2 + real(hilbert(H)).^2);
% figure;
% plot(t,y);
% hold on 
% plot(t,x,'r');
%% Phase representation

yh = hilbert(x);
amp = abs(yh);
sigphase = angle(yh);
% amplitude=abs(yh);
% phase=zscore(unwrap(angle(yh)));
figure,
plot(t,x);
hold on
plot(t,sigphase,'r');
plot(t,amp,'g');
figure, plot(real(yh)), hold on, plot(imag(yh),'r'), plot(abs(yh),'g'),  plot(x,'c')
figure, plot(real(yh),imag(yh))
%%
nt = length(t);  %length of t
n=2;                     % Number of cycles to be plotted
freq = 50;               % Frequency in Hz
T=1/freq;                % Time Period in aec
tt = 0:T/1200:n*T;                                           % Time for five cycle of selected frequency
w = 2*pi*freq;           % angluar velocity (rad/s)
circle = (cos(w*tt) + i*sin(w*tt));

r = real(yh);
I = imag(yh);
% Plot the reference circle
figure(2);
subplot 211
% plot(circle,'k','LineWidth',2.0);
plot(real(yh),imag(yh),'color',[0.5,0.5,0.5],'LineWidth',0.5)
hold on;
h0 = plot(real(yh(1)),imag(yh(1)),'ro','MarkerFaceColor','r');
h1 = line('XData',[0 real(yh(1))], ...
    'YData',[0 imag(yh(1))], ...
    'Color','k','LineStyle','-');
figure(2)
subplot 212
plot(t,x,'k','LineWidth',2.0);
hold on
h3 = plot(t(1),x(1),'ro','MarkerFaceColor','r'); axis off;
plot(t,sigphase,'r');

for i=2:50:nt
    figure(2)
    subplot 211
    axis off, axis equal
    set(h1,'XData',[0 r(i)]);
    set(h1,'YData',[0 I(i)]);
    delete(h0);
    h0 = plot(r(i),I(i),'ro','MarkerFaceColor','r');
    figure(2)
    subplot 212
    delete(h3);
    h3 = plot(t(i),x(i),'ro','MarkerFaceColor','r');
%     title(['i = ',num2str(i)]);
    pause(0.3)
end
