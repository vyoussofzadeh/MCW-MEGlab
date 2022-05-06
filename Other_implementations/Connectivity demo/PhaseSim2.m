clc;                     % clear screen
clear;                   % clear all the varibles
close all;               % Close all the open figures
%% Time specifications:
Fs = 500;                      % samples per second
dt = 1/Fs;                     % seconds per sample
StopTime = 1;                  % seconds
t = (0:dt:StopTime-dt)';
N = size(t,1);
%% Sine wave:
Fc = 8;                       % hertz
x = cos(2*pi*Fc*t);
%% Fourier Transform:
X = fftshift(fft(x));
% %% Frequency specifications:
% dF = Fs/N;                      % hertz
% f = -Fs/2:dF:Fs/2-dF;           % hertz
%% Fourier Transform:
N = size(t,1);
dt = t(2) - t(1);
fs = 1/dt;                     % FS
% Frequency specifications:
df = fs/N;                      % hertz
f = -fs/2:df:fs/2-df;           % hertz
% Plot the spectrum:
figure;
plot(f,abs(X)/N);
xlabel('Frequency (in hertz)');
title('Magnitude Response');
hold on
plot(f,unwrap(angle(X))/N);

%%
figure;
plot(t,x);
hold on

yh = hilbert(x);
amp = abs(yh);
sigphase = angle(yh);
plot(t,sigphase,'r');
plot(t,amp,'g');

figure, plot(real(yh),imag(yh))
title('phase')


nt = length(t);  %length of t
n=2;                     % Number of cycles to be plotted
freq = 50;               % Frequency in Hz
T=1/freq;                % Time Period in aec
tt = 0:T/1200:n*T;                                           % Time for five cycle of selected frequency
w = 2*pi*freq;           % angluar velocity (rad/s)
circle = (cos(w*tt) + i*sin(w*tt));

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

for i=2:nt
    figure(2)
    subplot 211
    axis off, axis equal
    set(h1,'XData',[0 real(yh(i))]);
    set(h1,'YData',[0 imag(yh(i))]);
    delete(h0);
    h0 = plot(real(yh(i)),imag(yh(i)),'ro','MarkerFaceColor','r');
    figure(2)
    subplot 212
    delete(h3);
    h3 = plot(t(i),x(i),'ro','MarkerFaceColor','r');
    pause(0.1)
end