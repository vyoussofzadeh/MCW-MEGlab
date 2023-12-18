clc;clear;close all;

% @author: Vahab Youssof zadeh

global Ae Ai Ae2 Ai2 a2 a1 b2 b1 DD ae ai c1 c2 c3 c4 c5 c6 c7 c8 c9 th P S a Neq A sel

% figure,imshow(imread('Picture1.jpg'));title('Extended Jansen Model 1995');

%% Neural mass params
Ae =  7.7094; % Average synaptic gain constants for excitatory populations
Ai = 35.5523;   % Average synaptic gain constants for inhibitory populations
ae = 100.8710;   % Average synaptic time constants for excitatory populations
ai = 50.9063;   % Average synaptic time constants for inhibitory populations

% %% Neural mass params
% Ae = 8.92; % Average synaptic gain constants for excitatory populations
% Ai = 39.22;   % Average synaptic gain constants for inhibitory populations
% ae = 107.32;  % Average synaptic time constants for excitatory populations
% ai = 52.85;   % Average synaptic time constants for inhibitory populations

%% Prepare simulation (timn scale)
h = 0.005;
Nwin=6;
Nt=round(2.56.*Nwin./h);
t=(1:Nt).*h;

DD = 0;

%% connection stregths
c=128;
c1 = 1*c;    c2 = 1*c;
c3 = 0.5*c; c4 = 0.5*c;

c=128;
c1 = 1*c;    c2 = 0.8*c;
c3 = 0.25*c; c4 = 0.25*c;

c5 = 0.01*c;
c6 = 0.01*c;
c7 = 0.01*c;

% c1 = 133.5201;
% c2 = 134.6967;
% c3 = 69.0351;
% c4 = 68.5286;
% c5 = 4.5135;
% c6 = 0.0481;
% c7 = -0.5387;

c8 = 1;
c9 = 1;
th = 1;

a1 = 55;
a2 = 140;
b1 = 27.5;
b2 = 55;

Ae2 = 1.65;
Ai2 = 32;

A = [32 16 4]; % extrinsic rates (forward, backward, lateral)

disp('without feedback                     : 1');
disp('feedback in In-In                    : 2');
disp('feedback in Ex-py                    : 3');
disp('feedback in In-In& Ex-py             : 4');
disp('feedback in In-In& Ex-py& Ex-In ver1 : 5');
disp('feedback in In-In& Ex-py& Ex-In ver2 : 6');

sel = input('please select model?');
P = input('Enter stimulus?');

if sel==1
    Neq=7;    % number of ode equations in defined model
    r1 = 2;
    r2 = 1;
    m1 = 1; m2 = 2; % EEG = v1-v2 = v7
    M = 7; % EEG = v7
elseif sel==2
    Neq = 11; % number of ode equations in defined model
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; % EEG = v2-v3
    M = 11; % EEG = v7
elseif sel==3
    Neq = 11; % number of ode equations in defined model
    r1 = 4;
    r2 = 2;
    m1 = 2; m2 = 3; % EEG = v2-v3
    M = 11; % EEG = v7
elseif sel==4
    Neq = 13; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; % EEG = v2-v3
    M = 13; % EEG = v7
elseif sel==5
    Neq = 15; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; m3 =7; % EEG = v2-v3-v7
    M = 15; % EEG = v7
elseif sel==6
    Neq = 15; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; m3 =7; % EEG = v2-v3-v7
    M = 15; % EEG = v7
elseif sel==9
    Neq = 22; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; m3 =7; % EEG = v2-v3-v7
    M = 15; % EEG = v7
end
S = @(x) 1./(1+exp(-r1*(x-r2)))-1/(1+exp(r1*r2));


%% Sigmiodal function -------------------------

% idx = ['f1';'f2';'f3']
t1=(-3:0.1:7);
figure,
plot(t1,S(t1)./max(S(t1)),'k','LineWidth',2);
% title('Sigmoidal function');
set(gca, 'Box', 'off');
% axis([-4 6 -0.2 1.2]);

hold on
a = 512e-3;
plot(t1,S(t1-a)./max(S(t1-a)),'r','LineWidth',2);
idx = ['Without adoptation effect';...
    'With adoptation effect   ';...
    ];
legend(idx)

%% External input
N = length(t);
U = 150 + 5e-3*randn(N,1);

ns = N;
dt = diff(t(1:2));
Ask = input('Add noise? yes(1),No(0), or Bump (3): ');

if Ask==1
    U(:,:) = P.*rand(N,1);
    
elseif Ask==0
    U(:,:) = P.*ones(N,1);
    
elseif Ask==3
    onset = 0.0724;
    amp = 6.0784;
    U = amp*exp(-(t - onset).^2/(45*dt^2));
    
end

ask2 = input('Euler or ODE45? yes(1),No(0)');
if ask2==1
    %% Numerical solution
    % initial values of input and output
    x = zeros(Neq,1);
    y = zeros(Neq,1);
    for i = 2:length(U) %construct a 'for' loop
        P = U(i-1);
        tmp1 = model_Jansen(t,x);
        y(:,i) = y(:,i-1)+ (dt.*tmp1);
        x = y(:,i);
    end
else
    x_ini = zeros(Neq,1);
    
    AbsTol = 1e-9.*ones(Neq,1);
    %     options = odeset('RelTol',1e-4,'InitialStep',0.1,'MaxStep',0.005);
    options = odeset('RelTol',1e-6,'AbsTol',AbsTol);
    [t,y]=ode45(@model_Jansen13,t,x_ini,options);
    y = y';
end

%% burn-in effect
L = 100;
y1=[];
y1(1,:) = y(m1,L:end);
y1(2,:) = y(m2,L:end);

%% EEG extraction
EEG = y(M,L:end);
idx = ['u=',num2str(P),'{,\it C_{5}=}',num2str(c5)];
% idx = ['u=',num2str(P)];


%% time and freq visualization
figure
subplot 311
plot(t(L:1000),EEG(L:1000),'k','LineWidth',2),set(gca,'color','none'),set(gca, 'Box', 'off');
title(idx,'FontSize', 20);xlabel('Time (sec)');ylabel('Amp');
fs = 1/h;
Y = fft(EEG);
ff=linspace(0,fs./2,length(Y)./2);
Pyy =(abs(Y));Pyy(1)=0;
Pyy = detrend(Pyy);
subplot 312
plot(ff,Pyy(1:floor(length(Y)./2)),'k','LineWidth',2);set(gca,'color','none'),legend boxoff,set(gca, 'Box', 'off');

[SS1,SS2] = max(Pyy(1:floor(length(Y)./2)));
Freq = ff(SS2)

xlabel('Frequency (Hz)');ylabel('Amp');
subplot 313
% plot3(y(2,L:end),y(1,L:end),y(3,L:end),'k','LineWidth',2);grid off;box off,set(gca,'color','none'),
plot3(y(3,L:end),y(1,L:end),y(2,L:end),'k','LineWidth',2);grid off;box off,set(gca,'color','none'),
xlabel('{\it V_{Inhib}} '); ylabel('{\it V_{Stellate}} '); zlabel('{\it  V_{Pyramidal}} ');

figure
plot(t(L:1000),EEG(L:1000),'k','LineWidth',2),set(gca,'color','none'),set(gca, 'Box', 'off');
title(idx,'FontSize', 20);xlabel('Time (sec)');ylabel('Amp');
