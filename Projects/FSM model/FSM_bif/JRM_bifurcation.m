clc;clear;close all;

% @author: Vahab Youssof zadeh


global Ae Ai ae ai c1 c2 c3 c4 c5 c6 c7 P S DD r1 r2 A sel Neq

%% Neural mass params
Ae = 3.25; % Average synaptic gain constants for excitatory populations
Ai = 22;   % Average synaptic gain constants for inhibitory populations
ae = 100;  % Average synaptic time constants for excitatory populations
ai = 50;   % Average synaptic time constants for inhibitory populations

% Ae =  8.03; % Average synaptic gain constants for excitatory populations
% Ai = 26.36;   % Average synaptic gain constants for inhibitory populations
% ae = 99;   % Average synaptic time constants for excitatory populations
% ai = 29;   % Average synaptic time constants for inhibitory populations

DD  = 0;
A = [32 16 4]; % extrinsic rates (forward, backward, lateral)

disp('without feedback                : 1');
disp('feedback in In-In               : 2');
disp('feedback in Ex-py               : 3');
disp('feedback in In-In& Ex-py        : 4');
disp('feedback in In-In& Ex-py& Ex-In : 5');
disp('Two Sources                     : 6');

sel = input('please select model?');

if sel==1
    Neq = 7;    % number of ode equations in defined model
    r1 =  2;
    r2 =  1;
    m1 =  1; m2 = 2; % EEG = v1-v2 = v7
    MM =  7; % EEG = v7
elseif sel==2
    Neq = 11; % number of ode equations in defined model
    r1  = 2;
    r2  = 1;
    m1  = 2; m2 = 3; % EEG = v2-v3
    MM  = 11; % EEG = v7
elseif sel==3
    Neq = 11; % number of ode equations in defined model
    r1 = 4;
    r2 = 2;
    m1 = 2; m2 = 3; % EEG = v2-v3
    MM = 11; % EEG = v7
elseif sel==4
    Neq = 13; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; % EEG = v2-v3
    MM = 13; % EEG = v7
elseif sel==5
    Neq = 15; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; m3 =7; % EEG = v2-v3-v7
    MM = 15; % EEG = v7
elseif sel==6
    Neq = 18; % number of ode equations in defined m
    r1 = 2;
    r2 = 1;
    m1 = 2; m2 = 3; m3 =7; % EEG = v2-v3-v7
    MM = 18; % EEG = v7
end
S = @(x) 1./(1+exp(-r1*(x-r2)))-1/(1+exp(r1*r2));

%% Jansen (1995) priors
% c = 135;
% c1 = 1*c;    c2 = 0.9*c;
% c3 = 0.25*c; c4 = 0.25*c;

%% Moran (2007) priors
% c=128;
% c1 = 1*c;    c2 = 0.8*c;
% c3 = 0.5*c; c4 = 0.25*c;

%% FSM priors
c=128;
c1 = 1*c;    c2 = 0.8*c;
c3 = 0.25*c; c4 = 0.25*c;

% c1 = 129;    c2 = 102;
% c3 = 90;     c4 = 61;

%% Prepare simulation
h = 0.005;
Nwin=6;
Nt=round(2.56.*Nwin./h);
t=(1:Nt).*h;
dt = diff(t(1:2));

N = length(t);ns = N;
U = 150 + 5e-3*randn(N,1);

ask2 = input('Euler or ODE45? yes(1),No(0)');

%% Numerical solution
x = ones(Neq,1);
y = ones(Neq,1);

disp('c5 : 1');
disp('c6 : 2');
disp('c7 : 3');
disp('P  : 4');
ask3 = input('Select?');

if ask3==1
    c5_min = 0*c; c5_max = 2*c; c5_step= 0.02*c;   % Input variations set: U = P_range in line 88
    c_range = [c5_min:c5_step:c5_max];
    c6 = 0.0001*c;
    c7 = 0.0001*c;
    P = 100;
elseif ask3==2
    c6_min = 0*c; c6_max = 2*c; c6_step= 0.04*c;   % Input variations set: U = P_range in line 88
    c_range = [c6_min:c6_step:c6_max,c6_max:-c6_step:c6_min];
    c5 = 0.1*c;
    c7 = 0.1*c;
    P = 100;
elseif ask3==3
    c7_min = 0*c; c7_max = 1*c; c7_step= 0.02*c;   % Input variations set: U = P_range in line 88
    c_range = [c7_min:c7_step:c7_max,c7_max:-c7_step:c7_min];
    c5 = 5*c;
    c6 = 5*c;
    P = 100;
elseif ask3==4
    p_min = -100; p_max = 700; p_step= 5;   % Input variations set: U = P_range in line 88
    c_range = [p_min:p_step:p_max];
end


c1 = 128.5201
c2 = 102.6967
c3 = 67.0351
c4 = 65.5286
c5 = 5.46
c6 = 3.11;
c7 = 3.33;

% 3-5)Predefine some matrices
max_num=100;
max_save1=zeros(length(c_range),max_num);
min_save1=zeros(length(c_range),max_num);

disp(' ')
disp('Total number of loops: '),
disp(length(c_range)),
disp('Current loop: ')

nframes = length(c_range);
Frames = moviein(nframes);
L = 100; % burn-in effect

i=1;
for P=c_range
    % for c7=c_range
    if ask2==1
        %% Numerical solution
        for j = 2:length(U) %construct a 'for' loop
            %         P = U(j-1);
            tmp1 = model_Jansen(t,x);
            y(:,j) = y(:,j-1)+ (dt.*tmp1);
            x = y(:,j);
        end
    else
        x_ini = zeros(Neq,1);
        AbsTol = 1e-9.*ones(Neq,1);
        %     options = odeset('RelTol',1e-4,'InitialStep',0.1,'MaxStep',0.005);
        options = odeset('RelTol',1e-6,'AbsTol',AbsTol);
        [t,y]=ode45(@model_Jansen13,t,x_ini,options);
        y = y';
    end
    x=y';
    plot(y(m1,L:end),y(m2,L:end),'k');xlabel('PY');ylabel('IN');title('phase space' )
    Frames(:,i) = getframe(gcf);
    pause(0.01)
    series=x(end-floor(3*length(x)/4):end,:);
    tmp = series(:,MM); % EEG (Y1-Y2)
    [max1, posmax1] = findpeaks(tmp);
    [min1, posmin1] = findpeaks(-tmp);
    m_diff = max([length(max1) length(min1)]);
    if m_diff > max_num
        max_save1=[max_save1 zeros(length(c_range),m_diff-max_num)];
        min_save1=[min_save1 zeros(length(c_range),m_diff-max_num)];
        max_num=size(min_save1,2);
    end
    if min([length(max1) length(min1)])<1
        max_save1(i,1) = tmp(end);
        min_save1(i,1) = tmp(end);
    else
        max_save1(i,1:length(max1)) = max_save1(i,1:length(max1))+max1';
        min_save1(i,1:length(min1)) = min_save1(i,1:length(min1))-min1';
    end
    M.x = x(end,:);
    disp(i)
    i=i+1;
end

z1=[];z2=[];
for i=1: length(c_range)
    z1(i)=nnz(max_save1(i,:));
    z2(i)=nnz(min_save1(i,:));
end

figure,
plot(c_range, z1,'b');
hold on;
plot(c_range, z2,'r');
legend('maxima','Minima');xlabel('Constant {\it u_{t}} ');ylabel('Number');
set(gca,'color','none')

max_save1(max_save1==0)=NaN;
min_save1(min_save1==0)=NaN;

figure,
plot(c_range,min_save1,'*b','MarkerSize',4);
% title('Bifurcation diagram','FontSize', 20);
hold on;
plot(c_range,max_save1,'dr','MarkerSize',4);set(gca, 'Box', 'off');
set(gca,'color','none')


disp('Completed ')

