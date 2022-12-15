% AGC, from Siljac book
% Marcello Farina, 19/12/2019
clc
close all
clear all
A1=[-0.2 0 0 0 0 0 0 -4
    4.75 -5 0 0 0 0 0 0
    0 0.1667 -0.1667 0 0 0 0 0
    0 0 2 -2 0 0 0 0
    0 -0.08 -0.0747 -0.112 -3.994 10 -0.928 -9.1011
    0 0 0 0 0.2 -0.5 0 0
    0 0 0 0 1.3194 0 -1.3889 -0.2778
    0 0.01 0.0093 0.014 -0.0632 0 0.1160 -0.1124];
A2=[-0.2 0 0 0 0 0 0 -4
    4.75 -5 0 0 0 0 0 0
    0 0.1667 -0.1667 0 0 0 0 0
    0 0 2 -2 0 0 0 0
    0 -0.1 -0.0933 -0.14 -4.096 10 -0.7442 -9.1079
    0 0 0 0 0.2 -0.5 0 0
    0 0 0 0 1.3194 0 -1.3889 -0.2778
    0 0.0125 0.0117 0.0175 -0.0506 0 0.0928 -0.1115];
at=[0 0 0 0 0.6667 0 0 -0.0833]';
d=[0 0 0 0 0 0 0 10]';
m=[0 0 0 0 0 0 0 22.2144]';
b1=[1.6 0 0 0 6 0 0 0]';
b2=[2 0 0 0 5 0 0 0]';

Atot =[A1 zeros(8,1) at zeros(8,9)
    d' 0 1 zeros(1,9)
    m' 0 0 0 -m'
    zeros(1,9) -1 0 d'
    zeros(8,9) -at zeros(8,1) A2];
         
Bdec{1}=[b1' zeros(1,11)]';
Bdec{2}=[zeros(1,11) b2']';
Btot = [Bdec{1} Bdec{2}];
    
Ctot=eye(19);
Cdec{1}=Ctot(1:10,:);
Cdec{2}=Ctot(11:end,:);
D = 0;
sys = ss(Atot,Btot,Ctot,D);


% figure
% step(sys(7:12,:))
% eig(Atot)
% 
% figure
% impulse(sys(7:12,:))
% % Riga 9 input 1 e riga 11 input 2
Ts = 0.1;
 sys_d = c2d(sys,Ts);

Bw = eye(19);
% 
% figure
% step(sys_d(7:12,:))
 [F,G,H,I] = ssdata(sys_d);

Gdec{1} = G(:,1);
Gdec{2} = G(:,2);

Hdec{1} = H(1:10,:);
Hdec{2}= H(11:end,:);


%% 1) Compute the eigenvalues
close all
lambda_cont = eig(Atot)
subplot(2,1,1)
impulse(sys(9,1)), grid;
title("9th row impulse response")
subplot(2,1,2)
impulse(sys(11,2)), grid;
title("11th row impulse response")

lambda_discr = eig(F)

%il sistema è semplicemente stabile ma non asintoticamente stabile nè
%instabile


%% CONTINUOS and DISCRETE TIME FIXED MODES

N=2; %subsystems
rounding_n = 10; %Fondamentale il rounding

%Centralized
ContStruc_c = [1 1;
                1 1];

[Difm_c_cont]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc_c,rounding_n)
[Difm_c_discr]=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_c,rounding_n);

%Decentralized
ContStruc_de = [1 0;
                0 1];

[Difm_de_cont]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc_de,rounding_n)
[Difm_de_discr]=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_de,rounding_n)

%Different Distributed
ContStruc_di1 = [1 1;
                0 1];

[Difm_di1_cont]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc_di1,rounding_n)
[Difm_di1_discr]=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_di1,rounding_n)

%Different Distributed
ContStruc_di2 = [1 0;
                1 1];

[Difm_di2_cont]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc_di2,rounding_n)
[Difm_di2_discr]=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_di2,rounding_n)

%% SIMULATION DATA
Tfinal=50;
T=[0:0.01:Tfinal];
%x0=[randn(19,1)];
x0 = [0.226677271849582;
0.824237059236596;
-0.616205044427732;
-0.845098448674646;
1.03458755998388;
0.810418118641796;
-0.766406942550575;
0.305396532258002;
-0.646292039547889;
-0.312411053862938,
-0.946728360430056;
0.121495551068268;
0.104374052679890;
0.448354953825058;
-0.514083120076390;
0.195948822329623;
0.771023707694623;
0.0741922629207770;
0.176145961870921];

%u_n=[randn(19,1)];
u_n = [-0.00625115632404493;
0.0789740600856282;
1.43266864767143;
-0.862688401742682;
0.442720184396889;
0.829725439521861;
1.11117461425725;
-1.08826254292861;
0.772688021097315;
-0.123388524620109;
-0.379705684663708;
0.669221526907101;
-0.108273901938890;
0.665466925189865;
-0.0273814139843615;
-1.53538016700313;
-0.766501678980656;
-1.91175040169875;
0.0592554202694067];

sys2=ss(Atot,eye(19), eye(19),0);
sys2d=c2d(sys2,Ts);
[Nw,Gw,Tw]=ssdata(sys2d);

%% CENTRALIZED
% -Continuos-
% LMI continuos for performance
[K1c_perf,rho1c_perf,feas1c_perf]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_c,1,0.345); %0.3 limite max
lambda_K1c_perf = eig(Atot+Btot*K1c_perf);
% LMI continuos region angle
[K1c_ang,rho1c_ang,feas1c_ang]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_c,2,20); %0.3 limite max
labda_K1c_ang = eig(Atot+Btot*K1c_ang);
% LMI limitation of control effort
[K1c_eff,rho1c_eff,feas1c_eff]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_c,3,0.08);
labda_K1c_eff = eig(Atot+Btot*K1c_eff);
% LMI minimization of H2 norm
[K1c_h2,rho1c_h2,feas1c_h2]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_c,4,999999);
labda_K1c_h2 = eig(Atot+Btot*K1c_h2);
sys_c_n1 = ss(Atot+Btot*K1c_h2,Bw,eye(19),D)
[x_1c_h2]=lsim(sys_c_n1,rand(19,length(T))*0.1,T,x0)
x_1c_h2 = x_1c_h2'

k=0;
for t=T
    k=k+1;
    x_1c_perf(:,k)=expm((Atot+Btot*K1c_perf)*t)*x0;
    x_1c_ang(:,k)=expm((Atot+Btot*K1c_ang)*t)*x0;
    x_1c_eff(:,k)=expm((Atot+Btot*K1c_eff)*t)*x0;
end
u_1c_eff = K1c_eff*x_1c_eff;

figure
subplot(4,1,1)
plot(T,x_1c_perf(10,:))
title('Centralized, continuos for performance')
subplot(4,1,2)
plot(T,x_1c_ang(10,:))
title('Centralized, continuos angle')
subplot(4,1,3)
plot(T,x_1c_eff(10,:),T,u_1c_eff)
title('Centralized, continuos effort limitation')
subplot(4,1,4)
plot(T,x_1c_h2(10,:))
title('Centralized, continuos H2')

% -Discrete- 
% LMI discrete for performance 
[K1d_perf,rho1d_perf,feas1d_perf]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_c,1,0.98);
lambda_K1d_perf = eig(F+G*K1d_perf);
% LMI discre time circle 
[K1d_circ,rho1c_circ,feas1c_circ]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_c,2,[0.2 -0.78]);
lambda_K1d_circ = eig(F+G*K1d_circ);
% LMI limitation of control effort
[K1d_eff,rho1d_eff,feas1d_eff]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_c,3,0.982);
lambda_K1d_eff = eig(F+G*K1d_eff);
% LMI minimization of H2 norm
[K1d_h2,rho1d_h2,feas1d_h2]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_c,3,Gw);
lambda_K1d_h2 = eig(F+G*K1d_h2);
sys_d_n1 = ss(F+G*K1d_h2,Gw,eye(19),D,Ts)
[x_1d_h2]=lsim(sys_d_n1,rand(19,Tfinal/Ts)*0.1,[Ts:Ts:Tfinal],x0)
x_1d_h2 = x_1d_h2'


for k=1:Tfinal/Ts
    x_1d_perf(:,k)=((F+G*K1d_perf)^k)*x0;
    x_1d_circ(:,k)=((F+G*K1d_circ)^k)*x0;
    x_1d_eff(:,k)=((F+G*K1d_eff)^k)*x0;
end
u_1d_eff = K1d_eff*x_1d_eff;

figure
subplot(4,1,1)
plot([Ts:Ts:Tfinal],x_1d_perf(10,:))
title('Centralized, discrete for performance')
subplot(4,1,2)
plot([Ts:Ts:Tfinal],x_1d_circ(10,:))
title('Centralized, discrete circle')
subplot(4,1,3)
plot([Ts:Ts:Tfinal],x_1d_eff(10,:),[Ts:Ts:Tfinal],u_1d_eff)
title('Centralized, discrete effort limitation')
subplot(4,1,4)
plot([Ts:Ts:Tfinal],x_1d_h2(10,:))
title('Centralized, discrete H2')

%% DECENTRALIZED
% -Continuos-
% LMI continuos for performance
[K2c_perf,rho2c_perf,feas2c_perf]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_de,1,0.086); %0.3 limite max
lambda_K2c_perf = eig(Atot+Btot*K2c_perf);
% LMI continuos region angle
[K2c_ang,rho2c_ang,feas2c_ang]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_de,2,56);
labda_K2c_ang = eig(Atot+Btot*K2c_ang);
% LMI limitation of control effort
[K2c_eff,rho2c_eff,feas2c_eff]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_de,3,0.03)
labda_K2c_eff = eig(Atot+Btot*K2c_eff);
% LMI minimization of H2 norm
[K2c_h2,rho2c_h2,feas2c_h2]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_de,4,999999)
labda_K2c_h2 = eig(Atot+Btot*K2c_h2);
sys_c_n2 = ss(Atot+Btot*K2c_h2,Bw,eye(19),D)
[x_2c_h2]=lsim(sys_c_n2,rand(19,length(T))*0.1,T,x0)
x_2c_h2 = x_2c_h2'

k = 0;
for t=T
    k=k+1;
    x_2c_perf(:,k)=expm((Atot+Btot*K2c_perf)*t)*x0;
    x_2c_ang(:,k)=expm((Atot+Btot*K2c_ang)*t)*x0;
    x_2c_eff(:,k)=expm((Atot+Btot*K2c_eff)*t)*x0;
end

figure
subplot(4,1,1)
plot(T,x_2c_perf(10,:))
title('Decentralized, continuos for performance')
subplot(4,1,2)
plot(T,x_2c_ang(10,:))
title('Decentralized, continuos angle')
subplot(4,1,3)
plot(T,x_2c_eff(10,:))
title('Decentralized, continuos effort limitation')
subplot(4,1,4)
plot(T,x_2c_h2(10,:))
title('Decentralized, continuos H2')

% -Discrete-
% LMI discrete for performance
[K2d_perf,rho2d_perf,feas2d_perf]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_de,1,0.993);
lambda_K2d_perf = eig(F+G*K2d_perf);
% LMI discre time circle 
[K2d_circ,rho2c_circ,feas2c_circ]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_de,2,[0.3 -0.693]);
lambda_K2d_circ = eig(F+G*K2d_circ);
% LMI limitation of control effort
[K2d_eff,rho2d_eff,feas2d_eff]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_de,3,0.993);
lambda_K2d_eff = eig(F+G*K2d_eff);
% LMI minimization of H2 norm
[K2d_h2,rho2d_h2,feas2d_h2]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_de,3,Gw);
lambda_K2d_h2 = eig(F+G*K2d_h2);
sys_d_n2 = ss(F+G*K2d_h2,Gw,eye(19),D,Ts)
[x_2d_h2]=lsim(sys_d_n2,rand(19,Tfinal/Ts)*0.1,[Ts:Ts:Tfinal],x0)
x_2d_h2 = x_2d_h2'

for k=1:Tfinal/Ts
    x_2d_perf(:,k)=((F+G*K2d_perf)^k)*x0;
    x_2d_circ(:,k)=((F+G*K2d_circ)^k)*x0;
    x_2d_eff(:,k)=((F+G*K2d_eff)^k)*x0;
end

figure
subplot(4,1,1)
plot([Ts:Ts:Tfinal],x_2d_perf(10,:))
title('Decentralized, discrete for performance')
subplot(4,1,2)
plot([Ts:Ts:Tfinal],x_2d_circ(10,:))
title('Decentralized, discrete circle')
subplot(4,1,3)
plot([Ts:Ts:Tfinal],x_2d_eff(10,:))
title('Decentralized, discrete effort limitation')
subplot(4,1,4)
plot([Ts:Ts:Tfinal],x_2d_h2(10,:))
title('Decentralized, discrete H2')
%% DISTRIBUTED 1
% -Continuos-
% LMI continuos for performance
[K3c_perf,rho3c_perf,feas3c_perf]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di1,1,0.086);
lambda_K3c_perf = eig(Atot+Btot*K3c_perf);
% LMI continuos region angle
[K3c_ang,rho3c_ang,feas3c_ang]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di1,2,59);
labda_K3c_ang = eig(Atot+Btot*K3c_ang);
% LMI limitation of control effort
[K3c_eff,rho3c_eff,feas3c_eff]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di1,3,0.03)
labda_K3c_eff = eig(Atot+Btot*K3c_eff);
% LMI minimization of H2 norm
[K3c_h2,rho3c_h2,feas3c_h2]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di1,4,999999)
labda_K3c_h2 = eig(Atot+Btot*K3c_h2);
sys_c_n3 = ss(Atot+Btot*K3c_h2,Bw,eye(19),D)
[x_3c_h2]=lsim(sys_c_n3,rand(19,length(T))*0.1,T,x0)
x_3c_h2 = x_3c_h2'

k = 0;
for t=T
    k=k+1;
    x_3c_perf(:,k)=expm((Atot+Btot*K3c_perf)*t)*x0;
    x_3c_ang(:,k)=expm((Atot+Btot*K3c_ang)*t)*x0;
    x_3c_eff(:,k)=expm((Atot+Btot*K3c_eff)*t)*x0;
end
% w = rand(19,1)
% Bw = eye(19)
figure
subplot(4,1,1)
plot(T,x_3c_perf(10,:))
title('Distributed1, continuos for performance')
subplot(4,1,2)
plot(T,x_3c_ang(10,:))
title('Distributed1, continuos angle')
subplot(4,1,3)
plot(T,x_3c_eff(10,:))
title('Distributed1, continuos effort limitation')
subplot(4,1,4)
plot(T,x_3c_h2(10,:))
title('Distributed1, continuos H2')

% -Discrete-
% LMI discrete for performance
[K3d_perf,rho3d_perf,feas3d_perf]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di1,1,0.992);
lambda_K3d_perf = eig(F+G*K3d_perf);
% LMI discre time circle 
[K3d_circ,rho3d_circ,feas3d_circ]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di1,2,[0.3 -0.693]);
lambda_K3c_circ = eig(F+G*K3d_circ);
% LMI limitation of control effort
[K3d_eff,rho3d_eff,feas3d_eff]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di1,3,0.992);
lambda_K3d_eff = eig(F+G*K3d_eff);
% LMI minimization of H2 norm
[K3d_h2,rho3d_h2,feas3d_h2]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di1,3,Gw);
lambda_K3d_h2 = eig(F+G*K3d_h2);
sys_d_n3 = ss(F+G*K3d_h2,Gw,eye(19),D,Ts)
[x_3d_h2]=lsim(sys_d_n3,rand(19,Tfinal/Ts)*0.1,[Ts:Ts:Tfinal],x0)
x_3d_h2 = x_3d_h2'

for k=1:Tfinal/Ts
    x_3d_perf(:,k)=((F+G*K3d_perf)^k)*x0;
    x_3d_circ(:,k)=((F+G*K3d_circ)^k)*x0;
    x_3d_eff(:,k)=((F+G*K3d_eff)^k)*x0;
end

figure
subplot(4,1,1)
plot([Ts:Ts:Tfinal],x_3d_perf(10,:))
title('Distributed1, discrete for performance')
subplot(4,1,2)
plot([Ts:Ts:Tfinal],x_3d_circ(10,:))
title('Distributed1, discrete circle')
subplot(4,1,3)
plot([Ts:Ts:Tfinal],x_3d_eff(10,:))
title('Distributed1, discrete effort limitation')
subplot(4,1,4)
plot([Ts:Ts:Tfinal],x_3d_h2(10,:))
title('Distributed1, discrete H2')
%% DISTRIBUTED 2
% -Continuos-
% LMI continuos for performance
[K4c_perf,rho4c_perf,feas4c_perf]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di2,1,0.087); %0.3 limite max
lambda_K4c_perf = eig(Atot+Btot*K4c_perf);
% LMI continuos region angle
[K4c_ang,rho4c_ang,feas4c_ang]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di2,2,62);
labda_K4c_ang = eig(Atot+Btot*K4c_ang);
% LMI limitation of control effort
[K4c_eff,rho4c_eff,feas4c_eff]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di2,3,0.03)
labda_K4c_eff = eig(Atot+Btot*K4c_eff);
% LMI minimization of H2 norm
[K4c_h2,rho4c_h2,feas4c_h2]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_di2,4,999999)
labda_K4c_h2 = eig(Atot+Btot*K4c_h2);
sys_c_n4 = ss(Atot+Btot*K4c_h2,Bw,eye(19),D)
[x_4c_h2]=lsim(sys_c_n4,rand(19,length(T))*0.1,T,x0)
x_4c_h2 = x_4c_h2'

k = 0;
for t=T
    k=k+1;
    x_4c_perf(:,k)=expm((Atot+Btot*K4c_perf)*t)*x0;
    x_4c_ang(:,k)=expm((Atot+Btot*K4c_ang)*t)*x0;
    x_4c_eff(:,k)=expm((Atot+Btot*K4c_eff)*t)*x0;
end

figure
subplot(4,1,1)
plot(T,x_4c_perf(10,:))
title('Distributed2, continuos for performance')
subplot(4,1,2)
plot(T,x_4c_ang(10,:))
title('Distributed2, continuos angle')
subplot(4,1,3)
plot(T,x_4c_eff(10,:))
title('Distributed2, continuos effort limitation')
subplot(4,1,4)
plot(T,x_4c_h2(10,:))
title('Distributed2, continuos H2')

% -Discrete-
% LMI discrete for performance
[K4d_perf,rho4d_perf,feas4d_perf]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di2,1,0.992);
lambda_K4d_perf = eig(F+G*K4d_perf);
% LMI discre time circle 
[K4d_circ,rho4d_circ,feas4d_circ]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di2,2,[0.3 -0.693]);
lambda_K4c_circ = eig(F+G*K4d_circ);
% LMI limitation of control effort
[K4d_eff,rho4d_eff,feas4d_eff]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di2,3,0.992);
lambda_K4d_eff = eig(F+G*K4d_eff);
% LMI minimization of H2 norm
[K4d_h2,rho4d_h2,feas4d_h2]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_di2,3,Gw);
lambda_K4d_h2 = eig(F+G*K4d_h2);
sys_d_n4 = ss(F+G*K4d_h2,Gw,eye(19),D,Ts)
[x_4d_h2]=lsim(sys_d_n4,rand(19,Tfinal/Ts)*0.1,[Ts:Ts:Tfinal],x0)
x_4d_h2 = x_4d_h2'

for k=1:Tfinal/Ts
    x_4d_perf(:,k)=((F+G*K4d_perf)^k)*x0;
    x_4d_circ(:,k)=((F+G*K4d_circ)^k)*x0;
    x_4d_eff(:,k)=((F+G*K4d_eff)^k)*x0;
end

figure
subplot(4,1,1)
plot([Ts:Ts:Tfinal],x_4d_perf(10,:))
title('Distributed2, discrete for performance')
subplot(4,1,2)
plot([Ts:Ts:Tfinal],x_4d_circ(10,:))
title('Distributed2, discrete circle')
subplot(4,1,3)
plot([Ts:Ts:Tfinal],x_4d_eff(10,:))
title('Distributed2, discrete effort limitation')
subplot(4,1,4)
plot([Ts:Ts:Tfinal],x_4d_h2(10,:))
title('Distributed2, discrete H2')
%% ::::::::COMPARISON PLOT:::::::::::
% Continuos
figure
subplot(4,1,1)
plot(T,x_1c_perf(10,:),T,x_2c_perf(10,:),T,x_3c_perf(10,:),T,x_4c_perf(10,:))
title('Continuos time for performance')
legend('Centr','Decentr','Distr1','Distr2')
subplot(4,1,2)
plot(T,x_1c_ang(10,:),T,x_2c_ang(10,:),T,x_3c_ang(10,:),T,x_4c_ang(10,:))
title('Continuos time with angle')
legend('Centr','Decentr','Distr1','Distr2')
subplot(4,1,3)
plot(T,x_1c_eff(10,:),T,x_2c_eff(10,:),T,x_3c_eff(10,:),T,x_4c_eff(10,:))
title('Continuos time effort limitation')
legend('Centr','Decentr','Distr1','Distr2')
subplot(4,1,4)
plot(T,x_1c_h2(10,:),T,x_2c_h2(10,:),T,x_3c_h2(10,:),T,x_4c_h2(10,:))
title('Continuos time h2')
legend('Centr','Decentr','Distr1','Distr2')
sgtitle('Continuos time comparison')
% Discrete
figure
subplot(4,1,1)
plot([Ts:Ts:Tfinal],x_1d_perf(10,:),[Ts:Ts:Tfinal],x_2d_perf(10,:),[Ts:Ts:Tfinal],x_3d_perf(10,:),[Ts:Ts:Tfinal],x_4d_perf(10,:))
title('Discrete time for performance')
legend('Centr','Decentr','Distr1','Distr2')
subplot(4,1,2)
plot([Ts:Ts:Tfinal],x_1d_circ(10,:),[Ts:Ts:Tfinal],x_2d_circ(10,:),[Ts:Ts:Tfinal],x_3d_circ(10,:),[Ts:Ts:Tfinal],x_4d_circ(10,:))
title('Discrete time with circle')
legend('Centr','Decentr','Distr1','Distr2')
subplot(4,1,3)
plot([Ts:Ts:Tfinal],x_1d_eff(10,:),[Ts:Ts:Tfinal],x_2d_eff(10,:),[Ts:Ts:Tfinal],x_3d_eff(10,:),[Ts:Ts:Tfinal],x_4d_eff(10,:))
title('Discrete time effort limitation')
legend('Centr','Decentr','Distr1','Distr2')
subplot(4,1,4)
plot([Ts:Ts:Tfinal],x_1d_h2(10,:),[Ts:Ts:Tfinal],x_2d_h2(10,:),[Ts:Ts:Tfinal],x_3d_h2(10,:),[Ts:Ts:Tfinal],x_4d_h2(10,:))
title('Discrete time h2')
legend('Centr','Decentr','Distr1','Distr2')
sgtitle('Discrete time comparison')

%% PLOT

sys_test = ss(F,Bw,eye(19),D,Ts)
[x_prova]=lsim(sys_test,rand(19,Tfinal/Ts)*0.1,[Ts:Ts:Tfinal],x0)
x_prova = x_prova'

figure
plot([Ts:Ts:Tfinal],x_1d_h2(10,:),[Ts:Ts:Tfinal],x_prova(10,:))


%% CONTROL EFFORT

figure
subplot(2,1,1)
plot(T,x_1c_eff(10,:),T,u_1c_eff)
legend('X','U_1','U_2')
title('More control effort')
grid on

[K1c_eff,rho1c_eff,feas1c_eff]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc_c,3,0.2);
labda_K1c_eff = eig(Atot+Btot*K1c_eff);

k=0;
for t=T
    k=k+1;
    x_1c_eff(:,k)=expm((Atot+Btot*K1c_eff)*t)*x0;
end
u_1c_eff_2 = K1c_eff*x_1c_eff;

subplot(2,1,2)
plot(T,x_1c_eff(10,:),T,u_1c_eff_2)
legend('X','U_1','U_2')
title('Less control effort')
grid on
sgtitle("Control effort comparison")

figure
plot(T,x_1c_perf(10,:),T,x_2c_perf(10,:),T,x_3c_perf(10,:),T,x_4c_perf(10,:))
legend('Centralized','Decentralized','Distributed 1','Distributed 2')
grid on