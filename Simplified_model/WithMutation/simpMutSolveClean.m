%Initial Conditions:
S_0 = 1000; % Susceptible
I_1_0=1; % Infected
I_2_0=0;
R_0=0; % Recovered
b=100; % birth rate into susceptible
D=0.1; % death rate (independent of disease)
N=1000;

detTime = 50;
endTime = 150;
T1 = 0:detTime;
T2 = detTime+1:endTime;
totalT=0:endTime;

mu=0.01; % Mutation rate
nu=0.2; % Recovery rate
beta=0.001; % Transmission rate
%% 0:50 - burn-in (no detection rate)
%%
% $R_{0}=\frac{\beta b}{\delta(\delta+\nu+\mu)}$

det=0;

[t, class]=ode45(@(t, class) simpMut(t, class, N, beta, nu, b, D, det, mu), T1,[S_0 I_1_0 I_2_0 R_0]);
S=class(:,1);
I_1=class(:,2);
I_2=class(:,3);
R=class(:,4);

figure(1)
p1=plot(t,S,'g','LineWidth',2); hold on
p2=plot(t,I_1,'r','LineWidth',2); hold on
p3=plot(t,I_2,'m','LineWidth',2); hold on
p4=plot(t,R,'b','LineWidth',2); hold on
legend("S","I_1","I_2","R")


clear all
clc