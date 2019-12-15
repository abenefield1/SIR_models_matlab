%Initial Conditions:
S_0 = 1000; % Susceptible
I_1_0=1; % Infected
I_2_0=0;
I_n_0=0;
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

a=1;

%[t, class]=ode45(@(t, class) simpMut(t, class, N, beta, nu, b, D, mu, a), T1,[S_0 I_1_0 I_2_0 R_0 I_n_0]);
[t, class]=ode45(@(t, class) simpMut(t, class, N, beta, nu, b, D, mu, a), totalT,[S_0 I_1_0 I_2_0 R_0 I_n_0]);
S=class(:,1);
I_1=class(:,2);
I_2=class(:,3);
I_n=class(:,4);
R=class(:,5);

green=[0.4660, 0.6740, 0.1880];
darkred=[0.6350, 0.0780, 0.1840];	
orange=[0.8500, 0.3250, 0.0980];
yellow=[0.9290, 0.6940, 0.1250];

figure(1);
p1=plot(t,S,'color', green,'LineWidth',2); hold on
p2=plot(t,I_1,'r','LineWidth',2); hold on
p3=plot(t,I_2,'color', darkred,'LineWidth',2); hold on
p4=plot(t,I_n,'color', yellow,'LineWidth',2); hold on
p5=plot(t,R, 'b', 'LineWidth',2); hold on
legend("Susceptible","I_0","I_1","I_2", "Recovered",'fontsize',14);
xlabel("Time",'fontsize',18);
ylabel("Incidence",'fontsize',18);

figure(2);
p1=plot(t,S,'color', green,'LineWidth',2); hold on
p2=plot(t,I_1,'r','LineWidth',2); hold on
p3=plot(t,I_2,'color', darkred,'LineWidth',2); hold on
p4=plot(t,I_n,'color', yellow,'LineWidth',2); hold on
p6=plot(t,I_1+I_2+I_n, 'k', 'LineWidth',2); hold on
p5=plot(t,R, 'b', 'LineWidth',2); hold on
legend("Susceptible","I_0","I_1","I_2", "I_0 + I_1 + I_2", "Recovered",'fontsize',14);
xlabel("Time",'fontsize',18);
ylabel("Incidence",'fontsize',18);


clear all