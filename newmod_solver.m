%Initial Conditions:
S_0 = 1000;
I_t_0= 25;
I_nt_0= 25;
%T_0= 0;
R_0= 0;

vaccTime = 100;
endTime = 150;
T1 = 0:40;
T2 = vaccTime:endTime;


%for t=1:T(150);

% deltaW=0.2;
% d=0.1;
% v=0.01;
% betaW=0.001;

N=1000;
beta= 0.001; % Rate of infection and tested by MDT
%gamma_t= 0.05; % Rate of recovery without treatment from Infected and tested by MDT class. Lower than recovery with treatment.
gamma_nt= 0.15; % Rate of recovery without treatment from Infected and NOT tested by MDT class. Lower than recovery with treatment. Same as gamma_t
%delta_t= 0.3; % Rate of treatment from Infected and tested with MDT class. Higher than treatment rate of not tested class
%delta_nt= 0.15; % Rate of treatment from Infected and NOT tested with MDT class. Lower than treatment rate of tested class.
%alpha= 0.4; % Rate of recovery from treatment. Higher than recovery without treatment
%% All classes plotted
gamma_tVec = [0.15, 0.2, 0.25, 0.5];
n = length(gamma_tVec);
figure
for i = 1:n
    gamma_t = gamma_tVec(i);

[t, state_variable]=ode45(@(t, state_variable) SIR_EVO_MD(t, state_variable, N, beta, gamma_t, gamma_nt), T1, [S_0 I_t_0 I_nt_0 R_0]);
S=state_variable(:,1);
I_t=state_variable(:,2);
I_nt=state_variable(:,3);
R=state_variable(:,4);


loyolagreen = 1/255*[50,205,50];
%orange = 1/255*[255,127,36];
orange = 1/255*[238,173,14];
subplot(2,2,i)
plot(t,S,'Color', loyolagreen,'LineWidth',2); hold on
plot(t,I_t,'m','LineWidth',3); hold on
plot(t,I_nt,'r','Linewidth',2); hold on
plot(t,R,'b','LineWidth',2);
axis([0 20 0 1200])
xlabel('Years')
ylabel('Disease Incidence')
h=legend('S', 'I_t','I_n_t','R','Location','northwest');
end

%h=legend('S', 'I_t','I_n_t','T','R','Location','northwest');
%% Only infected classes plotted:
%%
gamma_tVec = [0.15, 0.2, 0.25, 0.5];
n = length(gamma_tVec);
figure
for i = 1:n
    gamma_t = gamma_tVec(i);

[t, state_variable]=ode45(@(t, state_variable) SIR_EVO_MD(t, state_variable, N, beta, gamma_t, gamma_nt), T1, [S_0 I_t_0 I_nt_0 R_0]);
S=state_variable(:,1);
I_t=state_variable(:,2);
I_nt=state_variable(:,3);
R=state_variable(:,4);

subplot(2,2,i)
plot(t,I_nt+I_t,'k','LineWidth',3); hold on
plot(t,I_t,'m','LineWidth',2); hold on
plot(t,I_nt,'r','Linewidth',2); hold on
axis([0 30 0 1000])
%title('$\delta_{t}=$','Interpreter','latex');
end
xlabel('Years')
ylabel('Disease Incidence')
h=legend('I_0 + I_1', 'I_0','I_1','Location','northeast');