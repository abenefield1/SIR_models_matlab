%% Initial Conditions:
i_0_0=1;
i_k_0=0;
I_0=N^(-1);
Q_0=1-N^(-1);

%% Parameters
% Constant Parameters:
v=0.2;
n=60;
% Parameters that will vary eventually:
B= 0.24; % Rate of infection, beta
O=0.1; % herd immnity, theta
a=0.01; % immne escape per amino acid change, alpha
N=10^3; % host population size
m=0.001; % mtation rate, mu



%for i = 1:n
[t, state]=ode45(@(t, state) mutation_mod_Boni(t, state, N, B, m, O,T_l, T_k), T1, [i_0_0 i_k_0 I_0 Q_0]);
i_0=state(:,1);
i_k=state(:,2);
i_n=state(:,3);
Q=state(:,4);
I=state(:,5);

loyolagreen = 1/255*[50,205,50];
%orange = 1/255*[255,127,36];
orange = 1/255*[238,173,14];
%subplot(2,2,i)
plot(t,i_0,'Color', loyolagreen,'LineWidth',2); hold on
plot(t,i_k,'m','LineWidth',3); hold on
plot(t,i_n,'r','Linewidth',2); hold on
plot(t,Q,'b','LineWidth',2);
plot(t,I,'b','LineWidth',2);
%axis([0 8 0 1200])
xlabel('Years')
ylabel('Disease Incidence')
h=legend('S', 'I_0','I_k','I_n','Location','northwest');
%end