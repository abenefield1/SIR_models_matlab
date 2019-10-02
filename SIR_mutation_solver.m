%Test
%Initial Conditions:
S_0 = 1000;
I_0_0=25;
I_k_0=0;
I_n_0=0;
T1 = 0:40;

N=1025;
beta= 0.56; % Rate of infection
gamma=0.2;
mu=0.005;
%n=60

%for i = 1:n
[t, state]=ode45(@(t, state) SIR_mutation_mod(t, state, N, beta, gamma, mu), T1, [S_0 I_0_0 I_k_0 I_n_0]);
S=state(:,1);
I_0=state(:,2);
I_k=state(:,3);
I_n=state(:,4);

loyolagreen = 1/255*[50,205,50];
%orange = 1/255*[255,127,36];
orange = 1/255*[238,173,14];
%subplot(2,2,i)
plot(t,S,'Color', loyolagreen,'LineWidth',2); hold on
plot(t,I_0,'m','LineWidth',3); hold on
plot(t,I_k,'r','Linewidth',2); hold on
plot(t,I_n,'b','LineWidth',2);
axis([0 8 0 1200])
xlabel('Years')
ylabel('Disease Incidence')
h=legend('S', 'I_0','I_k','I_n','Location','northwest');
%end