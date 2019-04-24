function f=SIR_EVO_MD (t, state_variable, N, beta, gamma_t, gamma_nt);

%% 

S=state_variable(1);
I_t=state_variable(2);
I_nt=state_variable(3);
R=state_variable(4);

b = 100; % birth rate into susceptible
d=0.1; %death rate (independent of disease)

f(1,1)=b-beta*S*I_t - beta*S*I_nt; % Susceptible, S
f(2,1)= beta*S*I_t - gamma_t*I_t; % Infected treated, I_t
f(3,1)= beta*S*I_nt - gamma_nt*I_nt; % Infected not treated, I_nt
f(4,1)= gamma_t*I_t + gamma_nt*I_nt - d*R; % Recovered, R
end