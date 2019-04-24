function f=SIR_mutation (t, state, N, beta, gamma, mu);

S=state(1);
I_0=state(2);
I_k=state(3);
I_n=state(4);

f(1,1)=-beta*S*(I_0+I_k+I_n); % Susceptible, S
f(2,1)=beta*S*I_0-(gamma+mu)*I_0; % Infected strain 0, I_0
f(3,1)= beta*S*I_k-(gamma+mu)*I_k+mu*I_0; % Infected strain K, I_k
f(4,1)= beta*S*I_n-gamma*I_n+mu*I_k; % Infected strain n, I_n

end
