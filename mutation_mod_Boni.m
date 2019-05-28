function f=mutation_mod_Boni (t, state, N, B, m, O,T_l, T_k);

i_0=state(1);
i_k=state(2);
i_n=state(3);
Q=state(4);
I=state(5);

f(1,1)= B*Q*i_0*((1-O) - symsum(((1-O*T_l)*i_l),l,[0 n])) - m*i_0; % Eqn. 3.6: freq. of strain i_0
f(2,1)= B*Q*i_k*((1-O*T_k) - symsum(((1-O*T_l)*i_l),l,[0 n])) - m*i_k + m*i_k_1; % Eqn. 3.7: freq. of strain i_k. Need to deal with the k-1 bit, and also with the conditional: 0<k<n
f(3,1)= B*Q*i_n*((1-O*T_n) - symsum(((1-O*T_l)*i_l),l,[0 n])) + m*i_n_1; % Eqn. 3.8: freq. strain i_n. Need to deal with same issues here
f(4,1)= -B*Q*I; % Eqn. 3.9: freq. of susceptibles
f(5,1)= B*Q*I*(symsum((((1-O*T_k)*i_k)-v*I),k,[0 n])); % Eqn. 3.10: Total Infections
end