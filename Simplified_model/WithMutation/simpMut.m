function f=simpMut(t, class, N, beta, nu, b, D, det, mu);
%% 

S=class(1);
I_1=class(2);
I_2=class(3);
I_n=class(4);
R=class(5);

f(1,1)= -beta*S*(I_1 + I_2 + I_n) + b - D*S; %  Susceptible
f(2,1)= beta*S*I_1 - (nu + det + mu + D)*I_1; %  Strain I1
f(3,1)= beta*S*I_2 - (nu + det + mu + D)*I_2 + mu*I_1; %Strain I2
f(4,1)= beta*S*I_n - (nu + det + D)*I_n + mu*I_2; % Strain In
f(5,1)= (nu + det)*(I_1 + I_2 +I_n) - D*R; %  Recovered
end