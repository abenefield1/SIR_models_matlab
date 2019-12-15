function f=simpMut(t, class, N, beta, nu, b, D, mu, a);

S=class(1);
I_1=class(2);
I_2=class(3);
I_n=class(4);
R=class(5);

%dMax=1;
%det=dMax * exp(-k * a);
%k loop here

f(1,1)= -beta*S*(I_1 + I_2 + I_n) + b - D*S; %  Susceptible
f(2,1)= beta*S*I_1 - (nu + detection(0,a) + mu + D)*I_1; %  Strain I1
f(3,1)= beta*S*I_2 - (nu + detection(1,a) + mu + D)*I_2 + mu*I_1; %Strain I2
f(4,1)= beta*S*I_n - (nu + detection(2,a) + D)*I_n + mu*I_2; % Strain In
%f(5,1)= (nu + detection(3,a))*(I_1 + I_2 +I_n) - D*R; %  Recovered
f(5,1)= -D*R + ((nu + detection(0,a))*I_1) + ((nu + detection(1,a))*I_2) + ((nu + detection(2,a))*I_n);%  Recovered
%end k loop here
end