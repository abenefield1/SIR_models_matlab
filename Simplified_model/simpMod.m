function f=simpMod(t, class, N, beta, nu,b,D);
%% 

S=class(1);
I=class(2);
R=class(3);

%b = 100; % birth rate into susceptible
%d=0.1; %death rate (independent of disease)

%f(1,1)= $-\betaSI$;
%f(2,1)=$\beta SI-\nuI$

f(1,1)= -beta*S*I + b; %  Susceptible
f(2,1)= beta*S*I - (nu + D)*I; %  Iaccine-Infected
f(3,1)= nu*I - D*R; %  Recovered
end