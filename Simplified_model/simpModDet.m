function f=simpModDet(t, class, N, beta, nu, b, D, det);
%% 

S=class(1);
I=class(2);
R=class(3);

f(1,1)= -beta*S*I + b - D*S; %  Susceptible
f(2,1)= beta*S*I - (nu + det + D)*I; %  Iaccine-Infected
f(3,1)= (nu + det)*I - D*R; %  Recovered
end