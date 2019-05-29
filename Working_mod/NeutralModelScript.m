% neutral model call

S0 = 1000; % initial susceptibles
I0 = 1; % initial base strain/class infection
n = 60; % number of classes possible;
totalTime = 20;
a = 1; % immunity distance parameter


y0 = [S0; I0; zeros(n,1)]; % initial conditions as column vector
myBeta = 0.005;  % transmission parameter
mu = 0.01;  % strain mutation parameter
nu = 0.2;   %   

% dydt = NeutralModelFn(t, y, myBeta, nu, mu)

[time, abundances] = ode45( @(time, abundances) NeutralModelFn(time, abundances, myBeta, nu, mu, a), [0, totalTime], y0 );

plotCols = [1, 2, 10];
for i = 1:length(plotCols)
    subplot(1,length(plotCols),i);
    plot(time, abundances(:,plotCols(i)));
end