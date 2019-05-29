
% neutral model call

S0 = 1000; % initial susceptibles
I0 = 1; % initial base strain/class infection
n = 60; % number of classes possible;
totalTime = 120;
%a = 1; % immunity distance parameter


y0 = [S0; I0; zeros(n,1)]; % initial conditions as column vector
myBeta = 1.2;  % transmission parameter
mu = 0.01;  % strain mutation parameter
nu = 0.2;   %   

% dydt = NeutralModelFn(t, y, myBeta, nu, mu)

[time, abundances] = ode45( @(time, abundances) Neutral_wTest(time, abundances, myBeta, nu, mu, d), [0, totalTime], y0 );


% cols=[3:5:n];
% for i =1:length(cols)
%     plot(time, abundances(:,cols(i)));
% end


plotCols = [3, 10, 20, 30, 40, 50, 60];
for i = 1:length(plotCols)
    subplot(length(plotCols),1,i);
    plot(time, abundances(:,plotCols(i)));
end

