N=1020;
S0 = 1000; % initial susceptibles
I0 = 1; % initial base strain/class infection
I10=0;
I20=0;
I30=0;
n = 60; % number of classes possible;
totalTime = 100;
a = 0.06; % testing escape parameter - the higher a, the more rapid the escape
%b = 0.015; % births 0.015^-1/yr THIS IS THE CORRECT BIRTH RATE
b=0.1;

y0 = [S0; I0;I10;I20;I30 ;zeros(n-3,1)]; % initial conditions as column vector
myBeta = 0.0294;  % transmission parameter; based on calculations CORRECT
%myBeta = 0.003;
%mu = 0.9398;  % strain mutation parameter - definitely need to adjust, but currently in units of subs/site of the ompA C. trachomatis gene
mu = 0.01;  % strain mutation parameter - definitely need to adjust, but currently in units of subs/site of the ompA C. trachomatis gene. https://jb.asm.org/content/191/23/7182
death = 0.015; % death rate (natural) same as birth rate
gamma = 0.0384; % recovery rate
%nu = gamma + death;   % composite recovery rate/death rate parameter
nu=gamma;

dMax=0.5;


[time, abundances] = ode45( @(time, abundances) NeutFN(time, abundances, myBeta, nu, mu, b, a, N, dMax), [0, totalTime], y0 );


set(gca, 'FontName', 'Times New Roman')
plotCols = [3, 30, 60];
figure(1)
for i = 1:length(plotCols)
    subplot(length(plotCols),1,i);
    plot(time, abundances(:,plotCols(i)));
end
sgtitle('Abundance of Strain Variants')
xlabel('Years')
[ax,h2]=suplabel('Abundance of Infectious Classes','y'); 

figure(2)
Sum = sum(abundances, 2);
plot(time, Sum, 'linewidth', 2);
%xlim([0 15]);


dMax_Vec = [0, 0.3, 0.6, 1];
l = length(dMax_Vec);
figure
for i = 1:l
    dMax = dMax_Vec(i);
    
    [time, abundances] = ode45( @(time, abundances) NeutFN(time, abundances, myBeta, nu, mu, b, a, N, dMax), [0, totalTime], y0 );

    subplot(2,2,i)
    Sum = sum(abundances, 2);
    plot(time, Sum, 'linewidth', 2);
    axis([0 100 0 2000])
end

figure

for i = 1:l
    dMax = dMax_Vec(i);
    
    [time, abundances] = ode45( @(time, abundances) NeutFN(time, abundances, myBeta, nu, mu, b, a, N, dMax), [0, totalTime], y0 );

    subplot(2,2,i)
    plot(time, abundances(:,1), 'linewidth', 2); hold on
    plot(time, abundances(:,2), 'linewidth', 2); hold on
    plot(time, abundances(:,3), 'linewidth', 2); hold on
    plot(time, abundances(:,15), 'linewidth', 2); hold on
    plot(time, abundances(:,30), 'linewidth', 2); hold on
    plot(time, abundances(:,45), 'linewidth', 2); hold on
    plot(time, abundances(:,62), 'linewidth', 2); hold on
    h=legend('S','I_{0}','3', '15','30','45','62','Location','northeast');
    axis([0 100 0 1000])

end




X = [0 1 2; 3 4 5]
sum(X, 1)
sum(X, 2)