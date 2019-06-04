% Next steps:
% 1.) figure out if the 'd' function is working
    % if the code is working as expected, rethink the math or parameter values
% 2.) Figure out why code crashes for b*S (rather than b*N) about half the time - depending upon max time 
% 3.) Everything else - actually get at the selection pressure and pop-gen values


% neutral model call
N=1000;
S0 = 1000; % initial susceptibles
I0 = 25; % initial base strain/class infection
n = 60; % number of classes possible;
totalTime = 100;
a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
b = 0.7; % Annual population growth rate US


y0 = [S0; I0; zeros(n,1)]; % initial conditions as column vector
myBeta = 0.003;  % transmission parameter
%mu = 0.9398;  % strain mutation parameter - definitely need to adjust, but currently in units of subs/site of the ompA C. trachomatis gene
mu = 0.249;  % strain mutation parameter - definitely need to adjust, but currently in units of subs/site of the ompA C. trachomatis gene. https://jb.asm.org/content/191/23/7182
nu = 0.5;   % composite recovery and death rate parameter


[time, abundances] = ode45( @(time, abundances) Neutral_wTest(time, abundances, myBeta, nu, mu, b, a, N), [0, totalTime], y0 );


% cols=[3:5:n];
% for i =1:length(cols)
%     plot(time, abundances(:,cols(i)));
% end

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
