%%% Malaria Model incorporating R class but no vector population dynamics
S0 = 1000; % initial susceptibles
I0 = 20; % initial base strain/class infection
I10=0; %I10=15;
I20=0; %I20=10;
I30=0; %I30=5;
R0=0;
n=100; %n = 60; % number of classes possible;
ni=(0:1:100);
totalTime = 150;
a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
b = 0.015; % births 0.015^-1/yr THIS IS THE CORRECT BIRTH RATE AND IT CREATES GOOD MODEL
N=S0+I0+I10+I20+I30;


y0 = [S0; I0;I10;I20;I30 ;zeros(n-3,1); R0]; % initial conditions as column vector
myBeta = 0.5; %myBeta = 0.0294;  % transmission parameter; based on calculations CORRECT
mu=1.2;%mu = 0.9398; %mu = 0.249; %mu = 3; %mu = 0.002;% strain mutation parameter - definitely need to adjust, but currently in units of subs/site of the ompA C. trachomatis gene
  % strain mutation parameter - definitely need to adjust, but currently in units of subs/site of the ompA C. trachomatis gene. https://jb.asm.org/content/191/23/7182
death = 0.015; % death rate (natural) same as birth rate
gamma=0.06; % rate of loss of immunity
nu=0.0384; %recovery rate


[time, abundances] = ode45( @(time, abundances) malariaRec_FN(time, abundances, myBeta, nu, mu, b, a, N, gamma), [0, totalTime], y0 );
S=abundances(:,1);
I_0=abundances(:,2);
I_k=abundances(:,3:n+2);
R=abundances(:,n+3);

figure(7)
plot(time, S, 'linewidth', 2); hold on
plot(time, I_0, 'linewidth', 2); hold on
plot(time, sum(I_k, 2), 'linewidth', 2); hold on
plot(time, R, 'linewidth', 2); hold on
h=legend('S','I_{0}','I_{k}', 'R','Location','northeast');

% 
% figure(7)
% plot(time, S, 'linewidth', 2); hold on
% plot(time, I_0, 'linewidth', 2); hold on
% plot(time, sum(I_k, 2), 'linewidth', 2); hold on
% plot(time, R, 'linewidth', 2); hold on
% h=legend('S','I_{0}','I_{k}', 'R','Location','northeast');

% cols=[3:5:n];
% for i =1:length(cols)
%     plot(time, abundances(:,cols(i)));
% end

figure(1)
plotCols = [3, 30, 60];
for i = 1:length(plotCols)
    subplot(length(plotCols),1,i);
    plot(time, abundances(:,plotCols(i)), 'linewidth', 2);
    %plot(time, abundances(:,3), 'linewidth', 2);
    %plot(time, abundances(:,30), 'linewidth', 2);
    %plot(time, abundances(:,60), 'linewidth', 2);
    set(gca, 'FontName', 'Times New Roman')
    set(gca,'FontSize',30)
    %ylabel('Disease Prevalence')
end
sgtitle('Incidence of Strain Variants','FontName', 'Times New Roman','FontSize',30)
%title('Incidence of Strain Variants')
 %set(gca, 'FontName', 'Times New Roman')
    %set(gca,'FontSize',20)
xlabel('Time')
%[ax,h2]=suplabel('Abundance of Infectious Classes','y'); 

figure(2)
Sum = sum(abundances, 2);
f2 = plot(time, Sum, 'linewidth', 2);
set(gca,'FontSize',20)
set(gca, 'FontName', 'Times New Roman')
xlabel('Time')
ylabel('Force of Infection ($\sum_{k=0}^{n}I_{k}$)', 'Interpreter','latex')
title('Total Force of Infection with Mutating Pathogen Strains')


dMax_Vec = [0, 0.3, 0.6, 1];
l = length(dMax_Vec);
% figure(3)
% for i = 1:l
%     dMax = dMax_Vec(i);
% 
%     [time, abundances] = ode45( @(time, abundances) malariaRec_FN(time, abundances, myBeta, nu, mu, b, a, N, gamma), [0, totalTime], y0 );
% 
%      subplot(2,2,i)
%     Sum = sum(abundances(2:102), 2);
%     plot(time, Sum, 'linewidth', 2);
% end

figure(4)
for i = 1:l
    dMax = dMax_Vec(i);
    [time, abundances] = ode45( @(time, abundances) malariaRec_FN(time, abundances, myBeta, nu, mu, b, a, N, gamma), [0, totalTime], y0 );
    subplot(2,2,i)
    plot(time, abundances(:,1), 'linewidth', 2); hold on
    plot(time, abundances(:,2), 'linewidth', 2); hold on
    plot(time, abundances(:,3), 'linewidth', 2); hold on
    plot(time, abundances(:,15), 'linewidth', 2); hold on
    plot(time, abundances(:,30), 'linewidth', 2); hold on
    plot(time, abundances(:,45), 'linewidth', 2); hold on
    plot(time, abundances(:,62), 'linewidth', 2); hold on
    h=legend('S','I_{0}','3', '15','30','45','62','Location','northeast');
    %axis([0 100 0 1000])

end

figure(5)
dMax=0.5;
plot(time, abundances(:,1), 'linewidth', 2); hold on
plot(time, abundances(:,2), 'linewidth', 2); hold on
plot(time, abundances(:,3), 'linewidth', 2); hold on
plot(time, abundances(:,15), 'linewidth', 2); hold on
plot(time, abundances(:,30), 'linewidth', 2); hold on
plot(time, abundances(:,45), 'linewidth', 2); hold on
plot(time, abundances(:,60), 'linewidth', 2); hold on
plot(time, abundances(:,75), 'linewidth', 2); hold on
plot(time, abundances(:,90), 'linewidth', 2); hold on
plot(time, abundances(:,102), 'linewidth', 2); hold on
plot(time, abundances(:,103), 'linewidth', 2); hold on
h=legend('S','I_{0}','3', '15','30','45','62','75', '90', '102', 'R','Location','northeast')

figure(5)
dMax=0.5;
plot(time, abundances(:,2), 'linewidth', 2); hold on
plot(time, abundances(:,3), 'linewidth', 2); hold on
plot(time, abundances(:,15), 'linewidth', 2); hold on
plot(time, abundances(:,30), 'linewidth', 2); hold on
plot(time, abundances(:,45), 'linewidth', 2); hold on
plot(time, abundances(:,60), 'linewidth', 2); hold on
plot(time, abundances(:,75), 'linewidth', 2); hold on
plot(time, abundances(:,90), 'linewidth', 2); hold on
plot(time, abundances(:,102), 'linewidth', 2); hold on
plot(time, abundances(:,103), 'linewidth', 2); hold on
h=legend('I_{0}','3', '15','30','45','62','75', '90', '102', 'R','Location','northeast')

   
count=0;
dMax_Vec = [0, 0.3, 0.6, 1];
mu_Vec = [0, 0.3, 0.98, 1.02, 1.5];
Names = {'0'; '0.3'; '0.6'; '1'};
Names2 = {'0'; '0.3'; '0.98'; '1.02';'1.5';'0'; '0.3'; '0.98'; '1.02';'1.5';'0'; '0.3'; '0.98'; '1.02';'1.5';'0'; '0.3'; '0.98'; '1.02';'1.5'};

figure(6)
for i = 1:length(dMax_Vec)
    for j = 1:length(mu_Vec)
        dMax = dMax_Vec(i);
        mu = mu_Vec(j);
        [time, abundances] = ode45( @(time, abundances) malariaR_FN(time, abundances, myBeta, nu, mu, b, a, N), [0, totalTime], y0 );
        
        count = count + 1;
        subplot(4,5,count)
        plot(time, abundances(:,1), 'linewidth', 2); hold on
        plot(time, abundances(:,2), 'linewidth', 2); hold on
        plot(time, abundances(:,3), 'linewidth', 2); hold on
        plot(time, abundances(:,15), 'linewidth', 2); hold on
        plot(time, abundances(:,30), 'linewidth', 2); hold on
        plot(time, abundances(:,45), 'linewidth', 2); hold on
        plot(time, abundances(:,62), 'linewidth', 2); hold on
        caption = sprintf('d_{max}= %s',Names{i});
        ylabel(caption, 'FontSize', 12, 'FontName', 'Times New Roman');
        title(sprintf('$\\mu= %s$',Names2{j}),'Interpreter','latex', 'FontSize', 12, 'FontName', 'Times New Roman');
        %h=legend('S','I_{0}','3', '15','30','45','62','Location','northeast');
        grid on;
    end

end


% S=abundances(:,1);
% I_0=abundances(:,2);
% I_k=abundances(:,3:n+2);
% R=abundances(:,n+3);

% figure(7)
% plot(time, S, 'linewidth', 2); hold on
% plot(time, I_0, 'linewidth', 2); hold on
% plot(time, sum(I_k, 2), 'linewidth', 2); hold on
% plot(time, R, 'linewidth', 2); hold on
% h=legend('S','I_{0}','I_{k}', 'R','Location','northeast');
%     