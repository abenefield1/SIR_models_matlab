%Initial Conditions:
S_0 = 1000; % Susceptible
I_0=1; % Infected
R_0=0; % Recovered
b=100; % birth rate into susceptible
D=0.1; % death rate (independent of disease)
N=1000;

detTime = 50;
endTime = 150;
T1 = 0:detTime;
T2 = detTime+1:endTime;
totalT=0:endTime;

nu=0.2; % Recovery rate
beta=0.001; % Transmission rate
%% 0:50 - burn-in (no detection rate)
%%
% $R_{0}=\frac{\beta b}{\delta(\delta+\nu)}$

det=0;

[t, class]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T1,[S_0 I_0 R_0]);
S=class(:,1);
I=class(:,2);
R=class(:,3);


%% 50:150 - detection added after burn-in
%%
% $R_{0}=\frac{\beta b}{\delta(\delta+\nu+d_{k})}$

DetVec=[1, 0.5, 0.11, 0.09, 0.009, 0.0009,0.0001, 0];
Names=string(DetVec);
n = length(DetVec);
figure(1)
for i = 1:n
    det = DetVec(i);
    [t, class2]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T2, class(size(class,1),:) );
    S=class2(:,1);
    I=class2(:,2);
    R=class2(:,3);
    
    subplot(0.5*n,2,i)
    p1=plot(t,S,'g','LineWidth',2); hold on
    p2=plot(t,I,'r','LineWidth',2); hold on
    p3=plot(t,R,'b','LineWidth',2); hold on
    %axis([0 150 0 3000])
    ylabel('Incidence')
    title(sprintf('$d_{k}= %s$',Names{i}),'Interpreter','latex', 'FontSize', 12, 'FontName', 'Times New Roman');
    R_nought=(beta*b)/(D*(D + nu + det));
    text(100,max(S)*0.8,sprintf('$R_{0}= %.4f$',R_nought),'Interpreter','latex', 'FontSize', 12, 'FontName', 'Times New Roman')
    grid on
end
suplabel('Years');
hL = legend([p1,p2,p3],{'Susceptible (S)', 'Infected (I)','Recovered (R)'}, 'Orientation', 'horizontal');
newPosition = [0.4 0.87 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits, 'color','none','Box','off');


%%
DetVec=[1, 0.5, 0.11, 0.09, 0.009, 0.0009,0.0001, 0];
Names=string(DetVec);
n = length(DetVec);
figure(2)
for i = 1:n
    det = DetVec(i);
    [t, class2]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T2, class(size(class,1),:) );
    S=class2(:,1);
    I=class2(:,2);
    R=class2(:,3);
    
    %p1=plot(t,S,'g','LineWidth',2); hold on
    p2=plot(t,I,'r','LineWidth',2); hold on
    %p3=plot(t,R,'b','LineWidth',2); hold on
    %axis([0 150 0 3000])
    ylabel('Incidence')
    grid on
end
suplabel('Years');
hL = legend([p1,p2,p3],{'Susceptible (S)', 'Infected (I)','Recovered (R)'}, 'Orientation', 'horizontal');
newPosition = [0.4 0.87 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits, 'color','none','Box','off');


%% Plotting entire time span 0-150
% rbind class and class2
% figure(4)
% classes=vertcat(class,class2);
%     S=classes(:,1);
%     I=classes(:,2);
%     R=classes(:,3);
%     
%     p1=plot(totalT,S,'g','LineWidth',2); hold on
%     p2=plot(totalT,I,'r','LineWidth',2); hold on
%     p3=plot(totalT,R,'b','LineWidth',2); hold on
% 

%%
DetVec = [0, 0.09, 0.11, 1];
n = length(DetVec);
figure(3)
for i = 1:n
    det = DetVec(i);
    [t, class2]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T2, class(size(class,1),:) );
    S=class2(:,1);
    I=class2(:,2);
    R=class2(:,3);
    
    
    subplot(n,1,i)
    plot(t,S,'g','LineWidth',2); hold on
    plot(t,I,'r','LineWidth',2); hold on
    plot(t,R,'b','LineWidth',2); hold on
    %axis([0 50 0 500])
    ylabel('Incidence')
    h=legend('Susceptible (S)', 'Infected (I)','Recovered (R)','Location','northwest');
end
xlabel('Years')

%% 50:150 - post MDT: after burn-in
DetVec=[1, 0.5, 0.11, 0.09, 0.009, 0.0009,0.0001, 0];
Names=string(DetVec);
n = length(DetVec);
figure(4)
for i = 1:n
    det = DetVec(i);
    [t, class2]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T2, class(size(class,1),:) );
    classes=vertcat(class,class2);
    S=classes(:,1);
    I=classes(:,2);
    R=classes(:,3);
    
    subplot(0.5*n,2,i)   
    p1=plot(totalT,S,'g','LineWidth',2); hold on
    p2=plot(totalT,I,'r','LineWidth',2); hold on
    p3=plot(totalT,R,'b','LineWidth',2); hold on
    x1=50;
    xline(x1,'--');
    %axis([0 150 0 3000])
    ylabel('Incidence')
    title(sprintf('$d_{k}= %s$',Names{i}),'Interpreter','latex', 'FontSize', 12, 'FontName', 'Times New Roman');
    R_nought=(beta*b)/(D*(D + nu + det));
    text(100,max(S)*0.8,sprintf('$R_{0}= %.4f$',R_nought),'Interpreter','latex', 'FontSize', 12, 'FontName', 'Times New Roman')
    grid on

end
suplabel('Years');
hL = legend([p1,p2,p3],{'Susceptible (S)', 'Infected (I)','Recovered (R)'}, 'Orientation', 'horizontal');
newPosition = [0.4 0.87 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits, 'color','none','Box','off');

clear all
