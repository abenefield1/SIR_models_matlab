%Initial Conditions:
S_0 = 1000; % Susceptible
I_0=50; % Infected
R_0=0; % Recovered
N = S_0+ I_0 + R_0;
b=100; % birth rate into susceptible
D=0.1; % death rate (independent of disease)

%vaccTime = 100;
%endTime = 150;
%T1 = 0:vaccTime;
%T2 = vaccTime:endTime;
T2 = 0:100; % Time

nu=0.2; % Recovery rate
beta=0.0004; % Transmission rate
det=0.5
%% 0:100 - pre vaccination

[t, class]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T2,[S_0 I_0 R_0]);
S=class(:,1);
I=class(:,2);
R=class(:,3);


%% 100:150 - post vaccination

DetVec = [0, 0.0001, 0.09, 0.11, 0.5, 1];
n = length(DetVec);
figure
for i = 1:n
    det = DetVec(i);
    sigma=0.5;
    deltaI=0.2;
    [t, class2]=ode45(@(t, class) simpModDet(t, class, N, beta, nu, b, D, det), T2, class(size(class,1),:) );
    S=class2(:,1);
    I=class2(:,2);
    R=class2(:,3);
    
    
    subplot(0.5*n,2,i)
    plot(t,S,'g','LineWidth',2); hold on
    plot(t,I,'r','LineWidth',2); hold on
    plot(t,R,'b','LineWidth',2); hold on
    %axis([0 50 0 500])
    ylabel('Incidence')
    h=legend('Susceptible (S)', 'Infected (I)','Recovered (R)','Location','northwest');
end
xlabel('Years')

%%
DetVec = [0, 0.09, 0.11, 1];
n = length(DetVec);
figure
for i = 1:n
    det = DetVec(i);
    sigma=0.5;
    deltaI=0.2;
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