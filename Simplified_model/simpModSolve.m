%Initial Conditions:
S_0 = 1000;
I_0=50; 
R_0=0;
N = S_0+ I_0 + R_0;
b=100;
D=0.1;

%vaccTime = 100;
%endTime = 150;
%T1 = 0:vaccTime;
%T2 = vaccTime:endTime;
T2 = 0:150;

nu=0.2;
beta=0.0001;

%% 0:100 - pre vaccination

[t, class]=ode45(@(t, class) simpMod(t, class, N, beta, nu,b,D), T2,[S_0 I_0 R_0]);
S=class(:,1);
I=class(:,2);
R=class(:,3);


%% 100:150 - post vaccination

betavec = [0, 0.0003, 0.0004];
n = length(betavec);
figure
for i = 1:n
    beta = betavec(i);
    sigma=0.5;
    deltaI=0.2;
    [t, class2]=ode45(@(t, class) simpMod(t, class, N, beta, nu,b,D), T2, class(size(class,1),:) );
    S=class2(:,1);
    I=class2(:,2);
    R=class2(:,3);


    
    
    subplot(n,1,i)
    plot(t,S,'k','LineWidth',2); hold on
    plot(t,I,'r','LineWidth',2); hold on
    plot(t,R,'b','LineWidth',2); hold on
    %axis([0 50 0 500])
    ylabel('Incidence')
   % h=legend('R',-1);
    h=legend('Susceptible (S)', 'Infected (I)','Recovered (R)','Location','northwest');
end
xlabel('Years')