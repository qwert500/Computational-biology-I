%% 
clc,clear all,clf
%%%%%%%%%%%%%%%%%% task 1 b
%   Parameters   %
%%%%%%%%%%%%%%%%%%
K=10^3;
r=0.1;
b=1;
N0=[1 2 3 10];
tMax=100;
%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%
N=zeros(1,tMax);
%%%%%%%%%%%%%%%%%%
%      Main      %
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%   Eigenvalues  %
%%%%%%%%%%%%%%%%%%
lambda1=1-r*b/(r+1);
lambda2=r+1;
NSS=100;
N(1:4,1)=N0;
N2(1:4,1)=N0;
for t=2:tMax
  for i=1:size(N0,2)
    N(i,t)=((r+1)*N(i,t-1))/(1+(N(i,t-1)/K)^b);
    N2(i,t)=N0(i)*(lambda2)^t; %eigenvalue times population
  end
end

hold on
plot(log(1:tMax),log(N(1:4,1:tMax)))
plot(log(1:tMax),log(N2(1:4,1:tMax)),'k')


xlabel('time')
ylabel('Population')
title('Change of population over time')
legend('N0=1','N0=2','N0=3','N0=10')

%% 
clc,clear all,clf
%%%%%%%%%%%%%%%%%% task 1 c
%   Parameters   %
%%%%%%%%%%%%%%%%%%
K=10^3;
r=0.1;
b=1;
N0=[1 2 3 10];
tMax=100;
%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%
N=zeros(1,tMax);
%%%%%%%%%%%%%%%%%%
%      Main      %
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%   Eigenvalues  %
%%%%%%%%%%%%%%%%%%
lambda1=1-r*b/(r+1);
lambda2=r+1;
NSS=100;
N(1:4,1)=N0+NSS;
N2(1:4,1)=N0+NSS;
for t=2:tMax
  for i=1:size(N0,2)
    N(i,t)=((r+1)*N(i,t-1))/(1+(N(i,t-1)/K)^b);
    N2(i,t)=N0(i)*(lambda1)^t+NSS; %eigenvalue times population
  end
end

hold on
plot(log(1:tMax),log(N(1:4,1:tMax)))
plot(log(1:tMax),log(N2(1:4,1:tMax)),'k')


xlabel('time')
ylabel('Population')
title('Change of population over time')
legend('N0=1','N0=2','N0=3','N0=10')
%% task 2 a.) cannibalism

clc, clear all,clf
alpha=0.01;
eta0=900;
plotMatrix=[];
Rplot=[];
for R=1:0.1:30
  eta(1)=eta0;
  for t=1:300
    eta(t+1)=eta(t)*R*exp(-alpha*eta(t));
    RR(t)=R;
  end
  Rplot=[Rplot RR(201:300)];
  plotMatrix=[plotMatrix eta(201:300)];
end
%j=1;
%RR=0;
%for i=1:300
%  RR=RR+0.1;
%  for t=1:100
%    j=1+j;
%    Rplot(j)=RR;
%  end
%end
sz = 25;
c = linspace(1,30,length(Rplot));
scatter(Rplot,plotMatrix,sz,c,'filled')
xlabel('R')
ylabel('\eta')
title('Bifurcation Diagram')
%% task 2 b.) cannibalism

clc, clear all ,clf
alpha=0.01;
eta0=900;
plotMatrix=[];
R=[5 9 13];
eta(1,1)=eta0;
eta(2,1)=eta0;
eta(3,1)=eta0;
for i=1:3
  for t=1:300
    eta(i,t+1)=eta(i,t)*R(i)*exp(-alpha*eta(i,t));
  end
end
hold on

plot(1:301,eta(1:3,1:301),'.')
legend('R=5','R=9','R=13')
xlabel('t')
ylabel('\eta')
title('Population dynamics for different R')
%% task 3 a (just for fun)
clc, clear all, clf
%%%%%%%%%%%%%%%%%%
%   parameters   %
%%%%%%%%%%%%%%%%%%
%r=1;
%a=1;
aMax=40;
rMax=40;
nSteps=100;
%%%%%%%%%%%%%%%%%%
%      Main      %
%%%%%%%%%%%%%%%%%%
polynomialFunction=zeros(nSteps,nSteps);
lambda=-1;
for i=1:nSteps
  a=aMax/nSteps*i;
  for j=1:nSteps
    r=rMax/nSteps*j;
    N=r/(r+a);
    P=r/(r+a);
    polynomialFunction(i,j)=(1+r-2*r*N-a*P-lambda)*(1+a*N-2*a*P-lambda)+a*P*a*N;
  end
end
r=linspace(rMax/nSteps,rMax,nSteps);
a=linspace(aMax/nSteps,aMax,nSteps);

hold on
surf(a,r,polynomialFunction)
surf(a,r,zeros(nSteps,nSteps))
xlabel('a')
ylabel('r')
zlabel('P(lambda=1)')
%% Task 3 b 
clc,clear all, clf
%%%%%%%%%%%%%%%%%%%%%%
%     Parameters     %
%%%%%%%%%%%%%%%%%%%%%%
N0=0.25;
P0=0.15;
r=1.4;
a=2;
nTimeSteps=10^3;
%%%%%%%%%%%%%%%%%%%%%%
%       Main         %
%%%%%%%%%%%%%%%%%%%%%%

N=zeros(nTimeSteps,1);
P=zeros(nTimeSteps,1);
N(1)=N0;
P(1)=P0;

for t=2:nTimeSteps
  N(t)=N(t-1)+r*N(t-1)*(1-N(t-1))-a*N(t-1)*P(t-1);
  P(t)=P(t-1)+a*P(t-1)*(N(t-1)-P(t-1));
end
figure(1)
hold on
plot(1:nTimeSteps,N(:),'g')
plot(1:nTimeSteps,P(:),'r')
grid on
legend('Prey','Predator')
xlabel('time')
ylabel('population')
title('Population dynamics')

figure(2)
plot(N,P)
xlabel('population of prey')
ylabel('population of predators')
title('N-P diagram')
%% Task 3 b second part
clc,clear all, clf
%%%%%%%%%%%%%%%%%%%%%%
%     Parameters     %
%%%%%%%%%%%%%%%%%%%%%%
N0=0.25;
P0=0.15;
r=2.2;
a=2;
nTimeSteps=10^3;
%%%%%%%%%%%%%%%%%%%%%%
%       Main         %
%%%%%%%%%%%%%%%%%%%%%%

N=zeros(nTimeSteps,1);
P=zeros(nTimeSteps,1);
N(1)=N0;
P(1)=P0;

for t=2:nTimeSteps
  N(t)=N(t-1)+r*N(t-1)*(1-N(t-1))-a*N(t-1)*P(t-1);
  P(t)=P(t-1)+a*P(t-1)*(N(t-1)-P(t-1));
end
figure(1)
hold on

plot(1:nTimeSteps,N(:),'g.')
plot(1:nTimeSteps,P(:),'r.')

grid on
legend('Prey','Predator')
xlabel('time')
ylabel('population')
title('Population dynamics')

figure(2)
plot(N,P)
xlabel('population of prey')
ylabel('population of predators')
title('N-P diagram')
%% 3b stability check for r
clc, clear all
a=2;
lambda=1;
syms r
vpasolve(((1+r-2*r*(r/(r+a))-a*(r/(r+a))...
  )*(1+a*(r/(r+a))-2*a*(r/(r+a))...
  )+a*(r/(r+a))*a*(r/(r+a)))==1 , r)

%ans=2:)










