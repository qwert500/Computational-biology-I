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
%% task 1 c
clc,clear all,clf
%%%%%%%%%%%%%%%%%%
%   Parameters   %
%%%%%%%%%%%%%%%%%%
K=10^3;
r=0.1;
b=1;
N0=[1 2 3 10];
tMax=200;
deltaN0=0.1;
%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%
N=zeros(1,tMax);
%%%%%%%%%%%%%%%%%%
%      Main      %
%%%%%%%%%%%%%%%%%%
NSS=r^(1/b)*K
%%%%%%%%%%%%%%%%%%
%   Eigenvalues  %
%%%%%%%%%%%%%%%%%%
lambda1=1-r*b/(r+1);
lambda2=r+1;



N(1:4,1)=N0;
%N2(1:4,1)=N0;
for t=2:tMax
  for i=1:size(N0,2)
    N(i,t)=((r+1)*N(i,t-1))/(1+(N(i,t-1)/K)^b);
  end
end

N2(1)=NSS+deltaN0;
for t=2:tMax
  N2(t)=NSS+deltaN0*lambda1^t;
end
hold on
plot(log(1:tMax),log(N(1:4,1:tMax)))
plot(log(1:tMax),log(N2(1:tMax)),'k')


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

clc, clear all,clf
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
title('Bifurcation Diagram')


