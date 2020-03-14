%%%%%%%%%%%
%   PS3   %
%%%%%%%%%%%
clc;
clear all;
clf;

%%%%%%%%%%%%%%
% Constants  %
%%%%%%%%%%%%%%
alpha=2;
beta=1.5;
n0=10;
N=400;
numberOfTimeSteps=150;
numberOfAvrages=1;
SS=(1-beta/alpha)*N;

%%%%%%%%%%%%%%
%    Main    %
%%%%%%%%%%%%%%
n=n0;
nPlot=zeros(numberOfTimeSteps*2+1,2);
avrageInEachTime=ones(numberOfTimeSteps*2+1,1);
nPlot3=zeros(numberOfTimeSteps*2+1,2);
for i=1:numberOfAvrages
  t=0;
  n=n0;
  tt=0;
  tNumber=0;
  while t<numberOfTimeSteps
    
    b=alpha*(1-(n-1)/N)*(n-1);
    d=beta*(n+1);
    lambda=b+d;
    t=t+exprnd(1/lambda,1,1);
    probabilityForRecovery=d/lambda;
    randomNumber=rand();
    tt=0;
    
    if randomNumber<probabilityForRecovery
      n=n-1;
    elseif randomNumber>probabilityForRecovery
      n=n+1;
    end
    
    tNumber=tNumber+1;
    nPlot2(tNumber,1,i)=t;
    nPlot2(tNumber,2,i)=n;
    
%     for lowerLimit=0:0.5:numberOfTimeSteps
%       upperLimit=lowerLimit*2+0.5;
%       if lowerLimit<t && t<upperLimit
%         tt=lowerLimit*2+1;
%         avrageInEachTime(tt)=avrageInEachTime(tt)+1;
%         
%         nPlot3(tt,1)=lowerLimit;
%         nPlot3(tt,2)=n+nPlot3(tt,2);
%       end
%     end
    
    
  end
end

% for i=1:size(avrageInEachTime,1)
%   nPlot3(i,2)=nPlot3(i,2)/avrageInEachTime(i);
% end

hold on
for j=1:numberOfAvrages
plot(nPlot2(:,1,j),nPlot2(:,2,j))
end


x=linspace(0,numberOfTimeSteps);
plot(x,ones(size(x))*SS)
%Text=exp(N*(log(alpha/beta)-(1-beta/alpha)))









