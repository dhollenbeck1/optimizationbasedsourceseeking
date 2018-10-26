function f=ofun(x,c1,c2,ite,max)
% objective function (minimization)
of=(c1-c1/(10*max-ite))*(x(1)-50)^2+(c2-c2/(12*max-ite))*(x(1)+x(2)-100)^2;
% constraints (all constraints must be converted into <=0 type)
% if there is no constraints then comments all c0 lines below
c0=[];
% c0(1)=x(1)+x(2)+x(3)-5; % <=0 type constraints
% c0(2)=x(1)^2+2*x(2)-x(3); % <=0 type constraints
c0(1)=-1; % <=0 type constraints
c0(2)=-1; % <=0 type constraints
% defining penalty for each constraint
for i=1:length(c0)
if c0(i)>0
c(i)=1;
else
c(i)=0;
end
end
penalty=10000; % penalty on each constraint violation
f=of+penalty*sum(c)+100*rand(); % fitness function
end