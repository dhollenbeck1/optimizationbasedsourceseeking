tic
clc
clear;
close all
rng default
LB=[0 0]; %lower bounds of variables
UB=[100 100]; %upper bounds of variables
% pso parameters values
m=2; % number of variables
n=10; % population size
wmax=1; % inertia weight
wmin=1; % inertia weight
c1=2; % acceleration factor
c2=2; % acceleration factor
dis=2;
vmax = 2.5;
% pso main program----------------------------------------------------start
figure;
h=0.5;
ss = 100; noisemax = 1e-4;
p=0:h:ss;
y=0:h:ss;
len_p = length(p);
[yy,xx] = meshgrid(p,y);
Q = 0.5; hi=0; sig_x = 10; sig_y = 10; sig_z = 10; u = 2; z = 0;
% myfunc = @(x,y) (Q/(2*pi*sig_x*sig_y*u)).*exp(-(y-ss/2).^2./(2*sig_y^2)).*(exp(-(z-hi).^2./(2*sig_z^2))+exp(-(z+hi).^2./(2*sig_z^2)));

for i=1:100/h+1
    for j=1:100/h+1
        z(i,j)=1*(p(i)-50)^2+3*(p(i)+y(j)-100)^2+100*rand();
    end
end

c3=1;
c4=3;
maxite=500; % set maximum number of iteration

maxrun=1; % set maximum number of runs need to be
for run=1:maxrun
    
    % pso initialization----------------------------------------------start
    for i=1:n
        for j=1:m
            x0(i,j)=round(LB(j)+rand()*(UB(j)-LB(j)));
        end
    end
    x=x0; % initial population
    v=0.1*x0; % initial velocity
    
    myfunc = -10*(noisemax*randn(len_p,len_p) + GaussianPlume_example(len_p,u));
    for i=1:n
        % can change objective function
        %         f0(i,1)=ofun(x0(i,:),c3,c4,1,maxite);
        k1 = find(p == x0(i,1)); k2 = find(y == x0(i,2));
        if isempty(myfunc(k1,k2))
            f0(i,1) = 1;
        else
            f0(i,1) = myfunc(k1,k2);
        end
    end
    [fmin0,index0]=min(f0);
    pbest=x0; % initial pbest
    gbest=x0(index0,:); % initial gbest
    % pso initialization------------------------------------------------end
    % pso algorithm---------------------------------------------------start
    ite=1;
    tolerance=1;
    while ite<=maxite && tolerance>10^-12
        u = u*exp(-0.0001*ite);
        myfunc = -10*(noisemax*randn(len_p,len_p) + GaussianPlume_example(len_p,u));
        w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight
        % pso velocity updates
        for i=1:n
            for j=1:m
                v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))...
                    +c2*rand()*(gbest(1,j)-x(i,j));
                v(i,j)=min(vmax,norm(v(i,j)))*v(i,j)/norm(v(i,j));%./max(norm(v(i,j)),vmax);
            end
        end
        % pso position update
        for i=1:n
            for j=1:m
                x(i,j)=x(i,j)+v(i,j);
            end
        end
        x = round(x*10)/10;
        % handling boundary violations
        for i=1:n
            for j=1:m
                if x(i,j)<LB(j)
                    x(i,j)=LB(j);
                elseif x(i,j)>UB(j)
                    x(i,j)=UB(j);
                end
            end
        end
        % if round(ite/20)==ite/20
        % contour(p,y,z);
        % hold on;
        contour(xx,yy,myfunc,10); 
        hold on
        for i=1:n
          plot(x(i,1),x(i,2),'*');   
        end
        hold off;
        shg;
        xlim([0,100]);
        ylim([0,100])
%         pause(0.5);
        
        % end
        % evaluating fitness
        a(ite)=x(1,1);
        b(ite)=x(1,2);
        for i=1:n
            %             f(i,1)=ofun(x(i,:),c3,c4,ite,maxite);
            k1 = find(p == x(i,1)); k2 = find(y == x(i,2));
            if isempty(myfunc(k1,k2))
                f(i,1) = 1;
            else
                f(i,1) = myfunc(k1,k2);
            end
        end
        % updating pbest and fitness
        for i=1:n
            if f(i,1)<f0(i,1)
                pbest(i,:)=x(i,:);
                f0(i,1)=f(i,1);
            end
        end
        [fmin,index]=min(f0); % finding out the best particle
        ffmin(ite,run)=fmin; % storing best fitness
        ffite(run)=ite; % storing iteration count
        % updating gbest and best fitness
        if fmin<fmin0
            gbest=pbest(index,:);
            fmin0=fmin;
        end
        % calculating tolerance
        if ite>100
            tolerance=abs(ffmin(ite-100,run)-fmin0);
        end
        % displaying iterative results
        if ite==1
            fprintf('Iteration Best particle Objective fun\n');
        end
        fprintf('%8g %8g %8.4f\n',ite,index,fmin0);
        ite=ite+1;
    end
    % pso algorithm-----------------------------------------------------end
    gbest;
    fvalue=1*(gbest(1)-50)^2+3*(gbest(1)+gbest(2)-100)^2;
    fff(run)=fvalue;
    rgbest(run,:)=gbest;
    fprintf('--------------------------------------\n');
end
% pso main program------------------------------------------------------end
fprintf('\n\n');
fprintf('*********************************************************\n');
fprintf('Final Results-----------------------------\n');
[bestfun,bestrun]=min(fff)
best_variables=rgbest(bestrun,:)
fprintf('*********************************************************\n');
toc
% PSO convergence characteristic
figure;
plot(ffmin(1:ffite(bestrun),bestrun),'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('PSO convergence characteristic')

figure;
% h=0.1;
% p=0:h:100;
% y=0:h:100;
% for i=1:100/h+1
%     for j=1:100/h+1
%         z(i,j)=1*(p(i)-50)^2+3*(p(i)+y(j)-100)^2+100*rand();
%     end
% end
contour(xx,yy,myfunc); 
% contour(p,y,z);
hold on
plot(x(:,1),x(:,2),'*')
plot(0,50,'o')