function c = GaussianPlume_example(len_p,v)
% function GaussianPlume_example
% Gaussian plume model 
%    using MATLAB analytical solutions                   
%
%   $Ekkehard Holzbecher  $Date: 2006/08/21$
%-----------------------------------------------------------------------
Dy = 0.2197; Dz = 0.2197;           % diffusivities
% v = 2;                              % velocity
lambda = 0;                         % decay rate
Q = 0.5;                             % emission rate(s)
xstack = 0; ystack = 50;            % stack location(s)
xmin = 0; xmax = 100;              % x-axis interval
ymin = 0; ymax = 100;               % y-axis interval (used only for d>1)
H = 0.1;                            % effective stack height(s) 
z = 3;                              % height of observation (=0 for ground surface)
gplot = 1;                          % plot option (=1 yes; =0 no)
gcont = 1;                          % contour plot option (=2 filled; =1 yes; =0 none) 
res = len_p;

%----------------------------------execution-------------------------------
[x,y] = meshgrid (linspace(xmin,xmax,res),linspace(ymin,ymax,res));
c = zeros (size(x)); e = ones(size(x));
for i = 1:size(Q,2)
    xx = x - xstack(i); yy = y - ystack(i); 
    c = c + Q(i)*e./(4*pi*xx*sqrt(Dy*Dz)).*exp(-v*yy.*yy./(4*Dy*xx)).*... 
    (exp(-v*(z-H(i))*(z-H(i))*e./(4*Dz*xx))+exp(-v*(z+H(i))*(z+H(i))*e./(4*Dz*xx)))...
    .*exp(-lambda*xx/v);
end

%----------------------------------output----------------------------------
% if gplot
%     for i = 10:10:res
% 	    plot (c(:,i)); hold on;
%     end
% end
% if gcont
%     figure;
%     if gcont > 1
%         contourf (x,y,c); colorbar;
%     else
%         surf (x,y,c); 
%     end
% end

cmax = max(max(c));
cmin = min(min(c));
dc = cmax-cmin;
conc = dc*0.15+cmin;
% k = find(abs(c-conc)<1e-6);
% hold on
% plot3(x(k),y(k),conc*ones(length(k),1),'r*')
end