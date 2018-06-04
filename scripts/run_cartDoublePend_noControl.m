clc;clear all;close all

mc = 10;
m1 = 1;
m2 = 1;
g = -10;
L1 = 1;
L2 = 1;
b = 0;
c = 0;
d = 0;
% b = 0.01;
% c = 0.001;
% d = 0.001;
par.mc = mc; par.m1 = m1; par.m2 = m2; par.L1 = L1; par.L2 = L2; par.b = b; par.c=c; par.d=d; par.g=g;

y0 = [0, pi, pi*99/100,0,0,0];
u = 0;
dt = 0.01;
tInt = dt:dt:10;


%%
[t,yInt] = ode45(@(t,y) cartDoublePend(y, u, par ), tInt, y0);



%% plot figures 
x = yInt(:,1);
th1 = yInt(:,2);
th2= yInt(:,3);
xDot = yInt(:,4); 
th1Dot = yInt(:,5); 
th2Dot = yInt(:,6);
Ep = m1*L1*g*cos( th1 )   + m2*g*( L1*cos( th1 )+L2*cos( th2 ) );
Ek =  0.5*mc*xDot.^2 ...
    + 0.5*m1*( (xDot+L1.*cos(th1).*th1Dot).^2 + (L1.*sin(th1).*th1Dot).^2 ) ...
    + 0.5*m2*( (xDot + L1.*cos(th1).*th1Dot   + L2.*cos(th2).*th2Dot).^2 ...     
                     +(L1.*th1Dot.*sin(th1)   + L2.*th2Dot.*sin(th2) ).^2 ...
            );
% total energy 
figure( 'Position',[100,100,500,350]); hold on
    plot(tInt, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
    
% animate cart 
figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawCartDoublePend( yInt(j,:), tInt(j), par )
end

