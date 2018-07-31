clc;clear all;close all


run('config_cartSinglePend');

mc = 5;
m1 = .5;
m2 = .5;
m3 = .5;
g = -10;
L1 = 0.7;
L2 = 0.7;
L3 = 0.7;
b = 0.001;
c = 0.001;
par.mc = mc; par.m1 = m1; par.m2 = m2; par.m3 = m3; par.L1 = L1; par.L2 = L2; par.L3 = L3; par.b = b; par.c=c; par.g=g;

y0 = [0, pi, pi,pi*99/100 ,0,0,0,0];

u = 0;
dt = 0.002;
tInt = dt:dt:10;


%% simulate
[t,yInt] = ode45(@(t,y)cartTriplePend(y,u,par), tInt, y0);

%% check total energy  
x = yInt(:,1);
th1 = yInt(:,2);
th2= yInt(:,3);
th3= yInt(:,4);
xDot = yInt(:,5); 
th1Dot = yInt(:,6); 
th2Dot = yInt(:,7);
th3Dot = yInt(:,8);
Ep = (m1+m2+m3)*L1*g*cos(th1)   + (m2+m3)*g*(L2*cos(th2 ))+  m3*g*(L3*cos(th3 ));
Ek =  0.5*mc*xDot.^2 ...
    + 0.5*m1*( (xDot+L1.*cos(th1).*th1Dot).^2 + (L1.*sin(th1).*th1Dot).^2 ) ...
    + 0.5*m2*( (xDot + L1.*cos(th1).*th1Dot   + L2.*cos(th2).*th2Dot).^2 ...     
                     +(L1.*th1Dot.*sin(th1)   + L2.*th2Dot.*sin(th2) ).^2  )...
    + 0.5*m3*( (xDot + L1.*cos(th1).*th1Dot   + L2.*cos(th2).*th2Dot  + L3.*cos(th3).*th3Dot).^2 ...     
                     +(L1.*th1Dot.*sin(th1)   + L2.*th2Dot.*sin(th2)  + L3.*sin(th3).*th3Dot).^2 ...
            );
% total energy 
figure( 'Position',[100,100,500,350]); hold on
    plot(tInt, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
%     
% animate cart 
figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawCartTriplePend( yInt(j,:), tInt(j), par )
end

%% 
[U,S,V] = svd(yInt','econ');
figure( 'Position',[1100,350,800,700] );
    scatter3(V(:,4),V(:,2),V(:,3) )

    
%% 
 vidfile = VideoWriter( [rootPath filesep 'figs' filesep 'cartTriplePend_noControl.mp4'],'MPEG-4');
 vidfile.FrameRate = 30;
 open(vidfile);

figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawCartTriplePend( yInt(j,:), tInt(j), par )
    drawnow
    hold off 
    drawnow
    im = getframe(gcf); 
    writeVideo(vidfile, im);
end
close(vidfile)