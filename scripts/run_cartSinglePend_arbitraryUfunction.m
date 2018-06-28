clc;clear all;close all

% Initialize model parameters 
run('config_cartSinglePend');

% u = 0;
% %% simulate cart 
% [t,yInt] = ode45(@(t,y)cartSinglePend(y',u,par),tInt, y0);
% 
% % animate cart 
% figure( 'Position',[100,550,1000,400] );
% for j = 10:10:length(tInt)
%     drawCartSinglePend( yInt(j,:), tInt(j), par)
% end


%% simulate cart 

uPulse = [ 7,0.2,1   -8,0.2,3,   -20,0.1,4.5];
ut =@(t) uPulse(1)* gaussian(t, uPulse(2), uPulse(3) ) + ...
     uPulse(4)* gaussian(t, uPulse(5), uPulse(6) +...
     uPulse(7)* gaussian(t, uPulse(8), uPulse(9) ))  ;
ut =@(t) 0*t;

[t,yInt] = ode45(@(t,y)cartSinglePend(y',ut(t),par),tInt, y0);

% animate cart 
figure( 'Position',[100,550,1000,400] );
for j = 10:5:length(tInt)
    drawCartSinglePend( yInt(j,:), tInt(j), par)
end

%% 
figure();
hold on
scatter3( yInt(:,2), yInt(:,3), yInt(:,4), t)
xlabel('$\dot{x}$'); ylabel('$\theta$');zlabel('$\dot{\theta}$')
scatter3( 0,pi,0,'r','filled')
scatter3( 0,pi,0,40,'r','filled')
scatter3( 0,3*pi,0,40,'r','filled')
ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% ax.box ='off'
% ax.ZAxisLocation = 'origin';
% ax.oaxes
% oaxes