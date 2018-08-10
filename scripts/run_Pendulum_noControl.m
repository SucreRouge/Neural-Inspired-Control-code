clc;clear all;close all

% Initialize model parameters 
run('config_singlePend');
u = 0;

%% simulate cart 
[t,yInt] = ode45(@(t,y)singlePend(y',u,par),tInt, y0);

%% Plot and animate

thet = yInt(:,1);
thetDot = yInt(:,2);

V = -cos(thet)*L*mp*g;
T = 1/2* mp*thetDot.^2*L^2;

figure( 'Position',[100,100,500,350])
    hold on
    plot(tInt, V+T,'k')
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
    legend('Total')
    
%% animate cart 
figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawSinglePend( yInt(j,:), tInt(j), par)
    pause(0.02)
end
