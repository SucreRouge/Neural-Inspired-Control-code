clc;clear all;close all
mc = 5;
mp = 1;
g = -10;
L = 2;
b = 10;

par.mc = mc ; par.mp = mp; par.g = g; par.L = L; par.b = b; 
par.dist = 0.01;
F = 0;
uInt = [0,1,0,0]*F;

% y0 = [0,0,pi*119/120,0];
y0 = [1,0,pi*119/120,0];
% y0 = [0,0,0,0];
dt = 0.01; 
t0 = 0;
tLast = 7;
tInt = dt:dt:tLast;

K = [0,0,0,0];

load_system('simulink_cartpend_1')
set_param('simulink_cartpend_1', 'StopTime', num2str(tLast) )
sim('simulink_cartpend_1');
yout = simout.signals.values(:,1:4);
uout = simout.signals.values(:,5);
tout = simout.time;
yInt = interp1(tout,yout,tInt);
uInt = interp1(tout,uout,tInt);

% determine accelerations, for force on rod 
% u = -repmat(K,length(yInt),1) .* (yInt-repmat([1;0;pi;0]',length(yInt),1));
yDot = cartPendOwn(yInt,uInt',par);
x = yInt(:,1);
xDot = yInt(:,2); 
thet = yInt(:,3);
thetDot = yInt(:,4);
xDotDot = yDot(2,:)'; 
yDotDot = yDot(4,:)';
momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
Fp = 1./sin(thet) .* (par.mc*xDotDot + par.b*xDot - uInt(:) );
Ep = - par.mp* -cos( thet )* par.g;
Ek = 0.5*par.mc*xDot.^2 + 0.5*par.mp*( (xDot + par.L*cos(thet).*thetDot).^2 + (-sin(thet).*thetDot).^2  );

%% Metrics 
% total energy 
figure( 'Position',[100,100,500,350])
    plot(tInt, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
% momentum conservation
figure('Position',[600,100,500,350]);
    plot(tInt, momentum)
    xlabel('Time (s)')
    ylabel('Momentum in X (kg*m/s)')
% Rod tension 
figure( 'Position',[1100,100,500,350] );
    plot(tInt,Fp)
    xlabel('Time (s)')
    ylabel('Force in rod (N)')
% animate cart
figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawCartPendOwn( yInt(j,:), tInt(j), Fp(j), par)
end    

