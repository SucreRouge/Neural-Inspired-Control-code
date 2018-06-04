clc;clear all;close all

run('config_cartSinglePend');

%% Linearization
eqTheta = pi;  % which equilibrium? (theta = 0[down], or theta = pi[up] )
s = cos(eqTheta);
A = [ 0   1          0                   0;
      0   -bX/mc      -mp*g/mc            s*bT/L;
      0   0          0                   1;
      0   s*bX/(L*mc) s*(mc+mp)*g/(mc*L)  -bT*(mc+mp)/(mc*mp*L^2)];
B = [0 ;
    1/mc;
    0
    -s/(mc*L)];

%% set control gain  
K = place(A,B,p);

%% run simulink 
load_system('simulink_Cartpend_KControl')
set_param('simulink_Cartpend_KControl', 'StopTime', num2str(tLast) )
set_param('simulink_Cartpend_KControl', 'MaxStep', num2str( dt ) )
sim('simulink_Cartpend_KControl');

%% postprocessing
yout = simout.signals.values(:,1:4);
uout = simout.signals.values(:,5);
tout = simout.time;
yInt = interp1(tout,yout,tInt);
uInt = interp1(tout,uout,tInt);

% determine accelerations, for force on rod 
yDot = cartSinglePend_matrixOperation(yInt,uInt',par);
x = yInt(:,1);
xDot = yInt(:,2); 
thet = yInt(:,3);
thetDot = yInt(:,4);
xDotDot = yDot(2,:)'; 
yDotDot = yDot(4,:)';
momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
Fp = 1./sin(thet) .* (mc*xDotDot + bX*xDot - uInt(:) );
Ep = - (mp)* -cos( thet )*L* g;
Ek = 0.5*mc*xDot.^2 + 0.5*mp*( (xDot +L*cos(thet).*thetDot).^2+ ( L*sin(thet).*thetDot).^2  ) ;

%% Display results
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
    drawCartSinglePendForce( yInt(j,:), tInt(j), Fp(j), par)
end    

