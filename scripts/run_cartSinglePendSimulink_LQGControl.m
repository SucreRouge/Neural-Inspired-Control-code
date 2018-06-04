clc;clear all;close all

% Initialize model parameters 
run('config_cartSinglePend');


y0 = [-4; 0; 0; 0];
yGoal = y0 + [5;0;0;0];

%% Linearization
eqTheta = 0;  % which equilibrium? (theta = 0[down], or theta = pi[up] )
s = cos(eqTheta);
A = [ 0   1          0                   0;
      0   -bX/mc      -mp*g/mc            s*bT/L;
      0   0          0                   1;
      0   s*bX/(L*mc) s*(mc+mp)*g/(mc*L)  -bT*(mc+mp)/(mc*mp*L^2)];
B = [0 ;
    1/mc;
    0
    -s/(mc*L)];
% C = [0,1,0,0];
C = [1,0,1,0];
% C = [1,1,1,1];
D = zeros(size(C,1),size(B,2));
% obsv(A,C)
% det(obsv(A,C))

%% set control gain  
% Q = diag( [1,1,10,100]);
Q = diag( [1,1,1,10]);
R = 1e-3;
K = lqr( A,B,Q,R);


Vd = 0.1*eye(4); % disturbance covariance (state uncertainty) 
Vn = 1; %measurement noise 


[Lmat,P,E] = lqe(A,Vd,C,Vd,Vn);  % design Kalman filter
Kf = (lqr(A',C',Vd,Vn))';   % alternatively, possible to design using "LQR" code

sysKF = ss(A-Lmat*C,[B Lmat],eye(4),0*[B Lmat]);  % Kalman filter estimator
% eig(A-B*K)


% BF = [B Vd 0*B];
%% run simulink 
simFile = 'simulink_Cartpend_LQGControl';
load_system( simFile)
set_param( simFile , 'StopTime', num2str(tLast) )
set_param( simFile , 'MaxStep', num2str( dt ) )
sim( simFile );
% 
%% postprocessing
yout = simout.signals.values(:,1:4);
uout = simout.signals.values(:,5);
tout = simout.time;
yInt = interp1(tout,yout,tInt);
uInt = interp1(tout,uout,tInt);

% determine accelerations, for force on rod 
% yDot = cartPendMatrix(yInt,uInt',par);
% x = yInt(:,1);
% xDot = yInt(:,2); 
% thet = yInt(:,3);
% thetDot = yInt(:,4);
% xDotDot = yDot(2,:)'; 
% yDotDot = yDot(4,:)';
% momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
% Fp = 1./sin(thet) .* (mc*xDotDot + bX*xDot - uInt(:) );
% Ep = - (mp)* -cos( thet )*L* g;
% Ek = 0.5*mc*xDot.^2 + 0.5*mp*( (xDot +L*cos(thet).*thetDot).^2+ ( L*sin(thet).*thetDot).^2  ) ;

%% Display results
% total energy 
% figure( 'Position',[100,100,500,350])
%     plot(tInt, Ep+Ek)
%     xlabel('Time (s)')
%     ylabel('Total Energy T+V (J)')
% % momentum conservation
% figure('Position',[600,100,500,350]);
%     plot(tInt, momentum)
%     xlabel('Time (s)')
%     ylabel('Momentum in X (kg*m/s)')
% Rod tension 
figure( 'Position',[1100,100,500,350] );
%     plot(tInt,Fp)
    plot(tInt,uInt)
    xlabel('Time (s)')
    ylabel('Control input u')
% animate cart
figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawCartSinglePendForce( yInt(j,:), tInt(j),10, par)
end    

