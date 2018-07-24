clc;clear all;close all

% Initialize model parameters 
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
R = 1e3
K = lqr( A,B,Q,R);

par.dist = 0;
y0 = [0; 0; pi*12/12; 0];
y0 = [0; 0; pi*7/10; 0];
% y0 = [0; 0; pi*4/10; 0]; 
% y0 = [0; 0; pi*2/10; 0];

yGoal = [0; 0; pi; 0];
% tInt

dt = 0.01; 
tLast = 31;
% tSamp = 0.1;
tInt = dt:dt:tLast;
%% simulate cart 
% [t,yInt] = ode45(@(t,y)cartSinglePend(y, -K*(y-yGoal), par), tInt, y0);

[t,yInt] = ode45(@(t,y) cartSinglePendDist(t, y, hybrid_LQR_V(y,yGoal,K,par) , par),tInt, y0);
% [t,yInt] = ode45(@(t,y) cartSinglePend(y, 0 , par),tInt, y0);


%% postprocessing
% u = - (yInt-repmat(yGoal',length(yInt),1))*K';
% yDot = cartSinglePend_matrixOperation(yInt,u,par);
% 
% e = repmat(yGoal',[size(yInt,1),1]) - yInt;
% J = sum(sum( e.^2*Q' + u.^2*R ));
% display( ['Total cost J = ',num2str(sum(J)) ])
% 
x = yInt(:,1);
xDot = yInt(:,2); 
thet = yInt(:,3);
thetDot = yInt(:,4);
% xDotDot = yDot(2,:)'; 
% thetDotDot = yDot(4,:)';
% momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
Ep = - par.mp* -cos( thet )*L* par.g;
Ek = 0.5*par.mc*xDot.^2 + 0.5*par.mp*( (xDot + par.L*cos(thet).*thetDot).^2 + (-sin(thet).*L.*thetDot).^2  );
% Fp = 1./sin(thet) .* (par.mc*xDotDot + par.bX*xDot - u );

% total energy 
figure( 'Position',[100,100,500,350])
    plot(tInt, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
% % momentum conservation
% figure('Position',[600,100,500,350]);
%     plot(tInt, momentum)
%     xlabel('Time (s)')
%     ylabel('Momentum in X (kg*m/s)')
%     
% figure( 'Position',[1100,100,500,350] );
%     plot(tInt,u)
%     xlabel('Time (s)')
%     ylabel('Control input u')
% % animate cart 
figure( 'Position',[100,550,1000,400] );
% plotPar.Frange = [min(Fp),max(Fp)];
for j = 10:10:length(tInt)
    drawCartSinglePend( yInt(j,:), tInt(j), par)
end
%% 
for j = 1:length(yInt)
   u(j) = hybrid_LQR_V(yInt(j,:)', yGoal,K,par);  
end
% figure();plot(u)
figure();plot(u)
% %% 
% figure()
%     subplot(311)
%     plot(tInt,thet)
%     subplot(312)
%     plot(tInt,thetDot)
%     subplot(313)
%     plot(tInt,xDot)