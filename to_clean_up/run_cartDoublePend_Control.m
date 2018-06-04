clc;clear all;close all

mc = 10;
m1 = 1;
m2 = 1;
g = -10;
L1 = 1;
L2 = 1;
b = 0; c = 0; d = 0;
% b = 0.1; % c = 0.1; % d = 0.1;
% b = 1; % c = 1; % d = 1;
par.mc = mc; par.m1 = m1; par.m2 = m2; par.L1 = L1; par.L2 = L2; par.b = b; par.c=c; par.d=d; par.g=g;

% stable equilibrium
y0 = [0, 0, 0 ,0,0,0];
yGoal = y0;
theta1 = pi;
theta2 = pi;

% unstable equilibrium
y0 = [0, pi, 0 ,0,0,0];
yGoal = y0;
theta1 = 0;
theta2 = 0;

% y0 = [0, pi, pi,0,0,0];
% yGoal = y0;
% theta1 = pi;
% theta2 = pi;
% yGoal = [2, pi, pi ,0,0,0];
% y0 = [0, pi/6, pi/6 ,0,0,0];
% y0 = [0, pi*1.3, pi/2 ,0,0,0];
% y0 = [0, 0,0,0,pi/10,0];

u = 0;
dt = 0.01;
tInt = dt:dt:10;
%% Linearize
% theta1 = 0;
% theta2 = 0;
S1 = sin( theta1 );
C1 = cos( theta1 );
S2 = sin( theta2 );
C2 = cos( theta2 );
dTheta = theta1-theta2;

D = [ mc+m1+m2              (m1+m2)*L1*C1           m2*L2*C2    ; 
     (m1+m2)*L1*C1          (m1+m2)*L1^2            m2*L1*L2*cos(dTheta) ; 
      m2*L2*C2               m2*L1*L2*cos(dTheta)   m2*L2^2     ];
  
% D = [ mc+m1+m2              -(m1+m2)*L1*S1          -m2*L2*S2    ; 
%      -(m1+m2)*L1*S1          (m1+m2)*L1^2            -m2*L1*L2*sin(dTheta) ; 
%       m2*L2*S2               -m2*L1*L2*sin(dTheta)   m2*L2^2     ];

dGdx = [0   0                   0 ;
       0    -(m1+m2)*g*L1*C1    0;
       0    0                   -m2*g*L2*C2]; 

E = diag([b,c,d]);

H = [1;0;0];



A = [diag([0,0,0]) , eye(3); 
     inv(D)*dGdx   ,inv(D)*diag([b,c,d]) ];
 
 B = [zeros(3,1); inv(D)*H];
 
%  Q = diag([5,50,50,20,700,700]);
 Q = diag([1,1,1,1,1,1]);
 R = 1e3;
 
 K = lqr(A,B,Q,R);













%%

[t,yInt] = ode45(@(t,y) cartDoublePend(y, -K*(y-yGoal'), par ), tInt, y0);





%% check total energy  
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
    
    figure();
    plot(u)
%% animate cart 
figure( 'Position',[100,550,1000,400] );
for j = 10:10:length(tInt)
    drawCartDoublePend( yInt(j,:), tInt(j), par )
end

%% 
% [U,S,V] = svd(yInt','econ');
% figure( 'Position',[1100,350,800,700] );
%     scatter3(V(:,4),V(:,2),V(:,3) )
