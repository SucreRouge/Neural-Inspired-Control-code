clc;clear all;close all

mc = 5;
mp = 1;
g = -10;
L = 2;
b = 1;
par.mc = mc ; par.mp = mp; par.g = g; par.L = L; par.b = b; 
par.dist = 0.01;

t0 = 0;
tLast = 10;
% dt = .01;
t = 0:0.01:50;
% t = t0:0.02:tLast;
y0 = [-4,0,pi*9/12,0];

%% Linearized eqs pendulum on the cart 
eqTheta = 0;  % which equilibrium? (theta = 0[down], or theta = pi[up] )
s = cos(eqTheta);

A = [ 0  1         0                   0;
    0   -b/mc      -mp*g/mc             0;
    0   0          0                   1;
    0   s*b/(L*mc) s*(mc+mp)*g/(mc*L)  0];
B = [0 ;
    1/mc;
    0
    -s/(mc*L)];
C = [1 0 0 0];
D = zeros(size(C,1),size(B,1)); 

Vd = .1*eye(4);
Vn = 1;
BF = [B Vd, 0*B];

% obsv(A,C)
% det( obsv(A,C))

sysC = ss(A,BF,C,[0,0,0,0,0,Vn]);
sysFullOutput = ss(A,BF,eye(4),zeros(4,size(BF,2)) );

[Kf,P,E] = lqe(A,Vd,C,Vd,Vn);
Kf2 = (lqr(A',C',Vd,Vn))';
sysKF = ss(A-Kf*C,[B Kf],eye(4),0*[B Kf]);

%% 

uDIST = randn(4,size(t,2));
uNOISE = randn(size(t));
u = 0*t;
u(100:120) = 100;
u(1500:1520) = -100;

uAug = [u;Vd*Vd*uDIST; uNOISE];

[y,t] = lsim( sysC, uAug, t); 
plot(t,y,'b');
[xtrue,t] = lsim(sysFullOutput,uAug,t);

hold on
plot(t,xtrue(:,1),'r')

[x,t] = lsim(sysKF, [u;y'],t);
plot(t,x(:,1),'k--')



%% simulate 
% [t,yInt] = ode45(@(t,y)cartPendSingle(y, -K*(y-[1;0;pi;0]), par), tInt, y0);

yInt = xtrue;
% determine accelerations, for force on rod 
% u = -repmat(K,length(yInt),1) .* (yInt-repmat([1;0;pi;0]',length(yInt),1));
yDot = cartPendOwn(yInt,u,par);
x = yInt(:,1);
xDot = yInt(:,2); 
thet = yInt(:,3);
thetDot = yInt(:,4);
xDotDot = yDot(2,:)'; 
yDotDot = yDot(4,:)';
momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
Fp = 1./sin(thet) .* (par.mc*xDotDot + par.b*xDot - u(:,2) );
Ep = - par.mp* -cos( thet )* par.g;
Ek = 0.5*par.mc*xDot.^2 + 0.5*par.mp*( (xDot + par.L*cos(thet).*thetDot).^2 + (-sin(thet).*thetDot).^2  );

%% Metrics 
% total energy 
figure( 'Position',[100,100,500,350])
    plot(t, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
% momentum conservation
figure('Position',[600,100,500,350]);
    plot(t, momentum)
    xlabel('Time (s)')
    ylabel('Momentum in X (kg*m/s)')
% Rod tension 
figure( 'Position',[1100,100,500,350] );
    plot(t,Fp)
    xlabel('Time (s)')
    ylabel('Force in rod (N)')
% animate cart
figure( 'Position',[100,550,1000,400] );
for j = 1:10:length(t)
    drawCartPendOwn( yInt(j,:), t(j), Fp(j), par)
end    
