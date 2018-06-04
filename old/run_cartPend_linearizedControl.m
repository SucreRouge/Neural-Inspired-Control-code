clc;clear all;close all

mc = 5;
mp = 1;
g = -10;
L = 2;
b = 1;
par.mc = mc ; par.mp = mp; par.g = g; par.L = L; par.b = b; 
par.dist = 0.25;

t0 = 0;
tLast = 10;
tInt = t0:0.02:tLast;
y0 = [-4,0,pi*10/12,0];

%% Linearized eqs pendulum on the cart 
eqTheta = pi;  % which equilibrium? (theta = 0[down], or theta = pi[up] )
s = cos(eqTheta);

A = [ 0  1         0                   0;
    0   -b/mc      -mp*g/mc             0;
    0   0          0                   1;
    0   s*b/(L*mc) s*(mc+mp)*g/(mc*L)  0];
B = [0 ;
    1/mc;
    0
    -s/(mc*L)];

% p = [-3.1, -3.2, -3.3 -3.4];
% p = [-2.1, -2.2, -2.3 -2.4];
p = [-1.1, -1.2, -1.3 -1.4];
% p = [-0.5, -0.6, -0.7 -0.8];

K = place(A, B, p );
% eig(A)
% rank( ctrb(A,B))
% eig(A-B*K)

%% simulate 
[t,yInt] = ode45(@(t,y)cartPendSingle(y, -K*(y-[1;0;pi;0]), par), tInt, y0);

% determine accelerations, for force on rod 
u = -repmat(K,length(yInt),1) .* (yInt-repmat([1;0;pi;0]',length(yInt),1));
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
for j = 1:5:length(tInt)
    drawCartPendOwn( yInt(j,:), tInt(j), Fp(j), par)
end    
