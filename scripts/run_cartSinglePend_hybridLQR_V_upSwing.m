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
R = 1e3;
K = lqr( A,B,Q,R);

par.dist = 0;
% y0 = [0; 0; pi*11/12; 0];
y0 = [0; 0; pi*7/10; 0];
% y0 = [0; 0; pi*4/10; 0]; 
% y0 = [0; 0; pi*0; 0];

yGoal = [0; 0; pi; 0];


%% simulate cart 
% [t,yInt] = ode45(@(t,y)cartSinglePend(y, -K*(y-yGoal), par), tInt, y0);

[t,yInt] = ode45(@(t,y) cartSinglePend(y, hybrid_LQR_V_upswing(y,yGoal,K,par) , par),tInt, y0);
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
 %% 
figure(); hold on
%     plot(tInt,yInt)
%     subplot(311)
    plot(tInt,thet)
%     subplot(312)
    plot(tInt,thetDot)
    plot(tInt,u)
    xlabel('Time (s)')
    ylabel('rad, rad/s, N')
    legend('$\theta$','$\dot{\theta}$','u')
    axis([0,15,-10,40])
%     ylabel('Total Energy T+V (J)')
%     subplot(313)
%     plot(tInt,xDot)


%% create video
if false
     vidfile = VideoWriter( [rootPath filesep 'figs' filesep 'cartSinglePend_hybridLQR_V_upSwing.mp4'],'MPEG-4');
     vidfile.FrameRate = 15;
     open(vidfile);

    figure( 'Position',[100,550,1000,400] );
    for j = 10:10:length(tInt)
        drawCartSinglePend( yInt(j,:), tInt(j), par)
%         scatter( sin(psi_cg(j))*r_cg(j), -cos(psi_cg(j)).*r_cg(j),50,'or','filled')
        drawnow
        hold off 
        drawnow
        im = getframe(gcf); 
        writeVideo(vidfile, im);
    %     pause(0.02)
    end
    close(vidfile)
end

%% 
m = par.mp;
M = par.mc;
L = par.L;
% dimensions
W = 1*sqrt(M/5);  % cart width
H = .5*sqrt(M/5); % cart height
wr = .2; % wheel radius
mr = .3*sqrt(m); % mass radius
par.xl = [-4,4];
    figure( 'Position',[100,550,1000,400] );
for j = 10:5:length(tInt)/3% input wrangling 
x = yInt(j,1);
th = yInt(j,3);
% positions
y = wr/2+H/2; % cart vertical position
w1x = x-.9*W/2;
w1y = 0;
w2x = x+.9*W/2-wr;
w2y = 0;
px = x + L*sin(th);
py = y - L*cos(th);
col = [1,1,1]*0;
plot([-10 10],[0 0],'k','LineWidth',2)
hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[216,218,235]/255)
rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1,1,1]*0.5)
rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1,1,1]*0.5)
% plot([x px],[y py],'k','LineWidth',2)
plot([x px],[y py],'Color',col ,'LineWidth',2)

% j_swing = 
for k = 5:5:j
    x = yInt(k,1);
    th = yInt(k,3);
    px = x + L*sin(th);
    py = y - L*cos(th);
    col_r = [255,50,50,k/j*255/2]/255;
    rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',col_r,'EdgeColor',[0,0,0,(k)/j])
%     rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',col_r,'EdgeColor','none')
end

xl = par.xl;
yl = [-2.5, 2.5];

cartAx = { 'XTick', linspace(xl(1),xl(2),3), ...
            'YTick', linspace(yl(1),yl(2),3), ...
            'XLim', xl, ...
            'YLim', yl, ...
            'DataAspectRatio',[1 1 1],...
            'PlotBoxAspectRatio', [3 4 4],...
            };
set(gca(),cartAx{:})
title( ['Time = ', num2str( j*1e-2), ' (s)']);
    
drawnow
hold off
end
