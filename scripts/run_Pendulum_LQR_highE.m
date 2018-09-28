clc;clear all;close all

% Initialize model parameters 
run('config_singlePend');


L = 4*abs(g)/pi^2;
par.L = L

% eqTheta = pi;  % which equilibrium? (theta = 0[down], or theta = pi[up] )
% s = cos(eqTheta);
A = [ 0   1        ;
      0   -b/(mp*L^2)];     
B = [0 ;
    1/(mp*L^2)];

% time parameters
dt = 0.01; 
tLast = 40;
tSamp = 0.1;
tInt = dt:dt:tLast;
%% set control gain
Q = [1,0; 0,1];
R = 1e-3; 
K = lqr( A,B,Q,R);
par.b = 0;

% y0 = [-pi*0.90,0]';
% y0 = [-pi*0.8,3]';
% y0 = [-pi*0.5,0]';
% y0 = [pi*0.85,0]';
y0 = [pi*0.8,3]';
yGoal = [pi,0]';
%% simulate cart 
u = 0;
% [t,yInt] = ode45(@(t,y)singlePend(y',u,par),tInt, y0);
% [t,yInt] = ode45(@(t,y) singlePend(y,-K*(y-yGoal),par), tInt, y0);
[t,yInt] = ode23(@(t,y) singlePend(y, pendulum_hybrid_LQR(y,yGoal,K,par), par), tInt, y0);
%% Plot and animate

thet = yInt(:,1);
thetDot = yInt(:,2);
for j = 1:length(yInt)
    u(j) = pendulum_hybrid_LQR(yInt(j,:)',yGoal,K,par);
end


V = cos(thet)*L*mp*g;
T = 1/2* mp*thetDot.^2*L^2;
E = V+T; 
figure( 'Position',[100,100,500,350])
    hold on
    plot(tInt, V+T,'k')
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
    legend('Total')
    
figure( 'Position',[600,100,500,350])
    subplot(211)
    hold on
   
    plot(tInt, yInt)
    xlabel('Time (s)')
    ylabel('state')
    legend('$\theta$','$\dot{\theta}$')
    subplot(212)
    hold on
    plot(tInt, u)
    xlabel('Time (s)')
    ylabel('u')
%     legend('$\theta$','$\dot{\theta}$')
%% animate cart 
if false
    figure( 'Position',[100,550,1000,400] );
    for j = 10:10:length(tInt)
        drawSinglePend( yInt(j,:), tInt(j), par)
        pause(0.02)
    end
else
     vidfile = VideoWriter( [rootPath filesep 'figs' filesep 'singlePend_hybridLQR_highE.mp4'],'MPEG-4');
     vidfile.FrameRate = 15;
     open(vidfile);

    figure( 'Position',[100,550,1000,400] );
    for j = 10:10:length(tInt)
        subplot(121) 
        drawSinglePend( yInt(j,:), tInt(j), par)
%         scatter( sin(psi_cg(j))*r_cg(j), -cos(psi_cg(j)).*r_cg(j),50,'or','filled')
        subplot(222) 
            plot(tInt(1:j), E(1:j),'k')
    %         axis([0,10, 10,40])
            axis([0,15, 15,40])
            xlabel('Time (s)'); ylabel('Total Energy (J)')
        subplot(224) 
            plot(tInt(1:j), u(1:j),'k')
    %         axis([0,10, 10,40])
            axis([0,15, -5,5])
            xlabel('Time (s)'); ylabel('u')
        
        
        drawnow
        hold off 
        drawnow
        im = getframe(gcf); 
        writeVideo(vidfile, im);
    %     pause(0.02)
    end
    close(vidfile)
end