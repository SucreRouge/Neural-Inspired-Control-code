clc;clear all;close all

% Initialize model parameters 
run('config_cartSinglePend');
u = 0;

y0 = [0,0,pi*0.9,0];
%% simulate cart 
[t,yInt] = ode45(@(t,y)cartSinglePend(y',u,par),tInt, y0);
yDot = cartSinglePend_matrixOperation(yInt,u,par);

%% Plot and animate
x = yInt(:,1);
xDot = yInt(:,2); 
thet = yInt(:,3);
thetDot = yInt(:,4);
xDotDot = yDot(2,:)'; 
yDotDot = yDot(4,:)';
momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
Ep = - par.mp* -cos( thet )*L* par.g;
Ek = 0.5*par.mc*xDot.^2 + 0.5*par.mp*( (xDot + par.L*cos(thet).*thetDot).^2 + (-sin(thet).*L.*thetDot).^2  );
Fp = 1./sin(thet) .* (par.mc*xDotDot + par.bX*xDot - u );

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
plotPar.Frange = [min(Fp),max(Fp)];
for j = 10:10:length(tInt)
    drawCartSinglePendForce( yInt(j,:), tInt(j), Fp(j), par,plotPar)
%     pause(0.02)
end



%% 
if true
     vidfile = VideoWriter( [rootPath filesep 'figs' filesep 'cartSinglePend_noControl.mp4'],'MPEG-4');
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