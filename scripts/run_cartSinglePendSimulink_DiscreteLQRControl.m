clc;clear all;close all

run('config_cartSinglePend');
%% Initialize model parameters 
% 
tSamp = 0.1;
% tSTA = -40:1:0;
% STAw = [5];
% STAf = 1;
% STAd = 20;
% % STAFunc = @(t) cos( STAf*(t+ STAd)  ).*exp(-(t+STAd).^2 / STAw.^2);
% STAFunc = @(t) exp(-(t+STAd).^2 / STAw.^2);
% STA = STAFunc(tSTA);

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
K = lqr( A,B,Q,R);

%% run simulink 
load_system('simulink_CartpendDiscrete_KControl')
set_param('simulink_CartpendDiscrete_KControl', 'StopTime', num2str(tLast) )
set_param('simulink_CartpendDiscrete_KControl', 'MaxStep', num2str( dt ) )
sim('simulink_CartpendDiscrete_KControl');

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

% sensors out
yActual = sensorsOut.signal1.Data;
yActualT = sensorsOut.signal1.Time;
yActualDiscrete = sensorsOut.yActualDiscrete.Data;
yActualDiscreteT = sensorsOut.yActualDiscrete.Time ;
yMeasured = squeeze(sensorsOut.yMeasured.Data);
yMeasuredT = sensorsOut.yMeasured.Time;

distX = distAndNoise.signal1.Data;
tX = distAndNoise.signal1.Time;
distThet = distAndNoise.signal2.Data;
tThet = distAndNoise.signal2.Time;
noise = squeeze(distAndNoise.signal3.Data);
tNoise = distAndNoise.signal3.Time;
%% Display results
% total energy 
figure( 'Position',[100,100,500,350])
    plot(tInt, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
% control input
    figure('Position',[600,100,500,350]);
    plot(tInt,uInt)
    xlabel('Time (s)')
    ylabel('Control input u')
% state actual and sensed
figure('Position',[1100,100,500,350]);
    plot(yActualT,yActual,'LineWidth',2)
    hold on
    stairs(yActualDiscreteT,yActualDiscrete,'k')
    hold on
    stairs(yMeasuredT, yMeasured,'r')
    xlabel('Time (s)')
    ylabel('State')
    legend('x','$\dot{x}$','$\theta$','$\dot{\theta}$')
% look at noise and disturbance
figure('Position',[1600,100,500,350]);
    subplot(311)
    plot(tX,distX);
    subplot(312) 
    plot(tThet,distThet)
    subplot(313)
    plot(tNoise,noise)

% animate cart
figure( 'Position',[100,550,1000,400] );
frameStep = 1/dt/10;
for j = frameStep:frameStep:length(tInt)
    drawCartSinglePendForce( yInt(j,:), tInt(j), Fp(j), par)
end    


%% make video of  [y, x] = ndgrid(1:256);

% % % % % fig_anim = figure( 'Position',[100,550,1000,400] );
% % % % %  vidfile = VideoWriter('testmovie.mp4','MPEG-4');
% % % % %  open(vidfile);
% % % % % for j = frameStep:frameStep:length(tInt)
% % % % % %     z = sin(x*2*pi/ind)+cos(y*2*pi/ind);
% % % % % %     im = sc(z, 'hot'); 
% % % % %     drawCartPendForce( yInt(j,:), tInt(j), Fp(j), par)
% % % % %     frame = getframe(gcf);
% % % % %     writeVideo(vidfile, frame);
% % % % %  end
% % % % % close(vidfile)
% % % % % animation



% figure( 'Position',[100,550,1000,400] );
% frameStep = 1/dt/10;
% 
% % drawCartPendForce( yInt(1,:), tInt(1), Fp(1), par)
% for j = frameStep:frameStep:length(tInt)
%     drawCartPendForce( yInt(j,:), tInt(j), Fp(j), par)
% end    

% % axis tight
% % set(gca,'nextplot','replacechildren','visible','off')
% set(gca,'nextplot','replacechildren')
% f = getframe;
% % [im,map] = rgb2ind(f.cdata,256,'nodither');
% % im(1,1,1,20) = 0;
% for j = frameStep:frameStep:length(tInt)
%    drawCartPendForce( yInt(j,:), tInt(j), Fp(j), par)
% %   surf(cos(2*pi*k/20)*Z,Z)
%   f = getframe;
% %   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
% end
% % imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf) %g443800
% 
% 



%% Conv check

% hold on
% stairs(tInt,uInt,'k')
% %% 
% sig = squeeze(signal2.signals.values);
% tS = signal2.time;
% res = squeeze(signal1.signals.values);
% tR = signal1.time;
% % tC = convResult3.time;
% filt = squeeze(filter.signals.values);
% 
% figure();
% subplot(211)
% plot(tS,sig,'k')
% % figure();
% hold on
% plot(tR,res,'r')
% subplot(212)
% plot(filt(1,:))
% % legend('sig','sig1','filt')