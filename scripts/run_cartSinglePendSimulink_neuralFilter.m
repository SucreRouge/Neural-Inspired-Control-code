clc;clear all;close all

run('config_cartSinglePend');

% set initial conditions, goal conditions, and activation time
y0 = [0; 0; pi; 0];
tyGoal = 3;
yGoal = [2; 0; pi; 0];

% set time parameters
dt = 0.001; 
tLast = 30;
tInt = dt:dt:tLast;
tSamp = 0.001;

% Noise and filtering parameters
neuralEncodingOn = 1;
sensorNoiseVariance = [1 1 1 1]*0;%1e-4';
distXAccVariance =  0;%1e-4 ;
distThetaAccVariance = 0;%1e-4 ;

G = 3;
% neural encoding parameters 
tSTA = -0.5:tSamp:0;
lenSTA = length(tSTA);
nDelays = lenSTA -2; 
STAw = 0.01;   % width in seconds
STAd = 0;  
STAFunc = @(t) exp(- ((tSTA+STAd)./STAw).^2);
STA = STAFunc( tSTA ) / norm( STAFunc(tSTA), 1 )*G;
delayVal = find( cumsum(STA) >=0.5, 1);
% NLDs = 0;
% NLDg = 1;
NLDs = 0.5;
NLDg = 10;

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
r = 1e-3;
Q = diag(ones(4,1)); 
% Q = diag(zeros(4,1)); 
Q(1,1) = 10 ;
K = lqr( A,B,Q,R);

%% run simulink 
load_system('simulink_CartpendDiscrete_KControl_neuralSensing')
set_param('simulink_CartpendDiscrete_KControl_neuralSensing', 'StopTime', num2str(tLast) )
set_param('simulink_CartpendDiscrete_KControl_neuralSensing', 'MaxStep', num2str( dt ) )
sim('simulink_CartpendDiscrete_KControl_neuralSensing');

%% process simulink output
% state out 
yOut = stateOut.signals.values(:,1:4);
uOut = stateOut.signals.values(:,5);
tOut = stateOut.time;
yInt = interp1(tOut,yOut,tInt);
uInt = interp1(tOut,uOut,tInt);

% sensors out
yActual = sensorsOut.signal1.Data;
yActualT = sensorsOut.signal1.Time;
yActualDiscrete = sensorsOut.yActualDiscrete.Data;
yActualDiscreteT = sensorsOut.yActualDiscrete.Time ;
yMeasured = squeeze(sensorsOut.yMeasured.Data);
yMeasuredT = sensorsOut.yMeasured.Time;

% disturbance and noise out
distX = distNoiseOut.signal1.Data;
distThet = distNoiseOut.signal2.Data;
tDist = distNoiseOut.signal1.Time;
noise = squeeze(distNoiseOut.signal3.Data);
tNoise = distNoiseOut.signal3.Time;

% neural encoding out 
yDelta = squeeze(encodedOut.yDelta.Data);
tDelta = encodedOut.yDelta.Time;
encodedSTA = squeeze(encodedOut.ySTA.Data);
encodedNLD = squeeze(encodedOut.yNLD.Data);
tEncoded = squeeze( encodedOut.ySTA.time );

%% determine accelerations, for force on rod 
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
figure( 'Position',[100,100,500,350]); hold on
    plot(tInt, Ep+Ek)
    xlabel('Time (s)')
    ylabel('Total Energy T+V (J)')
% control input
figure('Position',[600,100,500,350]);
    plot(tInt,uInt)
    xlabel('Time (s)')
    ylabel('Control input u')
%% state actual and sensed
fig01 = figure('Position',[1100,100,500,350]); hold on
    stairs(yMeasuredT, yMeasured,'Color',[1,1,1]*0.7)
    stairs(yActualDiscreteT,yActualDiscrete,'k')
    legVec = plot(yActualT,yActual,'LineWidth',2);
    xlabel('Time (s)')
    ylabel('State')
    legend(legVec,'x','$\dot{x}$','$\theta$','$\dot{\theta}$')
    %
    
    width = 4;
    height = 3;
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% % % Here we preserve the size of the image when we save it.
set(fig01,'InvertHardcopy','on');
set(fig01,'PaperUnits', 'inches');
papersize = get(fig01, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig01, 'PaperPosition', myfiguresize);
% saving of image
print(fig01, ['Figure_neuralState'], '-dpng', '-r300');
% total hack, why does saving to svg scale image up???
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig01, 'PaperPosition', myfiguresize);

print(fig01, ['Figure_neuralState'], '-dsvg', '-r499');
    
    %%
% look at noise and disturbance
figure('Position',[1600,100,500,350]);
    subplot(311)
    plot(tDist,distX,'Color',[1,1,1]*0.7); ylabel('$\ddot{x}$ disturbance')
    subplot(312) 
    plot(tDist,distThet,'Color',[1,1,1]*0.7); ylabel('$\ddot{\theta}$ disturbance')
    subplot(313)
    plot(tNoise,noise,'Color',[1,1,1]*0.7); ylabel('Sensor noise')
% animate cart
figure( 'Position',[100,550,1000,400] );
frameStep = 1/dt/10;
for j = frameStep:frameStep:length(tInt)
    drawCartSinglePendForce( yInt(j,:), tInt(j), Fp(j), par)
end    


%% make video of  [y, x] = ndgrid(1:256);

% fig_anim = figure( 'Position',[100,550,1000,400] );
%  vidfile = VideoWriter('neuralControl1.mp4','MPEG-4');
%  open(vidfile);
% for j = frameStep:frameStep:length(tInt)
% %     z = sin(x*2*pi/ind)+cos(y*2*pi/ind);
% %     im = sc(z, 'hot'); 
%     drawCartPendForce( yInt(j,:), tInt(j), Fp(j), par)
%     frame = getframe(gcf);
%     writeVideo(vidfile, frame);
%  end
% close(vidfile)
% animation


%% Neural encoding effect 


figure('Position',[1500,550,400,400] )
    subplot(211)
    plot(tSTA,STA,'-o');
    hold on
    plot( tSTA(delayVal)*[1,1],[0,max(STA)],'--')
    xlabel('Time (s)'); ylabel('STA')
    subplot(212)
    s = -1:0.01:1;
    plot(s,NLDfun(s,NLDg,NLDs))
    xlabel('Similarity, $\xi$'); ylabel('NLD')
    
   
figure('Position',[1100,550,400,400] ); hold on; 
    p1 = plot(tDelta,yDelta,'Color',[1,1,1]*0.7);
%     plot(yActualDiscreteT,yActualDiscrete,'Color',[1,1,1]*0.7)
    p2 = plot(tEncoded,encodedSTA,'k');
    p3 = plot(tEncoded,encodedNLD,'r');
    legend([p1(1),p2(1),p3(1)],'Sensory input','STA filtered','STA + NLD filtered','Location','Best')
    
    
    
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
for j = 4000:500:length(tInt)% input wrangling 
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
for k = 4000:500:j
    x = yInt(k,1);
    th = yInt(k,3);
    px = x + L*sin(th);
    py = y - L*cos(th);
    col_r = [255,50,50,k/j*255/2]/255;
%     rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[216,218,235,25]/255,'EdgeColor',[0,0,0,0.1])
%     rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1,1,1,0.02]*0.5,'EdgeColor',[0,0,0,0.1])
%     rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1,1,1,0.02]*0.5,'EdgeColor',[0,0,0,0.1])
    % plot([x px],[y py],'k','LineWidth',2)
%     plot([x px],[y py],'Color',col ,'LineWidth',2)
%     plot([x px],[y py],'Color',[0,0,0,0.1] ,'LineWidth',2)
%     rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',col_r,'EdgeColor',[0,0,0,(k)/j])
    rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',col_r,'EdgeColor','none')
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
