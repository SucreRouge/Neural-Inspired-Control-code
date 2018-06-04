function drawCartDoublePend(y, tNow, pendPar)

% input wrangling 
x = y(1);
th1 = y(2);
th2 = y(3);
m2 = pendPar.m2;
m1 = pendPar.m1;
mc = pendPar.mc;
L1 = pendPar.L1;
L2 = pendPar.L2;

% dimensions
W = 1*sqrt(mc/5);  % cart width
H = .5*sqrt(mc/5); % cart height
wr = .2; % wheel radius
mr1 = .3*sqrt(m1); % mass radius
mr2 = .3*sqrt(m2); % mass radius

% positions
y = wr/2+H/2; % cart vertical position
w1x = x-.9*W/2;
w1y = 0;
w2x = x+.9*W/2-wr;
w2y = 0;

px1 = x + L1*sin(th1);
py1 = y - L1*cos(th1);
px2 = px1 + L2*sin(th2);
py2 = py1 - L2*cos(th2);

col = 'k';
plot([-10 10],[0 0],'k','LineWidth',2)
hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[216,218,235]/255)
rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1,1,1]*0.5)
rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1,1,1]*0.5)
% plot([x px],[y py],'k','LineWidth',2)
plot([x px1],[y py1],'Color',col ,'LineWidth',2)
plot([px1 px2],[py1 py2],'Color',col ,'LineWidth',2)
rectangle('Position',[px1-mr1/2,py1-mr1/2,mr1,mr1],'Curvature',1,'FaceColor',[216,218,235]/255)
rectangle('Position',[px2-mr2/2,py2-mr2/2,mr2,mr2],'Curvature',1,'FaceColor',[216,218,235]/255)

% figure makeup 
xl = [-5,5];
yl = [-2.5, 2.5];

% text(0, yl(2)-0.5, ['Force in rod = ',num2str( round(Fp,1) ), ' (N)'])
% text(3, yl(2)-0.5, ['Force eq. does not include', char(10),'torsional damping yet [TM]'])

cartAx = { 'XTick', linspace(xl(1),xl(2),3), ...
            'YTick', linspace(yl(1),yl(2),3), ...
            'XLim', xl, ...
            'YLim', yl, ...
            };
set(gca(),cartAx{:})
title( ['Time = ', num2str( tNow), ' (s)']);
    
drawnow
hold off