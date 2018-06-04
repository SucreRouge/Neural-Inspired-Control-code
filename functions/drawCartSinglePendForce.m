function drawCartPendForce(y, tNow, Fp, pendPar, plotPar)

% input wrangling 
x = y(1);
th = y(3);
m = pendPar.mp;
M = pendPar.mc;
L = pendPar.L;

% dimensions
W = 1*sqrt(M/5);  % cart width
H = .5*sqrt(M/5); % cart height
wr = .2; % wheel radius
mr = .3*sqrt(m); % mass radius

% positions
y = wr/2+H/2; % cart vertical position
w1x = x-.9*W/2;
w1y = 0;
w2x = x+.9*W/2-wr;
w2y = 0;
px = x + L*sin(th);
py = y - L*cos(th);

if Fp/abs(pendPar.g) > 0
    col = [1, 0, 0, abs(Fp/pendPar.g)] ;
elseif Fp/abs(pendPar.g) < 0
    col = [0, 0,1, abs(Fp/pendPar.g)];
else 
    col = [1,1,1,1];
end
if max(col(4))>1
    col(4) = col(4)/max(col(4));
end

plot([-10 10],[0 0],'k','LineWidth',2)
hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[216,218,235]/255)
rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1,1,1]*0.5)
rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1,1,1]*0.5)
% plot([x px],[y py],'k','LineWidth',2)
plot([x px],[y py],'Color',col ,'LineWidth',2)
rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[216,218,235]/255)

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