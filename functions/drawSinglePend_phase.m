function drawSinglePend_phase(y, tNow, pendPar)

% input wrangling 
% x = y(1);
th = y(1);
m = pendPar.mp;
% M = pendPar.mc;
L = pendPar.L;

% dimensions
mr = 1*sqrt(m); % mass radius
% % 
px = L*sin(th);
py = -L*cos(th);
x = 0;
y = 0;
% % col = [1,1,1]*0;

plot([-10 10],[0 0],'k','LineWidth',2)
hold on
plot([x px],[y py],'Color','k' ,'LineWidth',2)
rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[253,174,107]/255)

% figure makeup 
xl = [-5,5];
yl = [-5, 5];


cartAx = { 'XTick', linspace(xl(1),xl(2),3), ...
            'YTick', linspace(yl(1),yl(2),3), ...
            'XLim', xl, ...
            'YLim', yl, ...
            'DataAspectRatio',[1 1 1],...
            'PlotBoxAspectRatio', [3 4 4],...
            };
set(gca(),cartAx{:})
title( ['Time = ', num2str( tNow), ' (s)']);
    
drawnow
hold off