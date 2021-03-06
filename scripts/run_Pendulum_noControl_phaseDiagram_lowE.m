clc;clear all;close all

% Initialize model parameters 
run('config_singlePend');% physical parameters



% time parameters
dt = 0.01; 
tLast = 40;
tSamp = 0.1;
tInt = dt:dt:tLast;

% 
mp = 1;
g = -10;
% L = 2;
L = 4*abs(g)/pi^2;
b = 0;
% b = 0.0001;
par.mp = mp; par.g = g; par.L = L; par.b = b; 


% eqTheta = pi;  % which equilibrium? (theta = 0[down], or theta = pi[up] )
% s = cos(eqTheta);
A = [ 0   1        ;
      0   -b/(mp*L^2)];     
B = [0 ;
    1/(mp*L^2)];

% set control gain
Q = [10,0; 0,1];
R = 1e-3; 
K = lqr( A,B,Q,R);
par.b = 0;

% y0 = [-pi*0.90,0]';
% % y0 = [-pi*0.8,3]';
% y0 = [-pi*0.5,0]';
% y0 = [pi*0.1,0]';
y0 = [pi*0.01,0]';
% % y0 = [pi*0.85,0]';
% % y0 = [pi*0.8,3]';
yGoal = [pi,0]';

%% simulate cart 
u = 0;
% [t,yInt] = ode45(@(t,y)singlePend(y',u,par),tInt, y0);
% [t,yInt] = ode45(@(t,y) singlePend(y,-K*(y-yGoal),par), tInt, y0);
[t_sim,yInt_sim] = ode23(@(t,y) singlePend(y, pendulum_hybrid_LQR(y,yGoal,K,par), par), tInt, y0);

yInt_sim(:,1) = wrapToPi( yInt_sim(:,1));
for j = 1:length(yInt_sim)
    u_col(j) = pendulum_hybrid_LQR(yInt_sim(j,:),yGoal,K,par); 
end



%% Phase diagram

dt = 0.01; 
tLast = 15;
tSamp = 0.1;
tInt = dt:dt:tLast;

thetZoom = 0.5;

y0 = [  [-0.5;-1.1;-1.8;-2.4]+0.01 ,zeros( 4,1)+0.01 ;
         -pi,0.2;
         ones(4,1)*-pi ,  linspace(1.5,4,4)';
         ];
nT = [210,230,260,330,530,    320,250,200,150,150,100,100,100,50];

for j = 1:length(y0)
        [~,yInt(:,1:2,j) ] = ode45(@(t,y)singlePend(y',u,par),tInt, y0(j,:));
        yInt(1:end-1,3,j) = diff(yInt(:,2,j))/dt;
    
end
fig1 = figure();
    width = 8;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
    
    subplot(121)
cMap = [255,255,229
247,252,185
217,240,163
173,221,142
120,198,121
65,171,93
35,132,67
0,104,55
0,69,41] ;
cMapInterp = colorSchemeInterp( cMap/255, 50);

colormap(cMapInterp)
hold on
mod = [0,-2*pi,2*pi];

%
thetC = linspace(-pi*1.5,pi*1.5,25);
thetDC =  linspace(-pi*1.5,pi*1.5,25);
[X,Y] = meshgrid(thetC,thetDC); 
E = mp*g*(cos(X)-1)*L + 0.5*mp*(L*Y).^2;

surf(X,Y,E-500,'FaceAlpha',0.5)

fill3( [-1,1,1,-1]*thetZoom+pi,[-1,-1,1,1]*thetZoom ,[1,1,1,1]*-200,'w','FaceAlpha',0.4)
view(2)


%%
greyVal = 0.6;
for j = 1:length(y0)
    for k = 1:length(mod)
        
%         plot(yInt(1:nT(j),1,j) + mod(k),yInt(1:nT(j),2,j) ,'k')
%         plot( -yInt(1:nT(j),1,j) + mod(k),-yInt(1:nT(j),2,j) ,'k')
        plot(yInt(1:nT(j),1,j) + mod(k),yInt(1:nT(j),2,j) ,'Color',[1,1,1]*greyVal)
        plot( -yInt(1:nT(j),1,j) + mod(k),-yInt(1:nT(j),2,j) ,'Color',[1,1,1]*greyVal)
    end
end


% plot( yInt_sim(:,1),yInt_sim(:,2) ,'Color',[1,1,1]*0)
% plot( yInt_sim(:,1),yInt_sim(:,2) ,'Color','r','LineWidth',2)
%%
u_pos = u_col == 5;
u_neg = u_col == -5;
u_deg = logical( u_pos+u_neg );
u_0 = u_col == 0;
u_rest = ~(u_deg+u_0);

u_sig = zeros(length(u_pos),3);
u_sig(u_rest,:) = nonzeros(u_rest)*[1,0,0];
u_sig(u_0,:) = nonzeros(u_0)*[0,1,0];
u_sig(u_deg,:) = nonzeros(u_deg)*[0,0,1];
% u_sig = u_sig
% u_sig'*[1,1,1]

scatter( yInt_sim(:,1),yInt_sim(:,2),5,u_sig  ,'fill')


%%

for j = 1:length(y0)
    for k = 1:length(mod)
        for tj = 1%:25:length(tInt)
             pos =  [ yInt(tj,1,j) ,yInt(tj,2,j) ] ;
             phi = atan2( pos(2), pos(1) );
             er = (phi/ (pi/4) ); 
             der =  [  yInt(tj,2,j) ,yInt(tj,3,j) ];
             ang = atan2(der(2),der(1));
             R = [cos(ang),-sin(ang); sin(ang), cos(ang)];
             xyTri = [3, 0;  0, 1;  0, -1]*norm(der)*0.03;
             for jj = 1:3
                xyTurn(jj,:) = (  R * xyTri(jj,:)' )'; 
             end
             fill3( pos(1)+ xyTurn(:,1)+ mod(k) ,pos(2)+ xyTurn(:,2),[1,1,1]  ,'r');
        end
    end
   
end

% regular axes
xAx = [-10,10; 
       [1,1]*-pi; 
       [0,0]; 
       [1,1]*pi;
       ];
yAx = [ 0,0;
       -10,10;
       -10,10;
       -10,10;
       ];
plot(xAx',yAx','k')
axis square
axis( [[-pi,pi]*1.5, -6,6])
ax = gca();
phaseAx = {'XLim',[-1,1]*pi*1.4,'YLim',[-1,1]*pi*1.4,...
            'XTick', -pi:pi/2:pi ,'XTickLabel',{'$-\pi$','$-\frac{\pi}{2}$','0','$\frac{\pi}{2}$', '$\pi$'} ,...
            'YTick', -pi:pi/2:pi ,'YTickLabel',{'$-\pi$','$-\frac{\pi}{2}$','0','$\frac{\pi}{2}$', '$\pi$'} ,...
%             'CLim',[-205,-150],...
        };

xlabel('$\theta$','FontSize',14); 
ylabel('$\dot{ \theta}$','FontSize',14,'Rotation',0)
set(ax,phaseAx{:})











%% 
subplot(122) 
clear('yInt')

dt = 0.01; 
tLast = 5;
tSamp = 0.1;
tInt = dt:dt:tLast;

y0 = [linspace(0,thetZoom-0.05,8)'+pi,-ones(8,1)* (thetZoom-0.05);
        -linspace(0,thetZoom-0.05,8)'+pi,ones(8,1)* (thetZoom-0.05);
%         0.22328+pi,-(thetZoom-0.05);
%         -0.22328+pi,(thetZoom-0.05);
%         0.22+pi,-(thetZoom-0.05);
%         0.23+pi,-(thetZoom-0.05);
%         -0.22+pi,(thetZoom-0.05);
%         -0.23+pi,(thetZoom-0.05);
        ];

for j = 1:length(y0)
        [~,yInt(:,1:2,j) ] = ode45(@(t,y)singlePend(y',u,par),tInt, y0(j,:));
        yInt(1:end-1,3,j) = diff(yInt(:,2,j))/dt;
    
end
%
colormap(cMapInterp)
thetC = linspace(-thetZoom,thetZoom,15)+pi;
thetDC =  linspace(-thetZoom,thetZoom,15);
[X,Y] = meshgrid(thetC,thetDC); 
E = mp*g*(cos(X)-1)*L + 0.5*mp*(L*Y).^2;

surf(X,Y,E-max(E(:))*1.1,'FaceAlpha',0.4)
view(2)

hold on
mod = [0] ;
for j = 1:length(y0)
    for k = 1:length(mod)
%         plot(yInt(:,1,j) + mod(k),yInt(:,2,j) ,'k')
        plot(yInt(:,1,j) + mod(k),yInt(:,2,j) ,'Color',[1,1,1]*greyVal)
    end
   
end
for j = 1:length(y0)
    for k = 1:length(mod)
        for tj = 1%:25:length(tInt)
             pos =  [ yInt(tj,1,j) ,yInt(tj,2,j) ] ;
             phi = atan2( pos(2), pos(1) );
             er = (phi/ (pi/4) ); 
             der =  [  yInt(tj,2,j) ,yInt(tj,3,j) ];
             ang = atan2(der(2),der(1));
             R = [cos(ang),-sin(ang); sin(ang), cos(ang)];
             xyTri = [3, 0;  0, 1;  0, -1]*norm(der)*0.02;
             for jj = 1:3
                xyTurn(jj,:) = (  R * xyTri(jj,:)' )'; 
             end
%              fill( pos(1)+ xyTurn(:,1)+ mod(k) ,pos(2)+ xyTurn(:,2)  ,'r');
             fill3( pos(1)+ xyTurn(:,1)+ mod(k) ,pos(2)+ xyTurn(:,2) ,[1,1,1] ,'r');
        end
    end
   
end

% plot( yInt_sim(:,1),yInt_sim(:,2) ,'Color',[1,1,1]*0)
% plot( yInt_sim(:,1),yInt_sim(:,2) ,'Color','r','LineWidth',2)
scatter( yInt_sim(:,1),yInt_sim(:,2),5,u_sig  ,'fill')


% regular axes
xAx = [-10,10; 
       [1,1]*-pi; 
       [0,0]; 
       [1,1]*pi;
       ];
yAx = [ 0,0;
       -10,10;
       -10,10;
       -10,10;
       ];
plot(xAx',yAx','k')
axis square

axis square
ax = gca();
phaseAx = {'XLim',[-1,1]*thetZoom+pi,'YLim', [-1,1]*thetZoom ,...
            'XTick', [-1,0,1]*thetZoom+pi ,'XTickLabel',{'$\pi-0.5$','$\pi$','$\pi+0.5$'},...
            'YTick', [-1,0,1]*thetZoom,'YTickLabel',[-1,0,1]*thetZoom,...
        };
xlabel('$\theta$','FontSize',14); 
ylabel('$\dot{ \theta}$','FontSize',14,'Rotation',0)
set(ax,phaseAx{:})

%% 
if 0 
    set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
    % % % Here we preserve the size of the image when we save it.
    set(fig1,'InvertHardcopy','on');
    set(fig1,'PaperUnits', 'inches');
    papersize = get(fig1, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(fig1, 'PaperPosition', myfiguresize);
    % saving of image
    print(fig1, [rootPath filesep 'figs' filesep 'Figure_phasePortrait_lowE'], '-dpng', '-r300');
    % total hack, why does saving to svg scale image up???
    stupid_ratio = 15/16;
    myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
    set(fig1, 'PaperPosition', myfiguresize);

    print(fig1, [rootPath filesep 'figs' filesep 'Figure_phasePortrait_lowE'], '-dsvg', '-r499');
end

