clc;clear all;close all
run('config_singlePend');% physical parameters

% Set naming/saving parameters
save_figure = false;
play_animation = true; 
save_animation = false; 
animation_name = 'Animation_phasePortrain_highE';

% % % % low E 
% figure_name = 'Figure_phasePortrait_lowE';
y0 = [ -pi*0.1,0]';

% % % high E 
% figure_name = 'Figure_phasePortrait_highE';
% y0 = [-pi*0.8,3]';
% dist = 0.01;
dist = 0.002;
%% 
% time parameters for contours 
dt = 0.01; 
tLast = 25;
tInt = 0:dt:tLast;
thetZoom = 0.4;
%% % simulate cart 

load_system('simulink_singlePend_hybridLQR')
set_param('simulink_singlePend_hybridLQR', 'StopTime', num2str(tLast) )
set_param('simulink_singlePend_hybridLQR', 'MaxStep', num2str( dt ) )
sim('simulink_singlePend_hybridLQR');


% postprocessing
yout = simout.signals.values(:,1:2);
uout = simout.signals.values(:,3);
tout = simout.time;
yInt_sim = interp1(tout,yout,tInt);
u_Int_sim = interp1(tout,uout,tInt);

% wrap theta to -2pi to 2pi domain 
yInt_sim(:,1) = wrapToPi( yInt_sim(:,1) /2)*2;




%% 
% time parameters for contours 
dt = 0.01; 
tLast = 15;
tInt = dt:dt:tLast;

% determine boolean indices for control 
u_pos = u_Int_sim == 5;
u_neg = u_Int_sim == -5;
u_deg = logical( u_pos+u_neg );
u_0 = u_Int_sim == 0;
u_rest = ~(u_deg+u_0);

% determine colormap for control signal plot
controlCol = [
50,136,189
255,255,51
213,62,79
]/255;

% make color matrix for control vector 
u_sig = zeros(length(u_pos),3);
u_sig(u_deg,:) = nonzeros(u_deg) *controlCol(1,:);
u_sig(u_0,:) = nonzeros(u_0) *controlCol(2,:);
u_sig(u_rest,:) = nonzeros(u_rest) *controlCol(3,:) ;





% create colormap to use 
cMap = [255,255,229
    247,252,185
    217,240,163
    173,221,142
    120,198,121
    65,171,93
    35,132,67
    0,104,55
    0,69,41] ;

% create grid axes 
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
   
   
% parameters for phase diagram 
greyVal = 0.6;
edgC = 'none';
arrowCol = [1,1,1]*greyVal;

% square regions
LQR_region = [0.2,0.2];
thetZoom = 0.4;


%% make figure 
fig1 = figure();
    width = 8;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
    subplot(121)
    hold on
    
cMapInterp = colorSchemeInterp( cMap/255, 50);
colormap(cMapInterp)
mod = [0,-2*pi,2*pi];

% plot energy landscape 
thetC = linspace(-pi*1.5,pi*1.5,25);
thetDC =  linspace(-pi*1.5,pi*1.5,25);
[X,Y] = meshgrid(thetC,thetDC); 
E = mp*g*(cos(X)-1)*L + 0.5*mp*(L*Y).^2;
surf(X,Y,E-500,'FaceAlpha',0.5)

% plot zoomed in box 
fill3( [-1,1,1,-1]*thetZoom+pi,[-1,-1,1,1]*thetZoom ,[1,1,1,1]*-200,'w','FaceAlpha',0)
% plot LQR region 
fill3( [-1,1,1,-1]*LQR_region(1)+pi,[-1,-1,1,1]*LQR_region(1) ,[1,1,1,1]*-3,'w','FaceAlpha',1)

% draw gid at 0,pi
plot3(xAx',yAx',ones(size(xAx))*-2,'k')

% set starting positions for contour lines 
y0_contour = [  [-0.5;-1.1;-1.8;-2.4]+0.01 ,zeros( 4,1)+0.01 ;
         -pi,0.2;
         ones(4,1)*-pi ,  linspace(1.5,4,4)';
         ];
y0_contour  = [y0_contour; -y0_contour];
nT = [210,230,260,330,530,    320,250,200,150];
nT = [nT,nT];
% plot 
for j = 1:length(y0_contour)
    [~,yIntC1(:,1:2,j) ] = ode45(@(t,y)singlePend(y',0,par),tInt, y0_contour(j,:));
    yIntC1(1:end-1,3,j) = diff(yIntC1(:,2,j))/dt;
    for k = 1:length(mod)
        % plot contours
        plot3(yIntC1(1:nT(j),1,j) + mod(k),yIntC1(1:nT(j),2,j) , ones(size(yIntC1(1:nT(j),2,j)))*0,'Color',[1,1,1]*greyVal)
        % plot arrows 
         for tk = 1:500
             radial = ((abs(  abs(yIntC1(tk,1,j)) - abs(yIntC1(tk,2,j)) ) < 0.005 ) && (abs(yIntC1(tk,1,j)) < pi) );
             radial23 = ((abs(  abs(yIntC1(tk,1,j)*2) - abs(yIntC1(tk,2,j)) ) < 0.08 ) && (abs(yIntC1(tk,2,j)) > pi) );
             x0Line = (abs(yIntC1(tk,1,j))  < 0.008 );
             x0Line2 = (abs(yIntC1(tk,1,j)-pi)  < 0.008);
             y0Line = (abs(yIntC1(tk,2,j))  < 0.008 );
             if radial || x0Line || x0Line2 || y0Line ||radial23
                
                 pos =  [ yIntC1(tk,1,j) ,yIntC1(tk,2,j) ] ;
                 phi = atan2( pos(2), pos(1) );
                 er = (phi/ (pi/4) ); 
                 der =  [  yIntC1(tk,2,j) ,yIntC1(tk,3,j) ];
                 ang = atan2(der(2),der(1));
                 R = [cos(ang),-sin(ang); sin(ang), cos(ang)];
%                  xyTri = [3, 0;  0, 1;  0, -1]*norm(der)*0.03;
%                  xyTri = [3, 0;  0, 1;  0, -1]*0.1;
%                  xyTri = [0, 0;  -3, 1;  -3, -1]*0.1;
                  xyTri = [1.5, 0;  -1.5, 1;  -1.5, -1]*0.1;
                 for jj = 1:3
                    xyTurn(jj,:) = (  R * xyTri(jj,:)' )'; 
                 end
                 fill3( pos(1)+ xyTurn(:,1)+ mod(k) ,pos(2)+ xyTurn(:,2),-[1,1,1]  ,[1,1,1]*greyVal, 'EdgeColor',edgC);
             end
         end
    end
end

% plot trajectory 
scatter3( yInt_sim(:,1),yInt_sim(:,2),ones(size(yInt_sim(:,2)))*2,5,u_sig  ,'fill')
% plot initial/goal conditions 
scatter3( y0(1), y0(2) ,3,10, 'k'  ,'fill')
scatter3( yGoal(1), yGoal(2),3 ,10, 'k'  ,'fill')

% subplot axes limits, ticks 
view(2) % set view from top (surf makes 3d plot) 
axis square
ax = gca();
phaseAx = {'XLim',[-1,1]*pi*1.4,'YLim',[-1,1]*pi*1.5,...
            'XTick', -pi:pi/2:pi ,'XTickLabel',{'$-\pi$','$-\frac{\pi}{2}$','0','$\frac{\pi}{2}$', '$\pi$'} ,...
            'YTick', -pi:pi/2:pi ,'YTickLabel',{'$-\pi$','$-\frac{\pi}{2}$','0','$\frac{\pi}{2}$', '$\pi$'} ,...
        };
xlabel('$\theta$','FontSize',14); 
ylabel('$\dot{ \theta}$','FontSize',14,'Rotation',0)
set(ax,phaseAx{:})

%% subplot (122)
subplot(122) 
hold on

% set samp time small for arrows
dt = 0.01; 
tLast = 5;
tSamp = 0.1;
tInt = dt:dt:tLast;
nLines = 15;

% set starting positions of contours 
y0_contourP2 = [linspace(0,thetZoom*2,nLines )'+pi,-ones(nLines ,1)*thetZoom;
     - linspace(0,thetZoom*2,nLines )'+pi,ones(nLines ,1)*thetZoom;
        -0.262+pi,thetZoom;
        0.262+pi,-thetZoom;
        ];
nT = ones(size(y0_contourP2))*500;

% plot energy elevation
colormap(cMapInterp)
thetC = linspace(-thetZoom,thetZoom,15)+pi;
thetDC =  linspace(-thetZoom,thetZoom,15);
[X,Y] = meshgrid(thetC,thetDC); 
E = mp*g*(cos(X)-1)*L + 0.5*mp*(L*Y).^2;
surf(X,Y,E-max(E(:))*1.1,'FaceAlpha',0.4)

% plot LQR region
fill3( [-1,1,1,-1]*LQR_region(1)+pi,[-1,-1,1,1]*LQR_region(1) ,[1,1,1,1]*1,'w','FaceAlpha',0.4)

plot3(xAx',yAx',ones(size(xAx))*-2,'k')

% plot contours , arrows 
mod = [0] ;
for j = 1:length(y0_contourP2)
    [~,yIntC2(:,1:2,j) ] = ode45(@(t,y)singlePend(y',u,par),tInt, y0_contourP2(j,:));
    yIntC2(1:end-1,3,j) = diff(yIntC2(:,2,j))/dt;
    for k = 1:length(mod)
        % plot contour lines 
        plot(yIntC2(1:nT(j),1,j) + mod(k),yIntC2(1:nT(j),2,j) ,'Color',[1,1,1]*greyVal)
        
         for tk = 1:500
             radial1 = (abs(  abs(yIntC2(tk,1,j)-pi).^2 + abs(yIntC2(tk,2,j)).^2 - 0.2^2 ) < 0.001 );
             radial23 = (abs(  abs(yIntC2(tk,1,j)-pi).^2 + abs(yIntC2(tk,2,j)).^2 - 0.4^2) < 0.001 );
             x0Line = 0;%(abs(yIntC2(tk,1,j))  < 0.01 );
             x0Line2 = (abs(yIntC2(tk,1,j)-pi)  < 0.001);
             y0Line = (abs(yIntC2(tk,2,j))  < 0.001 );
             if radial1 || x0Line || x0Line2 || y0Line ||radial23
                
                 pos =  [ yIntC2(tk,1,j) ,yIntC2(tk,2,j) ] ;
                 phi = atan2( pos(2), pos(1) );
                 er = (phi/ (pi/4) ); 
                 der =  [  yIntC2(tk,2,j) ,yIntC2(tk,3,j) ];
                 ang = atan2(der(2),der(1));
                 R = [cos(ang),-sin(ang); sin(ang), cos(ang)];
%                  xyTri = [3, 0;  0, 1;  0, -1]*norm(der)*0.03;
%                  xyTri = [3, 0;  0, 1;  0, -1]*0.01;
%                  xyTri = [0, 0;  -3, 1;  -3, -1]*0.01;
                 xyTri = [1.5, 0;  -1.5, 1;  -1.5, -1]*0.01;
                 for jj = 1:3
                    xyTurn(jj,:) = (  R * xyTri(jj,:)' )'; 
                 end
                 fill3( pos(1)+ xyTurn(:,1)+ mod(k) ,pos(2)+ xyTurn(:,2),-[1,1,1]  ,[1,1,1]*greyVal, 'EdgeColor',edgC);
             end
         end
    end
end

% plot actual path 
scatter3( yInt_sim(:,1),yInt_sim(:,2),ones(size(yInt_sim(:,2)))*2,5,u_sig  ,'fill')

% plot initial/goal conditions 
scatter3( y0(1), y0(2),3 ,30, 'k'  ,'fill')
scatter3( yGoal(1), yGoal(2) ,3, 30, 'k'  ,'fill')

view(2)
axis square
ax = gca();
phaseAx = {'XLim',[-1,1]*thetZoom+pi,'YLim', [-1,1]*thetZoom ,...
            'XTick', [-1,0,1]*thetZoom+pi ,'XTickLabel',{'$\pi-0.4$','$\pi$','$\pi+0.4$'},...
            'YTick', [-1,0,1]*thetZoom,'YTickLabel',[-1,0,1]*thetZoom,...
        };
xlabel('$\theta$','FontSize',14); 
ylabel('$\dot{ \theta}$','FontSize',14,'Rotation',0)
set(ax,phaseAx{:})

%% 
% if true
if save_figure
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
    print(fig1, [rootPath filesep 'figs' filesep figure_name ], '-dpng', '-r300');
    % total hack, why does saving to svg scale image up???
    stupid_ratio = 15/16;
    myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
    set(fig1, 'PaperPosition', myfiguresize);

    print(fig1, [rootPath filesep 'figs' filesep figure_name ], '-dsvg', '-r499');
end
%% animate cart 
% if true
if play_animation
%     par.L = 2;
%     figure( 'Position',[100,550,1000,400] );
%     for j = 10:10:length(t_sim)
%         drawSinglePend( yInt_sim(j,:), t_sim(j), par)
%         pause(0.01)
%     end
end
% if save_animation
%     par.L = 2;
%     vidfile = VideoWriter( [rootPath filesep 'figs' filesep animation_name],'MPEG-4');
%      vidfile.FrameRate = 15;
%      open(vidfile);
% 
%     figure( 'Position',[100,550,1000,400] );
%     for j = 10:10:length(tInt)
%         subplot(121) 
%         drawSinglePend( yInt(j,:), tInt(j), par)
% %         scatter( sin(psi_cg(j))*r_cg(j), -cos(psi_cg(j)).*r_cg(j),50,'or','filled')
%         subplot(222) 
%             plot(tInt(1:j), E(1:j),'k')
%     %         axis([0,10, 10,40])
%             axis([0,15, 15,40])
%             xlabel('Time (s)'); ylabel('Total Energy (J)')
%         subplot(224) 
%             plot(tInt(1:j), u(1:j),'k')
%     %         axis([0,10, 10,40])
%             axis([0,15, -5,5])
%             xlabel('Time (s)'); ylabel('u')
%         
%         
%         drawnow
%         hold off 
%         drawnow
%         im = getframe(gcf); 
%         writeVideo(vidfile, im);
%     %     pause(0.02)
%     end
%     close(vidfile)
% end