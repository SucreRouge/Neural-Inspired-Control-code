clc;clear all;close all
run('config_singlePend');% physical parameters

% Set naming/saving parameters
save_figure = false;
play_animation = true; 
save_animation = true; 

% animation_name = 'Animation_phasePortrain_lowE';
% y0 = [ -pi*0.15,0]';

animation_name = 'Animation_phasePortrain_highE';
y0 = [-pi*0.8,3.6]';


anim_dT = 20; 

% dist = 0.005;
% dist = 0.008;
dist = 0.002;
%% 
% time parameters for contours 
dt = 0.01; 
tLast = 25;
t_Int_sim = 0:dt:tLast;
thetZoom = 0.4;
%% % simulate cart 

load_system('simulink_singlePend_hybridLQR')
set_param('simulink_singlePend_hybridLQR', 'StopTime', num2str(tLast) )
set_param('simulink_singlePend_hybridLQR', 'MaxStep', num2str( dt ) )
sim('simulink_singlePend_hybridLQR');

% postprocessing

yout = simout.signals.values(:,1:2);
uout = simout.signals.values(:,3);
ucat = simout.signals.values(:,4);
tout = simout.time;
yInt_sim = interp1(tout,yout,t_Int_sim);
u_Int_sim = interp1(tout,uout,t_Int_sim);
u_cat = interp1(tout,ucat,t_Int_sim,'nearest');

% wrap theta to -2pi to 2pi domain 
yInt_sim(:,1) = wrapToPi( yInt_sim(:,1) /2)*2;

%% 

% determine boolean indices for control 
u_bang = (u_cat == 1) | (u_cat == 2) ;
u_0 = (u_cat == 0);
u_LQR  = (u_cat == 3);

% determine colormap for control signal plot
controlCol = [
50,136,189
255,255,51
213,62,79
]/255;

% make color matrix for control vector 
u_col = zeros(length(u_bang),3);
u_col(u_bang,:) = nonzeros(u_bang) *controlCol(1,:);
u_col(u_0,:) = nonzeros(u_0) *controlCol(2,:);
u_col(u_LQR,:) = nonzeros(u_LQR) *controlCol(3,:) ;

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
    height = 7;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-height*50 width*100, height*100]); %<- Set size

% -------------------------------------------------------------------------
subplot(223);     hold on
    
cMapInterp = colorSchemeInterp( cMap/255, 50);
colormap(cMapInterp)
mod = [0,-2*pi,2*pi];

% plot energy landscape 
thetC = linspace(-pi*2,pi*2,25);
thetDC =  linspace(-pi*2,pi*2,25);
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
y0_contour = [  [-0.5;-1.1;-1.8;-2.4]+0.01 ,zeros( 4,1)-0.1 ;
         -pi,0.2;
         ones(4,1)*-pi ,  linspace(1.5,4,4)';
         ones(2,1)*-pi ,  linspace(5,6,2)';
         ];
% time parameters for contours 
t_cont = 0.01:0.01:15;
% do mirrored (negative) contours
y0_contour  = [y0_contour; -y0_contour];
nT = [210,230,260,330,530,    320,250,200,150,150,150];
nT = [nT,nT];


Xc = linspace(-2*pi,2*pi,50);
Yc =  zeros(1,50); 
Xc(2,:) = zeros(1,50);
Yc(2,:) = linspace(-2*pi,2*pi,50);
Xc(3,:) = Xc(2,:)-pi;
Yc(3,:) = Yc(2,:);
Xc(4,:) = Xc(2,:)+pi;
Yc(4,:) = Yc(2,:);
Xc(5,:) = linspace(-2*pi,2*pi,50) + pi;
Yc(5,:) = zeros(1,50);
Xc(6,:) = linspace(-1,1,50)*sqrt(2)/2*pi ;
Yc(6,:) = Xc(6,:) ;
Xc(7,:) = Xc(6,:) ;
Yc(7,:) = -Yc(6,:);
Xc(8,:) = linspace(1,3,50)*sin(0.4)*pi ;
Yc(8,:) = linspace(1,3,50)*cos(0.4)*pi ;
Xc(9,:) = Xc(8,:); 
Yc(9,:) = -Yc(8,:);
Xc(10,:) = -Xc(8,:); 
Yc(10,:) = Yc(8,:);
Xc(11,:) = -Xc(8,:); 
Yc(11,:) = -Yc(8,:);

% plot 
for j = 1:length(y0_contour)
    [~,yIntC1(:,1:2,j) ] = ode45(@(t,y)singlePend(y',0,par),t_cont, y0_contour(j,:));
    yIntC1(1:end-1,3,j) = diff(yIntC1(:,2,j))/dt;
    for k = 1:length(mod)
        % plot contours
        plot3(yIntC1(1:nT(j),1,j) + mod(k),yIntC1(1:nT(j),2,j) , ones(size(yIntC1(1:nT(j),2,j)))*0,'Color',[1,1,1]*greyVal)

        for jj = 2:size(Xc,1)
            [xout,yout,xper,yper] = intersections(Xc(jj,:),Yc(jj,:),yIntC1(:,1,j)' , yIntC1(:,2,j)',1);
            if ~isempty(xout)
                for kk = 1:length(xout)
                pos =  [xout(kk) + mod(k),yout(kk) ];
                der =  [  yIntC1(round(yper(kk)),2,j) ,yIntC1(round(yper(kk)),3,j) ];
                     ang = atan2(der(2),der(1));
                     R = [cos(ang),-sin(ang); sin(ang), cos(ang)];
                     xyTri = [1.5, 0;  -1.5, 1;  -1.5, -1]*0.1;
                     for jj = 1:3
                        xyTurn(jj,:) = (  R * xyTri(jj,:)' )'; 
                     end
                     fill3( pos(1)+ xyTurn(:,1) ,pos(2)+ xyTurn(:,2),-[1,1,1]  ,[1,1,1]*greyVal, 'EdgeColor',edgC);
                end
            end
        end        
    end
end

% plot initial/goal conditions 
scatter3( y0(1), y0(2) ,3,10, 'k'  ,'fill')
% subplot axes limits, ticks 
view(2) % set view from top (surf makes 3d plot) 
axis square
ax = gca();
phaseAx = {'XLim',[-1,1]*pi*2,'YLim',[-1,1]*pi*2,...
            'XTick', -pi:pi/2:pi ,'XTickLabel',{'$-\pi$','$-\frac{\pi}{2}$','0','$\frac{\pi}{2}$', '$\pi$'} ,...
            'YTick', -pi:pi/2:pi ,'YTickLabel',{'$-\pi$','$-\frac{\pi}{2}$','0','$\frac{\pi}{2}$', '$\pi$'} ,...
        };
xlabel('$\theta$','FontSize',14); 
ylabel('$\dot{ \theta}$','FontSize',14,'Rotation',0)
set(ax,phaseAx{:})

% -------------------------------------------------------------------------
subplot(224); hold on

% set samp time small for arrows
t_arr = 0.01:0.01:5;
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
Xc = linspace(-0.2,0.2,50) + pi;
Yc =  real( (0.2^2 - (Xc(1,:)-pi).^2 ).^0.5 ); 
Xc(2,:) = Xc(1,:);
Yc(2,:) = -Yc(1,:);
Xc(3,:) = linspace(-0.4,0.4,50) + pi;
Yc(3,:) = real( (0.4^2 - (Xc(3,:) -pi).^2 ).^0.5 );
Xc(4,:) = Xc(3,:);
Yc(4,:) = -Yc(3,:);
Xc(5,:) = linspace(-0.4,0.4,50) + pi;
Yc(5,:) = zeros(1,50);
Xc(6,:) = ones(1,50)*pi;
Yc(6,:) = linspace(-0.4,0.4,50);

for j = 1:length(y0_contourP2)
    [~,yIntC2(:,1:2,j) ] = ode45(@(t,y)singlePend(y',u,par),t_arr, y0_contourP2(j,:));
    % do the derivative to decide arrow direction
    yIntC2(1:end-1,3,j) = diff(yIntC2(:,2,j))/dt;
    
    % plot the contours
    plot(yIntC2(1:nT(j),1,j) ,yIntC2(1:nT(j),2,j) ,'Color',[1,1,1]*greyVal)

    for jj = 1:size(Xc,1)
        [xout,yout,xper,yper] = intersections(Xc(jj,:),Yc(jj,:),yIntC2(:,1,j)' , yIntC2(:,2,j)',1);
%         scatter(xout,yout,'or')
    %     scatter(  yIntC2( round(yper) ,1,j)' , yIntC2( round(yper) ,2,j)' ,'ok')
        if ~isempty(xout)
            for k = 1:length(xout)
            pos =  [xout(k),yout(k) ];
            der =  [  yIntC2(round(yper(k)),2,j) ,yIntC2(round(yper(k)),3,j) ];
                 ang = atan2(der(2),der(1));
                 R = [cos(ang),-sin(ang); sin(ang), cos(ang)];
                 xyTri = [1.5, 0;  -1.5, 1;  -1.5, -1]*0.01;
                 for jj = 1:3
                    xyTurn(jj,:) = (  R * xyTri(jj,:)' )'; 
                 end
                 fill3( pos(1)+ xyTurn(:,1) ,pos(2)+ xyTurn(:,2),-[1,1,1]  ,[1,1,1]*greyVal, 'EdgeColor',edgC);
            end
        end
    end
end
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

subplot(222); hold on 
	 t_vec = [1,2];
            a1= area( t_Int_sim( t_vec) , u_Int_sim(t_vec) ,'FaceColor', controlCol(1,:),'FaceAlpha',0.5, 'EdgeColor', 'none' );
            a2= area( t_Int_sim( t_vec) , u_Int_sim(t_vec) ,'FaceColor', controlCol(2,:),'FaceAlpha',0.5, 'EdgeColor', 'none' );
            a3= area( t_Int_sim( t_vec) , u_Int_sim(t_vec) ,'FaceColor', controlCol(3,:),'FaceAlpha',0.5, 'EdgeColor', 'none' );

            max_u = max(abs(u_Int_sim) ); 
            axis([0,t_Int_sim(end), [-1,1]*max_u*1.2 ])   
            xlabel('Time (s)'); ylabel('u (Nm)')
            legend([a1,a2,a3],{'BangBang control','Control off','LQR control'},'Location','SouthWest')

%% 
if save_animation
     vidfile = VideoWriter( [rootPath filesep 'figs' filesep animation_name],'MPEG-4');
     vidfile.FrameRate = 15;
     open(vidfile);
    for t_p = (anim_dT+1):anim_dT:length(t_Int_sim)
        t_vec = [(t_p-anim_dT),t_p];
        subplot(221)
            drawSinglePend_phase( yInt_sim(t_p,:), t_Int_sim(t_p), par)
        subplot(222); hold on 
            if u_cat(t_vec(1)) == u_cat(t_vec(2))
                area( t_Int_sim( t_vec) , u_Int_sim(t_vec) ,'FaceColor', u_col(t_vec(1),:),'FaceAlpha',0.5, 'EdgeColor', 'none' )
                plot( t_Int_sim( t_vec) , u_Int_sim(t_vec) , 'Color', u_col(t_vec(1),:) )
                  legend([a1,a2,a3],{'BangBang control','Control off','LQR control'},'Location','SouthWest')
            end
        subplot(223); hold on 
            if abs( diff( yInt_sim(  t_vec,1) )) < 5
                plot3( yInt_sim(  t_vec,1),yInt_sim(  t_vec,2), ones(size(yInt_sim( t_vec,2)))*2 ,'Color',u_col(t_p,:)  ,'LineWidth',3)
            end
        subplot(224) 
            plot3( yInt_sim(  t_vec,1),yInt_sim(  t_vec,2), ones(size(yInt_sim( t_vec,2)))*2 ,'Color',u_col(t_p,:)  ,'LineWidth',3)
            im = getframe(gcf); 
            writeVideo(vidfile, im);
    end
    close(vidfile)
end

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