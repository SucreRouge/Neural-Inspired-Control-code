% determine accelerations, for force on rod 
display('lala')
t = tout;
y = yout;
u = -repmat(K,length(yInt),1) .* (yInt-repmat([1;0;pi;0]',length(yInt),1));
yDot = cartPendOwn(yInt,u,par);
x = yInt(:,1);
xDot = yInt(:,2); 
thet = yInt(:,3);
thetDot = yInt(:,4);
xDotDot = yDot(2,:)'; 
yDotDot = yDot(4,:)';
momentum = xDot*mc + (xDot + L*cos(thet).*thetDot)*mp;
Fp = 1./sin(thet) .* (par.mc*xDotDot + par.b*xDot - u(:,2) );
Ep = - par.mp* -cos( thet )* par.g;
Ek = 0.5*par.mc*xDot.^2 + 0.5*par.mp*( (xDot + par.L*cos(thet).*thetDot).^2 + (-sin(thet)*L.*thetDot).^2  );

%% Metrics 
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
for j = 1:5:length(tInt)
    drawCartPendOwn( yInt(j,:), tInt(j), Fp(j), par)
end    