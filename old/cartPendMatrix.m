function [ dydt  ] = cartPendOwn( y, u, pendPar )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mc = pendPar.mc ;
mp = pendPar.mp;
g = pendPar.g;
L = pendPar.L;
b = pendPar.b; 
dist = pendPar.dist;

S = sin( y(:,3) );
C = cos( y(:,3) );
D = (mc+mp*S.^2);

% y
u
% u(:,2)
% u(:,4)
dydt(1,:) = y(:,2);
dydt(2,:) = 1./D .* ( -y(:,2)*b +y(:,4).^2.*mp.*L.*S - mp*g*S.*C + u(:,2)  );
dydt(3,:) =  y(:,4);
dydt(4,:) = 1./D .* ( y(:,2).*b.*C/L - y(:,4).^2.*mp.*S.*C + (mc+mp)*g*S/L + u(:,4)*C/L )    + dist*randn;


end

