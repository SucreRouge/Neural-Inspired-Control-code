function [ dydt  ] = cartSinglePend_matrixOperation( y, u, pendPar )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mc = pendPar.mc ;
mp = pendPar.mp;
g = pendPar.g;
L = pendPar.L;
bX = pendPar.bX; 
dist = pendPar.dist;

S = sin( y(:,3) );
C = cos( y(:,3) );
D = (mc+mp*S.^2);

dydt(1,:) = y(:,2);
dydt(2,:) = 1./D .* ( -y(:,2)*bX +y(:,4).^2.*mp.*L.*S - mp*g*S.*C + u  );
dydt(3,:) =  y(:,4);
dydt(4,:) = 1./D .* ( y(:,2).*bX.*C/L - y(:,4).^2.*mp.*S.*C + (mc+mp)*g*S/L + u.*C/L )    + dist*randn;


