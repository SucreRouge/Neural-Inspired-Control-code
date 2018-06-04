function [ dydt ] = cartSinglePend( y, u, par )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mc = par.mc ;
mp = par.mp ;
g = par.g;
L= par.L;
bT = par.bT; 
bX = par.bX; 
dist = par.dist;

S = sin( y(3) );
C = cos( y(3) );
D = (mc+mp*S.^2);

dydt(1,:) = y(2);
dydt(2,:) = 1/D * ( y(4)^2*mp*L*S + y(4)*bT/L*C - y(2)*bX - mp*g*S*C + u );
dydt(3,:) = y(4);
dydt(4,:) = 1/D * ( -y(4)^2*mp*S*C -y(4)*bT*(mc+mp)/(mp*L^2) +y(2)*bX*C/L    + (mc+mp)*g/L*S - u/L*C )  + dist*randn;