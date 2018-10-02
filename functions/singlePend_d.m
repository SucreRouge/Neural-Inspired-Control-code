function [ dydt ] = singlePend_d( y, u, par )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mp = par.mp ;
g = par.g;
L= par.L;
b = par.b; 

S = sin( y(1) );
% 
dydt(1,1) = y(2) + randn(1)*par.dist;
dydt(2,1) = -b/(mp*L^2)*y(2) + g/L*S  + 1/(mp*L^2)*u + randn(1)*par.dist;