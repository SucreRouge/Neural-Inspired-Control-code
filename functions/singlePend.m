function [ dydt ] = singlePend( y, u, par )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mp = par.mp ;
g = par.g;
L= par.L;
b = par.b; 

S = sin( y(1) );
% 
% size(y)
dydt(1,1) = y(2);
% % size(u)
% % size(y)
dydt(2,1) = -b/(mp*L^2)*y(2) + g/L*S  + 1/(mp*L^2)*u;