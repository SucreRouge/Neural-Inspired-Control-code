function [ dydt ] = cartDoublePend( y, u, par )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mc = par.mc ;
m1 = par.m1;
m2 = par.m2;
g = par.g;
L1= par.L1;
L2= par.L2;
b = par.b; 
c = par.c; 
d = par.d; 

S1 = sin( y(2) );
C1 = cos( y(2) );
S2 = sin( y(3) );
C2 = cos( y(3) );
dTheta = y(2)-y(3);

D = [ mc+m1+m2              (m1+m2)*L1*C1           m2*L2*C2    ; 
     (m1+m2)*L1*C1          (m1+m2)*L1^2            m2*L1*L2*cos(dTheta) ; 
      m2*L2*C2               m2*L1*L2*cos(dTheta)   m2*L2^2     ];
C = [0         -(m1+m2)*L1*S1           -m2*L2*S2               ;
     0          0                        m2*L1*L2*sin(dTheta)   ; 
     0          -m2*L1*L2*sin(dTheta)   0                       ];
F = [-u+b*y(4);
    -(m1+m2)*g*L1*S1+c*y(5)  ;
     -m2*g*L2*S2+d*y(6)];
 
Cstates = [y(4)^2; y(5)^2; y(6)^2; ];

dydt(1,:) = y(4);
dydt(2,:) = y(5);
dydt(3,:) = y(6);
dydt(4:6,:) = D^-1*( -C*Cstates -F ); 