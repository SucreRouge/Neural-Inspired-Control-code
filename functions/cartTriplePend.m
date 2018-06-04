function [ dydt ] = cartTriplePend( y, u, par )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mc = par.mc ;
m1 = par.m1;
m2 = par.m2;
m3 = par.m3;
g = par.g;
L1= par.L1;
L2= par.L2;
L3= par.L3;
b = par.b; 
c = par.c; 

S1 = sin( y(2) );
C1 = cos( y(2) );
S2 = sin( y(3) );
C2 = cos( y(3) );
S3 = sin( y(4) );
C3 = cos( y(4) );
d12 = y(2)-y(3);
d13 = y(2)-y(4);
d23 = y(3)-y(4);

D = [ mc+m1+m2+m3              (m1+m2+m3)*L1*C1           (m2+m3)*L2*C2             m3*L3*C3; 
     (m1+m2+m3)*L1*C1          (m1+m2+m3)*L1^2            (m2+m3)*L1*L2*cos(d12)    m3*L1*L3*cos(d13); 
     (m2+m3)*L2*C2             (m2+m3)*L1*L2*cos(d12)    (m2+m3)*L2^2              m3*L2*L3*cos(d23); 
     (m3)*L3*C3                 m3*L1*L3*cos(d13)          m3*L2*L3*cos(d23)       m3*L3^2 ];
  
C = [0         -(m1+m2+m3)*L1*S1           -(m2+m3)*L2*S2           -m3*L3*S3;
	 0         0                           (m2+m3)*L1*L2*sin(d12)   m3*L2*L3*sin(d13);
	 0         -(m2+m3)*L1*L2*sin(d12)     0                        m3*L2*L3*sin(d23);
	 0         -m3*L1*L3*sin(d13)          -m3*L2*L3*sin(d23)        0              ];

F = [-u+b*y(5);
    -(m1+m2+m3)*g*L1*S1+c*y(6)  ;
    -(m2+m3)*g*L2*S2+c*y(7)  ;
    -(m3)*g*L3*S3+c*y(8)  ];
 
Cstates = [y(5)^2; y(6)^2; y(7)^2; y(8)^2];

dydt(1,:) = y(5);
dydt(2,:) = y(6);
dydt(3,:) = y(7);
dydt(4,:) = y(8);
dydt(5:8,:) = inv(D)*( -C*Cstates -F ); 