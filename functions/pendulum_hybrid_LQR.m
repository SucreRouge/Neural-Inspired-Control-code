function u = hybrid_LQR_V( state, stateGoal, K, par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


mp = par.mp ;
g = par.g;
L= par.L;
b = par.b; 

thet = state(1);
thetDot = state(2);

V = cos(thet)*L*mp*g;
T = 1/2* mp*thetDot.^2*L^2;
E = V+T;

E0 = cos(stateGoal(1))*L*mp*g;

thet_threshold = pi*0.1;
E_error = ( cos(pi) - cos( thet_threshold - pi))*L*mp*g;


if   ( ( wrapToPi(abs(thet))-pi) < thet_threshold ) && ( abs(E-E0) < E_error )
    % LQR
    u = -K*(   [  sin(thet-pi); state(2) ]);  % why the sin? it wraps theta nicely
elseif ( E-E0 > 0)
    % remove energy 
    u = -sign(thetDot) ; 
elseif (E-E0 < 0 )
    % add energy 
    u = sign(thetDot) ; 
else
    % if energy levels are good, but not in linear control range 
    u = 0;
end

