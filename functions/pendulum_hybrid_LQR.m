function u = hybrid_LQR_V( state, stateGoal, K, par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


mp = par.mp ;
g = par.g;
L= par.L;
b = par.b; 

thet = state(1);
upDelta = wrapToPi(state(1)-pi);
thetDot = state(2);

V = cos(thet)*L*mp*g;
T = 1/2* mp*thetDot.^2*L^2;
E = V+T;

E0 = cos(stateGoal(1))*L*mp*g;

thet_threshold = 0.2;
E_error = ( cos(pi) - cos( thet_threshold - pi))*L*mp*g;



if ( E-E0 > 0.1) &&         ~(( abs( upDelta ) < thet_threshold ) &&  (abs(state(2)) < thet_threshold) )
    u = -sign( state(2) )*5 ; 
elseif (E-E0 < -0.1 ) &&    ~(( abs( upDelta ) < thet_threshold ) &&  (abs(state(2)) < thet_threshold) )
    u = sign( state(2) )*5 ; 
    
% elseif  ( abs( upDelta ) > thet_threshold ) && (abs(E-E0) < 5 )
elseif  ( abs( upDelta ) < thet_threshold ) &&  (abs(state(2)) < thet_threshold)
%     upDelta
%     thet
    u = -K*(   [  upDelta ; state(2) ]); 
else
    u = 0;
end







    
% 
