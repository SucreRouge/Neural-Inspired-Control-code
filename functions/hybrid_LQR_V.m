function u = hybrid_LQR_V( state, stateGoal, K, par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ep = - par.mp* -cos( state(3) )*par.L* par.g;
% Ek = 0.5*par.mc*state(2).^2 + 0.5*par.mp*( (state(2) + par.L*cos(state(3)).*state(4)).^2 + (-sin(state(3)).*par.L.*state(4)).^2  );
Ek = 0.5*par.mc*state(2).^2 + 0.5*par.mp*( (state(2) + par.L*cos(state(3)).*state(4)).^2 + (-sin(state(3)).*par.L.*state(4)).^2  );
% Ek_noCart = (0.5*par.mp*( par.L*cos(state(3)).*state(4)).^2 + (-sin(state(3)).*par.L.*state(4)).^2  );
V = Ep+Ek;

V0 = 20;
theta_error = pi*0.5;
E_error =4;
G=6;

if  ( abs(state(3)-pi ) <= theta_error ) && (abs(V-V0) <E_error)
    u = -K*(state-stateGoal);
elseif (abs(state(1)) > 5)% && (abs(state(2)) > 1)
    u  = -0.5*(state(1));
    
else
%     u = G*(V-V0)*state(4)*cos(state(3)) ;
%     u 
%     u = G*sin(state(4))*( cos(state(4)) );
%     u = 5* sin(state(4))*sign( cos(state(3)) );
    u = sign(V0-V)*G*sin(state(4));
%     u = G*cos(state(4));
%     u = G*sign(sin(state(4)));
end



