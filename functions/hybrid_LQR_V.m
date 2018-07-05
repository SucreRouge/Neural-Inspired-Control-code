function u = hybrid_LQR_V( state, stateGoal, K, par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ep = - par.mp* -cos( state(3) )*par.L* par.g;
Ek = 0.5*par.mc*state(2).^2 + 0.5*par.mp*( (state(2) + par.L*cos(state(3)).*state(4)).^2 + (-sin(state(3)).*par.L.*state(4)).^2  );
V = Ep+Ek;

V0 = 20;
theta_error = pi*1/12;
G=0.5;

if  ( abs(state(3)-pi ) <= theta_error )
    u = -K*(state-stateGoal);
else
    u = G*(V-V0)*state(4)*cos(state(3)) ;
end



