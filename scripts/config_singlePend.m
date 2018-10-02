% parameter config file
% contains parameterset to run inverted pendulum simulations. Add test description 
% Add paths required for main script
rootPath = fileparts(fileparts(mfilename('fullpath') ));
addpath([rootPath filesep 'functions'])
addpath([rootPath filesep 'scripts'])
addpath([rootPath filesep 'figs'])
addpath([rootPath filesep 'simulink'])

% physical parameters
mp = 1;
g = -10;
% L = 2;
L = 4*abs(g)/pi^2; % this L scales theta and theta dot to radius pi in phase portrait 
% b = 0.06;
b = 0;
dist = 0; 
yGoal = [pi; 0];
u = 0;

% linearized matrices
A = [ 0   1        ;
      0   -b/(mp*L^2)];     
B = [0 ;
    1/(mp*L^2)];

% set control gain
Q = [10,0; 0,1];
R = 1e-3; 
K = lqr( A,B,Q,R);

% time parameters for main simulation 
dt = 0.01; 
tLast = 20;
tSamp = 0.1;
tInt = dt:dt:tLast;

% convert into parameter struct 
par.mp = mp; par.g = g; par.L = L; par.b = b; par.dist = dist; 