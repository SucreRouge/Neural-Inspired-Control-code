% parameter config file
% contains parameterset to run inverted pendulum simulations. Add test description 
% Add paths required for main script
rootPath = fileparts(fileparts(mfilename('fullpath') ));
addpath([rootPath filesep 'functions'])
addpath([rootPath filesep 'scripts'])
addpath([rootPath filesep 'simulink'])

% physical parameters
mc = 5;
mp = 1;
g = -10;
L = 2;
bX = 0.5;
bT = 0.5;
dist = 0.1;
par.mc = mc ; par.mp = mp; par.g = g; par.L = L; par.bX = bX; par.bT = bT; par.dist = dist;


% xl = [-5,5];
par.xl = [-10,10];
% control penalties 
Q = diag( [1,1,1,10]);
R = 1e0;
p = [-1.1,-1.2,-1.3,-1.4];  % pole placements 

Vd = 0.1*eye(4); % disturbance covariance (state uncertainty) 
Vn = 1; %measurement noise 

% time parameters
dt = 0.01; 
tLast = 10;
tSamp = 0.1;
tInt = dt:dt:tLast;

% control IC and objective
% y0 = [-1; 0; pi*11/12; 0];
y0 = [-1; 0; pi; 0];
yGoal = [1; 0; pi; 0];

% noise parameters
sensorNoiseVariance = [1,1,1,1]'*1e-2;
% % sensorNoiseVariance = [0; 0; 0; 0 ];
distXAccVariance = 1e-1;
distThetaAccVariance = 1e-1;
