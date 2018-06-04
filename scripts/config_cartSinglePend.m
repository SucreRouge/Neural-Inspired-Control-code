% parameter config file

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
bX = 1;
bT = 1;
dist = 0;
par.mc = mc ; par.mp = mp; par.g = g; par.L = L; par.bX = bX; par.bT = bT; par.dist = dist;

% control penalties 
Q = diag( [1,100,1,1]);
R = 1;
p = [-1.1,-1.2,-1.3,-1.4];  % pole placements 

Vd = 0.1*eye(4); % disturbance covariance (state uncertainty) 
Vn = 1; %measurement noise 

% time parameters
dt = 0.01; 
tLast = 10;
tInt = dt:dt:tLast;

% control IC and objective
y0 = [-4; 0; pi*11/12; 0];
yGoal = [1; 0; pi; 0];

% noise parameters
sensorNoiseVariance = [1e-4 1e-4 1e-3 1e-4]';
% sensorNoiseVariance = [0; 0; 0; 0 ];
distXAccVariance = 1e-4;
distThetaAccVariance = 1e-4;
