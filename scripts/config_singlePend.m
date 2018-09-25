% parameter config file
% contains parameterset to run inverted pendulum simulations. Add test description 
% Add paths required for main script
rootPath = fileparts(fileparts(mfilename('fullpath') ));
addpath([rootPath filesep 'functions'])
addpath([rootPath filesep 'scripts'])
addpath([rootPath filesep 'figs'])
% addpath([rootPath filesep 'simulink'])

% physical parameters
mp = 1;
g = -10;
L = 2;
b = 0.06;
par.mp = mp; par.g = g; par.L = L; par.b = b; 
% time parameters
dt = 0.01; 
tLast = 15;
tSamp = 0.1;
tInt = dt:dt:tLast;

% % control IC and objective
y0 = [pi,0];