function run_sim
% Runs single simulation and analyzes results


%% General parameters

% Capture distance (i.e. distance where pred. captures prey)
param.d_capture     = 5e-3  ;    % m

% Proximity at which prey responds to predator
param.dist_thresh   = 1e-2;      % m

% Tank dimension
param.tank_radius   = 10e-2;     % m


%% Solver parameters

% Time span for simulation
param.t_span        = [0 80];     % s

% Relative tolerance (applies to all components of the solution vector)
param.rel_tol        = 1e-3;

% Absolute tolerance (applies to individual components of solution vector)
param.abs_tol       = 1e-5;


%% Prey parameters  

% Initial body position & orientation 
prey.x0             = -1e-2;      % m          
prey.y0             = -1e-2;      % m          
prey.theta0         = (180-30)/180*pi; % rad      

% Initial speed
prey.spd0           = 0.15e-2;        % m/s
% prey.spd0           = 0;        % m/s

% Escape parameters
prey.spdEscape      = 5e-2;       % m/s  
prey.omegaEscape    = 35;         % rad/s
prey.spdResp        = 2e-2;       % m/s
prey.latEscape      = 4e-3;       % s
prey.durEscape      = 40e-3;      % s

% Morphometrics
prey.bod_len        = 3.5e-3;           % m
prey.bod_width      = 0.5e-3;           % m
prey.COM_len        = prey.bod_len/4;   % m

% Length of head
prey.head_length        = 3.75e-6;    % m

% Interval between start of saccades
prey.saccade_interval   = 4;         % s

% Duration of a saccade
prey.saccade_period     = 0.25;       % s

% Peak rate of rotation during saccade
prey.saccade_omega      = 4;         % rad/s

% Peak rate of rotation encounter with wall
prey.wall_omega         = 20;         % rad/s

% Proportion region beyond head for responding to walls
prey.fieldSize       = prey.bod_width * 2;  % m
  

%% Predator parameters

% Initial body position & orientation 
pred.x0                 = 0e-2;         % m      
pred.y0                 = 0e-2;         % m      
pred.theta0             = 3*pi/4;      % rad  

% Initial speed
pred.spd0               = 2e-2;      % m/s 

% Interval between start of saccades
pred.saccade_interval   = 2;         % s

% Duration of a saccade
pred.saccade_period     = 0.5;       % s

% Peak rate of rotation during saccade
pred.saccade_omega      = 2;         % rad/s

% Peak rate of rotation encounter with wall
pred.wall_omega         = 10;         % rad/s

% Width of head
pred.bod_width          = 0.3e-2;    % m
pred.bod_len            = 2e-2;      % m

% Position of COM
pred.COM_len        = pred.bod_len/4; % m 

% Distance for sensing walls
pred.fieldSize       = pred.bod_width * 3;  % m

% Distance from prey for initiating a strike
pred.strike_thresh  = 1e-2;  % m

% Strike duration
pred.strike_dur     = 30e-3; % s

% Maximum distance of capture
pred.strike_reach   = 5e-3;  % m

% Angular spread of capture area
pred.strike_range   = pi/2;  % rad

% Receptive field of visual system (rad). Easter (1996) gives the functional 
% retinal field of 163 deg for 4 dpf. Though I give an additional 20 deg
% for saccades (based on vergenece angles of Patterson et al, 2013).
pred.vis_az     = 163/180*pi; % rad
%pred.vis_az     = (20)/180*pi; % rad

% Vergence angle: mean = 32 deg, max = 70 deg (value for larvae) (Patterson et al, 2013)
pred.verg_ang = 32/180*pi; % rad

% Deg per photorecptor pair (Easter, 1996) for 4 dpf
pred.rtnl_den = 3.5 * 2.9235/180*pi; % rad

% Sample rate of visual system
pred.vis_freq = 10; % 1/s

% Threshold number of retinal cells for noticing prey
pred.vis_thresh = 0*1;



%% Run simulation

R = simple_model(pred,prey,param);

% Store input parameters
R.param     = param;
R.pred      = pred;
R.prey      = prey;


%% Display results

% Animate results with even time intervals
%animate_sim(R)

% Plot orientation angles
%figure;
%vis_results('Turning data',R)

figure;
vis_results('Trajectories',R)

figure('DoubleBuffer','on')
vis_results('Animate',R)

% Display whether or not prey escaped
if isempty(R.tEnd)
    disp('prey successfully escaped')
else
    disp(['the prey was captured at time = ' num2str(R.tEnd)])
end


