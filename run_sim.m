function run_sim
% Runs single simulation and analyzes results


%% Prey parameters  

% Initial body position & orientation 
prey.x0             = 1e-2;      % m          
prey.y0             = 1e-2;      % m          
prey.theta0         = 45/180*pi; % rad      

prey.spdEscape      = 5e-2;      % m/s  
prey.spd0           = 0;         % m/s  

prey.bod_len        = 3.5e-3;    % m


%% Predator parameters

% Initial body position & orientation 
pred.x0                 = 0;         % m      
pred.y0                 = 0;         % m      
pred.theta0             = 0;         % rad  

% Initial speed
pred.spd0               = 2e-2;      % m/s 

% Interval between start of saccades
pred.saccade_interval   = 2;         % s

% Duration of a saccade
pred.saccade_period     = 0.5;       % s

% Peak rate of rotation during saccade
pred.saccade_omega      = 2;         % rad/s

% Peak rate of rotation encounter with wall
pred.wall_omega         = 40;         % rad/s

% Width of head
pred.body_width         = 0.3e-2;    % m

% Length of head
pred.head_length        = 0.4e-2;    % m 

% Proportion region beyond head for responding to walls
pred.fieldSize          = 2;        % Dimensionless

% Receptive field of visual system (rad). Easter (1996) gives the functional 
% retinal field of 163 deg for 4 dpf. Though I give an additional 20 deg
% for saccades (based on vergenece angles of Patterson et al, 2013).
pred.vis_az     = (163+20)/180*pi; % rad

% Vergence angle: mean = 32 deg, max = 70 deg (value for larvae) (Patterson et al, 2013)
pred.verg_ang = 32/180*pi; % rad

% Deg per photorecptor pair (Easter, 1996) for 4 dpf
pred.rtnl_den = 2.9235/180*pi; % rad


%% General parameters

% Capture distance (i.e. distance where pred. captures prey)
param.d_capture     = 5e-3  ;    % m

% Proximity at which prey responds to predator
param.dist_thresh   = 1e-2;      % m

% Tank dimension
param.tank_radius   = 10e-2;     % m


%% Solver parameters

% Time span for simulation
param.t_span        = [0 60];     % s

% Relative tolerance (applies to all components of the solution vector)
param.rel_tol        = 1e-3;

% Absolute tolerance (applies to individual components of solution vector)
param.abs_tol       = 1e-6;


%% Run simulation

R = simple_model(pred,prey,param);


%% Display results

% Animate results with even time intervals
%animate_sim(R)

% Plot orientation angles
figure;
vis_results(R,'Turning data')

figure;
vis_results(R,'Trajectories')

