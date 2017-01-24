function run_sim
% Runs single simulation and analyzes results


% Default parameter values
[p,d] = default_params;

p.param.sL = 1;
p.param.sT = 1;

p.param.t_span = [0 30];

% prey initial position & speed
p.prey.x0     = 3e-2;            % m
p.prey.y0     = 0;                  % m
p.prey.spd0   = 1.5e-2;                  % m/s
p.prey.theta0 = pi/4;                  % rad  
p.prey.thrshEscape    = 2e-2;       % m
p.prey.spdEscape      = 8e-2;       % m/s 

% Interval between start of saccades
p.prey.sccd_intvl       = 1.85;         % s

% Duration of a saccade
p.prey.sccd_prd         = 0.35;       % s

% Peak rate of rotation during saccade
p.prey.sccd_omega       = 5;         % rad/s

% prey escape duration (time to reach escape angle)
p.prey.durEscape = 25e-3;           % s

% predator initial position & speed
p.pred.x0 = -5e-2;                  % m
p.pred.y0 = 5e-2;                      % m
p.pred.spd0 = 4e-2;      		% m/s
p.pred.theta0 = -pi/3;                  % rad

% Duration of a saccade
p.pred.sccd_prd         = 0.25;       % s

p.pred.sccd_omega       = 1;         % rad/s

p.param.max_step = p.pred.sccd_prd/2;

p.pred.wall_omega       = 8;         % rad/s

%% Run simulation

R = simple_model(p,d);

R = reconstruct(R, p, d);


%% Report result


% % Display whether or not prey escaped
% if isempty(R.tEnd)
%     disp('prey successfully escaped')
% else
%     disp(['the prey was captured at time = ' num2str(R.tEnd)])
% end


%% Plot results

% Plot orientation angles
% vis_results('Turning data',R)


vis_results('Trajectories',R)


%% Animate simulation

%figure('DoubleBuffer','on')
%vis_results('Animate',R)
% TODO: Fix the animation code




