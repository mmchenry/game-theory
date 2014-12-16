function [max_min_dist, theta_min]= run_batch_sim
% RUN_BATCH_SIM runs batch simulations of 'simple_model' for Weihs
% situation

%% Define simulation parameters
% Default parameter values
[p,d] = default_params;

% Time span for simulation
p.param.t_span        = [0 4];      % s

p.param.tank_radius = 1;            % m  

% prey initial position & speed
p.prey.x0 = 1e-2;                   % m
p.prey.y0 = 0;                      % m
p.prey.spd0 = 0;                    % m/s
p.prey.theta0 = 0;                  % rad  
p.prey.thrshEscape    = 2e-2;       % m

% prey escape duration (time to reach escape angle)
p.prey.durEscape = 25e-3;           % s


% predator initial position & speed
p.pred.x0 = -1e-2;                  % m
p.pred.y0 = 0;                      % m
p.pred.spd0 = 1e-2;                 % m/s

% Values for K= predSpd/preySpd (gets coarse as values increase)
K = 0.1:0.1:0.9;                % unitless
K = [K, 1.25:0.25:6.25];        % unitless
K = [K, 6.5:0.5:15];            % unitless
K = [K, 16:1:40];               % unitless
% K = 2;                          % test value

% prey escape speed values 
prey_spdEscape = p.pred.spd0./K;    % m/s

M = length(prey_spdEscape);

% Angle of escape 
% prey_theta = 45;              % test value
prey_theta = 0:0.5:95;                  % deg
prey_theta_rad = prey_theta./180*pi;    % rad

% Rotational speed during escape
prey_rotSpd = prey_theta_rad ./ p.prey.durEscape;   % rad/s

% length of rotational speed vector
N = length(prey_rotSpd);

% preallocate minimum distance matrix
min_distance = zeros(N,M);
% t_at_min = zeros(N,M);

% preallocate max-min Distance, index where it occurs & corresponding theta
max_min_dist = zeros(M,1);
theta_index = zeros(M,1);
theta_min = zeros(M,1);

%% Run Simulation 
% (outer loop runs through speeds, inner loop runs through angles
for j=1:M
    % set current escape speed
    p.prey.spdEscape = prey_spdEscape(j);

    for k=1:N
        % set current speed of rotation 
        p.prey.rotSpdEscape = prey_rotSpd(k);
        
        % set current initial theta of prey (used for testing)
%         p.prey.theta0 = prey_theta_rad(k);

        
        % Run Simulation
        R = simple_model(p,d);
        
        % Store input parameters
        R.p     = p;
        R.d     = d;
        
        % find minimum distance of simulation and store output
        [min_dist,dist_vals] = distance_fun(R);
        min_distance(k,j) = min_dist;

        % Plot distance function
%         figure;
%         plot(R.t, dist_vals,'o','MarkerSize',10)
        
        % Plot orientation angles
%         figure;
%         vis_results('Turning data',R)
    end
    
    % Find the max_min distance over all angles; find corresponding angle 
    [max_min_d, theta_ind] = max(min_distance(:,j));
    max_min_dist(j) = max_min_d;
    theta_index(j) = theta_ind;
    theta_min(j) = prey_theta(theta_ind);
end

end
    
function [min_dist,dist_vals] = distance_fun(R)
dist_vals = (R.xPred-R.xPrey).^2 + (R.yPred-R.yPrey).^2;
min_dist = min(dist_vals);
end
