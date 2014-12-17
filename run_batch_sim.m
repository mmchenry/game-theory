function [max_min_dist, theta_min]= run_batch_sim
% RUN_BATCH_SIM runs batch simulations of 'simple_model' for Weihs
% situation


%% Parameter values

% Number of K values
num_K = 15;


%% Path definitions

% Identify Matt's computer
if isdir('/Users/mmchenry/Documents')
    root = '/Users/mmchenry/Documents/Projects/Game theory';
    
else
    %TO DO: Add path for saving data on Alberto's computer
    
end

% Path to data
dPath = [root filesep 'Batch weihs'];


%% Define simulation parameters
% Default parameter values
[p,d] = default_params;

% Time span for simulation
p.param.t_span        = [0 4];      % s

% Make tank way too big
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
p.pred.x0 = -2e-2;                  % m
p.pred.y0 = 0;                      % m
p.pred.spd0 = 1e-2;                 % m/s

% Values for K= predSpd/preySpd (gets coarse as values increase)
% K = 0.1:0.1:0.9;                % unitless
% K = [K, 1.25:0.25:6.25];        % unitless
% K = [K, 6.5:0.5:15];            % unitless
% K = [K, 16:1:40];               % unitless
% K = 2;                          % test value
K = 10.^linspace(0,1.6,num_K);

% Numbers of escape speeds to try
%M = length(prey_spdEscape);

% Angle of escape 
% prey_theta = 45;              % test value
prey_theta = [0:0.5:95] ./180*pi;    % rad

% Rotational speed during escape
%prey_rotSpd = prey_theta_rad ./ p.prey.durEscape;   % rad/s

% length of rotational speed vector
%N = length(prey_rotSpd);


%% Run Simulation 
% (outer loop runs through speeds, inner loop runs through angles

% Initialize timer, figure
tic

% Set up plot
f = figure('DoubleBuffer','on');
h = plot(K, 0.*K, 'o');
set(h,'Color',0.5.*[1 1 1])
hold on

% Loop thru escape speeds
for j = 1:length(K)
    
    % Store values for simulation
    B.K(j,1)  = K(j);
    
    % Set current escape speed
    p.prey.spdEscape = p.pred.spd0./K(j);
    
    % Loop thru angles of rotation
    for k = 1:length(prey_theta)
        
        % set current speed of rotation 
        p.prey.rotSpdEscape = prey_theta(k) ./ p.prey.durEscape;
        
        % set current initial theta of prey (used for testing)
%         p.prey.theta0 = prey_theta(k);

        % Run Simulation
        R   = simple_model(p,d);
        
        R = reconstruct(R, p, d);
        
        % Distance between predator and prey
        dist = hypot(R.xPred-R.xPrey,R.yPred-R.yPrey);
        
        % Store results
        B.prey_theta(j,k)   = R.thetaPrey(end);
        B.min_dist(j,k)     = min(dist);
        B.R{j,k}            = R;
        B.p{j,k}            = p;
        B.d{j,k}            = d;

        % Plot distance function
         figure;
         plot(R.t, dist,'-','MarkerSize',10)
        
        % Plot orientation angles
        vis_results('Turning data',R)
    end
  
    % Find the max_min distance over all angles; find corresponding angle 
    [max_min_d, theta_ind]  = max(B.min_dist(j,:));
    B.theta_maxmin(j,1)     = B.prey_theta(j,theta_ind);

    % Update plot
    figure(f)
    h = plot(B.K, B.theta_maxmin,'o-r');
    set(h,'MarkerFaceColor','r')
    pause(0.001)
    xlabel('K')
    ylabel('Theta maxmin')
    ylim([0 pi/2])
    
    % Save data so far
    save([dPath filesep date],'B')
    
    % Time so far
    tlapse = toc;
    
    % Time per loop
    time_per = tlapse/j;
    
    % Time remaining 
    time_left = round((length(K)-j)*time_per/60);
    
    % Update
    disp(['Completed ' num2str(j) ' of ' num2str(length(K)) ', '  ...
          num2str(time_left) ' min left']);
end

end
