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

% Relative tolerance
p.param.rel_tol  = 1e-9; 
 
% Absolute tolerance
p.param.abs_tol  = 1e-9;

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
p.pred.theta0 = 0;                  % rad

% Values for K= predSpd/preySpd (gets coarse as values increase)
K = 10.^linspace(0,2,num_K);

% Angle of escape 
% prey_theta = 45;              % test value
prey_theta = [0:0.5:120] ./180*pi;    % rad


%% Run Simulation 
% (outer loop runs through speeds, inner loop runs through angles

% Initialize timer, figure
tic
f = figure('DoubleBuffer','on');
h = semilogx(K, 0.*K, 'o');
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
        
        %R = reconstruct(R, p, d);
        
        % Distance between predator and prey
        dist = hypot(R.xPred-R.xPrey,R.yPred-R.yPrey);
        
        % Store results
        B.prey_theta(j,k)   = R.thetaPrey(end);
        B.min_dist(j,k)     = min(dist);
        B.R{j,k}            = R;
        B.p{j,k}            = p;
        B.d{j,k}            = d;

    end
  
    % Find the max_min distance over all angles; find corresponding angle 
    [max_min_d, theta_ind]  = max(B.min_dist(j,:));
    B.theta_maxmin(j,1)     = B.prey_theta(j,theta_ind);

    % Update plot
    figure(f)
    h = semilogx(B.K, B.theta_maxmin .*180/pi,'o-r');
    set(h,'MarkerFaceColor','r')
    pause(0.001)
    xlabel('K')
    ylabel('Theta maxmin')
    ylim([0 90])
    
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
