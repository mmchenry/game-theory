function [max_min_d, theta_min]= run_batch_sim
% RUN_BATCH_SIM runs batch simulations of 'simple_model' for Weihs
% situation


%% Parameter values

% Number of K values
num_K = 15;


%% Path definitions

% Identify Matt's computer
if isdir('/Users/mmchenry/Documents/Projects')
    root = '/Users/mmchenry/Documents/Projects/Game theory';
    
else
    root = '/Users/alberto/Dropbox/Review with Alberto/PredPrey Data';
    
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
p.prey.x0 = 0.25e-2;                   % m
p.prey.y0 = 0;                      % m
p.prey.spd0 = 0;                    % m/s
p.prey.theta0 = 0;                  % rad  
p.prey.thrshEscape    = 2e-2;       % m

% prey escape duration (time to reach escape angle)
p.prey.durEscape = 25e-3;           % s

% predator initial position & speed
p.pred.x0 = 0e-2;                  % m
p.pred.y0 = 0;                      % m
p.pred.spd0 = 1e-2;                 % m/s
p.pred.theta0 = 0;                  % rad

% Values for prey initial position (acts as response distance)
prey_x0 = linspace(1e-3,2e-2,10);

% Values for K= predSpd/preySpd (gets coarse as values increase)
K = 10.^linspace(-1,2,num_K);

% Angle of escape 
% prey_theta = 60 ./180*pi;              % test value
prey_theta = [0:0.5:120] ./180*pi;    % rad

% Set up cell array to store all output
S = cell(1,length(prey_x0));

% Store K values
B.K            = K;

% Store theta values
B.prey_theta   = prey_theta;

%% Run Simulation
% outermost loop runs through initial distances
% middle loop runs through speeds, inner loop runs through angles

% Initialize timer, figure
tic
% f = figure('DoubleBuffer','on');
% h = semilogx(K, 0.*K, 'o');
% set(h,'Color',0.5.*[1 1 1])
% hold on

% Loop through initial distance
for i = 1:length(prey_x0)
    
    % Store initial distance values for simulation
    B.prey_x0 = prey_x0(i);
    
    % set current initial distance
    p.prey.x0 = prey_x0(i);
    
    % Loop thru escape speeds
    for j = 1:length(K)
        
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
            B.min_dist(j,k)     = min(dist);
            B.R{j,k}            = R;
            B.p{j,k}            = p;
            B.d{j,k}            = d;
            
        end
        
        % Find the max_min distance over all angles; find corresponding angle
        [max_min_d, theta_ind]  = max(B.min_dist(j,:));
        B.theta_maxmin(j,1)     = B.prey_theta(theta_ind);
        
        % Update plot
%         figure(f)
%         h = semilogx(B.K, B.theta_maxmin .*180/pi,'o-r');
%         set(h,'MarkerFaceColor','r')
%         pause(0.001)
%         xlabel('K')
%         ylabel('Theta maxmin')
%         ylim([0 100])
        
    end
    
    % Time so far
    tlapse = toc;
    
    % Time per loop
    time_per = tlapse/i;
    
    % Time remaining
    time_left = round((length(prey_x0)-i)*time_per/60);
    
    % Update
    disp(['Completed ' num2str(i) ' of ' num2str(length(prey_x0)) ', '  ...
        num2str(time_left) ' min left']);
    
    % Save output from simulation
    S{1,i}         = B;
end
        % Date and Time string for file name
        % Example output: 18-Dec-2014-11h12m27s
        s = strcat(datestr(clock,'dd-mmm-yyyy-HH'),'h',...
            datestr(clock, 'MM'),'m',datestr(clock,'ss'),'s');
        
        % Save all data
        save([dPath filesep s],'S','-v7.3')


%% plot Weihs optimal angles and optimal angles from simulations

% optimal angles from Weihs (theta=0 for K<1, use eq. 42 for K>=1)
theta_weihs = [0.* K(K<1), acos(1./K(K>=1))*180./pi]';

F = figure;
semilogx(K,theta_weihs,'LineWidth',3)
hold on

% run through S and plot optimal angles
for l = 1:length(prey_x0)
    figure(F)
        semilogx(K, S{l}.theta_maxmin .*180/pi,'o-');
        xlabel('K')
        ylabel('Optimal Escape Angle')
        ylim([0 100])
end

% TO DO: Update colors and add labels/legend on plot.

end
