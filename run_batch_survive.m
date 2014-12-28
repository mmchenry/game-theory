function run_batch_survive(batch_type)
% RUN_BATCH_SURVIVE runs batch simulations of 'simple_model' to determine
% whether or not the prey survives given (x0,y0), k, and optimal angle. The
% optimal angle was computed previously and is contained in 'angles.mat'.
% Currently there's only one batch_type: 'Survive'

% Grid points used to find optimal angles: [X,Y] = meshgrid(0:0.005,0.305)
% K values used for optimal angles: K = 10.^linspace(0,2,18);

%% Batch parameters & load angles

% K values
num_K = 18;
min_K = 1;
max_K = 100;  

if (nargin<1) || strcmp(batch_type,'Survive')
    % step size for grid discretizing 
    step = 0.005;       % 5mm
    
    % initial position limits (m)
    min_x0 = 0;
    max_x0 = 0.305;
    min_y0 = 0;
    max_y0 = 0.305;
else
    error('Do not recognize batch_type')
end


% Define batch type
B.batch_type = batch_type;

% load optimal angles: angles.mat contains optimal angles at each (x0,y0)
% for various values of K. 'angles' is a 62 x 62 x 18 matrix. values are
% stored as (x0,y0,K)
angles = [];            % initialize matrix 
load 'angles.mat'


%% Path definitions

% Identify Matt's computer
if isdir('/Users/mmchenry/Dropbox/Literature')
    root = '/Users/mmchenry/Dropbox/Literature/Review with Alberto/predprey data';
    
% Or Alberto's
else    
%     root = '/Users/alberto/Dropbox/Review with Alberto/PredPrey Data';
    root = '/Users/A_Soto/Documents/MATLAB/McHenry_LabRotation/angle_data';
end

% Path to data
dPath = root;


%% Simulation parameters

% Default parameter values
[p,d] = default_params;


% Time span for simulation
p.param.t_span      = [0 7];        % s

% Make tank way too big
p.param.tank_radius = 0.5;          % m  

% prey initial position, speed & heading
p.prey.x0           = 0e-2;         % m
p.prey.y0           = 0;            % m
p.prey.spd0         = 0;            % m/s
p.prey.theta0       = 0;            % rad  

% distance threshold for prey escape response
p.prey.thrshEscape  = 2e-2;         % m

% prey escape duration (time to reach escape angle)
p.prey.durEscape    = 25e-3;        % s

% predator initial position, speed & heading
p.pred.x0           = 0e-2;         % m
p.pred.y0           = 0;            % m
p.pred.spd0         = 5e-2;         % m/s
p.pred.theta0       = 0;            % rad

% Sample rate of visual system
p.pred.vis_freq = 10; % 1/s

% Values for K= predSpd/preySpd (gets coarse as values increase)
B.K = 10.^linspace(log10(min_K), log10(max_K),num_K)';

% Store parameters
B.p = p;
B.d = d;

%% Parameters for batch type

if (nargin<1) ||strcmp(batch_type,'Survive')
    
    % Values for prey initial position 
    B.prey_x0         = min_x0:step:max_x0;
    B.prey_y0         = min_y0:step:max_y0;
    
    % Generic values to vary in a batch
    % column 1: x0
%     batch_vals(:,1)        =  B.prey_x0;
    batch_vals(:,1)        =  0.10;          % TEST value
    %column 2: y0
%     batch_vals(:,2)        =  B.prey_y0;
    batch_vals(:,2)        =  0.0;           % TEST value
    
    % Simulation type
    p.param.sim_type  = 'Weihs, survive';
    
    % Relative tolerance
    p.param.rel_tol  = 1e-9;
    
    % Absolute tolerance
    p.param.abs_tol  = 1e-9;
end

% preallocate escape matrix to store whether or not prey escapes at (x0,y0)
% for given K
B.escape = zeros(length(B.prey_x0),length(B.prey_y0),length(B.K));

% preallocate minimum distance vector
B.min_dist = zeros(length(B.prey_x0),length(B.prey_y0),length(B.K));

%% Run Simulation
% outermost loop runs through speed ratio K
% middle loop runs through x0, inner loop runs through y0

% Initialize timer
tStart = tic;

% Loop through initial escape speed (k values)
for i = 1:length(B.K)
    
        % Set current escape speed
        p.prey.spdEscape = p.pred.spd0./B.K(i);

    % Loop thru x0
    for j = 1:size(batch_vals,1)
        
     %----For Testing-----%        
        % find index of current x0 to extract angle
        x_ind = find(B.prey_x0==batch_vals(j,1));
     %---------------------%     
     
        % set current x0
        p.prey.x0 = batch_vals(j,1);
     
        for k = 1:size(batch_vals,1)
            
     %----For Testing-----%
            % find index of current x0 to extract angle
            y_ind = find(B.prey_y0==batch_vals(k,1));
            
            % set current optimal angle
            theta = angles(x_ind,y_ind,i);
     %---------------------%   
     
            % set current y0
            p.prey.y0 = batch_vals(k,2);    
            
            % set current optimal angle
%             theta = angles(j,k,i);

            % Set rotation speed according to optimal theta
            p.prey.rotSpdEscape = theta ./ p.prey.durEscape;
            
            % Run simulation
            R   = simple_model(p, d);
            
     %---For Testing--------%       
            % Display whether or not prey escaped
            if isempty(R.tEnd)
                disp('prey successfully escaped')
                
                % escape = 1 at point (x0,y0) if prey escapes
                B.escape(x_ind,y_ind,i) = 1;
            else
                disp(['the prey was captured at time = ' num2str(R.tEnd)])
                % escape = 0 at point (x0,y0) if prey captured
                B.escape(x_ind,y_ind,i) = 0;
            end
            
            % Reconstruct simulation data
            R = reconstruct(R, p, d);
            
            % Plot trajectories
            vis_results('Trajectories',R)
        end
     %--------------------%       
            
            % Distance between predator and prey with initial position of
            % prey at (x0,y0)
            min_dist = min(hypot(R.xPred-R.xPrey,R.yPred-R.yPrey));
            
            B.min_dist(x_ind,y_ind,i) = min_dist;
                        
        clear min_dist 
    end
    
    % Update status
    update_time(tStart, i, length(B.K), '');

end

%% Store data
% make escape matrix into logical matrix - to use color code quiver vectors
B.escape = logical(B.escape);

B.escape_angles = angles;

% Date and Time string for file name
% Example output: 18-Dec-2014-11h12m27s
s = strcat(datestr(clock,'dd-mmm-yyyy-HH'),'h',...
    datestr(clock, 'MM'),'m',datestr(clock,'ss'),'s');

if strcmp(batch_type,'Survive')
    s = ['Survive batch-' s]; 
end

% Save all data
save([dPath filesep s],'B','-v7.3')

%% Update status 
function update_time(tStart, idx, len, txt)
% Reports progress at the command line
    
    % Time so far
    tlapse = toc(tStart);

    % Time per loop
    time_per = tlapse/idx;

    % Time remaining
    time_left = ceil(((len-idx)*time_per/60));

    % Update
    disp([txt 'Done ' num2str(idx) ' of ' num2str(len) ', ~'  ...
        num2str(time_left) ' min left']);

end


end




