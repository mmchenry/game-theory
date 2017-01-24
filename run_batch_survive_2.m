function run_batch_survive_2(batch_type)
% RUN_BATCH_SURVIVE runs batch simulations of 'simple_model' to determine
% whether or not the prey survives given (x0,y0), k, and escape angle. The
% optimal escape angle is computed by calling the nested function
% 'optimTheta'. Currently there's only one batch_type: 'Survive'


%% Batch parameters & load angles

% K values
num_K = 1;
min_K = 1.001;
max_K = 1.197;  

% Values for K= predSpd/preySpd (gets coarse as values increase)
B.K = 10.^linspace(log10(min_K), log10(max_K),num_K)';
K_vals = B.K;

if (nargin<1) || strcmp(batch_type,'Survive')
    % number of points for x0 & y0
    num_x0 = 15;    
    num_y0 = 15;
    
    % initial position limits (m)
    min_x0 = 0;
    max_x0 = 0.065;
    min_y0 = 0;
    max_y0 = 0.05;
else
    error('Do not recognize batch_type')
end

% Define batch type
B.batch_type = batch_type;

%% Path definitions

% Identify Matt's computer
if isdir('/Users/mmchenry/Dropbox/Literature')
    root = '/Users/mmchenry/Dropbox/Literature/Review with Alberto/predprey data';
    
% Or Alberto's
else    
%     root = '/Users/alberto/Dropbox/Review with Alberto/PredPrey Data';
%     root = '/Users/alberto/Desktop/Temporary Files';
    root = '/Users/A_Soto/Documents/MATLAB/McHenry_LabRotation/angle_data';
end

% Path to data
dPath = root;


%% Simulation parameters

% Default parameter values
[p,d] = default_params;


% Time span for simulation
p.param.t_span      = [0 2];        % s

% Make tank way too big
p.param.tank_radius = 0.5;          % m  

% prey initial position, speed & heading
p.prey.x0           = 0e-2;         % m
p.prey.y0           = 0;            % m
p.prey.spd0         = 0;            % m/s
p.prey.theta0       = 0;            % rad  

% distance threshold for prey escape response
p.prey.thrshEscape  = 1.05e-2;         % m

% prey escape duration (time to reach escape angle)
p.prey.durEscape    = 25e-3;        % s

% set max step size
p.param.max_step = p.prey.durEscape/2;

% predator initial position, speed & heading
p.pred.x0           = 0e-2;         % m
p.pred.y0           = 0;            % m
p.pred.spd0         = 5e-2;         % m/s
p.pred.theta0       = 0;            % rad

% Angular spread of capture area
p.pred.strike_range   = 2*pi/3;  % rad

% Maximum distance of capture
p.pred.strike_reach   = 7.5e-3;  % m

% Store parameters
B.p = p;
B.d = d;

%% Parameters for batch type

if (nargin<1) ||strcmp(batch_type,'Survive')
    
    % Values for prey initial position 
    B.prey_x0         = linspace(min_x0, max_x0, num_x0)';
    B.prey_y0         = linspace(min_y0, max_y0, num_y0)';
    
    % Generic values to vary in a batch
    % vector 1: x0
    batch_vals1        =  B.prey_x0;
    
    % vector 2: y0
    batch_vals2        =  B.prey_y0;
    
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

% preallocate minimum distance matrix
B.min_dist = zeros(length(B.prey_x0),length(B.prey_y0),length(B.K));

% preallocate optimal angles matrix
escape_angles = zeros(length(B.prey_x0),length(B.prey_y0),length(B.K));

%% Compute optimal escape angles

for l = 1:length(K_vals)
    k_spd = K_vals(l);
    for m = 1:length(batch_vals1)
        x0 = batch_vals1(m);
        for n = 1:length(batch_vals2)
            y0 = batch_vals2(n);
            
            % compute current optimal escape angle
            theta = optimTheta(x0, y0, k_spd);
            
            % store optimal escape angle
            escape_angles(m,n,l) = theta;
        end
    end
end
disp('optimal angles computed')

B.escape_angles = escape_angles;


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
    for j = 1:length(batch_vals1)
     
        % set current x0
        p.prey.x0 = batch_vals1(j);
     
        for k = 1:length(batch_vals2)

            % set current y0
            p.prey.y0 = batch_vals2(k);

            % extract optimal escape angle
            escTheta = B.escape_angles(j,k,i);

            % Set rotation speed according to optimal theta
            p.prey.rotSpdEscape = escTheta ./ p.prey.durEscape;
            
            % Run simulation
            R   = simple_model(p, d);
                 
            % Distance between predator and prey with initial position of
            % prey at (x0,y0)
            min_dist = min(hypot(R.xPred-R.xPrey,R.yPred-R.yPrey));
            
            B.min_dist(j,k,i) = min_dist;
                
            % Display whether or not prey escaped & score simulation
            if isempty(R.tEnd)
%                 disp('prey successfully escaped')
                
                % escape = 1 at point (x0,y0,k) if prey escapes
                B.escape(j,k,i) = 1;
            else
%                 disp(['the prey was captured at time = ' num2str(R.tEnd)])
                
                % escape = 0 at point (x0,y0,k) if prey captured
                B.escape(j,k,i) = 0;
            end
            
            % Reconstruct simulation data
%             R = reconstruct(R, p, d);
            
            % Plot trajectories
%             vis_results('Trajectories',R)    

        end

    end
    
    % Update status
    update_time(tStart, i, length(B.K), '');

end

%% Store data
% make escape matrix into logical matrix - to color code quiver vectors
B.escape = logical(B.escape);

% Date and Time string for file name
% Example output: 18-Dec-2014-11h12m27s
s = strcat(datestr(clock,'dd-mmm-yyyy-HH'),'h',...
    datestr(clock, 'MM'),'m',datestr(clock,'ss'),'s');

if strcmp(batch_type,'Survive')
    s = ['SurviveBatch-' s]; 
end

% Save all data
% save([dPath filesep s],'B','-v7.3')

%% plot vectors
escape = permute(B.escape,[2 1 3]);
x0 = B.prey_x0;
y0 = B.prey_y0;
theta = permute(B.escape_angles,[2 1 3]);
[x,y]=meshgrid(x0,y0);
r = min(x0(2)-x0(1),y0(2)-y0(1)); % scaling factor for all vectors

for j = 1:length(K_vals)
    
    ind1 = find(escape(:,:,j)==1);
    ind2 = find(escape(:,:,j)==0);
    F = figure(j);
    title(['K = ' num2str(K_vals(j))])
    hold on
    thetaC = theta(:,:,j);
    
    % plot vectors of successful escape
    h1 = quiver(x(ind1),y(ind1),r*cos(thetaC(ind1)),r*sin(thetaC(ind1)),0);
    h1.Color = 'b';
    h1.LineWidth = 1;
    
    % plot vectors of capture
    h2 = quiver(x(ind2),y(ind2),r*cos(thetaC(ind2)),r*sin(thetaC(ind2)),0);
    h2.Color = 'r';
    h1.LineWidth = 1;
    
    axis image
    
    % save current figure
    fname = strcat(s, 'K=', num2str(K_vals(j)),'.eps');
    fpath = fullfile(dPath,fname);
%     print(F, '-depsc', fpath)

end



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
        num2str(time_left) ' min left'])

end


end




