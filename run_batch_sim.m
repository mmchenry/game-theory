function run_batch_sim(batch_type)
% RUN_BATCH_SIM runs batch simulations of 'simple_model'.  At each
% parameter value, the optimization attempts to find the escape angle that
% maximizes the minimum distance.


%% Code execution

% Re-rerun batch
run_batch = 1;

% Visualize simulations results during execution
vis_sims = 1;

% Default batch_type
if (nargin<1)
    batch_type = 'Distance';
end

opt_type = 'Distance';
%opt_type = 'Survival';


%% Batch parameters

% K values
num_K = 1;
min_K = 2;
max_K = 100;

% Angle values
min_ang = 0;
max_ang = 90;

if (nargin<1) || strcmp(batch_type,'Distance')
    % Response distances (m)
    num_dist = 5;
    min_dist = 1e-3;
    max_dist = 3e-2;
     
elseif strcmp(batch_type,'Acceleration')
    % Escape durations (s)
    num_dur = 5;
    min_dur = 5e-3;
    max_dur = 300e-3;
    
elseif strcmp(batch_type,'Strike thresh')
    % Strike threshold (m)
    num_thresh  = 3;
    min_thresh  = 0.5e-2;
    max_thresh  = 1.1e-2;
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
    root = '/Users/alberto/Dropbox/Review with Alberto/PredPrey Data'; 
end

% Path to data
dPath = root;


%% General parameters

% Default parameter values
[p,d] = default_params;

% Time span for simulation
p.param.t_span   = [0 4];      % s

% Make tank way too big
p.param.tank_radius = 1;            % m  

% prey initial position & speed
p.prey.x0     = 0.25e-2;            % m
p.prey.y0     = 0;                  % m
p.prey.spd0   = 0;                  % m/s
p.prey.theta0 = 0;                  % rad  
p.prey.thrshEscape    = 2e-2;       % m

% prey escape duration (time to reach escape angle)
p.prey.durEscape = 25e-3;           % s

% predator initial position & speed
p.pred.x0 = 0e-2;                  % m
p.pred.y0 = 0;                      % m
p.pred.spd0 = 1e-2;                 % m/s
p.pred.theta0 = 0;                  % rad

% Values for K= predSpd/preySpd (gets coarse as values increase)
B.K = 10.^linspace(log10(min_K), log10(max_K),num_K)';


%% Parameters for each batch type

if strcmp(batch_type,'Distance')
    % Values for prey initial position (acts as response distance)
    B.prey_x0         = linspace(min_dist, max_dist, num_dist)';
    
    % Generic values to vary in a batch
    batch_vals        =  B.prey_x0;
    
    % Simulation type
    p.param.sim_type  = 'Weihs';
    
    % Relative tolerance
    p.param.rel_tol  = 1e-9;
    
    % Absolute tolerance
    p.param.abs_tol  = 1e-9;
    
elseif strcmp(batch_type,'Acceleration')
    
    B.escape_dur     = linspace(min_dur, max_dur, num_dur); 
    batch_vals       = B.escape_dur;
    p.param.sim_type = 'Weihs, acceleration';
    p.prey.spd0      = 0;
    
    % Relative tolerance
    p.param.rel_tol  = 1e-11;
    
    % Absolute tolerance
    p.param.abs_tol  = 1e-11;
    
elseif strcmp(batch_type,'Strike thresh')
    
    % Strike parameters
    B.strike_thresh     = linspace(min_thresh, max_thresh, num_thresh);
    batch_vals          = B.strike_thresh;
    p.param.sim_type    = 'Simple strike';
    
    % Prey parameters 
    p.prey.spd0         = 0;
    p.prey.thrshEscape  = 0.5e-2; % m
    p.prey.x0           = 2e-2; % m
    p.prey.theta0       = 0;
    p.prey.lat          = 0;
    
    % Scaling parameters
    p.param.sL      = 1;   % m
    p.param.sT      = 1;      % s
    
    % Predator parameters
    p.pred.theta0       = 0;
    p.pred.spd0         = 10e-2; %m/s
    
    % Relative tolerance
    p.param.rel_tol  = 1e-9;
    
    % Absolute tolerance
    p.param.abs_tol  = 1e-9;
    
end

% Store parameters
B.p = p;
B.d = d;


%% Run Simulations
% outermost loop runs through initial distances
% middle loop runs through speeds, inner loop runs through angles

if run_batch

% Initialize timer
tStart = tic;

% Initialize plot
if vis_sims
    f = figure('DoubleBuffer','on');
    h = semilogx(B.K, 0.*B.K, 'o');
    set(h,'Color',0.5.*[1 1 1])
    hold on
    % Get color palette
    clrs = repmat(get(gca,'ColorOrder'), ceil(length(batch_vals)/7),1);
end

% Loop through initial distance
for i = 1:length(batch_vals)
    
    if strcmp(batch_type,'Distance')
        % set current initial distance
        p.prey.x0 = batch_vals(i);
    
    elseif strcmp(batch_type,'Acceleration')
        % Set current escape duration (to vary acceleration)
        p.prey.durEscape = batch_vals(i);
        p.param.max_step = p.prey.durEscape/2;
        
    elseif strcmp(batch_type,'Strike thresh')
        % Set strike threshold
        p.pred.strike_thresh = batch_vals(i);
        
    end  
    
    % Loop thru escape speeds
    for j = 1:length(B.K)
        
        % Set current approach speed
        p.prey.spdEscape = p.pred.spd0./B.K(j);
        %p.pred.spd0 = p.prey.spdEscape * B.K(j);
        
        % Check if min dist fails to vary with theta
        if simFunction(pi/8, p, d) == simFunction(3.3*pi/4, p, d)
            warning('Failed optimization: Theta has no effect on min distance');
            min_dist       = nan;
            theta_maxmin   = nan;
            capture        = nan;
        
        % If min dist does vary with theta . . .
        else          
            % Run requested optimization
            if strcmp(opt_type, 'Distance')
                % Optimization for finding the optimal theta
                theta_maxmin = fminbnd(@(theta) simFunction(theta, p, d), ...
                                       0,max_ang/180*pi);
            else
                error('Do not recognize requested opt_type')
            end
            
            % Set rotation speed according to optimal theta
            p.prey.rotSpdEscape = theta_maxmin ./ p.prey.durEscape;
            
            % Run simulation
            R   = simple_model(p, d);
            
            % Log capture
            capture = R.capture;
            
            % Reconstruct simulation data (for troubleshooting)
            %R = reconstruct(R, p, d);
            
            % Distance between predator and prey
            min_dist = min(hypot(R.xPred-R.xPrey,R.yPred-R.yPrey));
        end
        
        % Store optimal theta & minimum distance
        B.theta_maxmin(j,i)     = theta_maxmin;
        B.min_dist(j,i)         = min_dist;
        B.capture(j,i)          = capture;
            
        % Update plot
        if vis_sims
            figure(f)
            if (j > 1) , delete(h(j-1)); end
            h(j) = semilogx(B.K(1:j), B.theta_maxmin(1:j,i) .*180/pi,'or-');
            set(h(j),'MarkerFaceColor',clrs(i,:))
            set(h(j),'MarkerEdgeColor',clrs(i,:))
            set(h(j),'Color',clrs(i,:))
            pause(0.001)
            xlabel('K')
            ylabel('Theta maxmin')
            ylim([0 max_ang])
        end
        
        % Update status
        update_time(tStart, length(B.K)*(i-1) + j, ...
                    length(B.K)*length(batch_vals), ' ');
                
        clear k min_dist theta_ind max_min_d
    end
    
    % Update outer loop
    %update_time(tStart, i, length(B.prey_x0), '');

end

close(f)

% Date and Time string for file name
% Example output: 18-Dec-2014-11h12m27s
s = strcat(datestr(clock,'dd-mmm-yyyy-HH'),'h',...
    datestr(clock, 'MM'),'m',datestr(clock,'ss'),'s');

if strcmp(batch_type,'Distance')
    s = ['Dist batch ' s];
    
elseif strcmp(batch_type,'Acceleration')
    s = ['Accel batch ' s];
    
end

% Save all data
save([dPath filesep s],'B','-v7.3')

clear i h

else
    
    % Clear unneeded variables
    clear B p d
    
    % Select data file
    cd(dPath)
    [fName, dPath] = uigetfile('*.mat','Choose data file');
    
    % Load data
    load([dPath filesep fName]);
    
end

        
%% plot Weihs optimal angles and optimal angles from simulations


% Make smooth K values
K_ana = 10.^linspace(log10(min_K), log10(max_K), 500);

% optimal angles from Weihs (theta=0 for K<1, use eq. 42 for K>=1)
theta_weihs = [0.* K_ana(K_ana<1), acos(1./K_ana(K_ana>=1))*180./pi]';

F = figure;
subplot(2,1,2)
h(1) = semilogx(K_ana,theta_weihs,'LineWidth',1);

l_txt{1} = 'Weihs';

% Get color palette
%clrs = repmat(get(gca,'ColorOrder'), 1, ceil(length(B.prey_x0)/7));
hold on

% run through S and plot optimal angles
for i = 1:size(B.theta_maxmin, 2)
    figure(F)

    % Legend text (conditional)
    if strcmp(B.batch_type, 'Distance')
        l_txt{i+1} = [num2str(B.prey_x0(i)*100) ' cm'];
    elseif strcmp(B.batch_type, 'Acceleration')
        l_txt{i+1} = [num2str(B.escape_dur(i)*1000) ' ms'];
    end
    
    subplot(2,1,2)
    h(i+1) = semilogx(B.K, B.theta_maxmin(:,i).*180/pi,'o-');
    
    clr = get(h(i+1),'Color');
    set(h(i+1),'MarkerFaceColor',clr);
    
    subplot(2,1,1)
    hD(i+1) = semilogx(B.K, B.min_dist(:,i).*100,'o-');
    hold on
    set(hD(i+1),'MarkerFaceColor',clr);
    set(hD(i+1),'MarkerEdgeColor',clr);
    set(hD(i+1),'Color',clr);
end

subplot(2,1,1)
xlabel('Relative predator speed')
ylabel('Minimum distance (cm)')

subplot(2,1,2)
xlabel('Relative predator speed')
ylabel('Optimal Escape Angle')
ylim([min_ang max_ang])
legend(l_txt)




function y = simFunction(theta, p, d)
% Minimize this function to find the angle theta that maximizes the
% minimum distance
    % set current speed of rotation
    p.prey.rotSpdEscape = theta ./ p.prey.durEscape;
    
    % Run Simulation
    R   = simple_model(p, d);

    % Distance between predator and prey
    min_dist = min(hypot(R.xPred-R.xPrey,R.yPred-R.yPrey));

    % Take inverse of min dist b/c the optimization attempts to find the theta
    % that produces the smallest number 
    y = -min_dist;

end


function y = simFunction2(theta, p, d)
% Minimize this function to 

    % set current speed of rotation
    p.prey.rotSpdEscape = theta ./ p.prey.durEscape;
    
    % Run Simulation
    R   = simple_model(p, d);

    % Distance between predator and prey
    min_dist = min(hypot(R.xPred-R.xPrey,R.yPred-R.yPrey));

    % Take inverse of min dist b/c the optimization attempts to find the theta
    % that produces the smallest number 
    y = -min_dist;

end

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





