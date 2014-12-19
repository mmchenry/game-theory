function [max_min_d, theta_min]= run_batch_sim
% RUN_BATCH_SIM runs batch simulations of 'simple_model' for Weihs
% situation


%% Code execution

% Re-rerun batch
run_batch = 1;

% Visualize simulations results during execution
vis_sims = 1;


%% Batch parameters

% K values
num_K = 10;
min_K = 1;
max_K = 100;

% Angle values
num_ang = 100;
min_ang = 0;
max_ang = 120;

% Response distances
num_dist = 10;
min_dist = 1e-3;
max_dist = 3e-2;


%% Path definitions

% Identify Matt's computer
if isdir('/Users/mmchenry/Dropbox/Literature')
    root = '/Users/mmchenry/Dropbox/Literature/Review with Alberto/predprey data';
    
% Or Alberto's
else    
    root = '/Users/alberto/Dropbox/Review with Alberto/PredPrey Data'; 
end

% Path to data
dPath = [root];


%% Simulation parameters

% Default parameter values
[p,d] = default_params;

% Relative tolerance
p.param.rel_tol  = 1e-11; 
 
% Absolute tolerance
p.param.abs_tol  = 1e-11;

% Time span for simulation
p.param.t_span        = [0 2];      % s

% Make tank way too big
p.param.tank_radius = 1;            % m  

% prey initial position & speed
p.prey.x0 = 0.25e-2;                % m
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
B.prey_x0 = linspace(min_dist, max_dist, num_dist)';

% Values for K= predSpd/preySpd (gets coarse as values increase)
B.K = 10.^linspace(log10(min_K), log10(max_K),num_K)';

% Angle of escape 
% prey_theta = 60 ./180*pi;              % test value
B.prey_theta = linspace(min_ang, max_ang, num_ang)' ./180*pi;    % rad

% Store parameters
B.p = p;
B.d = d;


%% Run Simulation
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
    clrs = repmat(get(gca,'ColorOrder'), ceil(length(B.prey_x0)/7),1);
end

% Loop through initial distance
for i = 1:length(B.prey_x0)
    
    % set current initial distance
    p.prey.x0 = B.prey_x0(i);
    
    % Loop thru escape speeds
    for j = 1:length(B.K)
        
        % Set current escape speed
        p.prey.spdEscape = p.pred.spd0./B.K(j);
        
        
        theta_maxmin = fminbnd(@(theta) simFunction(theta, p, d),0,120/180*pi);
        
        
%         % Run simulations for each angle of rotation
%         for k = 1:length(B.prey_theta)
%             
%             % set current speed of rotation
%             p.prey.rotSpdEscape = B.prey_theta(k) ./ p.prey.durEscape;
%             
%             % Run Simulation
%             R   = simple_model(p,d);
%             
%             % Distance between predator and prey
%             dist = hypot(R.xPred-R.xPrey,R.yPred-R.yPrey);
%             
%             % Get minimum
%             min_dist(k,1) = min(dist);
%             
%             clear R dist
%         end
        
        % Find the max_min distance over all angles; find corresponding angle
        %[max_min_d, theta_ind]  = max(min_dist);
        B.theta_maxmin(j,i)     = theta_maxmin;

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
            ylim([0 120])
        end
        
        % Update status
        update_time(tStart, length(B.K)*(i-1) + j, ...
                    length(B.K)*length(B.prey_x0), ' ');
                
        clear k min_dist theta_ind max_min_d
    end
    
    % Update outer loop
    %update_time(tStart, i, length(B.prey_x0), '');

end

% Date and Time string for file name
% Example output: 18-Dec-2014-11h12m27s
s = strcat(datestr(clock,'dd-mmm-yyyy-HH'),'h',...
    datestr(clock, 'MM'),'m',datestr(clock,'ss'),'s');

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
h(1) = semilogx(K_ana,theta_weihs,'LineWidth',1);

l_txt{1} = 'Weihs';

% Get color palette
%clrs = repmat(get(gca,'ColorOrder'), 1, ceil(length(B.prey_x0)/7));
hold on

% run through S and plot optimal angles
for i = 1:size(B.theta_maxmin, 2)
    figure(F)
    h(i+1) = semilogx(B.K, B.theta_maxmin(:,i).*180/pi,'o-');
    l_txt{i+1} = [num2str(B.prey_x0(i)*100) ' cm'];
    clr = get(h(i+1),'Color');
    set(h(i+1),'MarkerFaceColor',clr);
end

xlabel('Relative predator speed')
ylabel('Optimal Escape Angle')
ylim([0 120])
legend(l_txt)




function y = simFunction(theta, p, d)

    % set current speed of rotation
    p.prey.rotSpdEscape = theta ./ p.prey.durEscape;

    % Run Simulation
    R   = simple_model(p,d);

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





