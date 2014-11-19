function simple_model
% This is a simple predator-prey interaction model


%% Parameters

% Scaling constants (helps the solver avoid crazy numbers)
sL = 2e-2;      % m
sT = 10e-3;     % s

% Prey parameters -----------------------------------------
prey.spd0           = 0;         % m/s          
prey.x0             = 1e-2;      % m          
prey.y0             = 1e-2;      % m          
prey.theta0         = 45/180*pi; % rad            
prey.spdEscape      = 5e-2;      % m/s  

% Predator parameters -------------------------------------
pred.spd0               = 2e-2;      % m/s  
pred.x0                 = 0;         % m      
pred.y0                 = 0;         % m      
pred.theta0             = 0;         % rad  

pred.saccade_interval   = 3;   
pred.saccade_period     = 0.5;       % s
pred.saccade_omega      = 2;         % rad/s

pred.body_width         = 0.3e-2;    % m
pred.head_length        = 0.4e-2;    % m 

pred.fieldSize          = 1.2;

% Capture distance (i.e. distance where pred. captures prey)
param.d_capture     = 5e-3  ;    % m

% Proximity at which prey responds to predator
param.dist_thresh   = 1e-2;      % m

% Time span for simulation
param.t_span        = [0 30];     % s

% Tank dimension
param.tank_radius   = 10e-2;     % m


%% Normalize input data

% Prey parameters -----------------------------------------
p.prey.spd0         = prey.spd0             ./sL .* sT;
p.prey.x0           = prey.x0               ./sL;
p.prey.y0           = prey.y0               ./sL;
p.prey.theta0       = prey.theta0;
p.prey.spdEscape    = prey.spdEscape        ./sL .* sT;

% Predator parameters -------------------------------------
p.pred.w            = pred.body_width       ./sL;
p.pred.L            = pred.head_length      ./sL;
p.pred.x0           = pred.x0               ./sL;
p.pred.y0           = pred.y0               ./sL;
p.pred.theta0       = pred.theta0;
p.pred.spd0         = pred.spd0             ./sL .* sT;
p.pred.sccd_prd     = pred.saccade_period   ./ sT;
p.pred.sccd_intvl   = pred.saccade_interval ./ sT;
p.pred.sccd_omega   = pred.saccade_omega    .* sT;
p.pred.fldSize      = pred.fieldSize;

% Interaction parameters (i.e. both play a role) ---------
p.both.dist_thresh  = param.dist_thresh     ./sL; 
p.both.d_capture    = param.d_capture       ./sL;

% General parameters -------------------------------------
p.t_span            = param.t_span          ./sT;
p.tank_rad          = param.tank_radius     ./sL;


%% Configure solver

% Solver options
options  = odeset('RelTol',1e-5,...
                  'Events',@capture_fnc, ...
                  'MaxStep',p.pred.sccd_prd/2);
              
                         
%% Solve & save results in SI units


% TODO: Adjust saccades to hug walls

% Run solver
[t,X] = solver(p,options);

% Results stored in SI units 
R.t = t                 .* sT;
R.xPrey = X(:,1)        .* sL;
R.yPrey = X(:,2)        .* sL;
R.thetaPrey = X(:,3);
R.xPred = X(:,4)        .* sL;
R.yPred = X(:,5)        .* sL;
R.thetaPred = X(:,6);

% Store input parameters
R.pred   = pred;
R.prey   = prey;
R.param  = param;

clear X t


%% Display results

% Animate results with even time intervals
%animate_sim(R)

% Plot orientation angles
figure;
vis_results(R,'Turning data')

figure;
vis_results(R,'Trajectories')

ttt=3

%% Solver

function [t,X] = solver(p,options)
% Numerical integration of the trajectory of predator and prey

X0(1,1) = p.prey.x0;
X0(2,1) = p.prey.y0;
X0(3,1) = p.prey.theta0;
X0(4,1) = p.pred.x0;
X0(5,1) = p.pred.y0;
X0(6,1) = p.pred.theta0;
%X0(7,1) = p.tank_rad - hypot(p.pred.x0,p.pred.y0);

% Initial predator saccade time
t_PredSccd = p.t_span(1);

% Initial saccade direction (randomly selected)
dir_PredSccd = rand(1);
dir_PredSccd = (dir_PredSccd-.5)./abs(dir_PredSccd-.5);

% Whether currently executing a saccade
on_PredSccd = 0;

% Whether currently bumping into a wall
bumpPred = 0;

% Points describing boundary of tank (global FOR)
theta = [linspace(0,2*pi,1000)]';
xTank = p.tank_rad .* cos(theta);
yTank = p.tank_rad .* sin(theta);
clear theta

% Points describing predator head (local FOR)
angHead    = [linspace(-pi/2, pi/2, 500)]';
xHeadPredL = p.pred.L    .* cos(angHead) - p.pred.L;
yHeadPredL = (p.pred.w/2) .* sin(angHead);

% Points describing predator head (local FOR)
xCollPredL = p.pred.fldSize*p.pred.L    .* cos(angHead) - p.pred.L;
yCollPredL = p.pred.fldSize*(p.pred.w/2) .* sin(angHead);


% Runge-Kutta solver
[t,X] = ode45(@gov_eqn,p.t_span,X0,options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % DEFINE INPUTS ----------------------
        % Prey position
        xPrey = X(1);
        yPrey = X(2);
        
        % Prey orientation
        thetaPrey  = X(3);
        
        % Predator position
        xPred = X(4);
        yPred = X(5);

        % Predator orientation
        thetaPred  = X(6);
        
        % DECISIONS ABOUT RATES OF CHANGE ----------------------
        
        %Distance from predator
        %dist = hypot(xPred-xPrey,yPred-yPrey);
        
        % Tank in predator FOR
        %[xTankPredL, yTankPredL] = coord_trans(thetaPred, [xPred yPred], ...
        %                    xTank, yTank, 'global to body');
                        
        % Transform head points into global FOR                 
        [xHead, yHead] = coord_trans(thetaPred, [xPred yPred], ...
                                xHeadPredL, yHeadPredL, 'body to global');                
        [xCollPred, yCollPred] = coord_trans(thetaPred, [xPred yPred], ...
                                xHeadPredL, yHeadPredL, 'body to global'); 
                            
        % Prevent movement outside of boundary
        %[xPred, yPred, xHead, yHead] = enforce_tank(xPred, yPred, ...
        %                      xHead, yHead, xTank, yTank, p.tank_rad);
                            
        % Input: 'Smart' Predator orientation: always toward Prey 
        %thetaPred = acos((xPrey-xPred)/dist);
        %thetaPred = X(6);
                 
        % Distance from predator
        %dist = hypot(xPred-xPrey,yPred-yPrey);
        
        % Predator speed
        spdPred = p.pred.spd0;
        
%         % Prey behavior -------------------------------------------
%         if dist < p.both.dist_thresh
%             spdPrey = p.prey.spdEscape;
%         else
            spdPrey = p.prey.spd0;
%        end
        
        % Predator behavior ----------------------------------------
        if 1 
            % TODO: set sensory criteria for strike
            [omegaPred, t_PredSccd, dir_PredSccd, on_PredSccd] = foraging_pred(...
                     xCollPred, yCollPred, thetaPred, p.pred, t_PredSccd, t, ...
                     dir_PredSccd, on_PredSccd, p.tank_rad, 'predator');
%         else
%             [spdPred,OmegaPred] = strike_pred(xPred,yPred,thetaPred,...
%                 xPrey,yPrey,thetaPrey);
%             %TODO: create this m-file
        end
        
        % Predator rate of rotation: 'dumb' predator
%         dX(6,1) = thetaPred;
        
        
        % Yaw rate of prey
        omegaPrey = 0;
        
        % Yaw rate of predator
        %OmegaPred = 0;
    
        % DEFINE OUTPUTS ---------------------- 
        
        % Prey velocity in x & y directions
        dX(1,1) = spdPrey * cos(thetaPrey);
        dX(2,1) = spdPrey * sin(thetaPrey);
        
        % Prey rate of rotation
        dX(3,1) = omegaPrey;
        
        % Predator velocity in x & y directions
        dX(4,1) = spdPred * cos(thetaPred);
        dX(5,1) = spdPred * sin(thetaPred);
        
        % Predator rate of rotation 
        dX(6,1) = omegaPred;
        
        %[t dX(6) onSccd]
    end
end


% capture function: use with 'Event' option in ode45
function [value, isterminal, direction] = capture_fnc(~,X)
    % the event occurs when distance is less than capture distance
    distance = hypot(X(4)-X(1),X(5)-X(2));
%     if distance < p.both.d_capture
%         value      = 0;
%         isterminal = 1;         % tells ode45 to stop integration
%         direction  = 0;
%     else
        value = 1;
        isterminal = 0;         
        direction  = 0;
%    end
end




end