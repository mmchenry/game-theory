function R = simple_model(pred,prey,param)
% ODE Predator-prey interaction model (no longer so "simple")


%% Parameters

% Scaling constants (helps the solver avoid crazy numbers)
sL = pred.bod_len;              % m
sT = pred.saccade_interval / 4; % s


%% Normalize input data

% PREY --------------------------------------------------------------------
% Initial body position & orientation 
p.prey.x0           = prey.x0               ./sL;
p.prey.y0           = prey.y0               ./sL;
p.prey.theta0       = prey.theta0;

% Prey escape parameters
p.prey.spdEscape    = prey.spdEscape        ./sL .* sT;
p.prey.rotSpdEscape = prey.omegaEscape           .* sT;
p.prey.spdResp      = prey.spdResp          ./sL .* sT;
p.prey.lat          = prey.latEscape        ./sT;
p.prey.durEscape    = prey.durEscape        ./sT;

% Morphometrics
p.prey.bod_width    = prey.bod_width        ./sL;
p.prey.bod_len      = prey.bod_len          ./sL;
p.prey.COM_len      = prey.COM_len          ./sL;
p.prey.numBodPts    = prey.numBodPts;

% Kinematic parameters
p.prey.spd0         = prey.spd0             ./sL .* sT;
p.prey.sccd_prd     = prey.saccade_period   ./ sT;
p.prey.sccd_intvl   = prey.saccade_interval ./ sT;
p.prey.sccd_omega   = prey.saccade_omega    .* sT;
p.prey.wall_omega   = prey.wall_omega       .* sT;

% Sensory parameters
p.prey.fieldSize      = prey.fieldSize       ./sL;


% PREDATOR ----------------------------------------------------------------
% Initial body position & orientation 
p.pred.x0           = pred.x0               ./sL;
p.pred.y0           = pred.y0               ./sL;
p.pred.theta0       = pred.theta0;

% Morphometrics
p.pred.bod_width       = pred.bod_width     ./sL;
p.pred.bod_len         = pred.bod_len       ./sL;
p.pred.COM_len         = pred.COM_len       ./sL;
p.pred.numBodPts       = pred.numBodPts;

% Kinematic parameters
p.pred.spd0                 = pred.spd0             ./sL .* sT;
p.pred.sccd_prd             = pred.saccade_period   ./ sT;
p.pred.sccd_intvl           = pred.saccade_interval ./ sT;
p.pred.sccd_omega           = pred.saccade_omega    .* sT;
p.pred.wall_omega           = pred.wall_omega       .* sT;

% Sensory parameters
p.pred.fieldSize    = pred.fieldSize  ./sL ;
p.pred.vis_az       = pred.vis_az;
p.pred.verg_ang     = pred.verg_ang;
p.pred.rtnl_den     = pred.rtnl_den;
p.pred.vis_freq     = pred.vis_freq   .*sT;
p.pred.vis_thresh   = pred.vis_thresh;

% Strike parameters
p.pred.strike_thresh   = pred.strike_thresh ./ sL;
p.pred.strike_dur      = pred.strike_dur    ./ sT;
p.pred.strike_reach    = pred.strike_reach  ./ sL;
p.pred.strike_range    = pred.strike_range;

% GENERAL -----------------------------------------------------------------
p.param.dist_thresh     = param.dist_thresh     ./sL; 
p.param.d_capture       = param.d_capture       ./sL;
p.param.tank_radius     = param.tank_radius     ./sL;


% SOLVER ------------------------------------------------------------------
p.t_span            = param.t_span          ./sT;
p.rel_tol           = param.rel_tol;
p.abs_tol           = param.abs_tol;

clear pred prey param


%% Configure solver

% Solver options
options  = odeset('RelTol',p.rel_tol,...
                  'AbsTol',p.abs_tol, ...
                  'Events',@capture_fnc, ...
                  'MaxStep',p.pred.sccd_prd/2);
              
                         
%% Solve & save results in SI units


% Run solver
[t,X,t_event] = solver(p,options);

% Results stored in SI units 
R.t = t .* sT;
R.tEnd = t_event .* sT;

R.xPrey         = X(:,1)    .* sL;
R.yPrey         = X(:,2)    .* sL;
R.thetaPrey     = X(:,3);
R.xPred         = X(:,4)    .* sL;
R.yPred         = X(:,5)    .* sL;
R.thetaPred     = X(:,6);

% Store scaling parameters
R.sT     = sT;
R.sL     = sL;


%% Solver

function [t,X,t_event] = solver(p,options)
% Numerical integration of the trajectory of predator and prey

X0(1,1) = p.prey.x0;
X0(2,1) = p.prey.y0;
X0(3,1) = p.prey.theta0;
X0(4,1) = p.pred.x0;
X0(5,1) = p.pred.y0;
X0(6,1) = p.pred.theta0;

% Intialize global status on a capture
global captured
captured = 0;

% Intitialize state variables
s = give_behavior('Initialize', p);


% RUN SOLVER --------------------------------------------------------------

% Runge-Kutta solver
[t,X,t_event] = ode45(@gov_eqn,p.t_span,X0,options);

    function dX = gov_eqn(t,X)
    % ODE for the dynamics of the system
    
        % Find body position of both fish
        s = give_behavior('Body positions', p, s, t, X);   
        
        % BEHAVIORAL CHANGES ----------------------------------------------

        % Set prey behavior
        s = give_behavior('Prey', p, s, t, X);
        
        % Set predator behavior
        s = give_behavior('Predator', p, s, t, X);
        
        % Determine whether prey is captured
        s = give_behavior('Capture', p, s, t, X);
        
        % Update global 'captured' variable
        captured = s.captured;

        
        % OUTPUTS --------------------------------------------------------- 
        
        % Prey velocity in x & y directions
        dX(1,1) = s.prey.spd * cos(s.prey.theta);
        dX(2,1) = s.prey.spd * sin(s.prey.theta);
        
        % Prey rate of rotation
        dX(3,1) = s.prey.omega;
       
        % Predator velocity in x & y directions
        dX(4,1) = s.pred.spd * cos(s.pred.theta);
        dX(5,1) = s.pred.spd * sin(s.pred.theta);
        
        % Predator rate of rotation 
        dX(6,1) = s.pred.omega;

    end
end


%% CAPTURE FUNCTION: use with 'Event' option in ode45 ----------------------
function [value, isterminal, direction] = capture_fnc(~,X)
    
    % Declare global variable
    global captured
    
    % the event occurs when distance is less than capture distance
    %distance = hypot(X(4)-X(1),X(5)-X(2));
    if captured
        %distance < p.param.d_capture
        value      = 0;
        isterminal = 1;         % tells ode45 to stop integration
        direction  = 0;
    else
        value = 1;
        isterminal = 0;         
        direction  = 0;
   end
end


end