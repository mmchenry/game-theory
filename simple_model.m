function R = simple_model(pred,prey,param)
% ODE Predator-prey interaction model


%% Parameters

% Scaling constants (helps the solver avoid crazy numbers)
sL = pred.head_length * 4;      % m
sT = pred.saccade_interval / 4; % s


%% Normalize input data

% PREY --------------------------------------------------------------------
p.prey.spd0         = prey.spd0             ./sL .* sT;
p.prey.x0           = prey.x0               ./sL;
p.prey.y0           = prey.y0               ./sL;
p.prey.theta0       = prey.theta0;
p.prey.spdEscape    = prey.spdEscape        ./sL .* sT;
p.prey.rotSpdEscape = prey.omegaEscape           .* sT;
p.prey.spdResp      = prey.spdResp          ./sL .* sT;
p.prey.lat          = prey.latEscape        ./sT;
p.prey.durEscape    = prey.durEscape        ./sT;

% Dimensions of head
p.prey.w            = prey.body_width       ./sL;
p.prey.L            = prey.head_length      ./sL;



% PREDATOR ----------------------------------------------------------------
% Initial body position & orientation 
p.pred.x0           = pred.x0               ./sL;
p.pred.y0           = pred.y0               ./sL;
p.pred.theta0       = pred.theta0;

% Dimensions of head
p.pred.w            = pred.body_width       ./sL;
p.pred.L            = pred.head_length      ./sL;

% Kinematic parameters
p.pred.spd0         = pred.spd0             ./sL .* sT;
p.pred.sccd_prd     = pred.saccade_period   ./ sT;
p.pred.sccd_intvl   = pred.saccade_interval ./ sT;
p.pred.sccd_omega   = pred.saccade_omega    .* sT;
p.pred.wall_omega   = pred.wall_omega       .* sT;

% Sensory parameters
p.pred.fldSize      = pred.fieldSize;
p.pred.vis_az       = pred.vis_az;
p.pred.verg_ang     = pred.verg_ang;


% INTERACTION -------------------------------------------------------------
p.both.dist_thresh  = param.dist_thresh     ./sL; 
p.both.d_capture    = param.d_capture       ./sL;


% SOLVER ------------------------------------------------------------------
p.t_span            = param.t_span          ./sT;
p.tank_rad          = param.tank_radius     ./sL;


%% Configure solver

% Solver options
options  = odeset('RelTol',1e-3,...
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

% Store input parameters
R.pred   = pred;
R.prey   = prey;
R.param  = param;

% clear X t t_event

%% Solver

function [t,X,t_event] = solver(p,options)
% Numerical integration of the trajectory of predator and prey

X0(1,1) = p.prey.x0;
X0(2,1) = p.prey.y0;
X0(3,1) = p.prey.theta0;
X0(4,1) = p.pred.x0;
X0(5,1) = p.pred.y0;
X0(6,1) = p.pred.theta0;

% Predator foraging parameters ---------------------------

% Initial predator saccade time
t_PredSccd = p.t_span(1);

% Initial saccade direction (randomly selected)
dir_PredSccd = rand(1);
dir_PredSccd = (dir_PredSccd-.5)./abs(dir_PredSccd-.5);

% Whether currently executing a saccade
on_PredSccd = 0;

% Whether predator is currently tracking a prey
seePrey = 0;

% Prey routine swimming parameters ----------------------

% Initial predator saccade time
t_PreySccd = p.t_span(1);

% Initial saccade direction (randomly selected)
dir_PreySccd = rand(1);
dir_PreySccd = (dir_PreySccd-.5)./abs(dir_PreySccd-.5);

% Whether currently executing a saccade
on_PreySccd = 0;

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

% Points describing prey head (local FOR)
angHead    = [linspace(-pi/2, pi/2, 500)]';
xHeadPreyL = p.prey.L    .* cos(angHead) - p.prey.L;
yHeadPreyL = (p.prey.w) .* sin(angHead);

% Logical that indicates whether the fish is presently escaping 
escapeOn = 0;

% Initiate stimTime variable
stimTime = nan;

% Set initial direction of escape 
dirEsc = 1;

% Runge-Kutta solver
[t,X,t_event] = ode45(@gov_eqn,p.t_span,X0,options);

    function dX = gov_eqn(t,X)
    % ODE for the dynamics of the system
        
        %  INPUTS ---------------------------------------------------------
        
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
        

        % Distance between prey and predator
        dist = hypot(xPred-xPrey,yPred-yPrey);

        % DECISIONS ABOUT RATES OF CHANGE ---------------------------------
                        
        % Transform head points of predator into global FOR                 
        [xHead, yHead] = coord_trans('body to global', ...
                      thetaPred, [xPred yPred], xHeadPredL, yHeadPredL);   
                            
        % Transform region outside of head points of predator into global FOR                    
        [xCollPred, yCollPred] = coord_trans('body to global', ...
                      thetaPred, [xPred yPred], xCollPredL, yCollPredL); 

        
        % Predator speed (fixed)
        spdPred = p.pred.spd0;
        
       % Prey speed         
        if dist < p.both.dist_thresh
            spdPrey = p.prey.spdResp;
        else
            spdPrey = p.prey.spd0;
       end
        
        % Prey behavior -------------------------------------------
        
        % Routine swimming 
        if ~escapeOn 
            [omegaPrey, t_PreySccd, dir_PreySccd, on_PreySccd] = prey_routine(...
                     xCollPred, yCollPred, thetaPred, p.prey, t_PreySccd, t, ...
                     dir_PreySccd, on_PreySccd, p.tank_rad, 'prey');
        end
        
        % Escape response
        
        % Check if in an escape response and within 4 x capture distance
        if ~escapeOn && (dist < 4*p.both.d_capture)

            % Indicate that we are in an escape response
            escapeOn = 1;
            
            % Note the time of its start
            stimTime = t;
            
            % Set escape response direction 
            dirEsc = rand(1);
            dirEsc = (dirEsc-.5)./abs(dirEsc-.5);

            % Determine speed of rotation of escape 
            % normrnd(m,s) returns rndm # from norm dist. w/ mean=m,stdev=s
            rotSpdEscape = normrnd(p.prey.rotSpdEscape/sT, 1.5);
            p.prey.rotSpdEscape = rotSpdEscape .* sT;
        end
        
        % If we are beyond the escape duration . . .
        if escapeOn && ((t-stimTime-p.prey.lat) > p.prey.durEscape)
            escapeOn = 0;
        end
        
        % If we are within the period of an escape . . .          
        if escapeOn 
            %stimTime
            [omegaPrey, spdPrey] = prey_escape(t, stimTime,...
                                               p.prey, 0, 0, dirEsc);
        else
            omegaPrey = 0;
        end 
        

        % Predator behavior ----------------------------------------
        
        % Determine whether prey is detected
        if ~seePrey
            %seePrey = find_prey(xPrey, yPrey, thetaPred, xPred, yPred);
        end

        
        if 1 
            [omegaPred, t_PredSccd, dir_PredSccd, on_PredSccd] = foraging(...
                     xCollPred, yCollPred, thetaPred, p.pred, t_PredSccd, t, ...
                     dir_PredSccd, on_PredSccd, p.tank_rad, 'predator');
%         else
%             [spdPred,OmegaPred] = strike(xPred,yPred,thetaPred,...
%                 xPrey,yPrey,thetaPrey);
%             %TODO: create this m-file
        end
        
        % Yaw rate of prey
%         omegaPrey = 0;
        
        
    
        % OUTPUTS --------------------------------------------------------- 
        
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
        
    end
end


% capture function: use with 'Event' option in ode45
function [value, isterminal, direction] = capture_fnc(~,X)
    % the event occurs when distance is less than capture distance
    distance = hypot(X(4)-X(1),X(5)-X(2));
    if distance < p.both.d_capture
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