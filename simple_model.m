function R = simple_model(pIn,d)
% ODE Predator-prey interaction model (no longer so "simple")


%% Parameters

% Scaling constants
sT = pIn.param.sT;
sL = pIn.param.sL;

% Convert parameter values, according to scaling constants
p = convert_params(pIn,d);


%% Configure solver

% Solver options
options  = odeset('RelTol',p.param.rel_tol,...
                  'AbsTol',p.param.abs_tol, ...
                  'Events',@capture_fnc, ...
                  'MaxStep',p.pred.sccd_prd/25);
                          
              
%% Solve & save results in SI units

% Run solver
[t,X,t_event] = solver(p,options);

% Results stored in SI units 
R.t             = t         .* sT;
R.X(:,1)        = X(:,1)    .* sL;
R.X(:,2)        = X(:,2)    .* sL;
R.X(:,3)        = X(:,3);
R.X(:,4)        = X(:,4)    .* sL;
R.X(:,5)        = X(:,5)    .* sL;
R.X(:,6)        = X(:,6);

% Same data, with field names
R.tEnd          = t_event   .* sT;
R.xPrey         = X(:,1)    .* sL;
R.yPrey         = X(:,2)    .* sL;
R.thetaPrey     = X(:,3);
R.xPred         = X(:,4)    .* sL;
R.yPred         = X(:,5)    .* sL;
R.thetaPred     = X(:,6);


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
[t,X,t_event] = ode45(@gov_eqn,p.param.t_span,X0,options);

    function dX = gov_eqn(t,X)
    % ODE for the dynamics of the system
    
        % Find body position of both fish
%         s = give_behavior('Body positions', p, s, t, X);   
        
        % BEHAVIORAL CHANGES ----------------------------------------------

        % Set prey behavior
%         s = give_behavior('Prey', p, s, t, X);
%         
%         % Set predator behavior
%         s = give_behavior('Predator', p, s, t, X);
%         
%         % Determine whether prey is captured
%         s = give_behavior('Capture', p, s, t, X);

        % Set predator and prey behavior for Weihs situation
 
        s = give_behavior('Weihs', p, s, t, X);
        
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