function simple_model
% This is a simple predator-prey interaction model


%% Parameters

% Scaling constants (helps the solver avoid crazy numbers)
sL = 1 %2e-2;      % m
sT = 1 %10e-3;     % s

% Initial conditions for the predator 
pred.spd0           = 6e-2;      % m/s  
pred.x0             = 0;         % m      
pred.y0             = 0;         % m      
pred.theta0         = 45/180*pi; % rad    

% Initial conditions and parameters for the prey
prey.spd0           = 0;         % m/s          
prey.x0             = 1e-2;      % m          
prey.y0             = 1e-2;      % m          
prey.theta0         = 45/180*pi; % rad
prey.spdEscape      = 9e-2;      % m/s  
prey.omegaEscape    = 10;       % rad/s
prey.spdResp        = 5e-2;      % m/s
prey.latEscape      = 4e-3;     % s
prey.durEscape      = 40e-3;     % s

% Capture distance (i.e. distance where predator captures prey)
param.d_capture     = 1e-3  ;    % m

% Proximity at which prey responds to predator
param.dist_thresh   = 1e-2;      % m

% Time span for simulation
param.t_span        = [0 5];     % s
              
              
%% Normalize input data

% Prey parameters
p.prey.spd0         = prey.spd0             ./sL .* sT;
p.prey.x0           = prey.x0               ./sL;
p.prey.y0           = prey.y0               ./sL;
p.prey.theta0       = prey.theta0;
p.prey.spdEscape    = prey.spdEscape        ./sL .* sT;
p.prey.rotSpdEscape = prey.omegaEscape           .* sT;
p.prey.spdResp      = prey.spdResp          ./sL .* sT;
p.prey.lat          = prey.latEscape        ./sT;
p.prey.durEscape    = prey.durEscape        ./sT;

% Predator parameters
p.pred.x0           = pred.x0               ./sL;
p.pred.y0           = pred.y0               ./sL;
p.pred.theta0       = pred.theta0;
p.pred.spd0         = pred.spd0             ./sL .* sT;

% Interaction parameters (i.e. both play a role)
p.both.dist_thresh  = param.dist_thresh     ./sL; 
p.both.d_capture    = param.d_capture       ./sL;
p.t_span            = param.t_span          ./sT;

%clear param prey pred 


%% Configure solver

% Solver options
options  = odeset('RelTol',1e-5,...
                  'Events',@capture_fnc, ...
                  'MaxStep',p.prey.durEscape/4);
              

%% Solve & post-processing

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

% Alberto -- Sorry, I broke your definition of R.tEnd.  Please put it back
% for me

% if isempty(R.tEnd)
%     disp('prey successfully escaped')
% else
%     disp(['the prey was captured at time = ' num2str(R.tEnd)])
% end


%% Display results

% Animate results with even time intervals
%animate_sim(R)

% Plot orientation angles
figure;
vis_results(R,'Turning data')

figure;
vis_results(R,'Trajectories')


%figure
%plot(R.xPrey,R.yPrey,'-b',R.xPrey(end),R.yPrey(end),'ob',...
%        R.xPred,R.yPred,'-r',R.xPred(end),R.yPred(end),'or')



%% Functions

function [t,X, teout] = solver(p,options)
% Numerical integration of the trajectory of predator and prey

X0(1,1) = p.prey.x0;
X0(2,1) = p.prey.y0;
X0(3,1) = p.prey.theta0;
X0(4,1) = p.pred.x0;
X0(5,1) = p.pred.y0;
X0(6,1) = p.pred.theta0;

% Logical that indicates whether the fish is presently escaping 
escapeOn = 0;

% Initiate stimTime variable
stimTime = nan;

[t,X, teout] = ode45(@gov_eqn,p.t_span,X0,options);


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
        dist = hypot(xPred-xPrey,yPred-yPrey);
        
        % Input: 'Smart' Predator orientation: always toward Prey 
        %thetaPred = acos((xPrey-xPred)/dist);
          
        % Predator speed
        spdPred = p.pred.spd0;
        
        % Prey speed
        if dist < p.both.dist_thresh
            spdPrey = p.prey.spdResp; 
        else
            spdPrey = p.prey.spd0;
        end
        
        
        % Alberto: Your control instructions need to act on the first
        % derivatives of your state variables. Otherwise, there is no
        % reason for the ODE solver.  Therefore the escape response
        % needs to determine omega, not theta. Also, this should occur
        % over some duration, rather than instantaneously.  I've taken a
        % stab at it to give you an idea of what I mean.
        
        % Initiate escape if current distance is twice capture distance
        if ~escapeOn && (dist < 4*p.both.d_capture)

            % Indicate that we are in an escape response
            escapeOn = 1;
            
            % Note the time of its start
            stimTime = t;
        end
        
        % If we are beyond the escape duration . . .
        if escapeOn && ((t-stimTime-p.prey.lat) > p.prey.durEscape)
            escapeOn = 0;
        end
        
        % If we are within the period of an escape . . .          
        if escapeOn 
            %stimTime
            [omegaPrey, spdPrey] = prey_escape(t, stimTime, p.prey, 0, 0);
        else
            omegaPrey = 0;
        end    
                %spdPrey     = p.prey.spdEscape;
        
        % Presently zero any changes in rotation for the predtaor 
        omegaPred = 0;         
       
                
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