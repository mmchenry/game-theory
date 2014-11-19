function simple_model
% This is a simple predator-prey interaction model


%% Parameters

% Scaling constants (helps the solver avoid crazy numbers)
sL = 2e-2;      % m
sT = 10e-3;     % s

% Initial conditions for the predator 
pred.spd0           = 6e-2;      % m/s  
pred.x0             = 0;         % m      
pred.y0             = 0;         % m      
pred.theta0         = 0;         % rad    

% Initial conditions and parameters for the prey
prey.spd0           = 0;         % m/s          
prey.x0             = 1e-2;      % m          
prey.y0             = 1e-2;      % m          
prey.theta0         = 45/180*pi; % rad
prey.spdEscape      = 9e-2;      % m/s     
prey.spdResp        = 5e-2;      % m/s

% Capture distance (i.e. distance where predator captures prey)
param.d_capture     = 1e-3  ;    % m

% Proximity at which prey responds to predator
param.dist_thresh   = 1e-2;      % m

% Time span for simulation
param.t_span        = [0 5];     % s

% Solver options
options             = odeset('RelTol',1e-3,'Events',@capture_fnc);


%% Normalize input data

% Prey parameters
p.prey.spd0         = prey.spd0             ./sL .* sT;
p.prey.x0           = prey.x0               ./sL;
p.prey.y0           = prey.y0               ./sL;
p.prey.theta0       = prey.theta0;
p.prey.spdEscape    = prey.spdEscape        ./sL .* sT;
p.prey.spdResp      = prey.spdResp          ./sL .* sT;

% Predator parameters
p.pred.x0           = pred.x0               ./sL;
p.pred.y0           = pred.y0               ./sL;
p.pred.theta0       = pred.theta0;
p.pred.spd0         = pred.spd0             ./sL .* sT;

% Interaction parameters (i.e. both play a role)
p.both.dist_thresh  = param.dist_thresh     ./sL; 
p.both.d_capture    = param.d_capture       ./sL;
p.t_span            = param.t_span          ./sT;

clear param prey pred 


%% Solve & post-processing

% Run solver
[t,X,teout] = solver(p,options);

R.t = t                 .* sT;
R.tEnd = teout          .* sT;
R.xPrey = X(:,1)        .* sL;
R.yPrey = X(:,2)        .* sL;
R.thetaPrey = X(:,3);
R.xPred = X(:,4)        .* sL;
R.yPred = X(:,5)        .* sL;
R.thetaPred = X(:,6);

clear X t teout


%% Display results

animate_sim(R)

if isempty(R.tEnd)
    disp('prey successfully escaped')
else
    disp(['the prey was captured at time = ' num2str(R.tEnd)])
end
%figure
%plot(R.xPrey,R.yPrey,'-b',R.xPrey(end),R.yPrey(end),'ob',...
%        R.xPred,R.yPred,'-r',R.xPred(end),R.yPred(end),'or')

function [t,X, teout] = solver(p,options)
% Numerical integration of the trajectory of predator and prey

X0(1,1) = p.prey.x0;
X0(2,1) = p.prey.y0;
X0(3,1) = p.prey.theta0;
X0(4,1) = p.pred.x0;
X0(5,1) = p.pred.y0;
X0(6,1) = p.pred.theta0;

[t,X, teout] = ode45(@gov_eqn,p.t_span,X0,options);


    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Input: Prey position
        xPrey = X(1);
        yPrey = X(2);
        
        % Input: Prey orientation
        thetaPrey  = X(3);     
        
        % Input: Predator position
        xPred = X(4);
        yPred = X(5);
        
        % Distance from predator
        dist = hypot(xPred-xPrey,yPred-yPrey);
     
        % Input: Predator orientation: 'dumb' predator
%         thetaPred  = X(6);

        % Input: 'Smart' Predator orientation: always toward Prey 
        thetaPred = acos((xPrey-xPred)/dist);
          
        % Predator speed
        spdPred = p.pred.spd0;
        
        % Prey speed
        if dist < p.both.dist_thresh
            spdPrey = p.prey.spdResp; 
        else
            spdPrey = p.prey.spd0;
        end
        
        % Initiate escape if current distance is twice capture distance
        if dist < 2*p.both.d_capture
            [thetaPrey] = prey_escape(thetaPrey);
            spdPrey     = p.prey.spdEscape;
        else
        end
        
        % Yaw rate  of prey
        OmegaPrey = 0;
        
        % Yaw rate of predator
        OmegaPred = 0;
    
        % Output: Prey velocity in x & y directions
        dX(1,1) = spdPrey * cos(thetaPrey);
        dX(2,1) = spdPrey * sin(thetaPrey);
        
        % Output: Prey rate of rotation
%         dX(3,1) = thetaPrey;
        dX(3,1) = 0;                
        
        % Output: Predator velocity in x & y directions
        dX(4,1) = spdPred * cos(thetaPred);
        dX(5,1) = spdPred * sin(thetaPred);

        % Output: Predator rate of rotation: 'dumb' predator
%         dX(6,1) = thetaPred;
        
        % Output: Predtator rate of rotation for 'smart' predator
        dX(6,1) = 0;
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