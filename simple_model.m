function simple_model
% This is a simple predator-prey interaction model


%% Parameters

% Initial conditions for the predator 
pred.spd0           = 1;        
pred.x0             = 0;        
pred.y0             = 0;        
pred.theta0         = 0;        

% Initial conditions and parameters for the prey
prey.spd0           = 0;            
prey.x0             = 1;            
prey.y0             = 0.5;          
prey.theta0         = 45/180*pi;    
prey.dist_thresh    = 1.5;          
prey.spdEscape      = 2;            

% Capture distance
param.d_capture     = 0.05;          % distance where pred. captures prey

% For the solver
param.t_span        = [0 5];

% Solver options
options             = odeset('RelTol',1e-3,'Events',@capture_fnc);


%fps = 30;

%% Solve & post-processing

% Run solver
[t,X] = solver(pred,prey,param,options);

R.t = t;
R.xPrey = X(:,1);
R.yPrey = X(:,2);
R.thetaPrey = X(:,3);
R.xPred = X(:,4);
R.yPred = X(:,5);
R.thetaPred = X(:,6);

clear X t

figure

plot(R.xPrey,R.yPrey,'-b',R.xPrey(end),R.yPrey(end),'ob',...
        R.xPred,R.yPred,'-r',R.xPred(end),R.yPred(end),'or')

function [t,X] = solver(pred,prey,param,options)
% Numerical integration of the trajectory of predator and prey

X0(1,1) = prey.x0;
X0(2,1) = prey.y0;
X0(3,1) = prey.theta0;
X0(4,1) = pred.x0;
X0(5,1) = pred.y0;
X0(6,1) = pred.theta0;

[t,X] = ode45(@gov_eqn,param.t_span,X0,options);


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
        spdPred = pred.spd0;
        
        % Prey speed
        if dist < prey.dist_thresh
            spdPrey = prey.spdEscape;
        else
            spdPrey = prey.spd0;
        end
        
        % Yaw rate  of prey
        OmegaPrey = 0;
        
        % Yaw rate of predator
        OmegaPred = 0;
    
        % Output: Prey velocity in x & y directions
        dX(1,1) = spdPrey * cos(thetaPrey);
        dX(2,1) = spdPrey * sin(thetaPrey);
        
        % Output: Prey rate of rotation
        dX(3,1) = thetaPrey;
        
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
    if distance < param.d_capture
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