function s = give_behavior(action, p, s, t, X)
% Determines the behavior of each fish for a single instant of time
%
% action    - String indicating which code to run (1xn)
% p         - Parameter strcuture for simulation
% s         - Structure of behavioral state variables 
% t         - Current time (1x1)
% X         - Vector of simulation results (6x1)


% Unpack simulation state variables, if provided
if nargin > 4
    
    % Prey position
    xPrey      = X(1);
    yPrey      = X(2);
    
    % Prey orientation
    thetaPrey  = X(3);
    
    % Predator position
    xPred      = X(4);
    yPred      = X(5);
    
    % Predator orientation
    thetaPred  = X(6);
    
    % Distance between prey and predator
    dist = hypot(xPred-xPrey,yPred-yPrey);
    
end


% Set default state for 'captured' field
s.captured = 0;


% Execute major code called for
switch action 

case 'Initialize'
    
    % Check inputs
    if nargin < 2
        error('Need to provide 2 inputs for this action')
    end
    
       
    % PREDATOR FORAGING PARAMETERS --------------------------------------------
    
    % Initial predator saccade time
    s.pred.tSccd = p.param.t_span(1);
    
    % Initial saccade direction (randomly selected)
    dir_PredSccd = rand(1);
    s.pred.dirSccd = (dir_PredSccd-.5)./abs(dir_PredSccd-.5);
    
    % Whether currently executing a saccade
    s.pred.onSccd = 0;
    
    % Time when prey enters predator's FOV
    s.pred.inFieldTime = nan;
    
    % Time when strike initiated (nan means strike not happening)
    s.pred.strikeTime = nan;
    
    % PREY ROUTINE SWIMMING PARAMETERS ----------------------------------------
    
    % Initial prey saccade time
    s.prey.tSccd = p.param.t_span(1);
    
    % Initial saccade direction (randomly selected)
    dir_PreySccd   = rand(1);
    s.prey.dirSccd = (dir_PreySccd-.5)./abs(dir_PreySccd-.5);
    
    % Whether currently executing a saccade
    s.prey.onSccd = 0;
    
    % PARAMETERS USEFUL FOR PREY ESCAPE ---------------------------------------
    
    % Logical that indicates whether prey is presently escaping
    s.prey.escapeOn = 0;
    
    % Initiate stimTime variable (start of clock for detecting predator)
    s.prey.stimTime = nan;
    
    % Set initial direction of escape
    s.prey.dirEsc = 1;

    % VISUAL SYSTEM PARAMETERS ------------------------------------------------
    
    % Angular location of target for predator (global FOR)
    s.pred.thetaTarget = nan;
    
    % Visual appearence of predator body
    [s.pred.xBodL, s.pred.yBodL] = give_coord('body coordinates',p.pred,p.pred.numBodPts);
    
    % Visual appearence of prey body
    [s.prey.xBodL, s.prey.yBodL] = give_coord('body coordinates',p.prey,p.prey.numBodPts);

    % Reset random number generator
    rng(1);
    
    
case 'Body positions' 
    
    % Check inputs
    if nargin < 5
        error('Need to provide 4 inputs for this action')
    end
    
    % Transform body points of predator into global FOR
    [s.pred.xBodG, s.pred.yBodG] = coord_trans('body to global', ...
        thetaPred, [xPred yPred], s.pred.xBodL, s.pred.yBodL);
    
    % Transform body points of prey into global FOR
    [s.prey.xBodG, s.prey.yBodG] = coord_trans('body to global', ...
        thetaPrey, [xPrey yPrey], s.prey.xBodL, s.prey.yBodL);
   
    % Store current orientation
    s.prey.theta = thetaPrey;
    s.pred.theta = thetaPred;
    
    
case 'Prey'
    
    % Check inputs
    if nargin < 5
        error('Need to provide 5 inputs for this action')
    end
    
    % Prey speed (fixed)      
    s.prey.spd = p.prey.spd0;
    
    % Max distance of body from origin
    body_dist_prey = max(hypot(xPrey, yPrey));
    
    % If wall is within sensory field  (trumps other behaviors) . . .
    if body_dist_prey > (p.param.tank_radius - p.prey.fieldSize)
        
        % Angular position of prey
        prey_ang = atan2(yPrey,xPrey);
        
        % Wall point
        xWallG = p.param.tank_radius * cos(prey_ang);
        yWallG = p.param.tank_radius * sin(prey_ang);
        
        % Local coordinates of wall
        [xWallL,yWallL] = coord_trans('global to body', thetaPrey, ...
            [xPrey yPrey], xWallG, yWallG);
        
        % Turn direction (away from wall point)
        turn_dir = - atan2(yWallL,xWallL) / norm(atan2(yWallL,xWallL));
        
        % Intensity of rotation
        turn_strength = (pi - abs(atan2(yWallL,xWallL)) + 0*pi*45/180) / pi;
        
        % Rate of turning, a function of angle to wall
        s.prey.omega = turn_dir * turn_strength * p.prey.wall_omega;
                
        % Reset saccade parameters
        s.prey.tSccd      = t;
        s.prey.onSccd     = 0;

    else
        % Escape response
        % If within an escape response AND within 4 x capture distance
        if ~s.prey.escapeOn && (dist < p.prey.thrshEscape)
            
            % Indicate that we are in an escape response
            s.prey.escapeOn = 1;
            
            % Note the time of its start
            s.prey.stimTime = t;
            
            % Set escape response direction
            dirEsc = rand(1);
            dirEsc = (dirEsc-.5)./abs(dirEsc-.5);
            
            % Determine speed of rotation of escape
            % normrnd(m,s) returns rndm # from norm dist. w/ mean=m,stdev=s
            rotSpdEscape = normrnd(p.prey.rotSpdEscape, 1.5);
            s.prey.omega = rotSpdEscape;
        end
        
        % If we are beyond the escape duration . . .
        if s.prey.escapeOn && ((t-s.prey.stimTime-p.prey.lat) > p.prey.durEscape)
            escapeOn = 0;
        end
        
        % If we are within the period of an escape . . .
        if s.prey.escapeOn
            %stimTime
            [s.prey.omega, s.prey.spd] = prey_escape(t, s.prey.stimTime,...
                p.prey, 0, 0, s.prey.dirEsc);
        else
            [s.prey.omega, s.prey.tSccd, s.prey.dirSccd, s.prey.onSccd] = ...
                foraging(p.prey, s.prey.tSccd, t, s.prey.dirSccd,...
                s.prey.onSccd);
        end
        
        % Note: if none of these conditions are met, then the prey
        % keeps moving along the same path
    end
    
      
case 'Predator'    
    
    % Check inputs
    if nargin < 5
        error('Need to provide 5 inputs for this action')
    end
    
    % Predator speed (fixed)      
    s.pred.spd = p.pred.spd0;
    
    % Determine whether prey is detected using 'see_fish'
    [s.pred.thetaTarget, s.pred.inFieldTime] = see_fish(t, [xPred yPred], ...
        thetaPred, s.prey.xBodG, s.prey.yBodG, p.pred, ...
        s.pred.inFieldTime, s.pred.thetaTarget);
    
    % Max distance of body from origin
    body_dist = max(hypot(xPred, yPred));
    
    % If wall is within sensory field  (trumps other behaviors) . . .
    if body_dist > (p.param.tank_radius - p.pred.fieldSize)
        
        % Angular position of predator
        pred_ang = atan2(yPred,xPred);
        
        % Wall point
        xWallG = p.param.tank_radius * cos(pred_ang);
        yWallG = p.param.tank_radius * sin(pred_ang);
        
        % Local coordinates of wall
        [xWallL,yWallL] = coord_trans('global to body', thetaPred, ...
            [xPred yPred], xWallG, yWallG);
        
        % Turn direction (away from wall point)
        turn_dir = - atan2(yWallL,xWallL) / norm(atan2(yWallL,xWallL));
        
        % Intensity of rotation
        turn_strength = (pi - abs(atan2(yWallL,xWallL)) + 0*pi*45/180) / pi;
        
        % Rate of turning, a function of angle to wall
        s.pred.omega = turn_dir * turn_strength * p.pred.wall_omega;
        
        % Disable targeting mode
        s.pred.thetaTarget = nan;
        s.pred.inFieldTime = nan;
        
        % Reset saccade parameters
        s.pred.tSccd      = t;
        s.pred.onSccd     = 0;
        
    % If prey visible (i.e. thetaTargetPred is not a nan) . . .
    elseif ~isnan(s.pred.thetaTarget)
        
        % Normalized deviation
        norm_dev = (s.pred.thetaTarget - thetaPred)/pi;
        
        % TARGETED SWIMMING: Adjust direction of predator, according
        % to position of prey
        s.pred.omega = norm_dev * p.pred.wall_omega;
        
        % Reset saccade parameters
        s.pred.tSccd      = t;
        s.pred.onSccd     = 0;
          
    % If prey not visible . . .
    else
        
        % Operate according to rules of foraging
        [s.pred.omega, s.pred.tSccd, s.pred.dirSccd, s.pred.onSccd] = foraging(...
            p.pred, s.pred.tSccd, t, s.pred.dirSccd, s.pred.onSccd);
        
    end

    
case 'Strike'    

% Default values
s.captured  = 0;
numSuc      = p.pred.numBodPts/2;
s.pred.xSuc = zeros(numSuc+2,1);
s.pred.ySuc = zeros(numSuc+2,1);
s.pred.rSuc = 0;

% If strike not started & within distance threshold . . .
if isnan(s.pred.strikeTime) &&  (dist <= p.pred.strike_thresh)
    % Set strike initiation
    s.pred.strikeTime = t;
end

% If strike started . . .
if ~isnan(s.pred.strikeTime)
    
    % If not finished . . .
    if t < (s.pred.strikeTime + p.pred.strike_dur)
        
        % Parameters for the each function
        t_c = t - s.pred.strikeTime;
        dur = p.pred.strike_dur;
        
        % Current reach
        s.pred.rSuc = p.pred.strike_reach .* (0.5.*(sin(2*pi*t_c./dur-pi/2)+1));
        
        % Plot local suction zone (local FOR)
        phi = linspace(-p.pred.strike_range/2, p.pred.strike_range/2, numSuc)';
        xSucL = [0; s.pred.rSuc.*cos(phi); 0];
        ySucL = [0; s.pred.rSuc.*sin(phi); 0];
        
        % Transform into global FOR
        [s.pred.xSuc, s.pred.ySuc] = coord_trans('body to global', ...
                                        thetaPred, [xPred yPred], ...
                                        xSucL, ySucL);
        
        % Coordinates of prey in pred FOR
        [xPreyPred,yPreyPred] = coord_trans('global to body', ...
                                    thetaPred, [xPred yPred], ...
                                    s.prey.xBodG, s.prey.yBodG);
        
        % Polar coordinates of prey body in pred FOR
        [phiPreyPred, rPreyPred] = cart2pol(xPreyPred, yPreyPred);
        
        % Index of points in capture zone
        idx = (rPreyPred < s.pred.rSuc) & ...
            (phiPreyPred >= -p.pred.strike_range/2) & ...
            (phiPreyPred <= p.pred.strike_range/2);
        
        % Captured, if more than half body in capture zone
        if sum(idx) > (length(idx)/2)
            s.captured = 1;
        end
        
    % If beyond duration . . .
    else
        % Turn off strike
        s.pred.strikeTime = nan;
    end
end

end








function [x,y] = give_coord(kind,p,num_pts)
% Returns coordinates requested, based on given parameter values


if strcmp(kind,'body coordinates')
% Coordinates describing the visual appearence of fish
    
    % Number of points to render each part of the body
    num1 = round(num_pts.*(p.COM_len/p.bod_len));
    num2 = num_pts - num1;
    
    % Angular position
    ang1        = linspace(-pi/2,pi/2,num1)';
    ang2        = linspace(pi/2,3*pi/2,num2)';
    
    % Trunk length
    trunkL      = p.bod_len - p.COM_len;
    
    % Coordinates
    x   = [p.COM_len.*cos(ang1); trunkL.*cos(ang2)]-p.COM_len;
    y   = [(p.bod_width/2).*sin(ang1); (p.bod_width/2).*sin(ang2)];
    
    
else
    error('"kind" input not recognized');
end