function R = simple_model(pred,prey,param)
% ODE Predator-prey interaction model (no longer so "simple")


%% Parameters

% Scaling constants (helps the solver avoid crazy numbers)
sL = pred.bod_len;              % m
sT = pred.saccade_interval / 4; % s

% Number of points used to render the periphery of fish bodies
num_bod_pts = 200;


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

% Kinematic parameters
p.pred.spd0                 = pred.spd0             ./sL .* sT;
p.pred.sccd_prd             = pred.saccade_period   ./ sT;
p.pred.sccd_intvl           = pred.saccade_interval ./ sT;
p.pred.sccd_omega           = pred.saccade_omega    .* sT;
p.pred.wall_omega           = pred.wall_omega       .* sT;

% Sensory parameters
p.pred.fieldSize    = pred.fieldSize ./sL ;
p.pred.vis_az       = pred.vis_az;
p.pred.verg_ang     = pred.verg_ang;
p.pred.rtnl_den     = pred.rtnl_den;
p.pred.vis_freq     = pred.vis_freq .* sT;
p.pred.vis_thresh   = pred.vis_thresh;


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

% PREDATOR FORAGING PARAMETERS --------------------------------------------

% Initial predator saccade time
t_PredSccd = p.t_span(1);

% Initial saccade direction (randomly selected)
dir_PredSccd = rand(1);
dir_PredSccd = (dir_PredSccd-.5)./abs(dir_PredSccd-.5);

% Whether currently executing a saccade
on_PredSccd = 0;

% Time when prey enters predator's FOV
inFieldTimePred = nan;

% PREY ROUTINE SWIMMING PARAMETERS ----------------------------------------

% Initial predator saccade time
t_PreySccd = p.t_span(1);

% Initial saccade direction (randomly selected)
dir_PreySccd = rand(1);
dir_PreySccd = (dir_PreySccd-.5)./abs(dir_PreySccd-.5);

% Whether currently executing a saccade
on_PreySccd = 0;

% POINTS FOR TANK, PREY & PREDATOR HEAD -----------------------------------

% % Points describing boundary of tank (global FOR)
% theta = [linspace(0,2*pi,1000)]';
% xTank = p.tank_rad .* cos(theta);
% yTank = p.tank_rad .* sin(theta);
% clear theta
% Angular location of target for predator (global FOR) 
thetaTargetPred = nan;

% Visual appearence of predator body
[xBodPredL, yBodPredL] = give_coord('body coordinates',p.pred,num_bod_pts);

% Visual appearence of prey body
[xBodPreyL, yBodPreyL] = give_coord('body coordinates',p.prey,num_bod_pts);

% % Points describing prey head (local FOR)
% angHead    = [linspace(-pi/2, pi/2, 500)]';
% xCollPreyL = p.prey.fldSize*p.prey.L    .* cos(angHead) - p.prey.L;
% yCollPreyL = p.prey.fldSize*(p.prey.w/2) .* sin(angHead);

% PARAMETERS USEFUL FOR PREY ESCAPE ---------------------------------------

% Logical that indicates whether prey is presently escaping 
escapeOn = 0;

% Initiate stimTime variable
stimTime = nan;

% Set initial direction of escape 
dirEsc = 1;


% Visual check
%plot(xBodPredL,yBodPredL,'.-',xFieldPredL,yFieldPredL,'.-');axis equal
%plot(xBodPreyL,yBodPreyL,'.-');axis equal

% Runge-Kutta solver
[t,X,t_event] = ode45(@gov_eqn,p.t_span,X0,options);

    function dX = gov_eqn(t,X)
    % ODE for the dynamics of the system
        
        %  INPUTS ---------------------------------------------------------
        
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

        % DECISIONS ABOUT RATES OF CHANGE ---------------------------------
                         
                            
%        % Transform region outside of head points of predator into global FOR                    
%        [xCollPred, yCollPred] = coord_trans('body to global', ...
%                      thetaPred, [xPred yPred], xCollPredL, yCollPredL); 
%        
%        % Transform region outside of head points of prey into global FOR                    
%        [xCollPrey, yCollPrey] = coord_trans('body to global', ...
%                      thetaPrey, [xPrey yPrey], xCollPreyL, yCollPreyL);                  

        % Transform field points of predator into global FOR                    
        %[xFieldPredG, yFieldPredG] = coord_trans('body to global', ...
        %               thetaPred, [xPred yPred], xFieldPredL, yFieldPredL); 
                  
        % Transform body points of predator into global FOR             
        [xBodPredG, yBodPredG] = coord_trans('body to global', ...
                       thetaPred, [xPred yPred], xBodPredL, yBodPredL); 
        
        % Transform body points of prey into global FOR             
        [xBodPreyG, yBodPreyG] = coord_trans('body to global', ...
                       thetaPrey, [xPrey yPrey], xBodPreyL, yBodPreyL); 
        
        % Predator speed (fixed)
        spdPred = p.pred.spd0;
        
       % Prey speed         
%         if dist < p.param.dist_thresh
%             spdPrey = p.prey.spd0;
%         else
            spdPrey = p.prey.spd0;
%         end
        
        % PREY BEHAVIOR ---------------------------------------------------
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
            omegaPrey = turn_dir * turn_strength * p.prey.wall_omega;
            
            % Disable targeting mode
            thetaTargetPred = nan;
            inFieldTimePred = nan;
            
            % Reset saccade parameters
            t_PreySccd      = t;
            on_PreySccd     = 0;
            
            clear xWallG yWallG xWallL yWallL turn_dir turn_strength 
     
        else
            % Escape response
            % If in an escape response AND within 4 x capture distance
            if ~escapeOn && (dist < 4*p.param.d_capture)
                
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
                [omegaPrey, t_PreySccd, dir_PreySccd, on_PreySccd] = ...
                    foraging(p.prey,t_PreySccd, t, dir_PreySccd,...
                    on_PreySccd);
            end
        end
        

        % PREDATOR BEHAVIOR -----------------------------------------------
        
        % Determine whether prey is detected using 'see_fish'
        [thetaTargetPred, inFieldTimePred] = see_fish(t, [xPred yPred], ...
                                thetaPred, xBodPreyG, yBodPreyG, p.pred, ...
                                inFieldTimePred, thetaTargetPred);
        
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
            omegaPred = turn_dir * turn_strength * p.pred.wall_omega;
            
            % Disable targeting mode
            thetaTargetPred = nan;
            inFieldTimePred = nan;
            
            % Reset saccade parameters
            t_PredSccd      = t;
            on_PredSccd     = 0;
            
            clear xWallG yWallG xWallL yWallL turn_dir turn_strength 
                            
            
        % If prey visible (i.e. thetaTargetPred is not a nan) . . .
        elseif ~isnan(thetaTargetPred)
            
            % Normalized deviation
            norm_dev = (thetaTargetPred-thetaPred)/pi;
            
            % TARGETED SWIMMING: Adjust direction of predator, according
            % to diretcion of prey
            omegaPred = norm_dev * p.pred.wall_omega;
            
            % TODO: Add strike code 
            
            % Reset saccade parameters
            t_PredSccd      = t;
            on_PredSccd     = 0;
            
            % Clear parameters
            clear norm_dev
          
            
        % If prey not visible . . .
        else
            
            % Operate according to rules of foraging
            [omegaPred, t_PredSccd, dir_PredSccd, on_PredSccd] = foraging(...
                     p.pred, t_PredSccd, t, dir_PredSccd, on_PredSccd);

        end
       
    
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


%% CAPTURE FUNCTION: use with 'Event' option in ode45 ----------------------
function [value, isterminal, direction] = capture_fnc(~,X)
    % the event occurs when distance is less than capture distance
    distance = hypot(X(4)-X(1),X(5)-X(2));
    if distance < p.param.d_capture
        value      = 0;
        isterminal = 1;         % tells ode45 to stop integration
        direction  = 0;
    else
        value = 1;
        isterminal = 0;         
        direction  = 0;
   end
end



% function [x,y] = coord_pred_field(p,num_pts)
% % Coordinates describing the sensory field of the predator
% 
% % Number of points to render each part of the body
% num1 = round(num_pts.*(p.COM_len/p.bod_len));
% 
% % Angular position
% ang1        = linspace(-pi/2,pi/2,num1)';
% 
% % Coordinates
% x   = p.fldSize.*([p.COM_len.*cos(ang1)])-p.COM_len;
% y   = p.fldSize.*([(p.bod_width/2).*sin(ang1)]);
% 
% end

end