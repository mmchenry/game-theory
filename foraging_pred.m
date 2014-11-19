function [omega, t_sccd, dirSccd, onSccd] = ...
               foraging_pred(xBod,yBod,theta, pred, t_sccd, t, dirSccd, ...
                             onSccd, tank_rad, player)
% Describes instantaneous routine foraging behavior 
%
% x       - x-position of body (nx1)
% y       - y-position of body (nx1)
% theta   - Body orientation of predator (1x1)
% pred    - Structure of predator (or prey) parameters
% t_sccd  - Time of start of last saccade (i.e. turn) (1x1)
% t       - Current time (1x1, can be vector for testing)
% dirSccd - Saccade direction (1 or -1)
% onSccd  - Currently executing a saccade (logical)
% xTank   - x coordinates of tank wall (nx1)
% yTank   - y coordinates of tank wall (nx1)
% player  - Either 'predator' or 'prey'
% 
% omega  - rate of body rotation

gain = 200;

% Max distance of body from origin
body_dist = max(hypot(xBod, yBod));

% Tank radius
%tank_rad = max(hypot(xTank, yTank));

% If bumped into to a wall . . .
if body_dist >= tank_rad
    
    % Index of points that exceed the boundary
    idx = hypot(xBod,yBod) > tank_rad;
    
    % Body points outside
    xOut    = xBod(idx);
    yOut    = yBod(idx);
    
    % Find distances for all head points to all tank points
    dist = hypot(xOut,yOut);
    
    % Find furthest distance
    iMax = find(dist==max(dist(:)),1,'first');
    
    % Angle of that furthest point
    ang_max = atan2(yOut(iMax),xOut(iMax));
    
    % Tangent angles
    tan_angs = [ang_max+pi/2 ang_max-pi/2];
    
    % Deviation angles
    dev_angs = abs(theta - tan_angs);
    
    % Normalized deviation of theta from normal angle
    dev_norm = (tan_angs(dev_angs==min(dev_angs)) - theta)/pi;
    
    % Make omega proportional to deviation
    omega = dev_norm * pred.sccd_omega * gain;
    
    % Vectorize coordinates for distance calculation
    %xTank = repmat(xTank',length(xBod),1);
    %yTank = repmat(yTank',length(yBod),1);
    %xBod  = repmat(xBod, 1, size(xTank,2));
    %yBod  = repmat(yBod, 1, size(xTank,2));
    
    % Find distances for all head points to all tank points
    %distBump = hypot(xTank-xBod,yTank-yBod);
    
    % Find closest distance
    %[r,c] = find(distBump==min(distBump(:)),1,'first');
    
    % Choose direction away from side closest to wall
    % (first half of head points assumed to be on right side)
%     if iMax < length(xBod)/2
%         dirSccd = 1;
%     else
%         dirSccd = -1;
%     end
    
    
    % Current angular position of collision surface
    %spot_ang = atan2(yBod(r,c),xBod(r,c));
    
    %coll_dist = hypot(xBod(r,c),yBod(r,c)
    
    % Current deviation from a normal angle
    %dev_ang = theta - tan_ang
    
    %omega =  pred.sccd_omega .*dev_ang .* gain;
    %omega = dirSccd*(max(dist(:))-tank_rad) .* gain;

    
    % Set saccade start to current time
    t_sccd = t;
    
    % Turn on off routine saccades
    onSccd = 0;
    
    
% If currently executing a saccade (and not bumping a wall) . . .
elseif onSccd
      
    % Calculated time within the saccade
    tc = t - t_sccd;
    
    % If this time has exceeded saccade period . . .
    if tc > pred.sccd_prd
        
        % Turn off logical
        onSccd = 0;
        
        % Zero out omega
        omega = 0;  
        
        % Determine next saccade direction 
        dirSccd = rand(1);dirSccd = (dirSccd-.5)./abs(dirSccd-.5);
        
    % Otherwise, continue with saccade function . . .
    else    
        
        % Define current rate of turn
        omega = saccade_func(tc,pred.sccd_omega,pred.sccd_prd,dirSccd);   
        
    end
       
% If not already executing a saccade (& not near a wall). . .
else
    
    % Next saccade time (useful for some code below)
    t_next = t_sccd + pred.sccd_intvl;

    % If current time exceeds the next scheduled saccade . . .
    if t >= t_next
        
        % Set new saccade start time
        t_sccd = t_next;
        
        % Turn on logical 
        onSccd = 1;
        
        % Define current rate of turn
        omega = saccade_func(t-t_next,pred.sccd_omega,pred.sccd_prd,dirSccd); 
    
    % Otherwise, we are between saccades . . .
    else            
        % . . . so zero rate of angular change
        omega = 0;         
    end
end




    
function omega = saccade_func(t, amp, dur, dirSccd)
% Curve describing the turning rate

omega = amp.*(0.5.*(sin(2*pi*t./dur-pi/2)+1));

% Set direction
omega = dirSccd .* omega;

