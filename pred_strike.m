function [strikeTime, captured] = pred_strike(t,strikeTime, origin, ...
                                    thetaPred, p, dist, xBodPreyG, yBodPreyG)
% Determine state of strike
     
% Default value
captured = 0;

% If strike not started & within distance threshold . . .
if isnan(strikeTime) &&  (dist <= p.strike_thresh)
    % Set strike initiation
    strikeTime = t;
end

% If strike started . . .
if ~isnan(strikeTime)
    
    % If not finished . . .
    if t < (strikeTime + p.strike_dur)
        
        % Parameters for the each function
        t_c = t - strikeTime;
        dur = p.strike_dur;
        
        % Current reach
        r = p.strike_range .* (0.5.*(sin(2*pi*t_c./dur-pi/2)+1));
        
        % Coordinates of prey in pred FOR
        [xPreyPred,yPreyPred] = coord_trans('global to body', ...
                                    thetaPred, origin, xBodPreyG, yBodPreyG);
        
        % Polar coordinates of prey body in pred FOR
        [phiPreyPred, rPreyPred] = cart2pol(xPreyPred, yPreyPred);
        
        % Index of points in capture zone
        idx = (rPreyPred < r) & ...
            (phiPreyPred >= -p.strike_range/2) & ...
            (phiPreyPred <= p.strike_range/2);
        
        % Captured, if more than half body in capture zone
        if sum(idx) > (length(idx)/2)
            captured = 1;
        end
        
    % If beyond duration . . .
    else
        % Turn off strike
        strikeTime = nan;
    end
end