%----PREY_ESCAPE() takes current direction of prey and outputs a new
% direction chosen randomly from a uniform distribution on the interval
% (pi/3, 2*pi/3). The direction is toward the left (L) or the right (R)
% with a bias toward the right, w.r.t. the midline of the prey

function [omega,spd] = prey_escape(t, stimTime, prey, omega0, spd0) 
% decide direction, right bias
%choice = round(1.1*rand);       
%if choice > 0
    %thetaPrey = unifrnd(pi/3,2*pi/3) - thetaPrey;
%else
    %thetaPrey = thetaPrey + unifrnd(pi/3,2*pi/3);
%end

% I'll leave it to you to enter a probabilatistic version

% If we are still within the latency period . . .
if (t-stimTime) < prey.lat 
    
    % Keep the rates the same
    omega = omega0;
    spd = spd0;

% Otherwise . . .
else
    
    % Time since start of escape
    escapeTime = t-stimTime-prey.lat;
    
    % Set omega and speed according to current time and change function
    omega = change_func(escapeTime, prey.rotSpdEscape, prey.durEscape);
    spd = change_func(escapeTime, prey.spdEscape, prey.durEscape);

end



function y = change_func(t, amp, dur)
% Curve describing the turning rate

y = amp.*(0.5.*(sin(2*pi*t./dur - pi/2)+1));
