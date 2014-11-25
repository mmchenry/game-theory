%----PREY_ESCAPE() returns omegaPrey and spdPrey based on change function
%below. The inputs direction and speed of rotation are chosen randomly at
%the onset of escape maneuver.  

function [omega,spd] = prey_escape(t, stimTime, prey, omega0, spd0, dirEsc) 

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
    omega = dirEsc .* omega;

end



function y = change_func(t, amp, dur)
% Curve describing the turning rate

y = amp.*(0.5.*(sin(2*pi*t./dur - pi/2)+1));
