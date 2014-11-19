%----PREY_ESCAPE() takes current direction of prey and outputs a new
% direction chosen randomly from a uniform distribution on the interval
% (pi/3, 2*pi/3). The direction is toward the left (L) or the right (R)
% with a bias toward the right, w.r.t. the midline of the prey

function [thetaPrey] = prey_escape(thetaPrey) 
% decide direction, right bias
choice = round(1.1*rand);       
if choice > 0
    thetaPrey = unifrnd(pi/3,2*pi/3) - thetaPrey;
else
    thetaPrey = thetaPrey + unifrnd(pi/3,2*pi/3);
end
end