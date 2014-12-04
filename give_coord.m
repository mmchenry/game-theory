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