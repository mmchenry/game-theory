%% Function to compute optimal escpae angle at (x0,y0,K)
function [theta,fval] =  optimTheta(x0,y0,k)

% Use 'fminbnd' by calling fminbnd(@fun,a,b)
[theta,fval] = fminbnd(@mindist,0,acos(1/k));

% Nested function that computes the objective function     
    function d = mindist(theta)
        
        % minimum distance function
        d = -1 * ((k*y0 - y0*cos(theta) + x0*sin(theta))^2 / ...
            (1+k^2-2*k*cos(theta)));
        
    end
end