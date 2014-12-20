function R = reconstruct(R, p, d)
% Reconstructs events in simulation

% Store input parameters
R.p     = p;
R.d     = d;

% Intitialize state variables
s = give_behavior('Initialize', R.p);

% Step thru time
for i = 1:length(R.t)

    % Find body position of both fish
    s = give_behavior('Body positions', R.p, s, R.t(i), R.X(i,:));
    
    % Set prey behavior
    s = give_behavior('Prey', R.p, s, R.t(i), R.X(i,:));
    
    % Set predator behavior
    s = give_behavior('Predator', R.p, s, R.t(i), R.X(i,:));
    
    % Determine whether prey is captured
    s = give_behavior('Strike', R.p, s, R.t(i), R.X(i,:));   
    
    % Store coordinates for suction
    R.pred.xSuc(:,i)    = s.pred.xSuc;
    R.pred.ySuc(:,i)    = s.pred.ySuc;
    
    % Store rates of change
    R.pred.omega(i,1)   = s.pred.omega;
    R.prey.omega(i,1)   = s.prey.omega;
    R.prey.spd(i,1)     = s.prey.spd;
    
    % Predator behavioral states
    R.pred.strike(i,1)  = ~isnan(s.pred.strikeTime);
    R.pred.seePrey(i,1) = ~isnan(s.pred.inFieldTime);
    R.pred.onSccd(i,1)  = s.pred.onSccd;
    R.pred.thetaTarget(i,1) = s.pred.thetaTarget;
    
    % Prey behavioral states
    R.prey.escape(i,1)         = s.prey.escapeOn;
    R.prey.onSccd(i,1)         = s.prey.onSccd;
    
    % TODO: Store other behaviors (escape, walls, saccades) & sensory
    % events
end