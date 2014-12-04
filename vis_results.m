function vis_results(graph_type,R)
% General function for plotting the results of an individual simulation
%
% R             - Structure that stores results from simulations
% graph_type    - String that designates the type of graph

% Get color palette
clrs = get(gca,'ColorOrder');


switch graph_type

case 'Turning data'
% Plot of body orientation and rate of change in orientaion

    % Plot orientaion angles
    subplot(2,1,1)
    h = plot(R.t,unwrap(R.thetaPred).*180/pi,'+-',...
        R.t,unwrap(R.thetaPrey).*180/pi,'+-');
    grid on
    set(h(2),'Color',clrs(1,:))
    set(h(1),'Color',clrs(2,:));
    xlabel('Time (s)')
    ylabel('Orientation (deg)')
    legend('Predator','Prey','Location','NorthWest')
    
    % Plot rates of change in orientation
    subplot(2,1,2)
    h = plot(R.t(2:end),unwrap(diff(R.thetaPred)).*180/pi,'-',...
        R.t(2:end),unwrap(diff(R.thetaPrey)).*180/pi,'-');
    grid on
    set(h(2),'Color',clrs(1,:))
    set(h(1),'Color',clrs(2,:));
    xlabel('Time (s)')
    ylabel('Rate of orientation change (deg/s)')
    
    
case 'Trajectories'
% Plot trajectories 

    % Number of points for the wall 
    num_pts = 200;
    
    % Radial coordinates for the wall
    phi = linspace(0,2*pi,500);

    % Extract wall radius 
    r = R.p.param.tank_radius;

    % Get limits
    xlims(1) = -r * 1.1;
    xlims(2) =  r * 1.1;
    ylims(1) = -r * 1.1;
    ylims(2) =  r * 1.1;

    % Time vector for rendering values
    t_render = linspace(0,R.t(end),num_pts);

    % Interpolate results for even timing
    xPred = interp1(R.t,R.xPred,t_render);
    yPred = interp1(R.t,R.yPred,t_render);
    xPrey = interp1(R.t,R.xPrey,t_render);
    yPrey = interp1(R.t,R.yPrey,t_render);

    % Render traces
    h = plot(R.xPred,R.yPred,'-',...
             xPred,yPred,'o',...           
             R.xPrey,R.yPrey,'-',...
             xPrey,yPrey,'o',...
             R.xPred(end),R.yPred(end),'o',...
             R.xPrey(end),R.yPrey(end),'o');

    hold on
    h2 = plot(r.*cos(phi),r.*sin(phi),'-k');
    hold off
    
    set(h(2),'MarkerSize',2)
    set(h(4),'MarkerSize',2)
    set(h(5),'MarkerSize',4)
    set(h(6),'MarkerSize',4)

    % Set axis
    axis square
    xlim(xlims)
    ylim(ylims)

    % Set colors
    set(h(1:2),'Color',clrs(2,:))
    set(h(3:4),'Color',clrs(1,:))
    set(h(2),'MarkerFaceColor',clrs(2,:))
    set(h(4),'MarkerFaceColor',clrs(1,:))
    set(h(5),'MarkerFaceColor',clrs(2,:))
    set(h(6),'MarkerFaceColor',clrs(1,:))
    set(h(5),'MarkerEdgeColor',clrs(2,:))
    set(h(6),'MarkerEdgeColor',clrs(1,:))

    % Removed axes
    set(gca,'Box','off')
    set(gca,'XColor','w')
    set(gca,'YColor','w')

    
case 'Animate'     
% Animate simulation

    % Duration for playing video (s)
    play_dur = 10;
    
    % Frame rate (Hz) 
    fps = 15;
    
     % Number of points for the wall 
    num_pts = 200;
    
    % Radial coordinates for the wall
    phi = linspace(0,2*pi,500);

    % Extract wall radius 
    r = R.param.tank_radius;
    
    % Get limits
    xlims(1) = -r * 1.1;
    xlims(2) =  r * 1.1;
    ylims(1) = -r * 1.1;
    ylims(2) =  r * 1.1;
    
    % Time vector for rendering values
    t_render = linspace(0,play_dur,fps*play_dur);
    
    % Normalize time of simulation results wrt animation duration
    t_norm = R.t/max(R.t)*play_dur;
      
    % Interpolate predator results for even timing
    xPred       = interp1(t_norm,R.xPred,t_render);
    yPred       = interp1(t_norm,R.yPred,t_render);
    thetaPred   = interp1(t_norm,R.thetaPred,t_render);
    
    % Interpolate prey results for even timing
    xPrey       = interp1(t_norm,R.xPrey,t_render);
    yPrey       = interp1(t_norm,R.yPrey,t_render);
    thetaPrey   = interp1(t_norm,R.thetaPrey,t_render);
    
    % Create figure window
    %figure('DoubleBuffer','on')
    
    % Start timer
    clock_start = tic;
    
    % Look through time
    for i = 1:length(t_render)
          
        % Draw walls
        h2 = plot(r.*cos(phi),r.*sin(phi),'-k');
        hold on
        
        % Render traces
        h = plot(xPred(1:i),yPred(1:i),'-',xPrey(1:i),yPrey(1:i),'-');
        
        % Render bodies of predtaor and prey
        hBod = vis_instant(xPred(i), yPred(i), thetaPred(i), xPrey(i), ...
                  yPrey(i), thetaPrey(i), R);
        
        % Set axis
        axis square
        xlim(xlims)
        ylim(ylims)
        
        % Get color palatte
        clrs = get(gca,'ColorOrder');
        
        % Set colors
        set(h(1),'Color',clrs(2,:))
        set(h(2),'Color',clrs(1,:))

        % Removed axes
        set(gca,'Box','off')
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        
        % Get elapsed time
        telaspsed = toc(clock_start);
        
        % Pause to adjust pacing
        if t_render(i)>telaspsed
            pause(t_render(i)-telaspsed);
        else
            pause(1e-6)
        end
        
        hold off
    end
    
    telaspsed = toc(clock_start);
    
    % Check duration
    if telaspsed > 1.2*play_dur
        warning(['Could not render animation fast enough to play at' ...
            ' requested duration']);
    end

    

end









