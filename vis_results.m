function vis_results(R,graph_type)
% General function for plotting the results of an individual simulation
%
% R             - Structure that stores results from simulations
% graph_type    - String that designates the type of graph

% Get color palette
clrs = get(gca,'ColorOrder');


switch graph_type
    
    % Plot of body orientation and rate of change in orientaion
    case 'Turning data'
        
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
        
    % Plot trajectories    
    case 'Trajectories'
        
        num_pts = 50;
        
        phi = linspace(0,2*pi,500);
        
        r = R.param.tank_radius;
        
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
end











