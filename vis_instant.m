function h = vis_instant(xPred, yPred, thetaPred, xPrey, yPrey, thetaPrey, p)
% Visualizes the state of a simulation at an instant of time

% Tank radius
r = p.param.tank_radius;

% Get color palette
clrs = get(gca,'ColorOrder');

num_bod_pts = 100;

% Get limits
xlims(1) = -r * 1.1;
xlims(2) =  r * 1.1;
ylims(1) = -r * 1.1;
ylims(2) =  r * 1.1;

phi = linspace(0,2*pi,500);

%% Resolve coordinates

% Visual appearence of predator body
[xBodPredL, yBodPredL] = give_coord('body coordinates',p.pred,num_bod_pts);

% Visual appearence of prey body
[xBodPreyL, yBodPreyL] = give_coord('body coordinates',p.prey,num_bod_pts);

% Transform body points of predator into global FOR
[xBodPredG, yBodPredG] = coord_trans('body to global', ...
    thetaPred, [xPred yPred], xBodPredL, yBodPredL);

% Transform body points of prey into global FOR
[xBodPreyG, yBodPreyG] = coord_trans('body to global', ...
    thetaPrey, [xPrey yPrey], xBodPreyL, yBodPreyL);

%% Plot

% Plot walls
plot(r.*cos(phi),r.*sin(phi),'-k');
hold on

% Render predator
h(1) = fill(xBodPredG,yBodPredG,clrs(2,:));
set(h(1),'EdgeColor','none')

% Render prey
h(2) = fill(xBodPreyG,yBodPreyG,clrs(1,:));
set(h(2),'EdgeColor','none')


%% Adjustments

% Set axis
axis square
xlim(xlims)
ylim(ylims)

% Removed axes
set(gca,'Box','off')
set(gca,'XColor','w')
set(gca,'YColor','w')

hold off