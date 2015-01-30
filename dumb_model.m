function dumb_model
% Find Dmin topography in the dumbest way possible

%% Code execution

% Dmin plot, as a fucntion of K and alpha
do_Dmin_topo = 0;

% Dmin plot, as a function of y0 and alpha (fixed K)
do_Dmin_topo_y0 = 1;

find_edge = 1;

%% Dmin topography (do_Dmin_topo)

if do_Dmin_topo
    
    % Paramaters
    alpha_max = 180;
    alpha_min = 0;
    t_max = 5;
    Upred = 1;
    x0 = 5;
    y0 = 0;
    
    
    % Number of alpha and K values inlcuded in plot
    num_alpha   = 65;
    num_K       = 65;
    num_t       = 10^4;
    
    % Creat vectors of alpha, K, and time
    alpha_vect = linspace(alpha_min,alpha_max,num_alpha)'.*pi/180;
    K_vect     = 10.^linspace(-2,2,num_K)';
    t = linspace(0,t_max,num_t)';
    
    Dmin0 = hypot(x0,y0);
    
    % Create matricies of alpha and K values
    [K,alfa] = meshgrid(K_vect,alpha_vect);
    
    % Create higher resolution version of K for the optimal curve
    K_lin = 10.^linspace(-2,2,num_K.*300)';

    % Solve high resolution alpha values and corresponding Dmins
    alpha_lin = [acos(K_lin(K_lin<=1)); acos(1./K_lin(K_lin>1))];
    Dmin_lin = eqn35(K_lin,alpha_lin);
    
    % Loop thru K and alpha values
    for i = 1:size(K,1)
        for j = 1:size(K,2)
            Uprey = Upred ./ K(i,j);
            Dmin(i,j)       = dumb_sim(t,alfa(i,j),Upred,Uprey,x0,y0)./Dmin0;
            %Dmin_f(i,j)     = dumb_follower(t,alfa(i,j),Upred,Uprey,x0,y0)./Dmin0;
        end
    end
    
    
    figure
    %subplot(2,1,1)
    
    h = mesh(K,alfa.*180/pi,Dmin);
    hold on
    %view([40 19])
    %view([67 22])
    %view([148 31])
    view([71 25])
    
    hL = plot3(K_lin,alpha_lin.*180/pi,Dmin_lin,'k-');
    hold off
    
    set(gca,'XScale','log')
    
    xlabel('K');
    ylabel('alpha (deg)')
    zlabel('D_m_i_n*')
    %set(hL,'LineWidth',2)
    
    caxis([0 1.2])
    set(gca,'TickLength',[0.03 0.03])
    set(gca,'YTick',0:30:alpha_max)
    set(gca,'ZTick',0:0.25:1)
    zlim([0 1])
    ylim([alpha_min alpha_max])
    colorbar('EastOutside')
    
end


%% Dmin topography wrt y0 (do_Dmin_topo_y0)

if do_Dmin_topo_y0
    
    % Paramaters
    alpha_max = 180;
    alpha_min = -180;
    
    theta = 80.6/180*pi;
    
    hyp = 20;
    y0_max = hyp*sin(theta);
    y0_min = -hyp*sin(theta);
    t_max = 5;
    Upred = 1;
    x0 = hyp*cos(theta);
    K = 0.9;
    
    step_size = 0.01;
    
    % Number of alpha and K values inlcuded in plot
    num_alpha   = 200;
    num_y0      = 200;
    num_t       = 10^4;
    
    % Creat vectors of alpha, K, and time
    alpha_vect = linspace(alpha_min,alpha_max,num_alpha)'.*pi/180;
    y0_vect    = linspace(y0_min,y0_max,num_y0)';
    t = linspace(0,t_max,num_t)';
    
    % Create matricies of alpha and K values
    [y0,alfa] = meshgrid(y0_vect,alpha_vect);

    % Loop thru K and alpha values
    for i = 1:size(y0,2)
        
        if find_edge
            % State of finding critical edge
            crit_state = 0;
        end
        
        for j = 1:size(y0,1)
            
            Dmin0 = hypot(x0,y0(j,i));
            
            Uprey = Upred ./ K;
            Dmin(j,i) = dumb_sim(t,alfa(j,i),Upred,Uprey,x0,y0(j,i))./Dmin0;
            
            if find_edge
                % It we just hit the plateau . . .
                if (crit_state==0) && (Dmin(j,i)==1)
                    
                    % Find critical alpha value
                    cAlfa = find_edge(t, x0, y0(j,i), Upred, Uprey, Dmin0, ...
                        alfa(j-1,i), alfa(j,i),'approaching',step_size);
                    
                    % Adjust critical state
                    crit_state = 1;
                    
                    % Store critical alpha value
                    crit(1:3,i) = [y0(j,i) cAlfa 1]';
                    
                    % If we are just leaving the plateau . . .
                elseif (crit_state==1) && (Dmin(j,i)~=1)
                    
                    % Find critical alpha value
                    cAlfa = find_edge(t, x0, y0(j,i), Upred, Uprey, Dmin0, ...
                        alfa(j-1,i), alfa(j,i),'leaving',step_size);
                    
                    % Adjust critical state
                    crit_state = nan;
                    
                    % Store critical alpha value
                    crit(4,i) = cAlfa;
                end
            end
        end
    end

    
    figure
    if find_edge
        subplot(2,1,1)
    end
    
    h = surf(y0,alfa.*180/pi,Dmin,'FaceColor','interp','FaceLighting','gouraud');
    set(h,'EdgeColor','none')
    
    %h = mesh(y0,alfa.*180/pi,Dmin);
    hold on
    %view([40 19])
    view([-49 20])
    %view([90 0])
    %view([0 90])
    
    if find_edge
        % Edges of plateau
        hL = plot3(crit(1,:), crit(2,:).*180/pi, crit(3,:),'k-',...
            crit(1,:), crit(4,:).*180/pi, crit(3,:),'k-');
        set(hL,'LineWidth',1.5)
    end
    
%     hL = plot3(crit(1,:), crit(2,:).*180/pi, crit(3,:).*0,'k-',...
%                crit(1,:), crit(4,:).*180/pi, crit(3,:).*0,'k-');        
%     set(hL,'LineWidth',2)
    
    %hL = plot3(K_lin,alpha_lin.*180/pi,Dmin_lin,'k-');
    hold off
    
    xlabel('y0');
    ylabel('alpha (deg)')
    zlabel('D_m_i_n*')
    title(['K = ' num2str(K)])
    %set(hL,'LineWidth',2)
    
    caxis([0 1.1])
    set(gca,'TickLength',[0.03 0.03])
    set(gca,'YTick',alpha_min:90:alpha_max)
    set(gca,'XTick',[y0_min 0 y0_max])
    set(gca,'ZTick',0:0.25:1)
    zlim([0 1])
    ylim([alpha_min alpha_max])
    xlim([y0_min y0_max])
    colorbar('EastOutside')
    axis square
    
    %figure
    if find_edge
        subplot(2,1,2)
        hL = plot(crit(1,:), crit(2,:).*180/pi,'k-',...
            crit(1,:), crit(4,:).*180/pi,'k-');
        set(gca,'YTick',alpha_min:45:alpha_max)
        ylim([alpha_min alpha_max])
        set(gca,'TickDir','out')
        xlabel('y0*');
        ylabel('alpha (deg)')
        grid on
        axis square
    end
end





function cAlfa = find_edge(t, x0, y0, Upred, Uprey, Dmin0, alfa1, ...
                           alfa2, status, step_size)

% Starting cut in the change in alpha
cutt = step_size;

% Check whether approaching or leaving the plateau
if strcmp(status, 'approaching')
    upwards = 1;
elseif strcmp(status, 'leaving')
    upwards = 0;
else
    error('status not recognized')
end

while 1==1
    
    % Change in alpha
    dA = alfa2 - alfa1;
    
    % Current alpha
    cAlfa = alfa1 + cutt.*dA;
    
    % Break, if we've run out of steps
    if (cutt<=0)
        break
    
    % If climbing to plateau and cAlfa == 1 . . .
    elseif upwards && (dumb_sim(t,cAlfa,Upred,Uprey,x0,y0)/Dmin0 == 1)
        
        % Break loop
        break
        
    % If decending from plateau and cAlfa ~= 1 . . .
    elseif ~upwards && (dumb_sim(t,cAlfa,Upred,Uprey,x0,y0)/Dmin0 ~= 1)
        
        % Break loop
        break
        
    % Otherwise, cut current alpha value
    else
        cutt = cutt + step_size;
    end
        
end


function Dmin = dumb_follower(t,alfa,Upred,Uprey,x0,y0)
% Predator follows prey

% Initial position values
xPrey = x0;
yPrey = y0;
xPred = 0;
yPred = 0;

for i = 2:length(t)
    
    % Time step
    dt = (t(i)-t(i-1));
    
    % Orientations
    angPred = atan2(yPrey(i-1)-yPred(i-1),xPrey(i-1)-xPred(i-1));
    angPrey = angPred + alfa;
    
    % Displacements
    dispPrey = Uprey * dt;
    dispPred = Upred * dt;
    
    % Prey position
    xPrey(i,1) = dispPrey.*cos(angPrey) + xPrey(i-1);
    yPrey(i,1) = dispPrey.*sin(angPrey) + yPrey(i-1);
    
    % Predator position
    xPred(i,1) = dispPred.*cos(angPred) + xPrey(i-1);
    yPred(i,1) = dispPred.*sin(angPred) + yPrey(i-1);
    
    % Distance between players
    D = hypot(xPrey(i)-xPred(i),yPrey-yPred(i));
    
end

if 0
    figure
%     for i = 1:length(t)
%         h = plot([xPrey(i) xPred(i)],[yPrey(i) yPred(i)],'-k');
%         h = set(h,'Color',.8.*[1 1 1]);
%         hold on
%     end
    
    plot(xPrey, yPrey, '-+', xPred, yPred,'-+')
    title(['alpha = ' num2str(alfa*180/pi) ', K = ' num2str(Upred/Uprey)])
    hold on
    xlabel('t')
    ylabel('D')
    axis equal
    grid on
end

        

Dmin = min(D);


function Dmin = dumb_sim(t,alfa,Upred,Uprey,x0,y0)
% Predator moves straight


% Predator position (starting at zero)
xPred = Upred.*t;

xPrey = Uprey.*t.*cos(alfa) + x0;
yPrey = Uprey.*t.*sin(alfa) + y0;

D = hypot(xPrey-xPred,yPrey);

if 0
    plot(t,D);
    title(['alpha = ' num2str(alfa*180/pi) ', K = ' num2str(Upred/Uprey)])
    hold on
    xlabel('t')
    ylabel('D')
    grid on
end

Dmin = min(D);


function Dmin = dumb_sim2(t,alfa,Upred,Uprey,x0,y0)
% Predator moves straight


% Predator position (starting at zero)
xPred = Upred.*t;

xPrey = Uprey.*t.*cos(alfa) + x0;
yPrey = Uprey.*t.*sin(alfa) + y0;

D = hypot(xPrey-xPred,yPrey);

if 0
    plot(t,D);
    title(['alpha = ' num2str(alfa*180/pi) ', K = ' num2str(Upred/Uprey)])
    hold on
    xlabel('t')
    ylabel('D')
    grid on
end

Dmin = min(D);

function Dmin = eqn35(K,alpha)
% Equation 35 from Weihs and Webb

Dmin = sqrt(sin(alpha).^2 ./ (K.^2 - 2.* K .* cos(alpha) +1 ) );
