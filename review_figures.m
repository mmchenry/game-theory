function review_figures

% Paramaters
alpha_max = 150;
alpha_min = 0;

% Number of alpha and K values inlcuded in plot
num_alpha   = 40;
num_K       = 40;

% Creat vectors of alpha and K
alpha_vect = linspace(alpha_min,alpha_max,num_alpha)'.*pi/180;
K_vect     = [10.^linspace(-1,0,round(num_K/3))'; ...
              10.^linspace(0,2,round(2*num_K/3))'];

% Create higher resolution version of K for the optimal curve
K_lin      = 10.^linspace(-1,2,num_K.*100)';

% Creat matricies of alpha and K values
[K,alpha] = meshgrid(K_vect,alpha_vect);

% Solve for minumum distance for meshgrid values
Dmin = eqn35(K,alpha);

% Solve high resolution alpha values and corresponding Dmins
alpha_lin = [acos(K_lin(K_lin<=1)); acos(1./K_lin(K_lin>1))];
Dmin_lin = eqn35(K_lin,alpha_lin);



figure
%subplot(2,1,1)

h = mesh(K,alpha.*180/pi,Dmin);
hold on
%view([40 19])
view([49 25])

hL = plot3(K_lin,alpha_lin.*180/pi,Dmin_lin,'k-');
hold off

set(gca,'XScale','log')

xlabel('K');
ylabel('alpha (deg)')
zlabel('D_m_i_n*')
set(hL,'LineWidth',2)

caxis([0 1.2])
set(gca,'TickLength',[0.03 0.03])
set(gca,'YTick',0:30:alpha_max)
set(gca,'ZTick',0:0.25:1)
zlim([0 1])
ylim([0 alpha_max])
colorbar('EastOutside')

%subplot(2,1,2)
%semilogx(K_lin,alpha_lin.*180/pi)


return


K_ex = 10.^linspace(-1,2,8)';

alpha_ex = [acos(K_ex(K_ex<=1)); acos(1./K_ex(K_ex>1))];

Dmin2 = eqn35(K_ex,alpha_ex);

for i = 1:length(K_ex)
    
    Dmin_ex(:,i) = eqn35(K_ex(i),alpha_lin);
    l_text{i} = num2str(K_ex(i));
end


figure
plot(alpha_lin.*180/pi,Dmin_ex);

xlabel('alpha (deg)')
ylabel('D_m_i_n*')
legend(l_text)

hold on
plot(alpha_ex.*180/pi,Dmin2,'ok');
hold off


function Dmin = eqn35(K,alpha)
% Equation 35 from Weihs and Webb

Dmin = sqrt(sin(alpha).^2 ./ (K.^2 - 2.* K .* cos(alpha) +1 ) );