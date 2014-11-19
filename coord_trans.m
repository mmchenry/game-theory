function [xRT,yRT] = coord_trans(theta,origin,x,y,type)
% Transforms points between coordinate systems 



switch type
    
    case 'global to body'
    
        % unit vector axes
        yaxis = [-sin(theta) cos(theta) 0]./norm([-sin(theta) cos(theta) 0]);
        xaxis = [cos(theta) sin(theta) 0]./norm([cos(theta) sin(theta) 0]);
        zaxis = [0 0 1];
        
        %Create rotation matrix (from inertial axes to local axes)
        R = [xaxis' yaxis' zaxis'];
        
        % Package points in a matrix
        pts = [x y z]';
        
        % Rotate points
        ptsT = [R * pts]';
        
        % Extract points
        xRT = ptsT(:,1);
        yRT = ptsT(:,2);        
    
        
   case 'body to global'
    
        % unit vector axes
        yaxis = [-sin(theta) cos(theta) 0]./norm([-sin(theta) cos(theta) 0]);
        xaxis = [cos(theta) sin(theta) 0]./norm([cos(theta) sin(theta) 0]);
        zaxis = [0 0 1];
        
        %Create rotation matrix (from inertial axes to local axes)
        R = [xaxis' yaxis' zaxis'];
        
        % Package points in a matrix
        pts = [x y y.*0]';
        
        % Rotate points
        ptsT = [R * pts]';
        
        % Extract points
        xRT = ptsT(:,1) + origin(1);
        yRT = ptsT(:,2) + origin(2);       
        

    case 'global to eye'
    % TODO: This needs work (will be based on 'model_vision.m')    
        
        % Define coordinate systems for the two eyes
        [rS,lS] = eye_systems(head,tail,verg_ang,vis_az);
        
        psi = pi/2 + theta;
        
        % unit vector axes
        yaxis = [-sin(psi) cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
        xaxis = [cos(psi) sin(psi) 0]./norm([cos(psi) sin(psi) 0]);
        zaxis = [0 0 1];
        
        %Create rotation matrix (from inertial axes to local axes)
        R = [xaxis' yaxis' zaxis'];
        
        % Package points in a matrix
        pts = [x y z]';
        
        % Rotate points
        ptsT = [R * pts]';
        
        % Extract points
        xRT = ptsT(:,1);
        yRT = ptsT(:,2);
        zRT = ptsT(:,3);
        
        % Visual test
        if 1
            figure
            plot3(xRT,yRT,zRT,'.'); %axis equal
            xlabel('X'); ylabel('Y');zlabel('Z')
            view(2)
        end
        

end






% % Vergence angle: mean = 32 deg, max = 70 deg (during predation) (Patterson et al, 2013)
% verg_ang = 32/180*pi;
% 
% % Receptive field of visual system (m, rad). Easter (1996) give the functional 
% % retinal field of 163 deg for 4 dpf. Though I give an additional 20 deg
% % for saccades (based on vergenece angles of Patterson et al, 2013)
% vis_depth  = 0.15;
% vis_az     = (163+20)/180*pi;
% vis_el     = (163)/180*pi;
% off_el     = 0;
% vis_num    = 100;




function [sR,sL] = eye_systems(head,tail,verg_ang,fov)

% Check dimensions
if size(head,1)~=1 || size(head,2)~=3 || size(tail,1)~=1 || size(tail,2)~=3 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - head(1);
xaxis(1,2) = tail(2) - head(2);
xaxis(1,3) = tail(3) - head(3);

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

% Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

% Create rotation matrix (from inertial axes to local axes)
sH = [xaxis; yaxis; zaxis];

% psi - forward tilt angle of an eye relative to the body
psi = pi/2 + verg_ang - fov/2;

% Unit vector axes for R eye (in body system)
yaxisR = [-sin(psi) cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
xaxisR = [cos(psi) sin(psi) 0]./norm([cos(psi) sin(psi) 0]);
zaxisR = [0 0 1];

% Create rotation matrix (from inertial axes to local axes)
sR = [xaxisR; yaxisR; zaxisR];

% Rotate the local system according to the eye position
sR = [sR * sH];

% Unit vector axes for L eye (in body system)
yaxisL = [-sin(psi) -cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
xaxisL = [-cos(psi) sin(psi) 0]./norm([-cos(psi) sin(psi) 0]);
zaxisL = [0 0 1];

% Create rotation matrix (from inertial axes to local axes)
sL = [xaxisL; yaxisL; zaxisL];

% Rotate the local system according to the eye position
sL = [sL * sH];



function R = prey_system(head,tail)

% Check dimensions
if size(head,1)~=1 || size(head,2)~=3 || size(tail,1)~=1 || size(tail,2)~=3 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - head(1);
xaxis(1,2) = tail(2) - head(2);
xaxis(1,3) = tail(3) - head(3);

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; zaxis];

% Take inverse
%R = inv(R);


function [pXg,pYg,pZg] = prey_local(head,R,pX,pY,pZ)
% Transforms polygon points from global to local FOR

tmp1 = global_to_local(head,R,[pX(1,:)' pY(1,:)' pZ(1,:)']);
tmp2 = global_to_local(head,R,[pX(2,:)' pY(2,:)' pZ(2,:)']);
tmp3 = global_to_local(head,R,[pX(3,:)' pY(3,:)' pZ(3,:)']);
tmp4 = global_to_local(head,R,[pX(4,:)' pY(4,:)' pZ(4,:)']);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];

function [pXg,pYg,pZg] = prey_global(head,R,pX,pY,pZ)
% Transforms polygon points from local to global FOR

tmp1 = local_to_global(head,R,[pX(1,:)' pY(1,:)' pZ(1,:)']);
tmp2 = local_to_global(head,R,[pX(2,:)' pY(2,:)' pZ(2,:)']);
tmp3 = local_to_global(head,R,[pX(3,:)' pY(3,:)' pZ(3,:)']);
tmp4 = local_to_global(head,R,[pX(4,:)' pY(4,:)' pZ(4,:)']);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];


function ptsT = global_to_local(head,R,pts)

% Translate
pts(:,1) = pts(:,1) - head(1);
pts(:,2) = pts(:,2) - head(2);
pts(:,3) = pts(:,3) - head(3);

% Rotate points
ptsT = [R * pts']';



function ptsT = local_to_global(head,R,pts)

% Rotate points
ptsT = [inv(R) * pts']';

% Translate global coordinates wrt head
ptsT(:,1) = ptsT(:,1) + head(1);
ptsT(:,2) = ptsT(:,2) + head(2);
ptsT(:,3) = ptsT(:,3) + head(3);


