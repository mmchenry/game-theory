function run_sim
% Runs single simulation and analyzes results


% Default parameter values
[p,d] = default_params;


%% Run simulation

R = simple_model(p,d);

% Store input parameters
R.p     = p;
R.d     = d;


%% Display results

% Animate results with even time intervals
%animate_sim(R)

% Plot orientation angles
%figure;
%vis_results('Turning data',R)

figure;
vis_results('Trajectories',R)

%figure('DoubleBuffer','on')
%vis_results('Animate',R)
% TODO: Fix the animation code

% Display whether or not prey escaped
if isempty(R.tEnd)
    disp('prey successfully escaped')
else
    disp(['the prey was captured at time = ' num2str(R.tEnd)])
end


