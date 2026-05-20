
close all
clear
clc

%addpath(genpath('./functions'));
%rng(5)
%rng(11)

%% 1. Smoothing
addpath('C:\Users\little_long\Desktop\SSY345\HA4')

% Generate the trajectory.
[X,T,Tvec] = generateTrueTrack();

% Set the intial prior
x0 = [0; 0; 0; -pi/2; 0];
P0 = diag([10^2, 10^2, 10^2, ((5*pi)/180)^2, ((5*pi)/180)^2 ]);

% Sensor position
s1 = [160; -300];
s2 = [420; -300];



% Measurement noise
sigma_phi_1 = 0.5*pi/180;
sigma_phi_2 = 0.5*pi/180;
R = diag([(sigma_phi_1).^2, (sigma_phi_2).^2]).*T;

% Process noise covariance 
sigma_v = 0.1;
sigma_w = pi/180;
Q = diag([0 ,0 ,(sigma_v)^2, 0 ,(sigma_w)^2]).*T;


% Generate measurements
h = @(x,T) dualBearingMeasurement(x,s1,s2);
Y = genNonLinearMeasurementSequence(X, h, R);
[xpos, ypos] = getPosFromMeasurement(Y, s1, s2);

% Dynamics model
proc_f = @(x,T) coordinatedTurnMotion(x, T);
sigmaPoints = @sigmaPoints;

% Dummy variable for sensor pos
S = zeros(2, size(Y, 2));

% Filter the position with non-linear Kalman filter
type = 'EKF';
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x0, P0, proc_f, T, Q, S, h, R, sigmaPoints, type);


% Plotting
plotFilterSmoother(X,xs, Ps, xf, Pf, xp, Pp, xpos, ypos, s1, s2, type, Tvec)



% b) Modify a measurement at k=400


Y(:,400) = Y(:,400) + rand();
[xpos, ypos] = getPosFromMeasurement(Y, s1, s2);

% Filter the position with non-linear Kalman filter
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x0, P0, proc_f, T, Q, S, h, R, sigmaPoints, type);

% Plot everything
plotFilterSmoother(X,xs, Ps, xf, Pf, xp, Pp, xpos, ypos, s1, s2, type, Tvec);





%% 2. Particle filters for linear/Gaussian systems
close all
clear
clc
rng(5)
addpath('C:\Users\little_long\Desktop\SSY345\HA4')

% Process noise
proc_Q = 1.5;
% Measurement noise
meas_R = 0.5;

% Initial prior and cov 2A
x0 = 2;
P0 = 8;

% Timesteps
K = 30;

% Dynamics model & Measurement model
proc_f = @(x) x;
meas_h = @(x) x;

% Generate state and measurement sequences
X = genLinearStateSequence(x0,P0,proc_f,proc_Q,K);
Y = genLinearMeasurementSequence(X, meas_h, meas_R);

% Incorrect prior UNCOMMENT FOR 2b
x0 = -20;
P0 = 2;

% Filter with kalman filter
[xf, Pf] = linearKalmanFilter(Y, x0, P0, 1, proc_Q, 1, meas_R);

% Filter with particle filter
N = 1000;
[xfp, Pfp, Xp, Wp] = pfFilter(x0, P0, Y, proc_f, proc_Q, meas_h, meas_R, N, false);
[xfp_resam, Pfp_resam, Xp_resam, Wp_resam] = pfFilter(x0, P0, Y, proc_f, proc_Q, meas_h, meas_R, N, true);


% Plot everything
plotParticleFilter(X,Y,K,xf,Pf,xfp,Pfp,xfp_resam,Pfp_resam);

% Plot posterior approx
%figure('Name', 'Task 2a: Posterior Approximations');
time_instances = [1, 15, 30];
plotPosteriorApproximation(xf, Pf, xfp, Pfp, xfp_resam, Pfp_resam, time_instances);



% Calculate RMSE 2a
err_KF = sqrt(immse( X(2:end) , xf ))
err_PF = sqrt(immse( X(2:end) , xfp))
err_PFresam = sqrt(immse( X(2:end) , xfp_resam))

% Plot errobar
%figure('Name', 'Task 2a: Trajectories & Errorbars');
% figure('Name', 'Task 2b: Catching up with Incorrect Prior');
plotParticleFilterErrorbar(X,Y,K,xf,Pf,xfp,Pfp,xfp_resam,Pfp_resam)



% Filter with particle filter
x0_2c = -20;
P0_2c = 2;
N_2c = 10;

bResampl = true;
% Without Resampling
[xfp_no_resam, Pfp_no_resam, Xp_no_resam, Wp_no_resam] = pfFilter(x0_2c, P0_2c, Y, proc_f, proc_Q, meas_h, meas_R, N_2c, false);
figure('Name', 'Task 2c: Trajectories Without Resampling');
plotParticleTrajectory(Xp_no_resam, Wp_no_resam, K, N_2c, X, xfp_no_resam, Pfp_no_resam, false);
title('Task 2c: N=10 Without Resampling');

% With Resampling
[xfp_with_resam, Pfp_with_resam, Xp_with_resam, Wp_with_resam] = pfFilter(x0_2c, P0_2c, Y, proc_f, proc_Q, meas_h, meas_R, N_2c, true);
figure('Name', 'Task 2c: Trajectories With Resampling');
plotParticleTrajectory(Xp_with_resam, Wp_with_resam, K, N_2c, X, xfp_with_resam, Pfp_with_resam, true);
title('Task 2c: N=10 With Resampling');







%% 3. Bicycle tracking in a village
close all;
clear;
clc;
rng(5)

% 3a  Generate bicycle trajectory interactively
%
%XP = MapProblemGetPoint() [2 x K] 
load('Xk.mat');   
XP = Xk;         

%  3b Generate velocity measurements
sigma_r = 0;   % noise sigma
R3 = sigma_r^2 * eye(2);

K3 = size(XP, 2);   

% true velocity: xv_k = xp_k - xp_{k-1}
XV = diff(XP, 1, 2);   % [2 x K-1]

% velocity measurement with noise
Y3 = XV + sigma_r * randn(size(XV));

% Task 3d - PF with known initial position
N_particles = 2000;
sigma_q = 0;

x_init = XP(:, 1);  % initial position
v_init = XP(:, 2) - XP(:, 1);  % initialization using the first step velocity

particles = zeros(4, N_particles);
particles(1:2, :) = repmat(x_init, 1, N_particles) + 0.01*randn(2, N_particles);
particles(3:4, :) = repmat(v_init, 1, N_particles) + 0.01*randn(2, N_particles);

weights = ones(1, N_particles) / N_particles;
x_est = zeros(4, K3);

figure(2);

for k = 1:K3
    
    particles(1:2, :) = particles(1:2, :) + particles(3:4, :) + ...
                        sigma_q * randn(2, N_particles);
    particles(3:4, :) = particles(3:4, :) + ...
                        sigma_q * randn(2, N_particles);
    
   
    on_road = isOnRoad(particles(1,:)', particles(2,:)')';
    weights = weights .* on_road;
    
  
    innovation = Y3(:, k) - particles(3:4, :);
    likelihoods = exp(-0.5 * sum(innovation .* (R3 \ innovation), 1));
    weights = weights .* likelihoods;
    
    if sum(weights) < 1e-10
        weights = ones(1, N_particles) / N_particles;
    else
        weights = weights / sum(weights);
    end
    
   
    x_est(:, k) = sum(particles .* weights, 2);
    
  
    [particles, weights, ~] = resampl(particles, weights);
    
   
    clf; hold on;
    plot([1+i 1+9*i 5+9*i])
    plot([7+9*i 11+9*i 11+i 7+i]); plot([5+i 1+i])
    plot([2+5.2*i 2+8.3*i 4+8.3*i 4+5.2*i 2+5.2*i])
    plot([2+3.7*i 2+4.4*i 4+4.4*i 4+3.7*i 2+3.7*i])
    plot([2+2*i 2+3.2*i 4+3.2*i 4+2*i 2+2*i])
    plot([5+i 5+2.2*i 7+2.2*i 7+i])
    plot([5+2.8*i 5+5.5*i 7+5.5*i 7+2.8*i 5+2.8*i])
    plot([5+6.2*i 5+9*i]); plot([7+9*i 7+6.2*i 5+6.2*i])
    plot([8+4.6*i 8+8.4*i 10+8.4*i 10+4.6*i 8+4.6*i])
    plot([8+2.4*i 8+4*i 10+4*i 10+2.4*i 8+2.4*i])
    plot([8+1.7*i 8+1.8*i 10+1.8*i 10+1.7*i 8+1.7*i])
    axis([0.8 11.2 0.8 9.2])
    scatter(particles(1,:), particles(2,:), 5, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    plot(x_est(1,1:k), x_est(2,1:k), 'r-', 'LineWidth', 2);
    plot(XP(1,1:k+1), XP(2,1:k+1), 'k-', 'LineWidth', 2);
    legend('','','','','','','','','','','','Particles','Estimated','True','Location','best')
    title(sprintf('PF tracking (known init) — step k=%d / %d', k, K3));
    axis([0.8 11.2 0.8 9.2])
    drawnow;
end

% Task 3e - PF without known initial position
N_particles_e = 5000;
sigma_q_e = 0.05;
sigma_r_e = 0.01;
R3_e = sigma_r_e^2 * eye(2);

n_try = N_particles_e * 20;
candidates_x = 1 + 10*rand(1, n_try);
candidates_y = 1 +  8*rand(1, n_try);


on_road_init = isOnRoad(candidates_x', candidates_y');
candidates_x = candidates_x(logical(on_road_init));
candidates_y = candidates_y(logical(on_road_init));

candidates_x = candidates_x(1:N_particles_e);
candidates_y = candidates_y(1:N_particles_e);

particles_e = zeros(4, N_particles_e);
particles_e(1, :) = candidates_x;
particles_e(2, :) = candidates_y;
particles_e(3:4, :) = 0.005 * randn(2, N_particles_e);


weights_e = ones(1, N_particles_e) / N_particles_e;
x_est_e = zeros(4, K3);

figure(3);

for k = 1:K3
   
    particles_e(1:2, :) = particles_e(1:2, :) + particles_e(3:4, :) + ...
                          sigma_q_e * randn(2, N_particles_e);
    particles_e(3:4, :) = particles_e(3:4, :) + ...
                          sigma_q_e * randn(2, N_particles_e);
    
   
    on_road_e = isOnRoad(particles_e(1,:)', particles_e(2,:)');
    weights_e = weights_e .* on_road_e';
    
   
    innovation_e = Y3(:, k) - particles_e(3:4, :);
    likelihoods_e = exp(-0.5 * sum(innovation_e .* (R3_e \ innovation_e), 1));
    weights_e = weights_e .* likelihoods_e;
    
    
   if sum(weights_e) < 1e-10
        disp(['Warning: weight degeneracy at k = ', num2str(k)])
        weights_e = ones(1, N_particles_e) / N_particles_e;
    else
        weights_e = weights_e / sum(weights_e);
    end
    
   x_est_e(:, k) = sum(particles_e .* weights_e, 2);
    
   
    [particles_e, weights_e, ~] = resampl(particles_e, weights_e);
    
   
clf; hold on;

plot([1+i 1+9*i 5+9*i])
plot([7+9*i 11+9*i 11+i 7+i]); plot([5+i 1+i])
plot([2+5.2*i 2+8.3*i 4+8.3*i 4+5.2*i 2+5.2*i])
plot([2+3.7*i 2+4.4*i 4+4.4*i 4+3.7*i 2+3.7*i])
plot([2+2*i 2+3.2*i 4+3.2*i 4+2*i 2+2*i])
plot([5+i 5+2.2*i 7+2.2*i 7+i])
plot([5+2.8*i 5+5.5*i 7+5.5*i 7+2.8*i 5+2.8*i])
plot([5+6.2*i 5+9*i]); plot([7+9*i 7+6.2*i 5+6.2*i])
plot([8+4.6*i 8+8.4*i 10+8.4*i 10+4.6*i 8+4.6*i])
plot([8+2.4*i 8+4*i 10+4*i 10+2.4*i 8+2.4*i])
plot([8+1.7*i 8+1.8*i 10+1.8*i 10+1.7*i 8+1.7*i])


h_particles = scatter(particles_e(1,:), particles_e(2,:), 5, 'b', ...
                      'filled', 'MarkerFaceAlpha', 0.3);
h_est  = plot(x_est_e(1,1:k), x_est_e(2,1:k), 'r-', 'LineWidth', 2);
h_true = plot(XP(1,1:k+1),    XP(2,1:k+1),    'k-', 'LineWidth', 2);

legend([h_particles, h_est, h_true], ...
       {'Particles', 'Estimated', 'True'}, 'Location', 'best')
n_alive = sum(weights_e > 1/N_particles_e * 0.01);
title(sprintf('PF (unknown init) — k=%d/%d,  active particles≈%d', k, K3, n_alive));
axis([0.8 11.2 0.8 9.2])
drawnow;
end