
close all
clear all
clc

% Parameters
mu = [2; 10];
Sigma = [1, 0; 0, 8];
level = 3;
npoints = 100;

A = [1, -0.5; 0, 1];

% 
[ xy_q ] = sigmaEllipse2D( mu, Sigma, level, npoints );

% Affine transformation
[mu_y, Sigma_y] = affineGaussianTransform(mu, Sigma, A, [0; 0])
[ xy_z ] = sigmaEllipse2D( mu_y, Sigma_y, level, npoints );


% Plot
figure();
plot(xy_q(1,:),xy_q(2,:),'linewidth',2)
axis equal
hold on
grid on
plot(mu(1), mu(2), '*')

plot(xy_z(1,:),xy_z(2,:),'linewidth',2)
axis equal
hold on
grid on
plot(mu_y(1), mu_y(2), '*')
legend('3\sigma - curve for q', 'mean q','3\sigma - curve for z','mean z')
xlabel('x_1'); ylabel('x_2')

%% 2 a)
close all
clear all
clc

% Define gaussian variable x
mu_x = 0;
sigma2_x = 1;
N = 50000;

x_samples = sqrt(sigma2_x) * randn(N, 1) + mu_x;
A = 5; b = 3;
mu_z_anal = A * mu_x + b;          % 3
sigma2_z_anal = A^2 * sigma2_x;    %  25
% Nummerically calculate z
z_samples = 5 * x_samples + 3;
[mu_z_approx, sigma2_z_approx] = approxGaussianTransform(mu_x, sigma2_x, @(x) 5*x + 3, N);


figure();
axi = linspace(mu_z_anal - 4*sqrt(sigma2_z_anal), mu_z_anal + 4*sqrt(sigma2_z_anal), 200);
histogram(z_samples, 100, 'Normalization', 'pdf', 'DisplayName', 'Transformed Samples');
hold on;

plot(axi, normpdf(axi, mu_z_anal, sqrt(sigma2_z_anal)), 'r', 'LineWidth', 2, 'DisplayName', 'Analytical PDF');

% 绘制数值近似 PDF
plot(axi, normpdf(axi, mu_z_approx, sqrt(sigma2_z_approx)), 'g--', 'LineWidth', 2, 'DisplayName', 'Approx Gaussian PDF');

title('Linear Transformation: z = 5x + 3');
xlabel('z'); ylabel('p(z)');
legend('Location', 'northeast');
grid on;


%% 2 b)
close all
clear all
clc

% Define gaussian variable x
mu_x = 0;
sigma2_x = 1;
N = 50000;

x_samples = sqrt(sigma2_x) * randn(N, 1) + mu_x;
z_samples = exp(x_samples);
% Nummerically calculate z
[mu_z_approx, sigma2_z_approx] = approxGaussianTransform(mu_x, sigma2_x, @(x) exp(x), N);


% plot
figure();

axi = linspace(0, 10, 500); 


histogram(z_samples, 'BinLimits', [0, 10], 'Normalization', 'pdf', 'DisplayName', 'True Samples (Log-normal)');
hold on;


plot(axi, normpdf(axi, mu_z_approx, sqrt(sigma2_z_approx)), 'r', 'linewidth', 2, 'DisplayName', 'Gaussian Approximation');

grid on;
title('Non-linear Transformation: z = e^x');
xlabel('z'); ylabel('p(z)');
legend('Location', 'northeast');




%% 3a
clear; clc; close all;

N = 100000;
sigma_r = 0.5;
r = sigma_r * randn(N, 1);

x = 10 * rand(N, 1) - 5; 

h_lin = @(x) 2*x;
h_nonlin = @(x) x.^2;

% p(y)
y_lin = h_lin(x) + r;
y_nonlin = h_nonlin(x) + r;

figure(1);
subplot(2,1,1);
histogram(y_lin, 100, 'Normalization', 'pdf');
title('Marginal p(y) with Linear h(x) and Uniform x');
grid on;

%  3b: p(y|x) 
%  x = 2
target_x = 2;
idx = abs(x - target_x) < 0.02; 
y_cond = y_nonlin(idx);

subplot(2,1,2);
histogram(y_cond, 50, 'Normalization', 'pdf'); hold on;
ax = linspace(min(y_cond), max(y_cond), 100);
plot(ax, normpdf(ax, h_nonlin(target_x), sigma_r), 'r', 'LineWidth', 2);
title(['Conditional p(y|x=', num2str(target_x), ') for Non-linear h(x)']);
legend('Samples', 'Analytical Gaussian');
grid on;





%% 4 a)
clear; clc; close all;
N = 100000;
sigma_w = 0.5;
theta = ones(N, 1);
theta(rand(N, 1) < 0.5) = -1;
w = sigma_w * randn(N, 1);
y = theta + w;

figure();
histogram(y, 100, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8]);
hold on;

% Gaussian Mixture
axi = linspace(-4, 4, 1000);
p_y = 0.5 * normpdf(axi, -1, sigma_w) + 0.5 * normpdf(axi, 1, sigma_w);
plot(axi, p_y, 'r', 'LineWidth', 2);

title('Histogram of Observation y = \theta + w (\sigma_w = 0.5)');
xlabel('y'); ylabel('p(y)');
legend('Samples', 'Analytical GMM');
grid on;
%%
copyfile('mainExample.m', 'mainExample.txt');




