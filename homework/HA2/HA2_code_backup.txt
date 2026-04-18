%% ASSIGNMENT 2

close all
clear all
clc


% Add path to functions folder
addpath('./functions');



%% Scenario 1 - A First Kalman filter and its properties
% a) Generate state sequence and measurements sequence
close all
clc
rng(2)

% Sequence length
N = 35;

% initial prior
x_0 = 2;
P_0 = 8;

% Measurments noise covariance and process disturbance covariance
Q = 1.5;
R = 3;

% Dynamics matrix
A = 1;
H = 1;


% Generate state sequence
x_seq = genLinearStateSequence(x_0, P_0, A, Q, N);

% Generate state measurements
y_seq = genLinearMeasurementSequence(x_seq, H, R);

% Plot the results
figure();
hold on;
grid on;
plot(0:N, x_seq, 'LineWidth', 2, 'DisplayName', 'State sequence x');
plot(1:N, y_seq, 'LineWidth', 2, 'DisplayName', 'Measurements y'); 
legend('Location', 'best');
title('Scenario 1: Linear Gaussian State Space Model');
xlabel('Time step k');
ylabel('Value');



% b) Filter the measurements through the Kalman filter
[x_kalmf, P] = kalmanFilter(y_seq, x_0, P_0, A, Q, H, R);

std = sqrt(reshape(P,1,35));


% Plot the results
figure();
hold on;
grid on;
plot(0:N, x_seq,'Linewidth', 2);
plot(1:N, y_seq,'linewidth', 2);
plot(1:N, x_kalmf,'linewidth', 2);
plot(1:N, x_kalmf+(3*std),':','linewidth', 2,'color','k');
plot(1:N, x_kalmf-(3*std),':','linewidth', 2,'color','k');
legend('state sequence x','state measurements y','Kalman filtered measurements' ...
    ,'kalmf + 3\sigma','kalmf - 3\sigma');
title('Linear Gaussian state space model');
xlabel('N');


% Plot error densities
% figure()
% subplot(4,1,1)
% counter = 0;
% 
% for i = [1, 2, 4, 30]
%     counter = counter + 1;
%     err = (x_seq(i+1) - x_kalmf(i))^2;
% 
% 
%     a = mvnrnd(0,err,10000);
%     subplot(4,1,counter)
%     histogram(a,75,'Normalization','pdf')
%     hold on
% end








% c) Filter the measurements through the Kalman filter (x0 = 12)
[x_kalmf_err, P] = kalmanFilter(y_seq, 12, P_0, A, Q, H, R);


% Plot the results
figure()
hold on
grid on
plot(1:N, x_kalmf,'linewidth', 2)
plot(1:N, x_kalmf_err,':','linewidth', 2)
legend('correct Kalman','Incorrect Kalman')
title('Kalman filters for right/wrong assumption on x_0')
xlabel('N');



% d) Plot yada ydada...
% We choose k = 23

k = 23;


figure()
ticks = -15:0.01:0

% Prior
one = normpdf(ticks,x_kalmf(k-1),std(k-1));
plot(ticks,one,'linewidth',2)
hold on
grid on

% Prediction
[x_pred, P_pred] = linearPrediction(x_kalmf(k-1), std(k-1)^2, A, Q);
two = normpdf(ticks,x_pred,sqrt(P_pred));
plot(ticks,two,'linewidth',2)

% Measurement
three = normpdf(ticks,y_seq(k),R);
plot(ticks,three,'linewidth',2)

% Update
[x_upd, P_upd] = linearUpdate(x_pred, P_pred, y_seq(k), H, R);
four = normpdf(ticks,x_upd,sqrt(P_upd));
plot(ticks,four,'linewidth',2)


% Plot details
legend('prior: P(x_{k-1} | y_{1:k-1})','prediction: P(x_{k} | y_{1:k-1})','measurements: y','update: P(x_{k} | y_{1:k})')
title('PDF for the different Kalman steps')


% e)
N_long = 1000; % 明确定义总步数为 1000 [cite: 44, 46]

% 1. 生成状态序列 (包含 x0 到 x1000，共 1001 个点)
x_true = genLinearStateSequence(x_0, P_0, A, Q, N_long); 

% 2. 生成测量序列 (对应 x1 到 x1000，共 1000 个点) 
% 注意：必须传入从索引 2 开始的状态 (即 x1)
y_meas = genLinearMeasurementSequence(x_true(:, 1:end), H, R); 

% 验证长度：此时 size(y_meas, 2) 应该等于 1000
if size(y_meas, 2) < N_long
    error('测量序列长度不足，请检查数据生成逻辑');
end

% 3. 初始化存储变量
x_est = zeros(1, N_long + 1);
P_est = zeros(1, N_long + 1);
innovations = zeros(1, N_long); % 存储 v_k 用于后续分析 [cite: 49, 51]

x_est(1) = x_0; % 初始先验均值 [cite: 22]
P_est(1) = P_0; % 初始先验协方差 [cite: 22]

% 4. 滤波循环
for k = 1:N_long
    % 预测步骤 (Prediction Step) [cite: 19, 39]
    x_pred = A * x_est(k);
    P_pred = A * P_est(k) * A' + Q;
    
    % 更新步骤 (Update Step) [cite: 19, 39, 51]
    v_k = y_meas(k) - H * x_pred; % 计算创新值 (Innovation) [cite: 49, 51]
    S_k = H * P_pred * H' + R;
    K_k = P_pred * H' / S_k;
    
    x_est(k+1) = x_pred + K_k * v_k;
    P_est(k+1) = (1 - K_k * H) * P_pred;
    
    innovations(k) = v_k; % 保存用于 Task e 的统计分析 [cite: 49]
end

% 5. 结果分析与绘图 (Task e 要求)
% 计算估计误差 (只分析达到稳态后的部分) [cite: 46]
err = x_true(2:end) - x_est(2:end); 

figure;
% A) 误差直方图与理论 PDF 对比 [cite: 46]
subplot(2,1,1);
histogram(err, 'Normalization', 'pdf', 'DisplayName', 'Error Histogram');
hold on;
% 绘制理论 PDF: N(0, P_k|k) [cite: 46]
P_inf = P_est(end); % 稳态协方差 [cite: 42]
range = linspace(min(err), max(err), 100);
plot(range, normpdf(range, 0, sqrt(P_inf)), 'r', 'LineWidth', 2);
title('Estimation Error Histogram vs Steady-state PDF');
legend('Empirical Error', 'Theoretical PDF');

% B) 创新值自相关分析 
subplot(2,1,2);
autocorr(innovations, 'NumLags', 20); % 
title('Innovation Process Autocorrelation');
%% Scenario 2 - Tuning a Kalman filter
close all
clear all
clc

% Load the data
trainMeasurements = load('SensorMeasurements.mat');

% Extract
pureNoise = trainMeasurements.CalibrationSequenceVelocity_v0;
trainSpeed10 = trainMeasurements.CalibrationSequenceVelocity_v10;
trainSpeed20 = trainMeasurements.CalibrationSequenceVelocity_v20;

% Determine C
C1 = pureNoise / pureNoise;
C2 = trainSpeed10 / (10 + pureNoise);
C3 = trainSpeed20 / (20 + pureNoise);


C2 = (mean(trainSpeed10) - mean(pureNoise)) / 10;
C2 = (mean(trainSpeed20) - mean(pureNoise)) / 20;
C = mean([C2 C3])


% Calculate variance
var2 = var(trainSpeed10 / C);
var3 = var(trainSpeed20 / C);

% Try to replicate the noisy behaviour
fakeTrainSpeed10 = C2*(10 + mvnrnd(mean(pureNoise), var2, 2000)');
fakeTrainSpeed20 = C3*(20 + mvnrnd(mean(pureNoise), var3, 2000)');

figure();
subplot(2,1,1);
plot(trainSpeed10);
hold on;
plot(fakeTrainSpeed10);
legend('Actual readings','generated readings');

subplot(2,1,2);
plot(trainSpeed20);
hold on;
plot(fakeTrainSpeed20);
legend('Actual readings','generated readings');


% =========================================================================
% Call Kalman filter (CONSTANT VEL)
h = 2*0.1;
alpha = 0.0001;
Y = Generate_y_seq;
x_0 = [0;0];
P_0 = [0.5 0.5; 
       0.5 0.5];

A_cv = [1  h; 
        0  1];

% Process disturbances
Q = alpha*[0   0
     0   1];

H = [1  0;
     0  1];

% Measurement noise
R = [1   0
    0  mean([var2, var3])];

% Remove every other measurement. Let the Kalman filter work on the
% instances where there is data available from both sensors. (double h)
Y = Y(:,1:2:end);
Y(2,:) = Y(2,:)./C;
time = h:h:200;


[X_cv, P] = kalmanFilter(Y, x_0, P_0, A_cv, Q, H, R);

figure();
subplot(2,1,1);
plot(time, X_cv(1,:),'r','linewidth',3);
title('Motion model CV');
grid on;
hold on;
plot(time, Y(1,:),':k','linewidth',2);
legend('x position of train CV','Positional measurements y','location','southeast');
xlabel('time [s]'); ylabel('Position of train [m]');

subplot(2,1,2);
plot(time, X_cv(2,:),'b','linewidth',3);
grid on;
hold on;
plot(time, Y(2,:),':k','linewidth',2);
legend('x Velocity of train CV','Velocity measurements y');
xlabel('time [s]'); ylabel('Velocity of train [m/s]');







% =========================================================================
% Call Kalman filter (CONSTANT ACC)
h = 2*0.1;
Y = Generate_y_seq;
x_0 = [0;0;0];
P_0 = eye(3);

A_ca = [1   h   h*h/2; 
        0   1   h;
        0   0   1];

% Process disturbances
Q = alpha*[0   0   0
     0   0   0
     0   0   1];

H = [1  0  0;
     0  1  0];

% Measurement noise
R = [1   0
    0  mean([var2, var3])];


Y = Y(:,1:2:end);
Y(2,:) = Y(2,:)./C;
time = h:h:200;
[X_ca, P] = kalmanFilter(Y, x_0, P_0, A_ca, Q, H, R);

figure();
subplot(2,1,1);
plot(time, X_ca(1,:),'r','linewidth',3);
title('Motion model CA');
grid on;
hold on;
plot(time, Y(1,:),':k','linewidth',2);
legend('x position of train CA','Positional measurements y','location','southeast');
xlabel('time [s]'); ylabel('Position of train [m]');

subplot(2,1,2);
plot(time, X_ca(2,:),'b','linewidth',3);
grid on;
hold on;
plot(time, Y(2,:),':k','linewidth',2);
legend('x Velocity of train CA','Velocity measurements y');
xlabel('time [s]'); ylabel('Velocity of train [m/s]');

