function [y_mean, P_y, y_s] = approxGaussianTransform(x, P, h, R, N)
% APPROXGAUSSIANTRANSFORM 通过蒙特卡罗采样近似非线性变换后的均值和协方差
%
% 输入:
%   x - [n x 1] 先验均值
%   P - [n x n] 先验协方差
%   h - 测量模型函数句柄
%   R - [m x m] 测量噪声协方差
%   N - 采样数量 (建议使用 10000)

    % 1. 从先验分布 p(x) 中抽取 N 个样本
    % 使用 mvnrnd (Multivariate Normal Random numbers)
    x_samples = mvnrnd(x, P, N)'; % 得到 [n x N] 的矩阵

    % 2. 将每个样本通过非线性函数 h(x)，并加上噪声 r
    % 这里的 h 是双轴承模型
    m = size(R, 1);
    y_s = zeros(m, N);
    
    % 生成测量噪声样本 r ~ N(0, R)
    r_samples = mvnrnd(zeros(m, 1), R, N)';

    for i = 1:N
        % 传播样本并加噪
        y_s(:, i) = h(x_samples(:, i)) + r_samples(:, i);
    end

    % 3. 计算样本均值 (Sample Mean)
    y_mean = mean(y_s, 2);

    % 4. 计算样本协方差 (Sample Covariance)
    % 注意：可以使用内置的 cov 函数，但要注意维度
    P_y = cov(y_s'); 
end