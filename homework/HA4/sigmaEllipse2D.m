function [xy] = sigmaEllipse2D(mu, Sigma, level, npoints)
% SIGMAELLIPSE2D 生成二维协方差椭圆的坐标点
%
% 输入:
%   mu      - [2 x 1] 均值向量
%   Sigma   - [2 x 2] 协方差矩阵
%   level   - 标准差倍数 (例如 3 代表 3-sigma)
%   npoints - 构成椭圆的点数

    % 生成单位圆上的点
    phi = linspace(0, 2*pi, npoints);
    unit_circle = [cos(phi); sin(phi)];

    % 通过 Cholesky 分解将单位圆缩放并旋转以匹配协方差矩阵
    % 如果 Sigma 不是严格正定，使用 eig 方法更稳健
    [V, D] = eig(Sigma);
    % 变换矩阵 = 旋转(特征向量) * 缩放(特征值的平方根) * 倍数
    A = V * sqrt(D) * level;

    % 变换点并平移到均值位置
    xy = A * unit_circle + mu;
end