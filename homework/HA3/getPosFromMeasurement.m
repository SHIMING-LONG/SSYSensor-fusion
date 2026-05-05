function [x, y] = getPosFromMeasurement(y1, y2, s1, s2)
% GETPOSFROMMEASUREMENT 通过两个传感器的位置和方位角测量值计算交点坐标（三角定位）
%
% INPUT: 
%   y1: 传感器 1 的方位角测量值
%   y2: 传感器 2 的方位角测量值
%   s1: 传感器 1 的 2D 坐标
%   s2: 传感器 2 的 2D 坐标
%
% OUTPUT:
%   x, y: 交点的横纵坐标

% 根据公式求解两个线性方程的交点：
% (y - s1(2)) = (x - s1(1)) * tan(y1)
% (y - s2(2)) = (x - s2(1)) * tan(y2)

x = (s2(2) - s1(2) + tan(y1).*s1(1) - tan(y2).*s2(1)) ./ (tan(y1) - tan(y2));
y = s1(2) + tan(y1) .* (x - s1(1));

end