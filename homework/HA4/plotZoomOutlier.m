function plotZoomOutlier(X, xs, Ps, xf, Pf, xpos, ypos, s1, s2, type, k_outlier)
%PLOTZOOMOUTLIER
%   Zoomed-in plot around the outlier measurement at k = k_outlier
%   Shows filter and smoother trajectories with 3-sigma ellipses

%% 颜色定义
c_true   = [0.10, 0.10, 0.10];
c_meas   = [0.00, 0.80, 0.80];   % 青色虚线，与你图中一致
c_filt   = [0.93, 0.69, 0.13];   % 黄橙 — 滤波
c_smooth = [0.49, 0.18, 0.56];   % 紫色 — 平滑
c_ell_f  = [0.10, 0.10, 0.10];   % 黑色椭圆 — 滤波
c_ell_s  = [0.65, 0.10, 0.20];   % 暗红椭圆 — 平滑

%% 确定缩放范围（以 k_outlier 对应的真实位置为中心）
% 取 outlier 前后各 ~60 步的 x 范围
k_range  = max(1, k_outlier-60) : min(size(xf,2), k_outlier+60);
x_center = X(1, k_outlier+1);    % +1 因为 X 从 k=0 开始
zoom_hw  = 30;                    % 水平半宽 [m]，可调
x_lim    = [x_center - zoom_hw,  x_center + zoom_hw];
% y 范围自动从数据决定，留一点余量
y_data   = [xf(2, k_range), xs(2, k_range), ypos(k_range)];
y_lim    = [min(y_data)-5, max(y_data)+5];

%% 绘图
figure('Position', [200, 200, 1100, 600]);
ax = gca;
hold(ax, 'on'); grid(ax, 'on');
ax.FontSize = 13;

% 真实轨迹
h_true   = plot(ax, X(1, k_range+1), X(2, k_range+1), '-', ...
                'Color', c_true, 'LineWidth', 2.5);
% 测量（Cartesian 投影，即 xpos/ypos）
h_meas   = plot(ax, xpos(k_range), ypos(k_range), ':', ...
                'Color', c_meas, 'LineWidth', 2.0, 'MarkerSize', 4);
% 滤波轨迹
h_filt   = plot(ax, xf(1, k_range), xf(2, k_range), '-', ...
                'Color', c_filt, 'LineWidth', 2.2);
% 平滑轨迹
h_smooth = plot(ax, xs(1, k_range), xs(2, k_range), '-', ...
                'Color', c_smooth, 'LineWidth', 2.2);

% 3σ 椭圆（每 10 步一个）
for i = k_range(1) : 10 : k_range(end)
    if i > size(xf,2); break; end
    xy_f = sigmaEllipse2D([xf(1,i); xf(2,i)], Pf(1:2,1:2,i), 3, 100);
    h_ell_f = plot(ax, xy_f(1,:), xy_f(2,:), '-', ...
                   'Color', c_ell_f, 'LineWidth', 0.9);
    xy_s = sigmaEllipse2D([xs(1,i); xs(2,i)], Ps(1:2,1:2,i), 3, 100);
    h_ell_s = plot(ax, xy_s(1,:), xy_s(2,:), '-', ...
                   'Color', c_ell_s, 'LineWidth', 0.9);
end

% 标记 outlier 测量点
h_out = scatter(ax, xpos(k_outlier), ypos(k_outlier), 150, ...
                'r', 'filled', 'pentagram');

% 坐标轴范围
xlim(ax, x_lim);
ylim(ax, y_lim);

legend(ax, [h_true, h_meas, h_filt, h_smooth, h_ell_f, h_ell_s, h_out], ...
    {'True position', 'Measured position', 'Filtered position', ...
     'Smoothed position', '3\sigma filter', '3\sigma smoother', 'Outlier'}, ...
    'Location', 'best', 'FontSize', 12);

title(ax, ['Zoomed view around outlier (k=', num2str(k_outlier), ') — ', type], ...
      'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax, 'x  [m]', 'FontSize', 13);
ylabel(ax, 'y  [m]', 'FontSize', 13);

end