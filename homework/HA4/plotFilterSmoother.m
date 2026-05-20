function plotFilterSmoother(X,xs, Ps, xf, Pf, xp, Pp, xpos, ypos, s1, s2, type, Tvec)
    %PLOTFILTERSMOOTHER 
    %   Helper function for plotting of the filter coordinated turn motion model

  
    c_true  = '#03396C'; 
    c_meas  = '#FF7F50'; 
    c_filt  = '#2E8B57';
    c_smth  = '#DC143C'; 
    c_cam   = '#555555'; 
    c_ell_f = '#8FBC8F'; 
    c_ell_s = '#FFB6C1'; 
    % --------------------

    %% 1. 二维轨迹对比图
    figure()
    grid on; hold on
    axis equal; % 确保 x 轴和 y 轴比例一致

    % 保存绘图句柄，以便准确生成图例
    h_true = plot(X(1,:), X(2,:), 'Color', c_true, 'LineWidth', 2);
    h_meas = plot(xpos, ypos, ':', 'color', c_meas, 'LineWidth', 2);
    h_filt = plot(xf(1,:), xf(2,:), 'Color', c_filt, 'LineWidth', 2);
    h_smooth = plot(xs(1,:), xs(2,:), 'Color', c_smth, 'LineWidth', 2);
    
    % 绘制传感器位置 (Camera/Sensor)
    h_cam1 = scatter(s1(1), s1(2), 50, 'MarkerFaceColor', c_cam, 'MarkerEdgeColor', c_cam, 'Marker', 'square');
    h_cam2 = scatter(s2(1), s2(2), 50, 'MarkerFaceColor', c_cam, 'MarkerEdgeColor', c_cam, 'Marker', 'square');

    % 绘制 Filter 的 3-sigma 椭圆 (要求每 10 个时刻画一次)
    for i = 10:10:size(xf, 2)
        P_i_f = Pf(1:2,1:2,i);
        [xy_f] = sigmaEllipse2D([xf(1,i); xf(2,i)], P_i_f, 3, 100);
        h_ell_f = plot(xy_f(1,:), xy_f(2,:), 'Color', c_ell_f); 
    end

    % 绘制 Smoother 的 3-sigma 椭圆 (要求每 10 个时刻画一次)
    for i = 10:10:size(xs, 2)
        P_i_s = Ps(1:2,1:2,i);
        [xy_s] = sigmaEllipse2D([xs(1,i); xs(2,i)], P_i_s, 3, 100);
        h_ell_s = plot(xy_s(1,:), xy_s(2,:), 'Color', c_ell_s); 
    end
    
    % 使用句柄精确设置图例
    legend([h_true, h_meas, h_filt, h_smooth, h_cam1, h_cam2, h_ell_f, h_ell_s], ...
        {'True position', 'Measured position', 'Filtered position', ...
         'Smoothed position', 'Sensor 1 pos', 'Sensor 2 pos', ...
         '3\sigma-contour filter', '3\sigma-contour smoother'}, ...
         'Location', 'best')
         
    title(['Coordinated turn motion model filtered with ', type])
    xlabel('x pos [m]'); ylabel('y pos [m]')
    
    %% 2. 状态分量对比图
    figure()
    
    % 速度 (Velocity)
    subplot(1,3,1)
    grid on; hold on
    plot(Tvec, X(3,:), 'Color', c_true, 'LineWidth', 2)
    plot(Tvec(2:end), xf(3,:), 'Color', c_filt)
    plot(Tvec(2:end), xs(3,:), 'Color', c_smth, 'LineWidth', 2)
    ylabel('Velocity [m/s]'); xlabel('time [s]')
    legend('True state', 'Filtered state', 'Smoothed state', 'Location', 'best')
    
    % 航向角 (Heading)
    subplot(1,3,2)
    grid on; hold on
    plot(Tvec, X(4,:), 'Color', c_true, 'LineWidth', 2)
    plot(Tvec(2:end), xf(4,:), 'Color', c_filt)
    plot(Tvec(2:end), xs(4,:), 'Color', c_smth, 'LineWidth', 2)
    ylabel('Turn (heading) [rad]'); xlabel('time [s]')
    legend('True state', 'Filtered state', 'Smoothed state', 'Location', 'best')
    
    % 偏航率 (Turn rate)
    subplot(1,3,3)
    grid on; hold on
    plot(Tvec, X(5,:), 'Color', c_true, 'LineWidth', 2)
    plot(Tvec(2:end), xf(5,:), 'Color', c_filt)
    plot(Tvec(2:end), xs(5,:), 'Color', c_smth, 'LineWidth', 2)
    ylabel('Turn rate [rad/s]'); xlabel('time [s]')
    legend('True state', 'Filtered state', 'Smoothed state', 'Location', 'best')

end