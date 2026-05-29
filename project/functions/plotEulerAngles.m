function plotEulerAngles(xhat, meas)
%PLOTEULERANGLES Plot Euler angles comparison

    % Convert quaternions to Euler angles
    euler_own    = q2euler(xhat.x(1:4,:));
    euler_google = q2euler(meas.orient(1:4,:));

    figure('Name', 'Euler Angles Comparison');

    % YAW
    subplot(3,1,1);
    plot(xhat.t, euler_own(1,:) * 180/pi, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(xhat.t, euler_google(1,:) * 180/pi, 'k--', 'LineWidth', 1.5);
    title('YAW (z-axis)');
    xlabel('Time [s]');
    ylabel('Angle [deg]');
    legend('Our estimate', 'MATLAB estimate');
    grid on;
    hold off;

    % ROLL
    subplot(3,1,2);
    plot(xhat.t, euler_own(3,:) * 180/pi, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(xhat.t, euler_google(3,:) * 180/pi, 'k--', 'LineWidth', 1.5);
    title('ROLL (y-axis)');
    xlabel('Time [s]');
    ylabel('Angle [deg]');
    legend('Our estimate', 'MATLAB estimate');
    grid on;
    hold off;

    % PITCH
    subplot(3,1,3);
    plot(xhat.t, euler_own(2,:) * 180/pi, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(xhat.t, euler_google(2,:) * 180/pi, 'k--', 'LineWidth', 1.5);
    title('PITCH (x-axis)');
    xlabel('Time [s]');
    ylabel('Angle [deg]');
    legend('Our estimate', 'MATLAB estimate');
    grid on;
    hold off;

end