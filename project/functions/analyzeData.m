function analyzeData(calAcc, calGyr, calMag)
    if size(calAcc,1) == 3, calAcc = calAcc'; end
    if size(calGyr,1) == 3, calGyr = calGyr'; end
    if size(calMag,1) == 3, calMag = calMag'; end

    %% mean and convariance
    fprintf('===== Mean =====\n');
    fprintf('Acc  [m/s^2]: x=%.4f  y=%.4f  z=%.4f\n', mean(calAcc));
    fprintf('Gyr  [rad/s]: x=%.4f  y=%.4f  z=%.4f\n', mean(calGyr));
    fprintf('Mag  [uT]:    x=%.4f  y=%.4f  z=%.4f\n', mean(calMag));

    fprintf('\n===== Covariance =====\n');
    fprintf('Acc:\n'); disp(cov(calAcc));
    fprintf('Gyr:\n'); disp(cov(calGyr));
    fprintf('Mag:\n'); disp(cov(calMag));

    %% figure
    labels   = {'x', 'y', 'z'};
    datasets = {calAcc, calGyr, calMag};
    names    = {'Accelerometer', 'Gyroscope', 'Magnetometer'};
    units    = {'m/s^2', 'rad/s', '\muT'};

    for d = 1:3      % sensor
        for i = 1:3  % roll,pitch,yaw
            figure('Name', [names{d} ' - ' labels{i}]);

            col = datasets{d}(:,i);
            col_clean = col(~isnan(col));

          
            subplot(2,1,1);
            plot(col);
            title([names{d} ' - ' labels{i} ' axis']);
            xlabel('Sample');
            ylabel(units{d});
            grid on;

           
            subplot(2,1,2);
            histogram(col_clean, 30, 'Normalization', 'pdf', ...
                      'FaceColor', [0.3 0.6 0.9], ...
                      'EdgeColor', 'white');
            hold on;

            mu_i    = mean(col_clean);
            std_i   = std(col_clean);
            x_range = linspace(mu_i - 4*std_i, mu_i + 4*std_i, 200);
            y_gauss = (1/(std_i*sqrt(2*pi))) * ...
                      exp(-0.5*((x_range - mu_i)/std_i).^2);
            plot(x_range, y_gauss, 'r-', 'LineWidth', 2);
            hold off;

            title(['Histogram - ' names{d} ' ' labels{i}]);
            xlabel(units{d});
            ylabel('Probability Density');
            legend('Data', 'Gaussian fit');
            grid on;
        end
    end
end