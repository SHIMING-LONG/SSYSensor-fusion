function meas = collectData(m, duration)

    fprintf('Recording for %d seconds...\n', duration);
    fprintf('Keep the phone still on the table.\n');

    % Start logging
    m.Logging = 1;
    pause(duration);
    m.Logging = 0;

    fprintf('Recording complete! Processing data...\n');

    % Read logged data
    [acc,    t_acc]    = accellog(m);
    [gyr,    t_gyr]    = angvellog(m);
    [mag,    t_mag]    = magfieldlog(m);
    [orient, t_orient] = orientlog(m);

    % Use accelerometer time as common reference
    t0       = t_acc(1);
    t_common = t_acc - t0;

    % Interpolate all sensors onto common time axis
    gyr_interp    = interp1(t_gyr    - t0, gyr,    t_common, 'linear', 'extrap');
    mag_interp    = interp1(t_mag    - t0, mag,    t_common, 'linear', 'extrap');
    orient_interp = interp1(t_orient - t0, orient, t_common, 'linear', 'extrap');

    % Pack into meas struct (3×N, same format as filterTemplate)
    meas.t      = t_common';
    meas.acc    = acc';
    meas.gyr    = gyr_interp';
    meas.mag    = mag_interp';
    meas.orient = orient_interp';

    fprintf('Done. Collected %d samples at %.1f Hz.\n', ...
            length(t_common), length(t_common)/duration);
end