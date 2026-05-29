function [t, acc, gyr, mag, dataOrient] = getiOSData(m)
    %GETIOSDATA This functions uses the connected iOS device
    % and gathers data from its sensors at the time instance when the function is
    % called.
   
    acc = m.Acceleration';
    gyr = m.AngularVelocity';
    mag = m.MagneticField';
    [dataOrient, timestamp] = orientlog(m);
   dataOrient = dataOrient(end,:); 
    t = timestamp(end);
end