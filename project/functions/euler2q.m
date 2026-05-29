function q = euler2q(roll,pitch,yaw)
    
    roll = deg2rad(roll);
    pitch = deg2rad(pitch);
    yaw = deg2rad(yaw);

    % Compute the corresponding orientation as quaternion
    q = [cos(roll/2)*cos(pitch/2)*cos(yaw/2)  +  sin(roll/2)*sin(pitch/2)*sin(yaw/2);
         sin(roll/2)*cos(pitch/2)*cos(yaw/2)  -  cos(roll/2)*sin(pitch/2)*sin(yaw/2);
         cos(roll/2)*sin(pitch/2)*cos(yaw/2)  +  sin(roll/2)*cos(pitch/2)*sin(yaw/2);
         cos(roll/2)*cos(pitch/2)*sin(yaw/2)  -  sin(roll/2)*sin(pitch/2)*cos(yaw/2)];


    % Normalize the quaternion
   [q, ~] = mu_normalizeQ(q, q);

    % return
    
end
