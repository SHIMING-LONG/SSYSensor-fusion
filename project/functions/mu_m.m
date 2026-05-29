function [x, P] = mu_m(x, P, mag, m0, Rm)
%MU_M Magnetometer measurement update
%
% Input:
%   x       [4x1]   State vector (quaternion)
%   P       [4x4]   State covariance
%   mag     [3x1]   Magnetometer measurement
%   m0      [3x1]   Nominal magnetic field in world frame
%   Rm      [3x3]   Measurement noise covariance
%
% Output:
%   x       [4x1]   Updated state
%   P       [4x4]   Updated covariance

    % Make sure m0 is column vector
    mag = mag / norm(mag);
m0  = m0  / norm(m0);
   % m0 = reshape(m0, [3,1]);

    % Predicted measurement
    y_hat = Qq(x)' * m0;

    % Compute Jacobian H = d(Q^T * m0)/dq
    [Q0, Q1, Q2, Q3] = dQqdq(x);
    H = [Q0'*m0, Q1'*m0, Q2'*m0, Q3'*m0];

    % Innovation
    innovation = mag - y_hat;

    % Innovation covariance
    S = H * P * H' + Rm;

    % Kalman gain
   K = P * H' * pinv(S);

    % Update state and covariance
    x = x + K * innovation;
    P = P - K * S * K';

    % Normalize quaternion
    [x, P] = mu_normalizeQ(x, P);

end