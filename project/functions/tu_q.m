function [x, P] = tu_q(x, P, T, Rw)
%TU_Q Time update without gyroscope measurement
%
% Input:
%   x       [4x1]   Current state (quaternion)
%   P       [4x4]   Current covariance
%   T               Time step [s]
%   Rw      [3x3]   Process noise covariance
%
% Output:
%   x       [4x1]   Predicted state
%   P       [4x4]   Predicted covariance

    % No angular rate, set omega = 0
    omega = zeros(3, 1);

    % F reduces to identity when omega = 0
    F = eye(4) + (T/2) * Somega(omega);  % = eye(4)

    % G matrix
    G = (T/2) * Sq(x);

    % Predict state (x unchanged since F = I)
    x = F * x;

    % Predict covariance (uncertainty grows over time)
    P = F * P * F' + G * Rw * G';

    % Normalize quaternion
    [x, P] = mu_normalizeQ(x, P);

end