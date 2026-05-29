function [x, P] = tu_qw(x, P, omega, T, Rw)
%
% Input:
%   x       [4x1]   Current state (quaternion)
%   P       [4x4]   Current covariance
%   omega   [3x1]   Angular rate measurement [rad/s]
%   T               Time step [s]
%   Rw      [3x3]   Process noise covariance
%
% Output:
%   x       [4x1]   Predicted state
%   P       [4x4]   Predicted covarian
% F matrix from Task 3
    F = eye(4) + (T/2) * Somega(omega);
    x = F * x;  % Predict the new state
    G = (T/2) * Sq(x);
    P = F * P * F' + G*Rw*G';  % Update the covariance
    [x, P] = mu_normalizeQ(x, P);
end
