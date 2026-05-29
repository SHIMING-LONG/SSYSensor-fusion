function [x, P] = mu_g(x, P, yacc, Ra, g0)
    g0 = reshape(g0, [3,1]);

    % Predicted measurement
    y_hat = Qq(x)' * g0;

    % Jacobian H [3x4]
    [Q0, Q1, Q2, Q3] = dQqdq(x);
    H = [Q0'*g0, Q1'*g0, Q2'*g0, Q3'*g0];

    % Innovation
    innovation = yacc - y_hat;

    % Innovation covariance
    S = H * P * H' + Ra;

    % Kalman gain
    K = P * H' * pinv(S);

    % Update
    x = x + K * innovation;
    P = P - K*S*K';

    [x, P] = mu_normalizeQ(x, P);
end