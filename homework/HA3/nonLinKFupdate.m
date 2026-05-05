function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%
n = length(x);
    switch type
        case 'EKF'
            
            % Your EKF update here
            [hx, Hx] = h(x);
            vk = y - hx;
            Sk = Hx * P * Hx' + R;
            K = P * Hx' / Sk;
            x = x + K * vk;
            P = P - K * Sk * K';
            
        case 'UKF'
    
            % Your UKF update here
            [SP, W] = sigmaPoints(x, P, 'UKF'); %[cite: 1]
            num_sp = size(SP, 2);
            m = size(y, 1);
            Y = zeros(m, num_sp);
            y_hat = zeros(m, 1);
            for i = 1:num_sp
                [Y(:,i), ~] = h(SP(:,i));
                y_hat = y_hat + W(i) * Y(:,i); 
            end
            Pyy = R;
            Pxy = zeros(n, m);
            for i = 1:num_sp
                y_diff = Y(:,i) - y_hat;
                x_diff = SP(:,i) - x;
                Pyy = Pyy + W(i) * (y_diff * y_diff'); 
                Pxy = Pxy + W(i) * (x_diff * y_diff'); 
            end
            K = Pxy / Pyy; 
            x = x + K * (y - y_hat);
            P = P - K * Pyy * K'; 
            if min(eig(P)) <= 0
                [v, e] = eig(P, 'vector');
                e(e < 0) = 1e-4; 
                P = v * diag(e) / v; 
            end
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
            
        case 'CKF'
    
            % Your CKF update here
            [SP, W] = sigmaPoints(x, P, 'CKF'); 
            num_sp = size(SP, 2);
            m = size(y, 1);
            Y = zeros(m, num_sp);
            y_hat = zeros(m, 1);
            for i = 1:num_sp
                [Y(:,i), ~] = h(SP(:,i));
                y_hat = y_hat + W(i) * Y(:,i); 
            end
            Pyy = R;
            Pxy = zeros(n, m);
            for i = 1:num_sp
                y_diff = Y(:,i) - y_hat;
                x_diff = SP(:,i) - x;
                Pyy = Pyy + W(i) * (y_diff * y_diff'); 
                Pxy = Pxy + W(i) * (x_diff * y_diff');
            end
            K = Pxy / Pyy; 
            x = x + K * (y - y_hat); 
            P = P - K * Pyy * K'; 
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

