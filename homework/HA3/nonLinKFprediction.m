function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
      n = length(x);

    switch type
        case 'EKF'
            
      
           [fx,Fx]=f(x); 
              
         
            x = fx;
            P = Fx*P*Fx' + Q;
        case 'UKF'
    
            % Your UKF code here
              % 1. Form a set of 2n+1 sigma points
            [SP, W] = sigmaPoints(x, P, 'UKF');
            num_sp = size(SP, 2);
            
            % 2. Compute the predicted moments
            SP_pred = zeros(n, num_sp);
            
            % Iterate through and sum
            x_pred = zeros(n, 1);
            for i = 1:num_sp
                
                fx = f(SP(:,i)); 
                SP_pred(:,i) = fx;
                x_pred = x_pred + W(i) * fx;
            end
            P_pred = Q; % 初始值设为过程噪声 Q
            for i = 1:num_sp
                % 注意：这里必须减去刚刚算出来的 x_pred
                diff = SP_pred(:,i) - x_pred; 
                P_pred = P_pred + W(i) * (diff * diff');
            end
            
            x = x_pred;
            P = P_pred;
            
           
            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
                
        case 'CKF'
            
            % Your CKF code here
             % 1. Form a set of 2n sigma points
            [SP, W] = sigmaPoints(x, P, 'CKF');
            num_points = size(SP, 2);;
            SP_transformed = zeros(n, num_points);
            % 2. Compute the predicted moments
           x_pred = zeros(n, 1);
            for i = 1:num_points
                fx = f(SP(:,i)); 
                SP_transformed(:,i) = fx;
                x_pred = x_pred + fx * W(i);
            end
            P_pred = Q; 
            for i = 1:num_points
                deviation = SP_transformed(:,i) - x_pred; 
                P_pred = P_pred + W(i) * (deviation * deviation');
            end
            % 3. Return prediction mean and covariance
            x = x_pred;
            P = P_pred;
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end