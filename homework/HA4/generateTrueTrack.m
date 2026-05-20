function [X,T,Tvec] = generateTrueTrack()
% True track
% Sampling period
T = 0.1;

% Length of time sequence
K = 800;

% Allocate memory
omega = zeros(1,K+1);

% Set turn-rate at turns
omega(200:350) = pi/301/T;
omega(450:600) = pi/301/T;

% Initial state
x0 = [0 0 20 -pi/2 0]';

% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;

% Create true track
for i=2:K+1
    % Simulate
    X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
    % Set turn−rate
    X(5,i) = omega(i);
end

Tvec = 0:T:(length(X)-1)*T;

end
