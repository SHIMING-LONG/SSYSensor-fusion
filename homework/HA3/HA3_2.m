%% HA3 Task 2 - Non-linear Kalman Filtering (Corrected)
% SSY345 - Sensor Fusion and Non-linear Filtering
%
% State: x = [px, py, v, phi, omega]
%   px, py : position (m)
%   v      : speed (m/s)
%   phi    : heading (rad), measured from East, counterclockwise
%   omega  : turn rate (rad/s)
%
% Key fix: true trajectory starts exactly at prior mean (no sampling of x0),
% and headings are chosen so target stays observable for all N=100 steps.

clear; clc; close all;

%% =========================================================
%  Parameters
%% =========================================================

sigma_v  = 0.3;
sigma_om = 0.3*pi/180;
Q = diag([0, 0, sigma_v^2, 0, sigma_om^2]);

sigma_phi = 0.1*pi/180;
R = sigma_phi^2 * eye(2);

s1 = [-50; 0];
s2 = [ 50; 0];

T = 1;       % sampling time [s]
N = 100;     % number of steps

P0 = diag([(5/3)^2, (5/3)^2, 2^2, (0.1*pi/180)^2, (0.1*pi/180)^2]);

% Case 1: [0, 200], v=2, heading East (phi=0)
%   -> moves along y=200, x goes 0..200 over 100 steps
%   -> always between the two sensors' "cone", good geometry
% Case 2: [200, 50], v=2, heading North (phi=pi/2)
%   -> moves from (200,50) upward to (200,250)
%   -> far to the right of both sensors -> poor lateral geometry
%   -> expected: elongated ellipses in x-direction
x0_case1 = [  0; 200; 2;    0; 0];   % heading East
x0_case2 = [200;  50; 2; pi/2; 0];   % heading North

cases      = {x0_case1, x0_case2};
case_names = {'Case 1: x_0=[0,200], heading East (good geometry)', ...
              'Case 2: x_0=[200,50], heading North (poor geometry)'};

%% =========================================================
%  Task 2a/b: Single sequence per case
%% =========================================================

for c = 1:2
    rng(7 + c);   % fixed seed for reproducibility

    x0 = cases{c};

    %--- True state sequence (start exactly at prior mean, add process noise) ---%
    X_true = zeros(5, N+1);
    X_true(:,1) = x0;          % <-- start at mean, NOT a sample from prior
    for k = 2:N+1
        w = [0; 0; sigma_v*randn; 0; sigma_om*randn];   % only v and omega noise
        X_true(:,k) = coordinatedTurn(X_true(:,k-1), T) + w;
    end

    %--- Measurements ---%
    Y_meas = zeros(2, N);
    for k = 1:N
        Y_meas(:,k) = dualBearing(X_true(:,k+1), s1, s2) + ...
                      sigma_phi * randn(2,1);
    end

    %--- EKF ---%
    [xf_ekf, Pf_ekf] = ekf(Y_meas, x0, P0, Q, R, T, s1, s2);

    %--- UKF ---%
    [xf_ukf, Pf_ukf] = ukf(Y_meas, x0, P0, Q, R, T, s1, s2);

    %--- Plot ---%
    figure('Name', sprintf('Case %d', c), 'Position', [50+c*20, 50, 1000, 750]);
    hold on; grid on; axis equal;
    title(sprintf('Non-linear Kalman Filter — %s', case_names{c}), 'FontSize', 11);
    xlabel('x position [m]');
    ylabel('y position [m]');

    % Measurements in Cartesian (intersection of two bearing lines)
    for k = 1:N
        [mx, my] = bearing2cart(Y_meas(1,k), Y_meas(2,k), s1, s2);
        plot(mx, my, '.', 'Color',[0.1 0.7 0.1], 'MarkerSize',5, 'HandleVisibility','off');
    end
    plot(nan,nan,'.','Color',[0.1 0.7 0.1],'MarkerSize',8,'DisplayName','Measurements');

    % True trajectory
    hTrue = plot(X_true(1,2:end), X_true(2,2:end), 'k-', 'LineWidth',2.5, ...
                 'DisplayName','True trajectory');
    plot(X_true(1,2), X_true(2,2), 'ks','MarkerSize',10,'MarkerFaceColor','k', ...
         'HandleVisibility','off');

    % EKF
    plot(xf_ekf(1,:), xf_ekf(2,:), 'b--', 'LineWidth',1.8, 'DisplayName','EKF estimate');
    for k = 10:10:N
        ellipse3sigma(xf_ekf(1:2,k), Pf_ekf(1:2,1:2,k), [0 0 0.85]);
    end
    patch(nan,nan,[0 0 0.85],'FaceAlpha',0,'EdgeColor',[0 0 0.85], ...
          'LineWidth',1.2,'DisplayName','EKF 3\sigma ellipse');

    % UKF
    plot(xf_ukf(1,:), xf_ukf(2,:), 'm-.', 'LineWidth',1.8, 'DisplayName','UKF estimate');
    for k = 10:10:N
        ellipse3sigma(xf_ukf(1:2,k), Pf_ukf(1:2,1:2,k), [0.75 0 0.75]);
    end
    patch(nan,nan,[0.75 0 0.75],'FaceAlpha',0,'EdgeColor',[0.75 0 0.75], ...
          'LineWidth',1.2,'DisplayName','UKF 3\sigma ellipse');

    % Sensors
    plot(s1(1), s1(2), '^', 'Color',[0 0 0.9], 'MarkerSize',14, ...
         'MarkerFaceColor',[0 0 0.9], 'DisplayName','Sensor 1');
    plot(s2(1), s2(2), '^', 'Color',[0.9 0 0], 'MarkerSize',14, ...
         'MarkerFaceColor',[0.9 0 0], 'DisplayName','Sensor 2');

    legend('Location','best','FontSize',9);

    % Add annotation about geometry
    if c == 1
        text(0.02, 0.05, 'Good geometry: target between sensors', ...
             'Units','normalized','FontSize',9,'Color',[0 0.5 0]);
    else
        text(0.02, 0.05, 'Poor geometry: target far to the side -> elongated ellipses', ...
             'Units','normalized','FontSize',9,'Color',[0.7 0 0]);
    end

    drawnow;
    fprintf('Case %d: true x range [%.1f, %.1f], y range [%.1f, %.1f]\n', c, ...
        min(X_true(1,:)), max(X_true(1,:)), min(X_true(2,:)), max(X_true(2,:)));
end

%% =========================================================
%  Task 2c: Monte Carlo (100 sequences each case)
%% =========================================================

MC = 100;

for c = 1:2
    x0 = cases{c};
    err_ekf = zeros(2, N*MC);
    err_ukf = zeros(2, N*MC);

    fprintf('\nMonte Carlo case %d ...\n', c);
    for imc = 1:MC
        X_true = zeros(5, N+1);
        X_true(:,1) = x0;
        for k = 2:N+1
            w = [0;0;sigma_v*randn;0;sigma_om*randn];
            X_true(:,k) = coordinatedTurn(X_true(:,k-1), T) + w;
        end
        Y_meas = zeros(2,N);
        for k = 1:N
            Y_meas(:,k) = dualBearing(X_true(:,k+1),s1,s2) + sigma_phi*randn(2,1);
        end

        [xfe,~] = ekf(Y_meas, x0, P0, Q, R, T, s1, s2);
        [xfu,~] = ukf(Y_meas, x0, P0, Q, R, T, s1, s2);

        idx = (imc-1)*N + (1:N);
        err_ekf(:,idx) = xfe(1:2,:) - X_true(1:2,2:end);
        err_ukf(:,idx) = xfu(1:2,:) - X_true(1:2,2:end);
    end

    figure('Name',sprintf('MC Case %d',c),'Position',[150 150 1000 520]);
    sgtitle(sprintf('MC Error Histograms — %s',case_names{c}),'FontSize',10);

    fdata  = {err_ekf, err_ukf};
    flabel = {'EKF','UKF'};
    dlabel = {'x error [m]','y error [m]'};
    fcolor = {[0.2 0.4 0.85],[0.7 0.1 0.7]};

    for fi = 1:2
        for di = 1:2
            subplot(2,2,(fi-1)*2+di);
            e = fdata{fi}(di,:);
            e = e(isfinite(e));               % drop any NaN/Inf
            e = e(abs(e) < 5*std(e));         % clip wild outliers
            histogram(e, 60, 'Normalization','pdf', ...
                      'FaceColor',fcolor{fi},'EdgeColor','none','FaceAlpha',0.75);
            hold on;
            mu_e=mean(e); s_e=std(e);
            xg=linspace(mu_e-4*s_e,mu_e+4*s_e,300);
            plot(xg,normpdf(xg,mu_e,s_e),'r-','LineWidth',2);
            xlabel(dlabel{di}); ylabel('pdf');
            title(sprintf('%s — %s',flabel{fi},dlabel{di}));
            grid on;
        end
    end
end
fprintf('\nDone.\n');

%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

function xn = coordinatedTurn(x, T)
    px=x(1);py=x(2);v=x(3);phi=x(4);om=x(5);
    if abs(om)<1e-10
        xn=[px+T*v*cos(phi); py+T*v*sin(phi); v; phi; om];
    else
        xn=[px+(v/om)*(sin(phi+om*T)-sin(phi));
            py+(v/om)*(-cos(phi+om*T)+cos(phi));
            v; phi+om*T; om];
    end
end

function y = dualBearing(x, s1, s2)
    y=[atan2(x(2)-s1(2), x(1)-s1(1));
       atan2(x(2)-s2(2), x(1)-s2(1))];
end

function [xc,yc] = bearing2cart(y1,y2,s1,s2)
    t1=tan(y1); t2=tan(y2);
    xc=(s2(2)-s1(2)+t1*s1(1)-t2*s2(1))/(t1-t2);
    yc=s1(2)+t1*(xc-s1(1));
end

%--- EKF ---%
function [xf,Pf] = ekf(Y,x0,P0,Q,R,T,s1,s2)
    N=size(Y,2); nx=5;
    xf=zeros(nx,N); Pf=zeros(nx,nx,N);
    xp=x0; Pp=P0;
    for k=1:N
        H=Hjac(xp,s1,s2);
        yp=dualBearing(xp,s1,s2);
        S=H*Pp*H'+R; K=Pp*H'/S;
        xu=xp+K*(Y(:,k)-yp);
        Pu=(eye(nx)-K*H)*Pp;
        xf(:,k)=xu; Pf(:,:,k)=Pu;
        F=Fjac(xu,T);
        xp=coordinatedTurn(xu,T);
        Pp=F*Pu*F'+Q;
    end
end

function F=Fjac(x,T)
    v=x(3);phi=x(4);om=x(5);
    F=eye(5);
    if abs(om)<1e-10
        F(1,3)=T*cos(phi);    F(1,4)=-T*v*sin(phi);
        F(2,3)=T*sin(phi);    F(2,4)= T*v*cos(phi);
        F(4,5)=T;
    else
        F(1,3)=(sin(phi+om*T)-sin(phi))/om;
        F(1,4)=v/om*(cos(phi+om*T)-cos(phi));
        F(1,5)=v/om^2*(om*T*cos(phi+om*T)-sin(phi+om*T)+sin(phi));
        F(2,3)=(-cos(phi+om*T)+cos(phi))/om;
        F(2,4)=v/om*(sin(phi+om*T)-sin(phi));
        F(2,5)=v/om^2*(om*T*sin(phi+om*T)+cos(phi+om*T)-cos(phi));
        F(4,5)=T;
    end
end

function H=Hjac(x,s1,s2)
    px=x(1);py=x(2);
    d1=(px-s1(1))^2+(py-s1(2))^2;
    d2=(px-s2(1))^2+(py-s2(2))^2;
    H=zeros(2,5);
    H(1,1)=-(py-s1(2))/d1; H(1,2)=(px-s1(1))/d1;
    H(2,1)=-(py-s2(2))/d2; H(2,2)=(px-s2(1))/d2;
end

%--- UKF ---%
function [xf,Pf] = ukf(Y,x0,P0,Q,R,T,s1,s2)
    N=size(Y,2); nx=5;
    alpha=1e-3;beta=2;kappa=0;
    lam=alpha^2*(nx+kappa)-nx;
    Wm=[lam/(nx+lam), repmat(0.5/(nx+lam),1,2*nx)];
    Wc=Wm; Wc(1)=Wc(1)+(1-alpha^2+beta);
    xf=zeros(nx,N); Pf=zeros(nx,nx,N);
    xp=x0; Pp=P0;
    for k=1:N
        sig=spSigma(xp,Pp,nx,lam);
        Ys=zeros(2,2*nx+1);
        for j=1:2*nx+1; Ys(:,j)=dualBearing(sig(:,j),s1,s2); end
        ym=Ys*Wm';
        S=R; Pxy=zeros(nx,2);
        for j=1:2*nx+1
            dy=Ys(:,j)-ym;
            S=S+Wc(j)*(dy*dy');
            Pxy=Pxy+Wc(j)*(sig(:,j)-xp)*dy';
        end
        K=Pxy/S;
        xu=xp+K*(Y(:,k)-ym);
        Pu=Pp-K*S*K';
        xf(:,k)=xu; Pf(:,:,k)=Pu;
        sig=spSigma(xu,Pu,nx,lam);
        xp=zeros(nx,1);
        for j=1:2*nx+1; xp=xp+Wm(j)*coordinatedTurn(sig(:,j),T); end
        Pp=Q;
        for j=1:2*nx+1
            d=coordinatedTurn(sig(:,j),T)-xp;
            Pp=Pp+Wc(j)*(d*d');
        end
    end
end

function sig=spSigma(x,P,nx,lam)
    A=chol((nx+lam)*P,'lower');
    sig=[x, x*ones(1,nx)+A, x*ones(1,nx)-A];
end

%--- 3-sigma ellipse ---%
function ellipse3sigma(mu,P,rgb)
    [V,D]=eig(P);
    t=linspace(0,2*pi,120);
    e=3*V*sqrt(abs(D))*[cos(t);sin(t)];
    plot(mu(1)+e(1,:),mu(2)+e(2,:),'-','Color',rgb,'LineWidth',1.0, ...
         'HandleVisibility','off');
end