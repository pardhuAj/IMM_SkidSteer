clc; clear; clf;
phi = [0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 0];
Gamma = [0 0;
    1 0;
    0 0;
    0 1];
dt = 1;
T = 600;
t = 0:dt:T;
u1 = zeros(1,length(t));
for i = 1:length(t)
    if t(i) < 400
        u1(i) = 0;
    else 
        u1(i) = 0.075;
    end
end

u2 = u1;
U = [u1; u2];
% x(k+1) = (I+phi*dt)*x(k) + Gamma*dt*u(k) States are [x xdot y ydot]'

phi = eye(4) +dt*phi;
Gamma = dt*Gamma;
x = zeros(4,length(t));
x(:,1) = [2000 0 10000 -15]';
for i = 1:length(t)-1
    x(:,i+1) = phi*x(:,i) + Gamma*U(:,i);
end

R1 = 1e4;
R2 = 1e4;
rng(0,'twister')
v1 = sqrt(R1)*randn(1,length(t));
v2 = sqrt(R2)*randn(1,length(t));

H = [1 0 0 0;
    0 0 1 0];

y = H*x+ [v1; v2];
figure(1)
plot(t, y(1,:),t,y(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Distance')
legend('X position', 'Y position')
grid minor


%% filter model 1 -- no process noise

phi1 = eye(6) + dt*[0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0];
Gamma1 =dt*[0 0;
    1 0;
    0 0;
    0 0;
    0 1;
    0 0];
upsilon = dt*[0 0;
    0 0;
    1 0;
    0 0;
    0 0;
    0 1];
H1 = [1 0 0 0 0 0;
    0 0 0 1 0 0];
w1 = 0*randn(1,length(t));
w2 = 0*randn(1,length(t));
W1 = [w1; w2];
R = diag([R1, R2]);
x1 = zeros(6, length(t));
x1(:,1) = [2000 0 0 10000 -15 0]';
for i = 1:length(t)-1
    x1(:,i+1) = phi1*x1(:,i) + upsilon*W1(:,i) + Gamma1*U(:,i);
end
figure(2)
plot(t, x1(1,:),t,x1(4,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Distance')
legend('X position', 'Y position')
grid minor

% Discrete kalman filter -- no process covariance

x1plus = zeros(6,length(t));
x1minus = zeros(6,length(t));
x1minus(:,1) = [2000 0 0 10000 -15 0]';
P1minus = zeros(6,6,length(t));
P1plus = zeros(6,6,length(t));
P1minus(:,:,1) = 1e-12*eye(6);
E1 = zeros(2,2,length(t));
for i = 1:length(t)-1
    % Kalman gain 
    W1 = P1minus(:,:,i)*H1'/(H1*P1minus(:,:,i)*H1'+R);
    % State update and covariance update
    x1plus(:,i) = x1minus(:,i) + W1*(y(:,i)-H1*x1minus(:,i));
    P1plus(:,:,i) = (eye(6)-W1*H1)*P1minus(:,:,i);
    % State and covariance propagation
    x1minus(:,i+1) = phi1*x1plus(:,i) + 0*Gamma1*U(:,i);
    P1minus(:,:,i+1)=phi1*P1plus(:,:,i)*phi1' ;
    E1(:,:,i) = H1*P1minus(:,:,i)*H1'+R;
end
W1 = P1minus(:,:,end)*H1'/(H1*P1minus(:,:,end)*H1'+R);
x1plus(:,end) = x1minus(:,end) + W1*(y(:,end)-H1*x1minus(:,end));
P1plus(:,:,end) = (eye(6)-W1*H1)*P1minus(:,:,end);
E1(:,:,end) = H1*P1minus(:,:,end)*H1'+R;
meas_res1 = y - H1*x1minus;


%% filter model 2 -- Process noise
phi1 = eye(6) + dt*[0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0];
upsilon = dt*[0 0;
    0 0;
    1 0;
    0 0;
    0 0;
    0 1];
Q = 1e-3*eye(2);
w3 = sqrt(1e-3)*randn(1,length(t));
w4 = sqrt(1e-3)*randn(1,length(t));
W2 = [w3; w4];
x2 = zeros(6, length(t));
x2(:,1) = [2000 0 0 10000 -15 0]';
for i = 1:length(t)-1
    x2(:,i+1) = phi1*x2(:,i) + upsilon*W2(:,i)+ 0*Gamma1*U(:,i);
end

figure(3)
plot(t, x2(1,:),t,x2(4,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Distance')
legend('X position', 'Y position')
grid minor


% Discrete kalman filter -- no process covariance

x2plus = zeros(6,length(t));
x2minus = zeros(6,length(t));
x2minus(:,1) = [2000 0 0 10000 -15 0]';
P2minus = zeros(6,6,length(t));
P2plus = zeros(6,6,length(t));
P2minus(:,:,1) = 1e-12*eye(6);
E2 = zeros(2,2,length(t));
for i = 1:length(t)-1
    % Kalman gain 
    W2 = P2minus(:,:,i)*H1'/(H1*P2minus(:,:,i)*H1'+R);
    % State update and covariance update
    x2plus(:,i) = x2minus(:,i) + W2*(y(:,i)-H1*x2minus(:,i));
    P2plus(:,:,i) = (eye(6)-W2*H1)*P2minus(:,:,i);
    % State and covariance propagation
    x2minus(:,i+1) = phi1*x2plus(:,i) + Gamma1*U(:,i);
    P2minus(:,:,i+1)=phi1*P2plus(:,:,i)*phi1' + upsilon*Q*upsilon' ;
    E2(:,:,i) = H1*P2minus(:,:,i)*H1'+R;
end
W2 = P2minus(:,:,end)*H1'/(H1*P2minus(:,:,end)*H1'+R);
x2plus(:,end) = x2minus(:,end) + W2*(y(:,end)-H1*x2minus(:,end));
P2plus(:,:,end) = (eye(6)-W2*H1)*P2minus(:,:,end);
E2(:,:,end) = H1*P2minus(:,:,end)*H1'+R;
meas_res2 = y - H1*x2minus;


%% Likelihood function
W = zeros(2,length(t));
W(:,1) = [0.5 0.5]';
for i = 2:length(t)
    W(1,i) = W(1,i-1)*exp(-0.5*(meas_res1(:,i))'*pinv(E1(:,:,i))*meas_res1(:,i))/(det(2*pi*E1(:,:,i)))^0.5;
    W(2,i) = W(2,i-1)*exp(-0.5*(meas_res2(:,i))'*pinv(E2(:,:,i))*meas_res2(:,i))/(det(2*pi*E2(:,:,i)))^0.5;
    k = W(1,i)+W(2,i);
    W(1,i) = W(1,i)/k;
    W(2,i) = W(2,i)/k;
end

% Updated state and error covariance
xplus = zeros(6,length(t));
Pplus = zeros(6,6,length(t));
P_est = zeros(2,length(t));
P_var = zeros(2,2,length(t));

for i = 1:length(t)
    xplus(:,i) = W(1,i)*x1plus(:,i) + W(2,i)*x2plus(:,i);
    Pplus(:,:,i) = W(1,i)*((x1plus(:,i)-xplus(:,i))*(x1plus(:,i)-xplus(:,i))'+P1plus(:,:,i)) ...
        + W(2,i)*((x2plus(:,i)-xplus(:,i))*(x2plus(:,i)-xplus(:,i))'+P2plus(:,:,i));
    P_est(:,i) = W(1,i)*[0 0]' + W(2,i)*[1e-3 1e-3]';
    P_var(:,:,i) = W(1,i)*(([0 0]'-P_est(:,i))*([0 0]'-P_est(:,i))') + W(2,i)*(([1e-3 1e-3]'-P_est(:,i))*([1e-3 1e-3]'-P_est(:,i))');
end
Pstd = Pplus.^0.5;
Pos_std = reshape(Pstd(1,1,:),1,[]);
figure(4)
plot(t,xplus(1,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Postion estimate')
title('MMAE estimate vs time')
grid minor
figure(5)
subplot(2,1,1)
plot(t, xplus(1,:)-x(1,:), t, 3*Pos_std, t, -3*Pos_std,'LineWidth',2)
xlabel('time in seconds')
ylabel('Position error')
title('MMAE Error')
grid minor
subplot(2,1,2)
plot(t, W(1,:), t, W(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('MMAE weights')
title('MMAE Weights vs time ')
grid minor
%% IMM filter
% filter 1 -----> initialization 
x3plus = zeros(6,length(t),2);
x3minus = zeros(6,length(t),2);
x3minus(:,1,1) = [2000 0 0 10000 -15 0]';
x3minus(:,1,2) = [2000 0 0 10000 -15 0]';
P3minus = zeros(6,6,length(t));
P3plus = zeros(6,6,length(t));
P3minus(:,:,1) = 1e-12*eye(6);
E3 = zeros(2,2,length(t));
% filter 2 ----> initialization
P4minus = zeros(6,6,length(t));
P4plus = zeros(6,6,length(t));
P4minus(:,:,1) = 1e-12*eye(6);
E4 = zeros(2,2,length(t));
% Initial conditions for each filter for mixing
x0plus = zeros(6,length(t),2);
P01plus = zeros(6,6,length(t));
P02plus = zeros(6,6,length(t));
% IMM states and covaraince
xplus_I = zeros(6,length(t));
Pplus_I = zeros(6,6,length(t));
% Weights
Wk = zeros(2,length(t));
Wk(:,1) = [0.5 0.5]';
% Conditonal weights
Wk_11 = zeros(1,length(t));
Wk_12 = zeros(1,length(t));
Wk_21 = zeros(1, length(t));
Wk_22 = zeros(1,length(t));
% transitional probabilities
p11 = 0.97; p12 = 0.03; p21 = 0.03; p22 = 0.97;
meas_res = zeros(2,length(t),2);
meas_res(:,1,1) = y(:,1) - H1*x3minus(:,1,1);
meas_res(:,1,2) = y(:,1) - H1*x3minus(:,1,2);
for i = 1:length(t)-1
    % Kalman gain for filter 1 and filter 2
    W3 = P3minus(:,:,i)*H1'/(H1*P3minus(:,:,i)*H1'+R);
    W4 = P4minus(:,:,i)*H1'/(H1*P4minus(:,:,i)*H1'+R);
    % State update and covariance update --- filter 1
    x3plus(:,i,1) = x3minus(:,i,1) + W3*(y(:,i)-H1*x3minus(:,i,1));
    P3plus(:,:,i) = (eye(6)-W3*H1)*P3minus(:,:,i);
    % filter 2
    x3plus(:,i,2) = x3minus(:,i,2) + W4*(y(:,i)-H1*x3minus(:,i,2));
    P4plus(:,:,i) = (eye(6)-W4*H1)*P4minus(:,:,i);
   % filter 1 residual error covariance
    E3(:,:,i) = H1*P3minus(:,:,i)*H1' + R;
      % filter 2residual error covariance
    E4(:,:,i) = H1*P4minus(:,:,i)*H1' + R;
  % Interaction --- mixing probabilities
    ck_1 = Wk(1,i)*p11 + Wk(2,i)*p21;
    ck_2 = Wk(1,i)*p12 + Wk(2,i)*p22;
    Wk_11(i) = Wk(1,i)*p11/(ck_1);
    Wk_12(i) = Wk(1,i)*p12/ck_2;
    Wk_21(i) = Wk(2,i)*p21/ck_1;
    Wk_22(i) = Wk(2,i)*p22/ck_2;

     % Mixed initial conditions ---> State
    x0plus(:,i,1) = Wk_11(i)*x3plus(:,i,1)  + Wk_21(i)*x3plus(:,i,2);
    x0plus(:,i,2) = Wk_12(i)*x3plus(:,i,1)  + Wk_22(i)*x3plus(:,i,2);

     % Mixed initial covariances 
    P01plus(:,:,i) = Wk_11(i)*(P3plus(:,:,i)+ (x3plus(:,i,1)-x0plus(:,i,1))*(x3plus(:,i,1)-x0plus(:,i,1))')...
         + Wk_21(i)*(P4plus(:,:,i)+ (x3plus(:,i,2)-x0plus(:,i,2))*(x3plus(:,i,2)-x0plus(:,i,2))');

    P02plus(:,:,i) = Wk_12(i)*(P3plus(:,:,i)+ (x3plus(:,i,1)-x0plus(:,i,1))*(x3plus(:,i,1)-x0plus(:,i,1))')...
         + Wk_22(i)*(P4plus(:,:,i)+ (x3plus(:,i,2)-x0plus(:,i,2))*(x3plus(:,i,2)-x0plus(:,i,2))');

    % State and covariance propagation -- filter 1
    x3minus(:,i+1,1) = phi1*x0plus(:,i,1) + 0*Gamma1*U(:,i);
    P3minus(:,:,i+1)=phi1*P01plus(:,:,i)*phi1' ;
    meas_res(:,i+1,1) = y(:,i+1) - H1*x3minus(:,i+1,1);
        % State and covariance propagation -- filter 2
    x3minus(:,i+1,2) = phi1*x0plus(:,i,2) + 0*Gamma1*U(:,i);
    P4minus(:,:,i+1)=phi1*P02plus(:,:,i)*phi1' + upsilon*Q*upsilon' ;
    meas_res(:,i+1,2) = y(:,i+1) - H1*x3minus(:,i+1,2);
    % Weight updates 
    Wk(1,i+1) = Wk(1,i)*exp(-0.5*(meas_res(:,i+1,1))'*pinv(E3(:,:,i))*meas_res(:,i+1,1))/(det(2*pi*E3(:,:,i)))^0.5;
    Wk(2,i+1) = Wk(2,i)*exp(-0.5*(meas_res(:,i+1,2))'*pinv(E4(:,:,i))*meas_res(:,i+1,2))/(det(2*pi*E4(:,:,i)))^0.5;

    k = Wk(1,i+1) + Wk(2,i+1);
    Wk(1,i+1) = Wk(1,i+1)/k;
    Wk(2,i+1) = Wk(2,i+1)/k;

    % IMM state and covariance update
    xplus_I(:,i) = Wk(1,i)*x3plus(:,i,1) + Wk(2,i)*x3plus(:,i,2);
    Pplus_I(:,:,i) = Wk(1,i)*((x3plus(:,i,1)-xplus_I(:,i))*(x3plus(:,i,1)-xplus_I(:,i))'+P3plus(:,:,i)) ...
        + Wk(2,i)*((x3plus(:,i,2)-xplus_I(:,i))*(x3plus(:,i,2)-xplus_I(:,i))'+P4plus(:,:,i));
end
i = i+1;
% Kalman gain for filter 1 and filter 2
W3 = P3minus(:,:,i)*H1'/(H1*P3minus(:,:,i)*H1'+R);
W4 = P4minus(:,:,i)*H1'/(H1*P4minus(:,:,i)*H1'+R);
% State update and covariance update --- filter 1
x3plus(:,i,1) = x3minus(:,i,1) + W3*(y(:,i)-H1*x3minus(:,i,1));
P3plus(:,:,i) = (eye(6)-W3*H1)*P3minus(:,:,i);
% filter 2
x3plus(:,i,2) = x3minus(:,i,2) + W4*(y(:,i)-H1*x3minus(:,i,2));
P4plus(:,:,i) = (eye(6)-W4*H1)*P4minus(:,:,i);
% filter 1 residual error covariance
E3(:,:,i) = H1*P3minus(:,:,i)*H1' + R;
  % filter 2residual error covariance
E4(:,:,i) = H1*P4minus(:,:,i)*H1' + R;
% Interaction --- mixing probabilities
ck_1 = Wk(1,i)*p11 + Wk(2,i)*p21;
ck_2 = Wk(1,i)*p12 + Wk(2,i)*p22;
Wk_11(i) = Wk(1,i)*p11/(ck_1);
Wk_12(i) = Wk(1,i)*p12/ck_2;
Wk_21(i) = Wk(2,i)*p21/ck_1;
Wk_22(i) = Wk(2,i)*p22/ck_2;

 % Mixed initial conditions ---> State
x0plus(:,i,1) = Wk_11(i)*x3plus(:,i,1)  + Wk_21(i)*x3plus(:,i,2);
x0plus(:,i,2) = Wk_12(i)*x3plus(:,i,1)  + Wk_22(i)*x3plus(:,i,2);

 % Mixed initial covariances 
P01plus(:,:,i) = Wk_11(i)*(P3plus(:,:,i)+ (x3plus(:,i,1)-x0plus(:,i,1))*(x3plus(:,i,1)-x0plus(:,i,1))')...
     + Wk_21(i)*(P4plus(:,:,i)+ (x3plus(:,i,2)-x0plus(:,i,2))*(x3plus(:,i,2)-x0plus(:,i,2))');

P02plus(:,:,i) = Wk_12(i)*(P3plus(:,:,i)+ (x3plus(:,i,1)-x0plus(:,i,1))*(x3plus(:,i,1)-x0plus(:,i,1))')...
     + Wk_22(i)*(P4plus(:,:,i)+ (x3plus(:,i,2)-x0plus(:,i,2))*(x3plus(:,i,2)-x0plus(:,i,2))');

% IMM state and covariance update
xplus_I(:,i) = Wk(1,i)*x3plus(:,i,1) + Wk(2,i)*x3plus(:,i,2);
Pplus_I(:,:,i) = Wk(1,i)*((x3plus(:,i,1)-xplus_I(:,i))*(x3plus(:,i,1)-xplus_I(:,i))'+P3plus(:,:,i)) ...
    + Wk(2,i)*((x3plus(:,i,2)-xplus_I(:,i))*(x3plus(:,i,2)-xplus_I(:,i))'+P4plus(:,:,i));

Pstd_IMM = Pplus_I.^0.5;
Pos_std_IMM = reshape(Pstd_IMM(1,1,:),1,[]);
figure(6)
plot(t,xplus_I(1,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Postion estimate')
title('IMM estimate vs time')
grid minor
figure(7)
subplot(2,1,1)
plot(t, xplus_I(1,:)-x(1,:), t, 3*Pos_std_IMM, t, -3*Pos_std_IMM,'LineWidth',2)
xlabel('time in seconds')
ylabel('Position error')
title('IMM Error')
grid minor
subplot(2,1,2)
plot(t, Wk(1,:), t, Wk(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('IMM weights')
title('IMM Weights vs time ')
grid minor