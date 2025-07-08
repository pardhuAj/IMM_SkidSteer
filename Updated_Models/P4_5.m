clc;clear;
phi = [0.9999 0.0099;
    -0.0296 0.9703];
upsilon = [0 0.01]';
H = [1 0];
q = 10;
r = 0.01;
rng(0,'twister');
dt = 0.01;
T = 10;
t = 0:dt:T;
v = sqrt(r)*randn(1,length(t));
w = sqrt(q)*randn(1,length(t));
x = zeros(2,length(t));
x(:,1) = [1 1]';

for i = 1:length(t)-1
    x(:,i+1) = phi*x(:,i) + upsilon*w(i);
end

y = H*x + v;

%%%% Discrete kalman filter

Q_0 = [0.1 0.5 1 10 20 100 1000 1e4 1e5];
R_0 = r;

% filter 1
x_est_minus1 = zeros(2,length(t));
x_est_plus1 = zeros(2,length(t));
P_minus1 = zeros(2,2,length(t));
P_plus1 = zeros(2,2,length(t));
x_est_minus1(:,1) = [0.999 0.999]';
P_minus1(:,:,1) = 0.001^2*eye(2);
E1 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W1 = P_minus1(:,:,i)*H'/(H*P_minus1(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus1(:,i) = x_est_minus1(:,i) + W1*(y(i) - H*x_est_minus1(:,i));
    P_plus1(:,:,i) = (eye(2) - W1*H)*P_minus1(:,:,i);
    % Propagate eqaution
    x_est_minus1(:,i+1) = phi*x_est_plus1(:,i);
    P_minus1(:,:,i+1) = phi*P_plus1(:,:,i)*phi'+upsilon*Q_0(1)*upsilon';
    E1(:,:,i) = H*P_minus1(:,:,i)*H'+R_0;
end
E1(:,:,end) = H*P_minus1(:,:,end)*H'+R_0;
meas_res1 = y - H*x_est_minus1;
% filter 2
x_est_minus2 = zeros(2,length(t));
x_est_plus2 = zeros(2,length(t));
P_minus2 = zeros(2,2,length(t));
P_plus2 = zeros(2,2,length(t));
x_est_minus2(:,1) = [0.999 0.999]';
P_minus2(:,:,1) = 0.001^2*eye(2);
E2 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W2 = P_minus2(:,:,i)*H'/(H*P_minus2(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus2(:,i) = x_est_minus2(:,i) + W2*(y(i) - H*x_est_minus2(:,i));
    P_plus2(:,:,i) = (eye(2) - W2*H)*P_minus2(:,:,i);
    % Propagate eqaution
    x_est_minus2(:,i+1) = phi*x_est_plus2(:,i);
    P_minus2(:,:,i+1) = phi*P_plus2(:,:,i)*phi'+upsilon*Q_0(2)*upsilon';
    E2(:,:,i) = H*P_minus2(:,:,i)*H'+R_0;
end
E2(:,:,end) = H*P_minus2(:,:,end)*H'+R_0;
meas_res2 = y - H*x_est_minus2;

% filter 3
x_est_minus3 = zeros(2,length(t));
x_est_plus3 = zeros(2,length(t));
P_minus3 = zeros(2,2,length(t));
P_plus3 = zeros(2,2,length(t));
x_est_minus3(:,1) = [0.999 0.999]';
P_minus3(:,:,1) = 0.001^2*eye(2);
E3 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W3 = P_minus3(:,:,i)*H'/(H*P_minus3(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus3(:,i) = x_est_minus3(:,i) + W3*(y(i) - H*x_est_minus3(:,i));
    P_plus3(:,:,i) = (eye(2) - W3*H)*P_minus3(:,:,i);
    % Propagate eqaution
    x_est_minus3(:,i+1) = phi*x_est_plus3(:,i);
    P_minus3(:,:,i+1) = phi*P_plus3(:,:,i)*phi'+upsilon*Q_0(3)*upsilon';
    E3(:,:,i) = H*P_minus3(:,:,i)*H'+R_0;
end
E3(:,:,end) = H*P_minus3(:,:,end)*H'+R_0;
meas_res3 = y - H*x_est_minus3;

% filter 4
x_est_minus4 = zeros(2,length(t));
x_est_plus4 = zeros(2,length(t));
P_minus4 = zeros(2,2,length(t));
P_plus4 = zeros(2,2,length(t));
x_est_minus4(:,1) = [0.999 0.999]';
P_minus4(:,:,1) = 0.001^2*eye(2);
E4 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W4 = P_minus4(:,:,i)*H'/(H*P_minus4(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus4(:,i) = x_est_minus4(:,i) + W4*(y(i) - H*x_est_minus4(:,i));
    P_plus4(:,:,i) = (eye(2) - W4*H)*P_minus4(:,:,i);
    % Propagate eqaution
    x_est_minus4(:,i+1) = phi*x_est_plus4(:,i);
    P_minus4(:,:,i+1) = phi*P_plus4(:,:,i)*phi'+upsilon*Q_0(4)*upsilon';
    E4(:,:,i) = H*P_minus4(:,:,i)*H'+R_0;
end
E4(:,:,end) = H*P_minus4(:,:,end)*H'+R_0;
meas_res4 = y - H*x_est_minus4;

% filter 5
x_est_minus5 = zeros(2,length(t));
x_est_plus5 = zeros(2,length(t));
P_minus5 = zeros(2,2,length(t));
P_plus5 = zeros(2,2,length(t));
x_est_minus5(:,1) = [0.999 0.999]';
P_minus5(:,:,1) = 0.001^2*eye(2);
E5 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W5 = P_minus5(:,:,i)*H'/(H*P_minus5(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus5(:,i) = x_est_minus5(:,i) + W5*(y(i) - H*x_est_minus5(:,i));
    P_plus5(:,:,i) = (eye(2) - W5*H)*P_minus5(:,:,i);
    % Propagate eqaution
    x_est_minus5(:,i+1) = phi*x_est_plus5(:,i);
    P_minus5(:,:,i+1) = phi*P_plus5(:,:,i)*phi'+upsilon*Q_0(5)*upsilon';
    E5(:,:,i) = H*P_minus5(:,:,i)*H'+R_0;
end
E5(:,:,end) = H*P_minus5(:,:,end)*H'+R_0;
meas_res5 = y - H*x_est_minus5;

% filter 6
x_est_minus6 = zeros(2,length(t));
x_est_plus6 = zeros(2,length(t));
P_minus6 = zeros(2,2,length(t));
P_plus6 = zeros(2,2,length(t));
x_est_minus6(:,1) = [0.999 0.999]';
P_minus6(:,:,1) = 0.001^2*eye(2);
E6 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W6 = P_minus6(:,:,i)*H'/(H*P_minus6(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus6(:,i) = x_est_minus6(:,i) + W6*(y(i) - H*x_est_minus6(:,i));
    P_plus6(:,:,i) = (eye(2) - W6*H)*P_minus6(:,:,i);
    % Propagate eqaution
    x_est_minus6(:,i+1) = phi*x_est_plus6(:,i);
    P_minus6(:,:,i+1) = phi*P_plus6(:,:,i)*phi'+upsilon*Q_0(6)*upsilon';
    E6(:,:,i) = H*P_minus6(:,:,i)*H'+R_0;
end
E6(:,:,end) = H*P_minus6(:,:,end)*H'+R_0;
meas_res6 = y - H*x_est_minus6;

% filter 7
x_est_minus7 = zeros(2,length(t));
x_est_plus7 = zeros(2,length(t));
P_minus7 = zeros(2,2,length(t));
P_plus7 = zeros(2,2,length(t));
x_est_minus7(:,1) = [0.999 0.999]';
P_minus7(:,:,1) = 0.001^2*eye(2);
E7 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W7 = P_minus7(:,:,i)*H'/(H*P_minus7(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus7(:,i) = x_est_minus7(:,i) + W7*(y(i) - H*x_est_minus7(:,i));
    P_plus7(:,:,i) = (eye(2) - W7*H)*P_minus7(:,:,i);
    % Propagate eqaution
    x_est_minus7(:,i+1) = phi*x_est_plus7(:,i);
    P_minus7(:,:,i+1) = phi*P_plus7(:,:,i)*phi'+upsilon*Q_0(7)*upsilon';
    E7(:,:,i) = H*P_minus7(:,:,i)*H'+R_0;
end
E7(:,:,end) = H*P_minus7(:,:,end)*H'+R_0;
meas_res7 = y - H*x_est_minus7;

% filter 8
x_est_minus8 = zeros(2,length(t));
x_est_plus8 = zeros(2,length(t));
P_minus8 = zeros(2,2,length(t));
P_plus8 = zeros(2,2,length(t));
x_est_minus8(:,1) = [0.999 0.999]';
P_minus8(:,:,1) = 0.001^2*eye(2);
E8 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W8 = P_minus8(:,:,i)*H'/(H*P_minus8(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus8(:,i) = x_est_minus8(:,i) + W8*(y(i) - H*x_est_minus8(:,i));
    P_plus8(:,:,i) = (eye(2) - W8*H)*P_minus8(:,:,i);
    % Propagate eqaution
    x_est_minus8(:,i+1) = phi*x_est_plus8(:,i);
    P_minus8(:,:,i+1) = phi*P_plus8(:,:,i)*phi'+upsilon*Q_0(8)*upsilon';
    E8(:,:,i) = H*P_minus8(:,:,i)*H'+R_0;
end
E8(:,:,end) = H*P_minus8(:,:,end)*H'+R_0;
meas_res8 = y - H*x_est_minus8;

% filter 9
x_est_minus9 = zeros(2,length(t));
x_est_plus9 = zeros(2,length(t));
P_minus9 = zeros(2,2,length(t));
P_plus9 = zeros(2,2,length(t));
x_est_minus9(:,1) = [0.999 0.999]';
P_minus9(:,:,1) = 0.001^2*eye(2);
E9 = zeros(1,1,length(t));
for i = 1:length(t)-1
    % Kalman gain calculation
    W9 = P_minus9(:,:,i)*H'/(H*P_minus9(:,:,i)*H' +R_0);
    % Update equation
    x_est_plus9(:,i) = x_est_minus9(:,i) + W9*(y(i) - H*x_est_minus9(:,i));
    P_plus9(:,:,i) = (eye(2) - W9*H)*P_minus9(:,:,i);
    % Propagate eqaution
    x_est_minus9(:,i+1) = phi*x_est_plus9(:,i);
    P_minus9(:,:,i+1) = phi*P_plus9(:,:,i)*phi'+upsilon*Q_0(9)*upsilon';
    E9(:,:,i) = H*P_minus9(:,:,i)*H'+R_0;
end
E9(:,:,end) = H*P_minus9(:,:,end)*H'+R_0;
meas_res9 = y - H*x_est_minus9;


%%%% likelihood 
W = zeros(length(t),9);
W(1,:) = 1/9;
for i = 2:length(t)

    W(i,1) = W(i-1,1)*((det(2*pi*E1(:,:,i)))^-0.5)*exp(-0.5*(meas_res1(:,i)'*pinv(E1(:,:,i)))*meas_res1(:,i)) ;
    W(i,2) = W(i-1,2)*((det(2*pi*E2(:,:,i)))^-0.5)*exp(-0.5*meas_res2(:,i)'*pinv(E2(:,:,i))*meas_res2(:,i)) ;
    W(i,3) = W(i-1,3)*((det(2*pi*E3(:,:,i)))^-0.5)*exp(-0.5*meas_res3(:,i)'*pinv(E3(:,:,i))*meas_res3(:,i)) ;
    W(i,4) = W(i-1,4)*((det(2*pi*E4(:,:,i)))^-0.5)*exp(-0.5*meas_res4(:,i)'*pinv(E4(:,:,i))*meas_res4(:,i)) ;
    W(i,5) = W(i-1,5)*((det(2*pi*E5(:,:,i)))^-0.5)*exp(-0.5*meas_res5(:,i)'*pinv(E5(:,:,i))*meas_res5(:,i)) ;
    W(i,6) = W(i-1,6)*((det(2*pi*E6(:,:,i)))^-0.5)*exp(-0.5*meas_res6(:,i)'*pinv(E6(:,:,i))*meas_res6(:,i)) ;
    W(i,7) = W(i-1,7)*((det(2*pi*E7(:,:,i)))^-0.5)*exp(-0.5*meas_res7(:,i)'*pinv(E7(:,:,i))*meas_res7(:,i)) ;
    W(i,8) = W(i-1,8)*((det(2*pi*E8(:,:,i)))^-0.5)*exp(-0.5*meas_res8(:,i)'*pinv(E8(:,:,i))*meas_res8(:,i)) ;
    W(i,9) = W(i-1,9)*((det(2*pi*E9(:,:,i)))^-0.5)*exp(-0.5*meas_res9(:,i)'*pinv(E9(:,:,i))*meas_res9(:,i)) ;
    k = W(i,1)+W(i,2) + W(i,3) + W(i,4) + W(i,5)+W(i,6) + W(i,7) + W(i,8) + W(i,9);
    W(i,1) = W(i,1)/k;
    W(i,2) = W(i,2)/k;
    W(i,3) = W(i,3)/k;
    W(i,4) = W(i,4)/k;
    W(i,5) = W(i,5)/k;
    W(i,6) = W(i,6)/k;
    W(i,7) = W(i,7)/k;
    W(i,8) = W(i,8)/k;
    W(i,9) = W(i,9)/k;
end


% Estimate and updated error covaiance -- state
% Estimate and updated error covaiance -- parameter

X_plus = zeros(2,length(t));
P_plus = zeros(2,2,length(t));
P_est = zeros(1,length(t));
Rho_est = zeros(1,length(t));

for i = 1:length(t)
    X_plus(:,i) = W(i,1)*x_est_plus1(:,i) + W(i,2)*x_est_plus2(:,i) + W(i,3)*x_est_plus3(:,i) + W(i,4)*x_est_plus4(:,i)...
         + W(i,5)*x_est_plus5(:,i) + W(i,6)*x_est_plus6(:,i) + W(i,7)*x_est_plus7(:,i) + W(i,8)*x_est_plus8(:,i) ...
          + W(i,9)*x_est_plus9(:,i);
    P_plus(:,:,i) = W(i,1)*(P_plus1(:,:,i)+(x_est_plus1(:,i)-X_plus(:,i))*(x_est_plus1(:,i)-X_plus(:,i))') ...
        + W(i,2)*(P_plus2(:,:,i)+(x_est_plus2(:,i)-X_plus(:,i))*(x_est_plus2(:,i)-X_plus(:,i))') ...
        + W(i,3)*(P_plus3(:,:,i)+(x_est_plus3(:,i)-X_plus(:,i))*(x_est_plus3(:,i)-X_plus(:,i))') ...
        + W(i,4)*(P_plus4(:,:,i)+(x_est_plus4(:,i)-X_plus(:,i))*(x_est_plus4(:,i)-X_plus(:,i))') ...
        + W(i,5)*(P_plus5(:,:,i)+(x_est_plus5(:,i)-X_plus(:,i))*(x_est_plus5(:,i)-X_plus(:,i))') ...
        + W(i,6)*(P_plus6(:,:,i)+(x_est_plus6(:,i)-X_plus(:,i))*(x_est_plus6(:,i)-X_plus(:,i))') ...
        + W(i,7)*(P_plus7(:,:,i)+(x_est_plus7(:,i)-X_plus(:,i))*(x_est_plus7(:,i)-X_plus(:,i))') ...
        + W(i,8)*(P_plus8(:,:,i)+(x_est_plus8(:,i)-X_plus(:,i))*(x_est_plus8(:,i)-X_plus(:,i))') ...
        + W(i,9)*(P_plus9(:,:,i)+(x_est_plus9(:,i)-X_plus(:,i))*(x_est_plus9(:,i)-X_plus(:,i))');
    P_est(i) = W(i,1)*Q_0(1) + W(i,2)*Q_0(2) + W(i,3)*Q_0(3) + W(i,4)*Q_0(4) + W(i,5)*Q_0(5) ...
        + W(i,6)*Q_0(6) + W(i,7)*Q_0(7) +W(i,8)*Q_0(8) + W(i,9)*Q_0(9) ;
    Rho_est(i) = W(i,1)*((Q_0(1)-P_est(i))*(Q_0(1)-P_est(i))')+ W(i,2)*((Q_0(2)-P_est(i))*(Q_0(2)-P_est(i))') ...
        + W(i,3)*((Q_0(3)-P_est(i))*(Q_0(3)-P_est(i))') + W(i,4)*((Q_0(4)-P_est(i))*(Q_0(4)-P_est(i))')...
        + W(i,5)*((Q_0(5)-P_est(i))*(Q_0(5)-P_est(i))') ...
        + W(i,6)*((Q_0(6)-P_est(i))*(Q_0(6)-P_est(i))') + W(i,7)*((Q_0(7)-P_est(i))*(Q_0(7)-P_est(i))')...
        +W(i,8)*((Q_0(8)-P_est(i))*(Q_0(8)-P_est(i))') + W(i,9)*((Q_0(9)-P_est(i))*(Q_0(9)-P_est(i))') ;
end

rho_std = Rho_est.^0.5;
figure(1)
plot(t,3*rho_std, t,-3*rho_std,t, P_est-q,'linewidth',2)
xlabel('time in seocnds')
ylabel('Parameter error')
title('MMAE Parameter estimate errors')
ylim([-100 100])
grid minor





