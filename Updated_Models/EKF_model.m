clc; clear;
generate_input_output_data
n = 5; % Five states -- Global X coordinate, Global Y coordinate, yaw angle, body velocity and angular velocity. we have a nonlinear system 
q = 1; % Process noise covariance with wheel angular velocities
r = 1; % measurement noise covariance
m = length(Output1(:,1)); % available measurements
np = 4; % filters numbers
dt = t(2)-t(1);
% Model 1 
K1 = 0.0763; M1 = 0.1134;
Pv = 0.012; Pw = 0.0025;
% Model 2 
K2 = 0.0784; M2 = 0.1257;
% Model 3 
K3 = 0.0776; M3 = 0.1196;
% Model 4 
K4 = 0.0774; M4 = 0.1216;
K = [K1 K2 K3 K4;
    M1 M2 M3 M4];
P = [Pv Pv Pv Pv;
    Pw Pw Pw Pw];
Model_num = 7;
if Model_num == 1
    U = Input1;
    K_M = K(:,1);
    Pvw = P(:,1);
    x_0 = Output1(:,1);
    y_gt = Output1;
elseif Model_num == 2
    U = Input2;
    K_M = K(:,1);
    Pvw = P(:,1);
    x_0 = Output2(:,1);
    y_gt = Output2;
elseif Model_num == 3
    U = Input3;
    K_M = K(:,2);
    t = t2;
    Pvw = P(:,2);
    x_0 = Output3(:,1);
    y_gt = Output3;
elseif Model_num == 4
    U = Input4;
    K_M = K(:,3);
    Pvw = P(:,2);
    x_0 = Output4(:,1);
    y_gt = Output4;
elseif Model_num == 5
    U = Input5;
    K_M = K(:,3);
    t = t2;
    Pvw = P(:,3);
    x_0 = Output5(:,1);
    y_gt = Output5;
elseif Model_num == 6
    U = Input6;
    K_M = K(:,3);
    Pvw = P(:,3);
    x_0 = Output6(:,1);
    y_gt = Output6;
elseif Model_num == 7
    U = Input7;
    K_M = K(:,4);
    t = t2;
    Pvw = P(:,4);
    x_0 = Output7(:,1);
    y_gt = Output7;
else
    U = Input8;
    K_M = K(:,4);
    Pvw = P(:,4);
    x_0 = Output8(:,1);
    y_gt = Output8;
end
B = [0 0;
 0 0;
 0 0;
 Pvw(1)*(1-K_M(1)) 0;
 0 Pvw(2)*(1-K_M(2))]; % Input matrix
upsilon = B;
H = eye(m);
Q = upsilon*q*dt*eye(2)*upsilon'; % process noise
R = eye(m)*r;
P0 = eye(n)*0.02; % Initial state error covariance
X_hat = zeros(n,length(t)); % states output for each filter through out the entire time -- column 1 - filter 1, column 2 - filter 2 ...
X_hat(:,1) = x_0;
P_plus = zeros(n,n,length(t));
P_plus(:,:,1) = P0; % 1st filter covariance 
ne = 1; % switch for effects of noise...
% 1 --> process noise is considered
% anyything else --> no effects of noise
for i = 1:length(t)-1
    U1 = U(:,i);
    y_op = y_gt(:,i+1) + sqrt(r)*randn(5,1);
    
   [X_hat(:,i+1), P_plus(:,:,i+1)] = EKF(K_M, Pvw, X_hat(:,i), y_op, P_plus(:,:,i), U1, dt, Q,R,H);
  
end
x = zeros(n,length(t));
x(:,1) = x_0;
for i = 1:length(t)-1
    xdot = dxdt(x(:,i),Input1(:,i),K_M,Pvw);
        if ne == 1
        upsilon = [0 0;
                   0 0;
                   0 0;
                   Pvw(1)*(1-K_M(1)) 0;
                   0 Pvw(2)*(1-K_M(2))];
        x(:,i+1) = x(:,i) + xdot*dt + upsilon*sqrt(q*dt)*randn(2,1);
        else 
        x(:,i+1) = x(:,i) + xdot*dt;
        end

end
[t1,x1]=ode45(@(t1,x1) odefun(t1,x1,t,U, dt, q,K_M,Pvw,ne),t,x_0);

figure
plot(x1(:,1),x1(:,2),X_hat(1,:),X_hat(2,:),'Linewidth',2)
xlabel('Longitudinal position')
ylabel('Lateral position')
title('trajectory')
legend('True','Estimate')
grid minor
% Plots compating Euler and ode45
figure
subplot(5,1,1)
plot(t1,x1(:,1),t,X_hat(1,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Global X position')
title('Global X vs time')
legend('True','Estimate')
grid minor
subplot(5,1,2)
plot(t1,x1(:,2),t,X_hat(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Global Y position')
title('Global Y vs time')
legend('True','Estimate')
grid minor
subplot(5,1,3)
plot(t1,x1(:,3),t ,X_hat(3,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Yaw angle')
title('Yaw angle vs time')
legend('True','Estimate')
grid minor
subplot(5,1,4)
plot(t1,x1(:,4),t,X_hat(4,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Body linear velocity')
title('v vs time')
legend('True','Estimate')
grid minor

subplot(5,1,5)
plot(t1,x1(:,5),t,X_hat(5,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Angular velocity')
title('W vs time')
legend('True','Estimate')
grid minor

%% function space
function [X_plus, P_plus, Likelihood] = EKF(K, P, X_mix_init, y, P_mix_init, U, dt, Q,R,H)
K1 = K(1);
M1 = K(2);
Pv = P(1);
Pw = P(2);
B = [0 0;
 0 0;
 0 0;
 Pv*(1-K1) 0;
 0 Pw*(1-M1)]; % input matrix
 F = [0 0 -X_mix_init(4)*sin(X_mix_init(3)) cos(X_mix_init(3)) 0;
     0 0 X_mix_init(4)*cos(X_mix_init(3)) sin(X_mix_init(4)) 0;
     0 0 0 0 1;
     0 0 0 0 0;
     0 0 0 0 0]; % Jacobain matrix
xdot = [X_mix_init(4)*cos(X_mix_init(3));
    X_mix_init(4)*sin(X_mix_init(3));
    X_mix_init(5);
    0;
    0] + B*U;
% Propagation
X_minus = X_mix_init + dt*xdot; 
F_discrete = eye(5)+dt*F;
P_minus = F_discrete*P_mix_init*F_discrete' + Q;
% Kalman gain
KG = P_minus*H'/(H*P_minus*H'+R);
% Update
X_plus = X_minus + KG*(y - H*X_minus);
P_plus = (eye(5)-KG*H)*P_minus;
% Likelihood fucntion
e = y - H*X_minus;
E = H*P_minus*H'+R;
Likelihood = ((det(2*pi*E))^0.5)*exp(-0.5*e'*(inv(E))*e);
end

function xdot = dxdt(state, Input,M,P)
% Model 1 
x = state;
K1 = M(1); M1 = M(2);
Pv = P(1); Pw = P(2);
x1dot = x(4)*cos(x(3));
x2dot = x(4)*sin(x(3));
x3dot = x(5);
x4dot =  Pv*(1-K1)*Input(1);
x5dot = Pw*(1-M1)*Input(2);

xdot = [x1dot x2dot x3dot x4dot x5dot]';
end

function xdot = odefun(t1,state,t,U, dt, q,M,P,ne)
% Model 1 
x = state;
K1 = M(1); M1 = M(2);
Pv = P(1); Pw = P(2);
upsilon = [0 0;
 0 0;
 0 0;
 Pv*(1-K1) 0;
 0 Pw*(1-M1)];
Input1 = interp1(t,U(1,:),t1);
Input2 = interp1(t,U(2,:),t1);

x1dot = x(4)*cos(x(3));
x2dot = x(4)*sin(x(3));
x3dot = x(5);
x4dot = Pv*(1-K1)*Input1;
x5dot = Pw*(1-M1)*Input2;
if ne == 1
    xdot = [x1dot x2dot x3dot x4dot x5dot]' + upsilon*sqrt(q*dt)*randn(2,1);
else 
    xdot = [x1dot x2dot x3dot x4dot x5dot]';
end
end
