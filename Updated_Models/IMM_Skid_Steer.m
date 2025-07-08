%% Interacting Mixing Model approach for off road navigation
%%% We have four different models to capture the uncertain parameters

%%% Wheel angular velocities are transformed to body velcoity and angular
%%% velocity
clc; clear; 
generate_input_output_data
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
Model_num = 1;
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
% Transitional probability matrix 
P_tran = [0.95 0.02 0.02 0.01;
                          0.05 0.85 0.05 0.05;
                          0.05 0.05 0.85 0.05;
                          0.05 0.05 0.05 0.85]; 
% noise information
n = 5; % Five states -- Global X coordinate, Global Y coordinate, yaw angle, body velocity and angular velocity. we have a nonlinear system 
q = 0.1; % Process noise covariance with wheel angular velocities
r = 1; % measurement noise covariance
m = length(Output1(:,1)); % available measurements
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
%%% True state information
ne = 1; % switch for effects of noise...
% 1 --> process noise is considered
% anyything else --> no effects of noise
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
% for  i = 1:length(t)
%     y_gt(:,i) = H*x1(i,:)';
% end
% Plots compating Euler and ode45
figure(1)
subplot(5,1,1)
plot(t,x(1,:),t1,x1(:,1),'LineWidth',2)
xlabel('time in seconds')
ylabel('Global X position')
title('Global X vs time')
legend('Euler','Ode45')
grid minor
subplot(5,1,2)
plot(t,x(2,:),t1,x1(:,2),'LineWidth',2)
xlabel('time in seconds')
ylabel('Global Y position')
title('Global Y vs time')
legend('Euler','Ode45')
grid minor
subplot(5,1,3)
plot(t,x(3,:),t1,x1(:,3),'LineWidth',2)
xlabel('time in seconds')
ylabel('Yaw angle')
title('Yaw angle vs time')
legend('Euler','Ode45')
grid minor
subplot(5,1,4)
plot(t,x(4,:),t1,x1(:,4),'LineWidth',2)
xlabel('time in seconds')
ylabel('Body linear velocity')
title('v vs time')
legend('Euler','Ode45')
grid minor
subplot(5,1,5)
plot(t,x(5,:),t1,x1(:,5),'LineWidth',2)
xlabel('time in seconds')
ylabel('Angular velocity')
title('W vs time')
legend('Euler','Ode45')
grid minor

%% IMM filters
np = 4; % 4 banks of EKFs in IMM filter
y = y_gt + 1*sqrt(r)*randn(m,length(t)); % measurement
X_hat = zeros(n,length(t)); % Combined state estimate
P_hat = zeros(n,n,length(t)); % combined state covariaNCE

% Inititalizations
weight = ones(np,1)/np; % Initial weight
Weight = zeros(np,length(t)); % Weight vector
X_mix_init = zeros(n,np,length(t)); % Mixed states for each filter through out the entire time
X_hat_bank = zeros(n,np,length(t)); % states output for each filter through out the entire time -- column 1 - filter 1, column 2 - filter 2 ...
P_plus_bank(:,:,1) = P0; % 1st filter covariance 
P_plus_bank(:,:,2) = P0; % 2nd filter covariance 
P_plus_bank(:,:,3) = P0; % 3rd filter covariance 
P_plus_bank(:,:,4) = P0; % 4th filter covariance 
likelihood = zeros(1,np,length(t)); % Likelihood of measurement based on filters
P_std = zeros(n, length(t)); % Std deviation
for i = 1:length(t)
    % finding the interacting weights and normalizing sum
    weight_inter = weight.*P_tran; % Interacting model weights

    normalizing_sum = sum(weight_inter);  % normalizing factor

    % Setting the condition in case the weights are too small
    for ii = 1:np
        for j = 1:np
            if normalizing_sum(j) > 1e-20
                weight_inter(ii,j) = weight_inter(ii,j) ./ normalizing_sum(j);
            else
                normalizing_sum(j) = 0;
                weight_inter(ii,j) = 0;
            end
        end
    end
    
    % Mixed initial states and covarinnces
   
    for j = 1:np
        mixed_state = zeros(n, 1);
        cov_mixed = zeros(n,n);
        for ii = 1:np
            % mixed_state = mixed_State + w_ij * xi
            mixed_state = mixed_state + weight_inter(ii,j) * X_hat_bank(:,ii,i);
        end
        X_mix_init(:,j,i) = mixed_state;
        for ii = 1:np
            % Mixed cov = cov_mixed + weight*(filter_cov + covaraince of mixed state)
            error = X_hat_bank(:,ii,i) - X_mix_init(:,j,i);  % Mixed state error 
            cov_mixed = cov_mixed + weight_inter(ii,j)*(P_plus_bank(:,:,ii) + error*error' );
        end
        P_mixed(:,:,j) = cov_mixed;
    end

    % Output 
    y_op = y(:,i);
 
    % Individual filters -- states and covariance calculation
    for j = 1:np
        [X_hat_bank(:,j,i), P_plus_bank(:,:,j), likelihood(:,j,i)] = EKF(K(:,j), P(:,j), X_mix_init(:,j,i), y_op, P_mixed(:,:,j), U(:,i), dt, Q,R,H);
    end

   % Model Probability update
   weight = weight.*likelihood(:,:,i)';
   weight = weight / (sum(weight));

   Weight(:,i) = weight; 

   % Conditional mean and variance
   k = zeros(n,1);
   P_plus_dum = zeros(n);
   for j = 1:np
        k = k + weight(j)*X_hat_bank(:,j,i); 
   end

%    X_hat(:,i) = weight(1)*X_hat_bank(:,1,i) + weight(2)*X_hat_bank(:,2,i) + weight(3)*X_hat_bank(:,3,i) + weight(4)*X_hat_bank(:,4,i) ;
   X_hat(:,i) = k; % Updated states 
   e1 = X_hat_bank(:,1,i)-X_hat(:,i);
   e2 = X_hat_bank(:,2,i)-X_hat(:,i);
   e3 = X_hat_bank(:,3,i)-X_hat(:,i);
   e4 = X_hat_bank(:,4,i)-X_hat(:,i);
   for j = 1:np
      P_plus_dum = P_plus_dum + weight(j)*((X_hat_bank(:,j,i)-X_hat(:,i))...
          *(X_hat_bank(:,j,i)-X_hat(:,i))'+ P_plus_bank(:,:,j));   
       
   end
   P_hat(:,:,i) = P_plus_dum; % Updated covariance
   P_std(:,i) = (diag(P_hat(:,:,i))).^0.5; % Std deviation
%    P_hat(:,:,i) = weight(1)*(e1*e1' + P_plus_bank(:,:,1)) + weight(2)*(e2*e2' + P_plus_bank(:,:,2)) + weight(3)*(e3*e3' + P_plus_bank(:,:,3)) + weight(4)*(e4*e4' + P_plus_bank(:,:,4));
end

figure(2)
subplot(5,1,1)
plot(t,-3*P_std(1,:), 'b', t, X_hat(1,:)-x1(:,1)','r', t, 3*P_std(1,:),'b','LineWidth',2)
xlabel('time in seconds')
ylabel('Error in horizontal position')
grid minor
subplot(5,1,2)
plot(t,-3*P_std(2,:), 'b', t, X_hat(2,:)-x1(:,2)','r', t, 3*P_std(2,:),'b','LineWidth',2)
xlabel('time in seconds')
ylabel('Error in vertical position')
grid minor
subplot(5,1,3)
plot(t,-3*P_std(3,:), 'b', t, X_hat(3,:)-x1(:,3)','r', t, 3*P_std(3,:),'b','LineWidth',2)
xlabel('time in seconds')
ylabel('Error in yaw angle')
grid minor
subplot(5,1,4)
plot(t,-3*P_std(4,:), 'b', t, X_hat(4,:)-x1(:,4)','r', t, 3*P_std(4,:),'b','LineWidth',2)
xlabel('time in seconds')
ylabel('Error in velocity')
grid minor
subplot(5,1,5)
plot(t,-3*P_std(5,:), 'b', t, X_hat(5,:)-x1(:,5)','r', t, 3*P_std(5,:),'b','LineWidth',2)
xlabel('time in seconds')
ylabel('Error in angular velocity')
grid minor

figure(3)
plot(t, Weight, 'linewidth',2)
xlabel('time in seconds')
ylabel('weights')
legend('Model1', 'Model 2', 'Model 3', 'Model 4')
grid minor
%%
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
KG = P_minus*H'/(H*P_minus*H' + R + 0*1e-6*eye(size(R))); % Add regularization to prevent singular matrices

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
x4dot = Pv*(1-K1)*Input(1);
x5dot = Pw*(1-M1)*Input(2);

xdot = [x1dot x2dot x3dot x4dot x5dot]';
end

function xdot = odefun(t1,state,t,U, dt, q, M, P, ne)
% Model 1 
x = state;
upsilon = [0 0;
    0 0;
    0 0;
    M(1) 0;
    0 M(2)];
Input1 = interp1(t,U(1,:),t1);
Input2 = interp1(t,U(2,:),t1);
K1 = M(1); M1 = M(2);
Pv = P(1); Pw = P(2);
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