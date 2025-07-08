generate_input_output_data
n = 2; % Five states -- Global X coordinate, Global Y coordinate, yaw angle, body velocity and angular velocity. we have a nonlinear system 
q = 0.01; % Process noise covariance with wheel angular velocities
r = 0.001; % measurement noise covariance
m = 2;
upsilon = [1 0;
    0 1];
dt = t(2)-t(1);
% Model 1 
K1 = 0.0763; M1 = 0.1134;
% Model 2 
K2 = 0.0784; M2 = 0.1257;
% Model 3 
K3 = 0.0776; M3 = 0.1196;
% Model 4 
K4 = 0.0774; M4 = 0.1216;
K = [K1 K2 K3 K4;
    M1 M2 M3 M4];
B = [1 0;
 0 1]; % Input matrix

Model_num = 1;
if Model_num == 1
    U = Input1;
    K_M = K(:,1);
elseif Model_num == 2
    U = Input2;
    K_M = K(:,1);
elseif Model_num == 3
    U = Input3;
    K_M = K(:,2);
    t = t2;
elseif Model_num == 4
    U = Input4;
    K_M = K(:,3);

elseif Model_num == 5
    U = Input5;
    K_M = K(:,3);
    t = t2;
elseif Model_num == 6
    U = Input6;
    K_M = K(:,3);
elseif Model_num == 7
    U = Input7;
    K_M = K(:,4);
    t = t2;
else
    U = Input8;
    K_M = K(:,4);
end

H = [ 1 0;
 0 1]; % Measurement matrix
Q = upsilon*q*dt*eye(2)*upsilon'; % process noise covariance
R = eye(m)*r; % Measurement noise covariance

ne = 1; % switch for effects of noise...
% 1 --> process noise is considered
% anyything else --> no effects of noise

x = zeros(n,length(t));
x(:,1) = zeros(n,1);
xdo = zeros(n,1);
for i = 1:length(t)-1
    xdot = dxdt(x(:,i),U(:,i),K_M);
    xdo(:,i) = xdot;
    if ne == 1
        upsilon = [K_M(1) 0;
            0 K_M(2)];
        x(:,i+1) = x(:,i) + xdot*dt + upsilon*sqrt(q*dt)*randn(2,1);
    else 
        x(:,i+1) = x(:,i) + xdot*dt;
    end
end
[t1,x1]=ode45(@(t1,x1) odefun(t1,x1,t,U, dt, q, K_M,ne),t,zeros(n,1));

figure
plot(x1(:,1),x1(:,2),'Linewidth',2)
xlabel('Longitudinal position')
ylabel('Lateral position')
title('trajectory')
grid minor
% Plots compating Euler and ode45
figure
subplot(1,2,1)
plot(t,x(1,:),t1,x1(:,1),'LineWidth',2)
xlabel('time in seconds')
ylabel('Global velocity')
title('Global vel vs time')
legend('Euler','Ode45')
grid minor
subplot(1,2,2)
plot(t,x(2,:),t1,x1(:,2),'LineWidth',2)
xlabel('time in seconds')
ylabel('Global Angular velocity')
title('Global Angulat velocity vs time')
legend('Euler','Ode45')
grid minor

function xdot = dxdt(state, Input,M)
% Model 1 
x = state;
K1 = M(1); M1 = M(2);
x4dot = -K1*x(1) + K1*Input(1);
x5dot = -M1*x(2) + M1*Input(2);

xdot = [x4dot x5dot]';
end

function xdot = odefun(t1,state,t,U, dt, q,M,ne)
% Model 1 
x = state;
Input1 = interp1(t,U(1,:),t1);
Input2 = interp1(t,U(2,:),t1);
K1 = M(1); M1 = M(2);
x4dot = -K1*x(1) + K1*Input1;
x5dot = -M1*x(2) + M1*Input2;
upsilon = [K1 0;
    0 M1];
if ne == 1
    xdot = [x4dot x5dot]' + upsilon*sqrt(q*dt)*randn(2,1);
else 
    xdot = [x4dot x5dot]';
end

end
