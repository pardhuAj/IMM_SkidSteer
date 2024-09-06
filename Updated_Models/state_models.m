generate_input_output_data
n = 5; % Five states -- Global X coordinate, Global Y coordinate, yaw angle, body velocity and angular velocity. we have a nonlinear system 
q = 0.01; % Process noise covariance with wheel angular velocities
r = 0.001; % measurement noise covariance
m = 2; % no of measurements 
dt = t(2)-t(1); % time step in secs
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
    x_0 =  Output1(:,1);
elseif Model_num == 2
    U = Input2;
    K_M = K(:,1);
    Pvw = P(:,1);
    x_0 = Output2(:,1);
elseif Model_num == 3
    U = Input3;
    K_M = K(:,2);
    t = t2;
    Pvw = P(:,2);
    x_0 = Output3(:,1);
elseif Model_num == 4
    U = Input4;
    K_M = K(:,3);
    Pvw = P(:,2);
    x_0 = Output4(:,1);
elseif Model_num == 5
    U = Input5;
    K_M = K(:,3);
    t = t2;
    Pvw = P(:,3);
    x_0 = Output5(:,1);
elseif Model_num == 6
    U = Input6;
    K_M = K(:,3);
    Pvw = P(:,3);
    x_0 = Output6(:,1);
elseif Model_num == 7
    U = Input7;
    K_M = K(:,4);
    t = t2;
    Pvw = P(:,4);
    x_0 = Output7(:,1);
else
    U = Input8;
    K_M = K(:,4);
    Pvw = P(:,4);
    x_0 = Output8(:,1);
end
B = [0 0;
 0 0;
 0 0;
 Pvw(1)*(1-K_M(1)) 0;
 0 Pvw(2)*(1-K_M(2))]; % Input matrix
upsilon = B;
H = [0 0 0 1 0;
 0 0 0 0 1]; % Measurement matrix
Q = upsilon*q*dt*eye(2)*upsilon'; % process noise covariance
R = eye(m)*r; % Measurement noise covariance

ne = 1; % switch for effects of noise...
% 1 --> process noise is considered
% anyything else --> no effects of noise

x = zeros(n,length(t));
x(:,1) = x_0;
for i = 1:length(t)-1
    xdot = dxdt(x(:,i),U(:,i),K_M,Pvw);
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
[t1,x1]=ode45(@(t1,x1) odefun(t1,x1,t,U, dt, q, K_M,Pvw,ne),t,x_0);
f=1;
figure(f)
plot(x1(:,1),x1(:,2),'Linewidth',2)
xlabel('Longitudinal position')
ylabel('Lateral position')
title('trajectory')
grid minor
% Plots compating Euler and ode45
f= f+1;
figure(f)
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
function xdot = dxdt(state,Input,M,P)
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

function xdot = odefun(t1,state,t,U,dt,q,M,P,ne)
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
