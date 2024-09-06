%% Plots of the input and output
generate_input_output_data
f=1;
figure(f)
subplot(2,2,1)

plot(t,Input1(1,:),t,Input1(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 1')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,2)

plot(t,Output1(1,:),t,Output1(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 1')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor
subplot(2,2,3)

plot(t,Input2(1,:),t,Input2(2,:),'LineWidth',2)

xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 1')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,4)

plot(t,Output2(1,:),t,Output2(2,:),'LineWidth',2)

xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 1')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor

%% Model 2 plots
f= f+1;
figure(f)
subplot(2,2,1)
plot(t2,Input3(1,:),t2,Input3(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 2')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,2)
plot(t2,Output3(1,:),t2,Output3(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 2')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor
subplot(2,2,3)
plot(t,Input4(1,:),t,Input4(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 2')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,4)
plot(t,Output4(1,:),t,Output4(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 2')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor

%% Model 3 Plots
f= f+1;
figure(f)
subplot(2,2,1)
plot(t2,Input5(1,:),t2,Input5(2,:),'LineWidth',2)

xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 3')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,2)
hold on
plot(t2,Output5(1,:),t2,Output5(2,:),'LineWidth',2)

xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 3')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor
subplot(2,2,3)

plot(t,Input6(1,:),t,Input6(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 3')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,4)

plot(t,Output6(1,:),t,Output6(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 3')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor

%% Model 4 Plots
f= f+1;
figure(f)
subplot(2,2,1)
plot(t2,Input7(1,:),t2,Input7(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Input signal - Model 4')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,2)
plot(t2,Output7(1,:),t2,Output7(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 4')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor
subplot(2,2,3)
plot(t,Input8(1,:),t,Input8(2,:),'LineWidth',2)
ylabel('Signal')
title('Input signal - Model 4')
legend('Commanded linear vel','Commanded Angular vel','Location','east')
grid minor
subplot(2,2,4)
plot(t,Output8(1,:),t,Output8(2,:),'LineWidth',2)
xlabel('time in seconds')
ylabel('Signal')
title('Output signal - Model 4')
legend('measured linear vel','measured Angular vel','Location','east')
grid minor

% trajectory and angular velocity plots
% model 1
f= f+1;
figure(f)
subplot(2,2,1)
plot(Output1(3,:),Output1(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 1')
grid minor
subplot(2,2,2)
plot(t,Output1(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 1')
grid minor
subplot(2,2,3)
plot(Output2(3,:),Output2(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 1')
grid minor
subplot(2,2,4)
plot(t,Output2(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 1')
grid minor
% Model 2
f= f+1;
figure(f)
subplot(2,2,1)
plot(Output3(3,:),Output3(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 2')
grid minor
subplot(2,2,2)
plot(t2,Output3(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 2')
grid minor
subplot(2,2,3)
plot(Output4(3,:),Output4(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 2')
grid minor
subplot(2,2,4)
plot(t,Output4(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 2')
grid minor

% Model 3 
f= f+1;
figure(f)
subplot(2,2,1)
plot(Output5(3,:),Output5(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 3')
grid minor
subplot(2,2,2)
plot(t2,Output5(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 3')
grid minor
subplot(2,2,3)
plot(Output6(3,:),Output6(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 3')
grid minor
subplot(2,2,4)
plot(t,Output6(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 3')
grid minor
% Model 4
f= f+1;
figure(f)
subplot(2,2,1)
plot(Output7(3,:),Output7(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 4')
grid minor
subplot(2,2,2)
plot(t2,Output7(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 4')
grid minor
subplot(2,2,3)
plot(Output8(3,:),Output8(4,:),'Linewidth',2)
xlabel('Global X coordinate')
ylabel('Global Y coordinate')
title('trajectory for model 4')
grid minor
subplot(2,2,4)
plot(t,Output8(5,:),'Linewidth',2)
xlabel('time in seconds')
ylabel('Yaw angle in rad')
title('Yaw angle for model 4')
grid minor