clc; clear;
generate_input_output_data

weight1 = IMM_Skid_Steer_func(Input1, Output1,t,dt);
figure
plot(t, weight1, 'linewidth',2)
xlabel('time in seconds')
ylabel('weights')
legend('Model1', 'Model 2', 'Model 3', 'Model 4')
grid minor