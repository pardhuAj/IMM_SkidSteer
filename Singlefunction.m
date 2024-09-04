clc; clear;
generate_input_output_data % Input and output data arrays
n = 5; % no of states
np = 4; % Models in IMM bank
m = length(Output8(:,1)); % No of measurements is a parameter which can be changed from 1 to 5
Weight = zeros(np,length(t)); % Mode probabalities for each filter
Initial_weight = ones(np,1)/np;
Weight(:,1) = Initial_weight;
dt = 0.1; % time step for discrete process
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
Output8(1,:) = movmean(Output8(1,:),5);
Output8(2,:) = movmean(Output8(2,:),5);
for i = 1:length(t)-1
    Weight(:,i+1) = IMM_SS_single_func(Input8(:,i),Output8(:,i),dt,m,3,Weight(:,i),K);
end
figure
hold on
plot(t, Weight, 'linewidth',2)
xlabel('time in seconds')
ylabel('weights')
legend('Model1', 'Model 2', 'Model 3', 'Model 4','Location','northwest')
grid minor