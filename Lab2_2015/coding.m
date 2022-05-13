clc
clear all

%% loading data
[data_h, header_h] = lab_read_edf('Healthy_control.edf');
[data_p, header_p] = lab_read_edf('Cardiomyopathy.edf');

[N_channels,T] = size(data_h);

%% plotting

fs = 1000;
T_axis = linspace(0,T/fs,T);

figure;
plot(T_axis,data_h(7,:),T_axis,data_p(7,:)); title('v1')
xlim([0 5]); legend('Healthy','Patient')
%% analysing data

PeakDetect2(data_h(7,:), fs)