% ME 8843
% Wearable Robotics
% Final Project data analysis
%
% This model generates data of musculotendon (MT) force
% based on exoskeleton spring stiffness and muscle atrophy
%% Clearing variables and figures
close all
clear
clc
tic
%% Variables

grav_range = 2:2:24;
f_max_range = 6000 .* [0.7, 0.85, 1, 1.15, 1.3];
v_max_range = -0.45 .* [0.5, 0.75, 1, 1.25, 1.5];


%% 3 a.

W_mtu_net = zeros(1, length(phase_range));
W_mtu_pos = zeros(1, length(phase_range));
W_mtu_neg = zeros(1, length(phase_range));

for i = 1:length(phase_range)
    
    phase = phase_range(i);
    file_id = strcat('phase_',num2str(phase),'_stiff_180000.mat');

    if exist(file_id, 'file') == 2          % Checking if file exists      

        d = load(file_id);          % Loads data from file
        t = d.data.time;        
        data = d.data.data;
        t0 = t>=1.4 & t<=1.8;
        
        P = data(t0,12);
        W_mtu_net(i) = sum(P);
        W_mtu_pos(i) = sum(P(P>0));
        W_mtu_neg(i) = sum(P(P<0));
       
    end
end
toc