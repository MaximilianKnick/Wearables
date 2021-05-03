% ME 8843
% Wearable Robotics
% Final Project data analysis
%
% This model generates data of musculotendon (MT) force
% based on exoskeleton spring stiffness and muscle atrophy
%% Clearing variables and figures
close all
clear all
clc
tic
%% Variables

exo_stiff_range = linspace(100,200000,5);
grav_range = 2:1:14;
fmax_range = 6000 * [0.7 0.85 1 1.15 1.3];
vmax_range = -0.45 * [0.5 0.75 1 1.25 1.5];
act_range = linspace(0.5,1,5);


%% 3 c.

grav_range   = 2:1:14;
exo_stiff_range = linspace(100,200000,5);

vmax = vmax_range(3);

W_mtu_net = zeros(5, length(grav_range));
W_mtu_pos = zeros(5, length(grav_range));
W_mtu_neg = zeros(5, length(grav_range));

for j = 1:length(grav_range)
    
    grav = grav_range(j);
    
    for i = 1:length(exo_stiff_range)
        
        exo_stiff = exo_stiff_range(i);
        
        fid = sprintf('exoData_grav_%s_vmax_%s_stiff_%s.mat', num2str(grav),num2str(vmax),num2str(exo_stiff));

        if exist(fid, 'file') == 2          % Checking if file exists      

            d = load(fid);          % Loads data from file
            t = d.data.time;        
            data = d.data.data;
            act = data(:,1);
            crossings = find(diff(act) > 0);
            start = crossings(end-1);
            last = crossings(end);

            P = data(start:last,20);
            W_mtu_net(i, j) = sum(P);
            W_mtu_pos(i, j) = sum(P(P>0));
            W_mtu_neg(i, j) = sum(P(P<0));

        end
    end
end

figure(3)
subplot(3,1,1);
bar(grav_range, W_mtu_net);
h = bar(grav_range, W_mtu_net);
% legend(h, mat2cell(exo_stiff_range), 'Location','southwest');
legend('boxoff')
ylabel('Net Mechanical Work (J)');

subplot(3,1,2);
bar(grav_range, W_mtu_pos);
ylabel('Positive Mechanical Work (J)');

subplot(3,1,3);
bar(grav_range, W_mtu_neg); 
xlabel ('Gravity [m/s^2]')
ylabel('Negative Mechanical Work (J)');

sgtitle('Q.3 c.')
toc