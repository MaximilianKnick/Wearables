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

exo_stiff_range = linspace(50000,200000,5);
grav_range = 2:1:14;
fmax_range = 6000 * [0.7 0.85 1 1.15 1.3];
vmax_range = -0.45 * [0.5 0.75 1 1.25 1.5];
act_range = linspace(0.5,1,5);


%% Work - v_max, f_max, act Nominal ---------- Nick

vmax = vmax_range(3); % -0.45

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
            u = data(:,1);
            crossings = find(diff(u) > 0);
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
legend(h, {'50 kN/m','87.5 kN/m','125 kN/m','162.5 kN/m','200 kN/m'},'Location','southwest');
legend('boxoff')
ylabel('Net Mechanical Work (J)');

subplot(3,1,2);
bar(grav_range, W_mtu_pos);
ylabel('Positive Mechanical Work (J)');

subplot(3,1,3);
bar(grav_range, W_mtu_neg); 
xlabel ('Gravity [m/s^2]')
ylabel('Negative Mechanical Work (J)');

sgtitle('Work vs. Gravity at Multiple Exo Stiffnesses')
toc

%% P_met - stiffness vs. vmax (fixed g = 10)  ------------ Nick
P_met = zeros(length(vmax_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
grav = 10;
f_max = 6000;

for i = 1:length(vmax_range)
    
    vmax = vmax_range(i);
    
    for j = 1:length(exo_stiff_range)
        
        exo_stiff = exo_stiff_range(j);
        
        fid = sprintf('exoData_grav_%s_vmax_%s_stiff_%s.mat', num2str(grav),num2str(vmax),num2str(exo_stiff));

        if exist(fid, 'file') == 2          % Checking if file exists      

            d = load(fid);          % Loads data from file
            t = d.data.time;        
            data = d.data.data;
            
            act = data(:,4);
            v_m = data(:,11);   % USED v_m - Check if v_mtu(10) works better
            phi = zeros(length(v_m), 1);
            vpos = v_m > 0;
            vneg = ~vpos;
            
            
            phi(vpos) = 0.01 - 0.11.*(v_m(vpos)./vmax) + 0.06*exp(23*v_m(vpos)./vmax);
            phi(vneg) = 0.23 - 0.16*exp(-8*v_m(vneg)./vmax);
            
            integrand = (f_max * act * abs(vmax).* phi)./ mass;

            P_met(j, i) = trapz(t, integrand)/T_stim;

        end
    end
end

[exo_stiff_range,vmax_range] = meshgrid(exo_stiff_range,vmax_range);

c = contourf(exo_stiff_range,vmax_range, P_met , 20,'edgecolor', 'none');
h = colorbar;
ylabel(h,'Contour Values')
xlabel('Exo Stiffness (N/m)')
ylabel('v_m_a_x (m/s^2)')
h.Label.String = 'P_m_e_t (W/kg)';
title('Average Metabolic Rate - v_m_a_x changing', 'FontSize', 20)
%set(gca,'FontSize',20)



%% P_met - stiffness vs. fmax (fixed g = 10) -------------- Rish

P_met = zeros(length(fmax_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
grav = 10;
vmax = -0.45;

for i = 1:length(fmax_range)
    
    fmax = fmax_range(i);
    
    for j = 1:length(exo_stiff_range)
        
        exo_stiff = exo_stiff_range(j);
        
        fid = sprintf('exoData_grav_%s_fmax_%s_stiff_%s.mat', num2str(grav),num2str(fmax),num2str(exo_stiff));

        if exist(fid, 'file') == 2          % Checking if file exists      

            d = load(fid);          % Loads data from file
            t = d.data.time;        
            data = d.data.data;
            
            act = data(:,4);
            v_m = data(:,11);   % USED v_m - Check if v_mtu(10) works better
            phi = zeros(length(v_m), 1);
            vpos = v_m > 0;
            vneg = ~vpos;
            
            
            phi(vpos) = 0.01 - 0.11.*(v_m(vpos)./vmax) + 0.06*exp(23*v_m(vpos)./vmax);
            phi(vneg) = 0.23 - 0.16*exp(-8*v_m(vneg)./vmax);
            
            integrand = (fmax * act * abs(vmax).* phi)./ mass;

            P_met(j, i) = trapz(t, integrand)/T_stim;

        end
    end
end

[exo_stiff_range,fmax_range] = meshgrid(exo_stiff_range,fmax_range);

c = contourf(exo_stiff_range,fmax_range, P_met , 20,'edgecolor', 'none');
h = colorbar;
ylabel(h,'Contour Values')
xlabel('Exo Stiffness (N/m)')
ylabel('f_m_a_x (N)')
h.Label.String = 'P_m_e_t (W/kg)';
title('Average Metabolic Rate - f_m_a_x changing')
% set(gca,'FontSize',20)


%% P_met - stiffness vs. act (fixed g = 10) ----------- Nolan

P_met = zeros(length(act_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
grav = 10;
vmax = -0.45;
fmax = 6000;

for i = 1:length(act_range)
    
    Act = act_range(i);
    
    for j = 1:length(exo_stiff_range)
        
        exo_stiff = exo_stiff_range(j);
        
        fid = sprintf('exoData_grav_%s_act_%s_stiff_%s.mat', num2str(grav),num2str(Act),num2str(exo_stiff));

        if exist(fid, 'file') == 2          % Checking if file exists      

            d = load(fid);          % Loads data from file
            t = d.data.time;        
            data = d.data.data;
            
            act = data(:,4);
            v_m = data(:,11);   % USED v_m - Check if v_mtu(10) works better
            phi = zeros(length(v_m), 1);
            vpos = v_m > 0;
            vneg = ~vpos;
            
            
            phi(vpos) = 0.01 - 0.11.*(v_m(vpos)./vmax) + 0.06*exp(23*v_m(vpos)./vmax);
            phi(vneg) = 0.23 - 0.16*exp(-8*v_m(vneg)./vmax);
            
            integrand = (fmax * act * abs(vmax).* phi)./ mass;

            P_met(j, i) = trapz(t, integrand)/T_stim;

        end
    end
end

[exo_stiff_range,act_range] = meshgrid(exo_stiff_range,act_range);

c = contourf(exo_stiff_range,act_range, P_met , 20,'edgecolor', 'none');
h = colorbar;
ylabel(h,'Contour Values')
xlabel('Exo Stiffness (N/m)')
ylabel('Activation')
h.Label.String = 'P_m_e_t (W/kg)';
title('Average Metabolic Rate - Activation changing', 'FontSize', 20)
%set(gca,'FontSize',20)

%% P_met - stiffness vs. gravity (nominal fmax/vmax/act)  -------- Nolan

% P_met = zeros(length(grav_range), length(exo_stiff_range));

P_met = [];
T_stim = 20;

mass = 35;
vmax = -0.45;
fmax = 6000;

for i = 1:length(grav_range)
    
    grav = grav_range(i);
    
    for j = 1:length(exo_stiff_range)
        
        exo_stiff = exo_stiff_range(j);
        
        fid = sprintf('exoData_grav_%s_act_%s_stiff_%s.mat', num2str(grav),num2str(1),num2str(exo_stiff));

        if exist(fid, 'file') == 2          % Checking if file exists      

            d = load(fid);          % Loads data from file
            t = d.data.time;        
            data = d.data.data;
            
            act = data(:,4);
            v_m = data(:,11);   % USED v_m - Check if v_mtu(10) works better
            phi = zeros(length(v_m), 1);
            vpos = v_m > 0;
            vneg = ~vpos;
            
            
            phi(vpos) = 0.01 - 0.11.*(v_m(vpos)./vmax) + 0.06*exp(23*v_m(vpos)./vmax);
            phi(vneg) = 0.23 - 0.16*exp(-8*v_m(vneg)./vmax);
            
            integrand = (fmax * act * abs(vmax).* phi)./ mass;

            P_met(j, i) = trapz(t, integrand)/T_stim;

        end
    end
end

[exo_stiff_range,grav_range] = meshgrid(exo_stiff_range,grav_range);

c = contourf(exo_stiff_range,grav_range, P_met' , 20,'edgecolor', 'none');
h = colorbar;
ylabel(h,'Contour Values')
xlabel('Exo Stiffness (N/m)')
ylabel('Gravity (m/s^2)')
h.Label.String = 'P_m_e_t (W/kg)';
title('Average Metabolic Rate - Gravity changing', 'FontSize', 20)
% set(gca,'FontSize',20)






