% ME 8843
% Wearable Robotics
% Final Project data analysis
%
% This model generates data of musculotendon (MT) force
% based on exoskeleton spring stiffness and muscle atrophy

% ------------ Sections To Run --------------

% Nick - 1, 2, 6
% Ishrat/Rish - 3, 7, 9, 10
% Nolan - 4, 5, 8

% ------------ *************** --------------

%% Clearing variables and figures
close all
clear
clc

%% Variables

exo_stiff_range = linspace(50000,200000,5);
grav_range = 2:1:14;
fmax_range = 6000 * [0.7 0.85 1 1.15 1.3];
vmax_range = -0.45 * [0.5 0.75 1 1.25 1.5];
act_range = linspace(0.5,1,5);


%% 1. 

% Work - v_max, f_max, act Nominal ---------- Nick

vmax = vmax_range(3); % -0.45

W_mtu_net = zeros(5, length(grav_range));
W_mtu_pos = zeros(5, length(grav_range));
W_mtu_neg = zeros(5, length(grav_range));

for j = 1:length(grav_range)
    
    grav = grav_range(j);
    
    for c = 1:length(exo_stiff_range)
        
        exo_stiff = exo_stiff_range(c);
        
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
            W_mtu_net(c, j) = sum(P);
            W_mtu_pos(c, j) = sum(P(P>0));
            W_mtu_neg(c, j) = sum(P(P<0));

        end
    end
end

figure(1)
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

%% 2. 

% P_met - stiffness vs. vmax (fixed g = 10)  ------------ Nick
P_met = zeros(length(vmax_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
grav = 10;
f_max = 6000;

for c = 1:length(vmax_range)
    
    vmax = vmax_range(c);
    
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

            P_met(j, c) = trapz(t, integrand)/T_stim;

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
title('Average Metabolic Rate - v_m_a_x changing')
%set(gca,'FontSize',20)



%% 3. 

% P_met - stiffness vs. fmax (fixed g = 10) -------------- Rish

P_met = zeros(length(fmax_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
grav = 10;
vmax = -0.45;

for c = 1:length(fmax_range)
    
    fmax = fmax_range(c);
    
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

            P_met(j, c) = trapz(t, integrand)/T_stim;

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


%% 4. 

% P_met - stiffness vs. act (fixed g = 10) ----------- Nolan

P_met = zeros(length(act_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
grav = 10;
vmax = -0.45;
fmax = 6000;

for c = 1:length(act_range)
    
    Act = act_range(c);
    
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

            P_met(j, c) = trapz(t, integrand)/T_stim;

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
title('Average Metabolic Rate - Activation changing')
%set(gca,'FontSize',20)

%% 5. 

% P_met - stiffness vs. gravity (nominal fmax/vmax/act)  -------- Nolan

P_met = [];
T_stim = 20;

mass = 35;
vmax = -0.45;
fmax = 6000;

for c = 1:length(grav_range)
    
    grav = grav_range(c);
    
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

            P_met(j, c) = trapz(t, integrand)/T_stim;

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
title('Average Metabolic Rate - Gravity changing')
% set(gca,'FontSize',20)

%% 6. 

% P_met - stiffness vs. vmax vs grav ----------- Nick 

P_met = zeros(length(vmax_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
fmax = 6000;

for k = 1:length(grav_range)
    grav = grav_range(k);
    
    for c = 1:length(vmax_range)
        
        vmax = vmax_range(c);
        
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
                
                integrand = (fmax * act * abs(vmax).* phi)./ mass;
                
                P_met(j, c, k) = trapz(t, integrand)/T_stim;
                
            end
        end
    end
    
end


[exo_stiff_range,vmax_range] = meshgrid(exo_stiff_range,vmax_range); % Define mesh

% Generate P_met Colormaps for various g's
for c = 1:length(grav_range)
    grav = grav_range(c);
    figure(c)
    c = contourf(exo_stiff_range,vmax_range, P_met(:,:,c) , 20,'edgecolor', 'none');
    h = colorbar;
    ylabel(h,'Contour Values')
    xlabel('Exo Stiffness (N/m)')
    ylabel('vmax (m/s)')
    h.Label.String = 'P_m_e_t (W/kg)';
    title(sprintf('Average Metabolic Rate - v_m_a_x vs. Exo Stiffness , g = %s m/s^2', num2str(grav) ))
    %set(gca,'FontSize',20)
    saveas(gcf,sprintf('Pmet_vmax_%s.png', num2str(c+1))) %save figure
end

%% 7. 

% P_met - stiffness vs. fmax vs grav ----------- Rish

P_met = zeros(length(fmax_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
vmax = -0.45;
Act = 1;

for k = 1:length(grav_range)
    grav = grav_range(k);
    
    for c = 1:length(fmax_range)
        
        fmax = fmax_range(c);
        
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
                
                P_met(j, c, k) = trapz(t, integrand)/T_stim;
                
            end
        end
    end
    
end


[exo_stiff_range,fmax_range] = meshgrid(exo_stiff_range,fmax_range); % Define mesh

% Generate P_met Colormaps for various g's
for c = 1:length(grav_range)
    grav = grav_range(c);
    figure(c)
    c = contourf(exo_stiff_range,fmax_range, P_met(:,:,c) , 20,'edgecolor', 'none');
    h = colorbar;
    ylabel(h,'Contour Values')
    xlabel('Exo Stiffness (N/m)')
    ylabel('fmax (N)')
    h.Label.String = 'P_m_e_t (W/kg)';
    title(sprintf('Average Metabolic Rate - f_m_a_x vs. Exo Stiffness , g = %s m/s^2', num2str(grav) ))
    %set(gca,'FontSize',20)
    saveas(gcf,sprintf('Pmet_fmax_%s.png', num2str(c+1))) %save figure
end

%% 8. 

% P_met - stiffness vs. act vs grav ----------- Nolan

P_met = zeros(length(act_range), length(exo_stiff_range));

T_stim = 20;

mass = 35;
vmax = -0.45;
fmax = 6000;

for k = 1:length(grav_range)
    grav = grav_range(k);
    
    for c = 1:length(act_range)
        
        Act = act_range(c);
        
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
                
                P_met(j, c, k) = trapz(t, integrand)/T_stim;
                
            end
        end
    end
    
end

[exo_stiff_range,act_range] = meshgrid(exo_stiff_range,act_range); % Define mesh

% Generate P_met Colormaps for various g's
for c = 1:length(grav_range)
    grav = grav_range(c);
    figure(c)
    c = contourf(exo_stiff_range,act_range, P_met(:,:,c) , 20,'edgecolor', 'none');
    h = colorbar;
    ylabel(h,'Contour Values')
    xlabel('Exo Stiffness (N/m)')
    ylabel('Activation')
    h.Label.String = 'P_m_e_t (W/kg)';
    title(sprintf('Average Metabolic Rate - Activation vs. Exo Stiffness, g = %s m/s^2', num2str(grav) ))
    %set(gca,'FontSize',20)
    saveas(gcf,sprintf('Pmet_act_%s.png', num2str(c+1))) %save figure
end


%% 9.

% Mars - Ideal Stiffness Plot ----------- Rish

ideal_stiff_jump = zeros(length(fmax_range), length(vmax_range));
ideal_stiff_land = zeros(length(fmax_range), length(vmax_range));
ideal_stiff_minWork = zeros(length(fmax_range), length(vmax_range));

for a = 1:length(fmax_range)
    
    fmax = fmax_range(a);
    
    for b = 1:length(vmax_range)
        
        vmax = vmax_range(b);
        
        W_mtu_net = zeros(1, length(exo_stiff_range));
        W_mtu_pos = zeros(1, length(exo_stiff_range));
        W_mtu_neg = zeros(1, length(exo_stiff_range));
        
        for c = 1:length(exo_stiff_range)
            
            exo_stiff = exo_stiff_range(c);
            
            fid = sprintf('exoData_Mars_fmax_%s_vmax_%s_stiff_%s.mat', num2str(fmax),num2str(vmax),num2str(exo_stiff));
            
            if exist(fid, 'file') == 2          % Checking if file exists
                
                d = load(fid);          % Loads data from file
                t = d.data.time;
                data = d.data.data;
                u = data(:,1);
                crossings = find(diff(u) > 0);
                start = crossings(end-1);
                last = crossings(end);
                
                P = data(start:last,20);
                W_mtu_net(c) = sum(P);
                W_mtu_pos(c) = sum(P(P>0));
                W_mtu_neg(c) = sum(P(P<0));
            end
        end
        
        [~, ind] = max(W_mtu_pos);
        ideal_stiff_jump(a,b) = exo_stiff_range(ind);
        
        [~, ind] = min(W_mtu_neg);
        ideal_stiff_land(a,b) = exo_stiff_range(ind);
        
        [~, ind] = min(abs(W_mtu_net));
        ideal_stiff_minWork(a,b) = exo_stiff_range(ind);
        
    end
end

X = fmax_range;
Y = vmax_range;

figure(1)
Z = ideal_stiff_jump;
s = surf(X,Y,Z,'FaceAlpha',0.5);
shading interp
s.EdgeColor = [0.1,0.1,0.1];
xlabel('F_m_a_x (N)');
ylabel('v_m_a_x (m/s)');
zlabel('Ideal Exo Stiffness (N/m)');
title('Ideal Exo Stiffness for Mars with Varied Atrophy [Jumping]')

figure(2)
Z = ideal_stiff_land;
s = surf(X,Y,Z,'FaceAlpha',0.5);
shading interp
s.EdgeColor = [0.1,0.1,0.1];
xlabel('F_m_a_x (N)');
ylabel('v_m_a_x (m/s)');
zlabel('Ideal Exo Stiffness (N/m)');
title('Ideal Exo Stiffness for Mars with Varied Atrophy [Landing]')

figure(3)
Z = ideal_stiff_minWork;
s = surf(X,Y,Z,'FaceAlpha',0.5);
shading interp
s.EdgeColor = [0.1,0.1,0.1];
xlabel('F_m_a_x (N)');
ylabel('v_m_a_x (m/s)');
zlabel('Ideal Exo Stiffness (N/m)');
title('Ideal Exo Stiffness for Mars with Varied Atrophy [Minimum Total Work]')


%% 10.

% Moon - Ideal Stiffness Plot ----------- Rish

ideal_stiff_jump = zeros(length(fmax_range), length(vmax_range));
ideal_stiff_land = zeros(length(fmax_range), length(vmax_range));
ideal_stiff_minWork = zeros(length(fmax_range), length(vmax_range));

for a = 1:length(fmax_range)
    
    fmax = fmax_range(a);
    
    for b = 1:length(vmax_range)
        
        vmax = vmax_range(b);
        
        W_mtu_net = zeros(1, length(exo_stiff_range));
        W_mtu_pos = zeros(1, length(exo_stiff_range));
        W_mtu_neg = zeros(1, length(exo_stiff_range));
        
        for c = 1:length(exo_stiff_range)
            
            exo_stiff = exo_stiff_range(c);
            
            fid = sprintf('exoData_Moon_fmax_%s_vmax_%s_stiff_%s.mat', num2str(fmax),num2str(vmax),num2str(exo_stiff));
            
            if exist(fid, 'file') == 2          % Checking if file exists
                
                d = load(fid);          % Loads data from file
                t = d.data.time;
                data = d.data.data;
                u = data(:,1);
                crossings = find(diff(u) > 0);
                start = crossings(end-1);
                last = crossings(end);
                
                P = data(start:last,20);
                W_mtu_net(c) = sum(P);
                W_mtu_pos(c) = sum(P(P>0));
                W_mtu_neg(c) = sum(P(P<0));
            end
        end
        
        [~, ind] = max(W_mtu_pos);
        ideal_stiff_jump(a,b) = exo_stiff_range(ind);
        
        [~, ind] = min(W_mtu_neg);
        ideal_stiff_land(a,b) = exo_stiff_range(ind);
        
        [~, ind] = min(abs(W_mtu_net));
        ideal_stiff_minWork(a,b) = exo_stiff_range(ind);
        
    end
end

X = fmax_range;
Y = vmax_range;

figure(1)
Z = ideal_stiff_jump;
s = surf(X,Y,Z,'FaceAlpha',0.5);
shading interp
s.EdgeColor = [0.1,0.1,0.1];
xlabel('F_m_a_x (N)');
ylabel('v_m_a_x (m/s)');
zlabel('Ideal Exo Stiffness (N/m)');
title('Ideal Exo Stiffness for the Moon with Varied Atrophy [Jumping]')

figure(2)
Z = ideal_stiff_land;
s = surf(X,Y,Z,'FaceAlpha',0.5);
shading interp
s.EdgeColor = [0.1,0.1,0.1];
xlabel('F_m_a_x (N)');
ylabel('v_m_a_x (m/s)');
zlabel('Ideal Exo Stiffness (N/m)');
title('Ideal Exo Stiffness for the Moon with Varied Atrophy [Landing]')

figure(3)
Z = ideal_stiff_minWork;
s = surf(X,Y,Z,'FaceAlpha',0.5);
shading interp
s.EdgeColor = [0.1,0.1,0.1];
xlabel('F_m_a_x (N)');
ylabel('v_m_a_x (m/s)');
zlabel('Ideal Exo Stiffness (N/m)');
title('Ideal Exo Stiffness for the Moon with Varied Atrophy [Minimum Total Work]')

