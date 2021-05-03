% ME 8843
% Wearable Robotics
% Final Project data generation
%
% This model generates data of musculotendon (MT) force
% based on exoskeleton spring stiffness and muscle atrophy
%% Clearing variables and figures ------- 1st
close all
clear
clc
%% Variables ------------ 2nd

% ranges of variables for traversing through the for loops
exo_stiff_range = linspace(100,200000,5);
grav_range = 2:2:24;
fmax_range = 6000 * [0.7 0.85 1 1.15 1.3];
vmax_range = -0.45 * [0.5 0.75 1 1.25 1.5];
act_range = linspace(0.5,1,5);

% % test values (comment out and use ^above^ definitions for actual data)
% exo_stiff_range = 180000 * 1;
% grav_range = 9.81 * 1;
% fmax_range = [6000];
% vmax_range = [-0.45];
% act_range = [1];

% values of variables for when they are being held constant over the
% iterations
% fmax = 6000;
% vmax = -0.45;
% act = 1;
% grav = 9.81;

file_name = 'FullHopper_passiveExo_PWM.slx';

%------------- Run F_MAX OR V_MAX OR ACT OR Fixed g ----------------------
%% F_MAX - Ishrat
parfor a = 1:length(grav_range)    
    grav = grav_range(a);
    for b = 1:length(fmax_range)
        fmax = fmax_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system(file_name);                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo_PWM/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo_PWM/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo_PWM/Fmax_mus (N)','Value',num2str(fmax));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_grav_%s_fmax_%s_stiff_%s.mat', num2str(grav),num2str(fmax),num2str(stiff));
                set_param('FullHopper_passiveExo_PWM/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim(file_name);                            %Simulates model
                close_system(file_name,0);                             %Closes model                
            catch err
                err.identifier
                'caught error 1'
            end
        end
    end
end

%% V_MAX - Nicholas
parfor a = 1:length(grav_range)    
    grav = grav_range(a);
    for b = 1:length(vmax_range)
        vmax = vmax_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system(file_name);                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo_PWM/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo_PWM/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo_PWM/Vmax_mus (m\s)','Value',num2str(vmax));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_grav_%s_vmax_%s_stiff_%s.mat', num2str(grav),num2str(vmax),num2str(stiff));
                set_param('FullHopper_passiveExo_PWM/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim(file_name);                            %Simulates model
                close_system(file_name,0);                             %Closes model                
            catch err
                err.identifier
                'caught error 2'
            end
        end
    end
end

%% ACTIVATION - Nolan
parfor a = 1:length(grav_range)    
    grav = grav_range(a);
    for b = 1:length(act_range)
        act = act_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system(file_name);                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo_PWM/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo_PWM/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo_PWM/act_gain','Value',num2str(act));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_grav_%s_act_%s_stiff_%s.mat', num2str(grav),num2str(act),num2str(stiff));
                set_param('FullHopper_passiveExo_PWM/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim(file_name);                            %Simulates model
                close_system(file_name,0);                             %Closes model                
            catch err
                err.identifier
                'caught error 3'
            end
        end
    end
end

%% Fixed g - F_MAX and V_MAX

% grav = 3.7;   % Mercury
% grav = 8.87;  % Venus
grav = 9.81;  % Earth
% grav = 1.62;  % Moon
% grav = 3.72;  % Mars
% grav = 24.8;  % Jupiter
% grav = 10.44; % Saturn
% grav = 8.8;   % Uranus
% grav = 11.15; % Neptune

parfor a = 1:length(fmax_range)    
    fmax = fmax_range(a);
    for b = 1:length(vmax_range)
        vmax = vmax_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system(file_name);                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo_PWM/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo_PWM/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo_PWM/Fmax_mus (N)','Value',num2str(fmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo_PWM/Vmax_mus (m\s)','Value',num2str(vmax));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_fmax_%s_vmax_%s_stiff_%s.mat', num2str(fmax),num2str(vmax),num2str(stiff));
                set_param('FullHopper_passiveExo_PWM/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim(file_name);                            %Simulates model
                close_system(file_name,0);                             %Closes model                
            catch err
                err.identifier
                'caught error 4'
            end
        end
    end
end