% ME 8843
% Wearable Robotics
% Final Project data generation
%
% This model generates data of musculotendon (MT) force
% based on exoskeleton spring stiffness and muscle atrophy
%% Clearing variables and figures
close all
clear
clc
tic
%% Variables

% *** need to find reasonable stiffness values and values for muscle atrophy ***
% exo_stiff_range = 180000 * 0.8:0.1:1.2;   % Range of tendon stiffness values
% grav_range = 9.81 * 0.5:0.1:2;
% atrophy = 0.8:0.1:1.2 * [6000; -0.45];
exo_stiff_range = 180000 * 1;   % Range of tendon stiffness values
grav_range = 9.81 * 1;
fmax_range = [6000];
vmax_range = [-0.45];
act_range = [1];
fmax = 6000;
vmax = -0.45;
act = 1;


%% F_MAX
for a = 1:length(grav_range)    
    grav = grav_range(a);
    for b = 1:length(fmax_range)
        fmax = fmax_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system('FullHopper_passiveExo.slx');                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo/Fmax_mus (N)','Value',num2str(fmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/Vmax_mus (m\s)','Value',num2str(vmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/pulse','Amplitude',num2str(act));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_grav_%s_fmax_%s_stiff_%s.mat', num2str(grav),num2str(vmax),num2str(stiff));
                set_param('FullHopper_passiveExo/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim('FullHopper_passiveExo.slx');                            %Simulates model
                close_system('FullHopper_passiveExo.slx',0);                             %Closes model                
            catch err
                err.identifier
                'caught error 1'
            end
        end
    end
end

%% V_MAX
parfor a = 1:length(grav_range)    
    grav = grav_range(a);
    for b = 1:length(vmax_range)
        vmax = vmax_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system('FullHopper_passiveExo.slx');                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo/Fmax_mus (N)','Value',num2str(fmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/Vmax_mus (m\s)','Value',num2str(vmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/u(t)','Amplitude',num2str(act));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_grav_%s_vmax_%s_stiff_%s.mat', num2str(grav),num2str(vmax),num2str(stiff));
                set_param('FullHopper_passiveExo/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim('FullHopper_passiveExo.slx');                            %Simulates model
                close_system('FullHopper_passiveExo.slx',0);                             %Closes model                
            catch err
                err.identifier
                'caught error 2'
            end
        end
    end
end

%% ACTIVATION
parfor a = 1:length(grav_range)    
    grav = grav_range(a);
    for b = 1:length(vmax_range)
        act = act_range(b);
        for c = 1:size(exo_stiff_range)
            stiff = exo_stiff_range(c);
            try
                load_system('FullHopper_passiveExo.slx');                                %Loading model                
                % CHANGE MODEL PARAMETERS
                set_param('FullHopper_passiveExo/stiffness','Value',num2str(stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo/Fmax_mus (N)','Value',num2str(fmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/Vmax_mus (m\s)','Value',num2str(vmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/u(t)','Amplitude',num2str(act));    %Setting muscle parameters in model
                % FILE NAME
                name = sprintf('exoData_grav_%s_act_%s_stiff_%s.mat', num2str(grav),num2str(act),num2str(stiff));
                set_param('FullHopper_passiveExo/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim('FullHopper_passiveExo.slx');                            %Simulates model
                close_system('FullHopper_passiveExo.slx',0);                             %Closes model                
            catch err
                err.identifier
                'caught error 3'
            end
        end
    end
end


toc