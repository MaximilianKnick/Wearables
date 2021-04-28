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
atrophy = 1 * [6000; -0.45];

%% Main Loop

for a = 1:length(exo_stiff_range)
    exo_stiff = exo_stiff_range(a);
    for b = 1:length(grav_range)
        grav = grav_range(b);
        for c = 1:size(atrophy, 2)
            fmax = atrophy(1,c);
            vmax = atrophy(2,c);
            try
                load_system('FullHopper_passiveExo.slx');                                %Loading model
                
                % change model parameters
                set_param('FullHopper_passiveExo/stiffness','Value',num2str(exo_stiff));    %Setting exo stiffness in model
                set_param('FullHopper_passiveExo/LoadDynamics/gravity','Value',num2str(grav));    %Setting gravity constant in model
                set_param('FullHopper_passiveExo/Fmax_mus (N)','Value',num2str(fmax));    %Setting muscle parameters in model
                set_param('FullHopper_passiveExo/Vmax_mus (m\s)','Value',num2str(vmax));    %Setting muscle parameters in model
                
%                 fileName = strcat('exoData_stiff_',num2str(stiff),'_grav_',num2str(grav), '_Fmax_',num2str(fmax),'_Vmax_',num2str(vmax), '.mat');       
                name = sprintf('exoData_stiff_%s__grav_%s_Fmax_%s_Vmax_%s.mat', num2str(exo_stiff),num2str(grav),num2str(fmax),num2str(vmax));
                set_param('FullHopper_passiveExo/To File','FileName', name);  %Setting File Name                
                
                simout{a} = sim('FullHopper_passiveExo.slx');                            %Simulates model
                close_system('FullHopper_passiveExo.slx',0);                             %Closes model
                
            catch err
                err.identifier
                'caught error'
            end
        end
    end
    
end

toc