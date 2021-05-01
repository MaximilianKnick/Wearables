function [work] = MaxJump_Obj(x)
% input x is stiffness
% output work is work
assignin('base','x',x); %https://www.reddit.com/r/matlab/comments/3t12b5/fmincon_and_simulink/cx26s83/
set_param('FullHopper_passiveExo/stiffness', 'Value', num2str(x));
simout = sim('FullHopper_passiveExo'); % simout is output of simulation
work = max(simout.yout.signals.values)*-1;
end

