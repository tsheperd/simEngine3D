%% Generic Driver Function
clear; close all; clc;

% addpath(genpath('/Functions'));
% addpath(genpath('/Functions/Constraints'));
addpath('./Functions');
addpath('./Functions/Constraints');

simulation = simEngine3D;

simulation.ReadInputDeck("inputDeck.mdl");
simulation.KinematicSolver(0,1,1,1e-3);

N_Cons = size(simulation.input.constraints,1)-2;

simulation.input.constraints;
simulation.input.constraints{1};
size(simulation.input.constraints,1);

simulation.t;
simulation.N_t;
simulation.f;
simulation.f_dot;
simulation.f_ddot;

simulation.q
simulation.q_dot
simulation.q_ddot
i=2;
tt = 1;
simulation.q(7*(i-1)+1+3:7*(i-1)+7,tt);
p_i = simulation.q(7*(i-1)+1+3:7*(i-1)+7,tt);
p_i'*p_i - 1;

simulation.Phi_G
simulation.nu_G
simulation.gamma_G
simulation.Jacobian_G