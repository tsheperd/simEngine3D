%% simEngine3D_A6P2 Driver Function
clear; close all; clc;


%% Initial state
theta_ini = pi/4*cos(2*0);
A_ini = [0,0,1;...
	 sin(theta_ini), cos(theta_ini), 0;...
	 -cos(theta_ini), sin(theta_ini), 0;];
e0_ini = ((trace(A_ini)+1)/4)^(1/2)
e1_ini = ((A_ini(1,1)-trace(A_ini)+1)/4)^(1/2)
e2_ini = ((A_ini(2,2)-trace(A_ini)+1)/4)^(1/2)
e3_ini = ((A_ini(3,3)-trace(A_ini)+1)/4)^(1/2)
r_ini = A_ini*[2,0,0]'


%% Add library of functions to path
addpath('./Functions');
addpath('./Functions/Constraints');


%% Solver
% Initialize class for the simEngine3D
simulation = simEngine3D;

% Read the input deck
simulation.ReadInputDeck("revJoint.mdl");

% Run the kinematic solver: (t_initial, dt, t_final, tolerance)
simulation.KinematicSolver(0, 1, 0, 1e-10);


%% Output final timestep information
disp("Phi");
simulation.Phi_G
disp("nu");
simulation.nu_G
disp("gamma");
simulation.gamma_G
disp("Jacobian");
simulation.Jacobian_G


