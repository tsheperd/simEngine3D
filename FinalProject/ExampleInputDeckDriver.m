%% simEngine3D_InputDeckTest Driver Function
clear; close all; clc;

% Timer
tic;

% Booleans
SAVE_PLOTS = 0;


%% Add library of functions to path
addpath('./Functions');
addpath('./Functions/Constraints');


%% Solver
% Initialize class for the simEngine3D
simulation = simEngine3D;

% Read the input deck
simulation.ReadInputDeck("TestInputDeck.mdl");

% Initialize the Solver
simulation.initializeSolver(0, 0.005, 10, 1e-3)

% Check the GCons information
simulation.N_GCons
simulation.N_CCons
simulation.N_CCons_GCons
simulation.N_CCons_GCons_tot


% End timer
toc;