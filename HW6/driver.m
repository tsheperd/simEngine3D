%% Generic Driver Function
clear; close all; clc;

% addpath(genpath('/Functions'));
% addpath(genpath('/Functions/Constraints'));
addpath('./Functions');
addpath('./Functions/Constraints');

simulation = simEngine3D;

simulation.ReadInputDeck("inputDeck.mdl");
simulation.KinematicSolver();

