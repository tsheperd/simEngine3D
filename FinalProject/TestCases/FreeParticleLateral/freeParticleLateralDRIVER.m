%% Free Particle Lateral Driver Function
clear; close all; clc;
% Profiler shows the timings, calls, etc.
%profile on

% Timer
tic;

% Booleans
SAVE_PLOTS = 1;


%% Initial state
rho = 1;
R = 1;
m = rho*4/3*pi*R^2

J_xx_bar = 2/5*m*R^2
J_yy_bar = J_xx_bar
J_zz_bar = J_xx_bar

z_0 = -4;
g = -9.81;


%% Add library of functions to path
addpath('./Functions');
addpath('./Functions/Constraints');


%% Solver
% Initialize class for the simEngine3D
simulation = simEngine3D;

% Read the input deck
simulation.ReadInputDeck("freeParticleLateralINPUT.mdl");


% Run the solver: (t_initial, dt, t_final, tolerance)
simulation.DynamicsSolver(0, 0.005, 10, 1e-4);


%% Output final timestep information
disp("Phi");
simulation.Phi_G
disp("nu");
simulation.nu_G
disp("gamma");
simulation.gamma_G
disp("Jacobian");
simulation.Jacobian_G


%% Origin Plots
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,simulation.q(1,:));
plot(simulation.t,simulation.q(2,:));
plot(simulation.t,simulation.q(3,:));
title("Free Particle Lateral: O Global Position");
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t,simulation.q_dot(1,:));
plot(simulation.t,simulation.q_dot(2,:));
plot(simulation.t,simulation.q_dot(3,:));
title("Free Particle Lateral: O Global Velocity");
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t,simulation.q_ddot(1,:));
plot(simulation.t,simulation.q_ddot(2,:));
plot(simulation.t,simulation.q_ddot(3,:));
title("Free Particle Lateral: O Global Acceleration");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'freeParticleLateral_O_Plot.png');
end

% Quaternian Plot
figure;
hold on;
plot(simulation.t,simulation.p(1,:));
plot(simulation.t,simulation.p(2,:));
plot(simulation.t,simulation.p(3,:),'-.');
plot(simulation.t,simulation.p(4,:),'--');
title("Free Particle Lateral: Quaternian Components");
xlabel("t (s)");
ylabel("Quaternian (-)");
legend('1','2','3','4');
ylim([-1,2]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'freeParticleLateral_Quaternian_Plot.png');
end


%% Analytical Solution
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,1/2*10*simulation.t.^2);
plot(simulation.t,z_0+1/2*g*simulation.t.^2);
title("Free Particle Lateral: O Global Position (Analytical)");
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,10*simulation.t);
plot(simulation.t,g*simulation.t);
title("Free Particle Lateral: O Global Velocity (Analytical)");
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,10*ones(size(simulation.t)));
plot(simulation.t,g*ones(size(simulation.t)));
title("Free Particle Lateral: O Global Acceleration (Analytical)");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'freeParticleLateral_O_Plot_Analytical.png');
end

% CONSTRAINED
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,0*simulation.t);
plot(simulation.t,z_0+1/2*g*simulation.t.^2);
title("Free Particle Lateral Constrained: O Global Position (Analytical)");
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,0*simulation.t);
plot(simulation.t,g*simulation.t);
title("Free Particle Lateral Constrained: O Global Velocity (Analytical)");
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,0*simulation.t);
plot(simulation.t,g*ones(size(simulation.t)));
title("Free Particle Lateral Constrained: O Global Acceleration (Analytical)");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'freeParticleLateralConstrained_O_Plot_Analytical.png');
end

% Quaternian Plot
figure;
hold on;
plot(simulation.t,simulation.t*0+1);
plot(simulation.t,simulation.t*0);
plot(simulation.t,simulation.t*0,'-.');
plot(simulation.t,simulation.t*0,'--');
title("Free Particle Lateral Constrained: Quaternian Components (Analytical)");
xlabel("t (s)");
ylabel("Quaternian (-)");
legend('1','2','3','4');
ylim([-1,2]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'freeParticleLateralConstrained_Quaternian_Plot_Analytical.png');
end

toc;
%profile viewer