%% simEngine3D_A6P3 Driver Function
clear; close all; clc;
% Profiler shows the timings, calls, etc.
%profile on

% Timer
tic;

% Booleans
SAVE_PLOTS = 0;

%% Initial state
theta_ini = pi/4*cos(2*0);
A_ini = [0,				0,					1;...
		sin(theta_ini), cos(theta_ini),		0;...
		-cos(theta_ini), sin(theta_ini),	0;];
e0_ini = ((trace(A_ini)+1)/4)^(1/2)
e1_ini = ((2*A_ini(1,1)-trace(A_ini)+1)/4)^(1/2)
e2_ini = ((2*A_ini(2,2)-trace(A_ini)+1)/4)^(1/2)
e3_ini = ((2*A_ini(3,3)-trace(A_ini)+1)/4)^(1/2)
p_ini = [e0_ini, e1_ini, e2_ini, e3_ini]
p_ini_norm = norm(p_ini)
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
simulation.KinematicSolver(0, 0.1, 10, 1e-6);


%% Output final timestep information
disp("Phi");
simulation.Phi_G
disp("nu");
simulation.nu_G
disp("gamma");
simulation.gamma_G
disp("Jacobian");
simulation.Jacobian_G


%% Calculate rho (vector to Q), rho_dot, rho_ddot for all times
for tt = 1:simulation.N_t
	% Local vector to Q
	a_i_bar = [-2, 0, 0]';
	% From the simulation
	r_i = simulation.q(1:3,tt);
	r_i_dot = simulation.q_dot(1:3,tt);
	r_i_ddot = simulation.q_ddot(1:3,tt);
	p_i = simulation.q(4:7,tt);
	p_i_dot = simulation.q_dot(4:7,tt);
	p_i_ddot = simulation.q_ddot(4:7,tt);
	
	% Calculate rho, rho_dot, rho_ddot
	rho(1:3,tt) = r_i + A(p_i)*a_i_bar;
	rho_dot(1:3,tt) = r_i_dot + B(p_i, a_i_bar)*p_i_dot;
	rho_ddot(1:3,tt) = r_i_ddot + B(p_i_dot, a_i_bar)*p_i_dot + B(p_i, a_i_bar)*p_i_ddot;
end


%% Q Plots
% Q Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,rho(1,:));
plot(simulation.t,rho(2,:));
plot(simulation.t,rho(3,:));
title("Q Global Position");
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t,rho_dot(1,:));
plot(simulation.t,rho_dot(2,:));
plot(simulation.t,rho_dot(3,:));
title("Q Global Velocity");
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t,rho_ddot(1,:));
plot(simulation.t,rho_ddot(2,:));
plot(simulation.t,rho_ddot(3,:));
title("Q Global Acceleration");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'Q_Plot.png');
end


%% Origin Plots
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,simulation.q(1,:));
plot(simulation.t,simulation.q(2,:));
plot(simulation.t,simulation.q(3,:));
title("O Global Position");
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
title("O Global Velocity");
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
title("O Global Acceleration");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'O_Plot.png');
end


%% Analytical Solution
L = 2;

theta = @(t)(pi/4*cos(2*t));
theta_dot = @(t)(-1/2*pi*sin(2*t));
theta_ddot = @(t)(-pi*cos(2*t));

%ff = @(t)cos(pi/4*cos(2*t));
%figure;plot(0:.01:10,ff(0:.01:10));

for tt = 1:simulation.N_t
	ttt = simulation.t(tt);
	r_ana_i = [0,L*sin(theta(ttt)),-L*cos(theta(ttt))]';
	omega_ana = [theta_dot(ttt), 0, 0]';
	omega_ana_dot = [theta_ddot(ttt), 0, 0]';
	r_ana(:,tt) = r_ana_i;
	r_ana_dot(:,tt) = ToTilde(omega_ana)*r_ana_i;
	r_ana_ddot(:,tt) = ToTilde(omega_ana)*ToTilde(omega_ana)*r_ana_i + ToTilde(omega_ana_dot)*r_ana_i;
end

% O Analytical Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t, r_ana(1,:));
plot(simulation.t, r_ana(2,:));
plot(simulation.t, r_ana(3,:));
title("O Global Position: Analytical");
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Analytical Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t, r_ana_dot(1,:));
plot(simulation.t, r_ana_dot(2,:));
plot(simulation.t, r_ana_dot(3,:));
title("O Global Velocity: Analytical");
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Analytical Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t, r_ana_ddot(1,:));
plot(simulation.t, r_ana_ddot(2,:));
plot(simulation.t, r_ana_ddot(3,:));
title("O Global Acceleration: Analytical");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'Q_Analytical_Plot.png');
end

%{
% Deviations from the analytical
dev_r_x = norm(simulation.q(1,:)-r_ana(1,:))
dev_r_y = norm(simulation.q(2,:)-r_ana(2,:))
dev_r_z = norm(simulation.q(3,:)-r_ana(3,:))

dev_r_dot_x = norm(simulation.q_dot(1,:)-r_ana_dot(1,:))
dev_r_dot_y = norm(simulation.q_dot(2,:)-r_ana_dot(2,:))
dev_r_dot_z = norm(simulation.q_dot(3,:)-r_ana_dot(3,:))

dev_r_ddot_x = norm(simulation.q_ddot(1,:)-r_ana_ddot(1,:))
dev_r_ddot_y = norm(simulation.q_ddot(2,:)-r_ana_ddot(2,:))
dev_r_ddot_z = norm(simulation.q_ddot(3,:)-r_ana_ddot(3,:))
%}

toc;
%profile viewer