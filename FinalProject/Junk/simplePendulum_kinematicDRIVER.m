%% Simple Pendulum Driver Function
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


L = 2;
rho = 7800;
a = 2*L;
b = 0.05;
c = b;
m = rho*a*b*c
J_xx_bar = 1/12*m*(b^2 + c^2)
J_yy_bar = 1/12*m*(a^2 + c^2)
J_zz_bar = 1/12*m*(a^2 + b^2)




%% Add library of functions to path
addpath('./Functions');
addpath('./Functions/Constraints');


%% Solver
% Initialize class for the simEngine3D
simulation = simEngine3D;

% Read the input deck
simulation.ReadInputDeck("simplePendulum_kinematicINPUT.mdl");

% Run the solver: (t_initial, dt, t_final, tolerance)
simulation.KinematicSolver(0, 0.01, 10, 1e-6);
%simulation.InverseDynamicsSolver(0, 0.01, 10, 1e-6);
%simulation.KinematicsSolver(0, 0.01, 10, 1e-4);


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
title("Simple Pendulum: Q Global Position");
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
title("Simple Pendulum: Q Global Velocity");
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
title("Simple Pendulum: Q Global Acceleration");
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
title("Simple Pendulum: O Global Position");
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
title("Simple Pendulum: O Global Velocity");
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
title("Simple Pendulum: O Global Acceleration");
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
	
	omega_ana_store(:,tt) = omega_ana;
end

% O Analytical Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t, r_ana(1,:));
plot(simulation.t, r_ana(2,:));
plot(simulation.t, r_ana(3,:));
title("Simple Pendulum: O Global Position (Analytical)");
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
title("Simple Pendulum: O Global Velocity (Analytical)");
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
title("Simple Pendulum: O Global Acceleration (Analytical)");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'simplePendulum_Q_Analytical_Plot.png');
end


% Analytical Omega plot
figure;
hold on;
plot(simulation.t, omega_ana_store(1,:));
plot(simulation.t, omega_ana_store(2,:));
plot(simulation.t, omega_ana_store(3,:));
title("Simple Pendulum: \omega (Analytical)");
xlabel("t (s)");
ylabel("\omega (rad/s)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'simplePendulum_omega_Analytical_Plot.png');
end


%% Calculate and Plot omega of each body for all times
for i = 1:simulation.N_Bodies
	% Calculate omega per body
	for tt = 1:simulation.N_t
		omega{i}(:,tt) = 2*E(simulation.p(4*(i-1)+1:4*(i-1)+4,tt))*simulation.p_dot(4*(i-1)+1:4*(i-1)+4,tt);
	end
	
	% Omega Plot per body
	figure;
	hold on;
	plot(simulation.t,omega{i}(1,:));
	plot(simulation.t,omega{i}(2,:));
	plot(simulation.t,omega{i}(3,:));
	title(['Simple Pendulum: \omega Body ' num2str(i)]);
	xlabel("t (s)");
	ylabel("\omega (rad/s)");
	legend('x','y','z');
	hold off;
end


toc;
%profile viewer