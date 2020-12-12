%% Simple Pendulum Driver Function
clear; close all; clc;
% Profiler shows the timings, calls, etc.
%profile on

% Timer
tic;

% Booleans
SAVE_PLOTS = 0;


%% Solution Type Selection
% Set this
case_id = 1;

if case_id == 1
	SolveType = "Kinematics";
	name = "kinematic_GCon";
	name_nice = " (Kinematic, GCon)";
elseif case_id == 2
	SolveType = "InverseDynamics";
	name = "inverseDynamics_GCon";
	name_nice = " (Inverse Dynamics, GCon)";
elseif case_id == 3
	SolveType = "Dynamics";
	name = "dynamics_GCon";
	name_nice = " (Dynamics, GCon)";
elseif case_id == 4
	SolveType = "Dynamics";
	name = "dynamicsFree_GCon";
	name_nice = " (Dynamics Free, GCon)";
elseif case_id == 5
	SolveType = "Kinematics";
	name = "kinematic_CCon";
	name_nice = " (Kinematic, CCon)";
elseif case_id == 6
	SolveType = "InverseDynamics";
	name = "inverseDynamics_CCon";
	name_nice = " (Inverse Dynamics, CCon)";
elseif case_id == 7
	SolveType = "Dynamics";
	name = "dynamics_CCon";
	name_nice = " (Dynamics, CCon)";
elseif case_id == 8
	SolveType = "Dynamics";
	name = "dynamicsFree_CCon";
	name_nice = " (Dynamics Free, CCon)";
elseif case_id == 9
	SolveType = "Dynamics";
	name = "dynamicsTorqued_CCon";
	name_nice = " (Dynamics Torqued, CCon)";
elseif case_id == 10
	SolveType = "Dynamics";
	name = "dynamicsTorqued2_CCon";
	name_nice = " (Dynamics Torqued 2, CCon)";
elseif case_id == 11
	SolveType = "Dynamics";
	name = "dynamicsDistance_GCon";
	name_nice = " (Dynamics GCon Distance, GCon)";
end


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
simulation.ReadInputDeck("simplePendulum_"+name+"_INPUT.mdl");

% Run the solver: (t_initial, dt, t_final, tolerance)
if SolveType == "Kinematics"
	simulation.KinematicSolver(0, 0.01, 10, 1e-6);
elseif SolveType == "InverseDynamics"
	simulation.InverseDynamicsSolver(0, 0.01, 10, 1e-6);
elseif SolveType == "Dynamics"
	simulation.DynamicsSolver(0, 0.005, 10, 1e-4);
end


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
title("Simple Pendulum: Q Global Position"+name_nice);
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% Q Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t,rho_dot(1,:));
plot(simulation.t,rho_dot(2,:));
plot(simulation.t,rho_dot(3,:));
title("Simple Pendulum: Q Global Velocity"+name_nice);
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% Q Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t,rho_ddot(1,:));
plot(simulation.t,rho_ddot(2,:));
plot(simulation.t,rho_ddot(3,:));
title("Simple Pendulum: Q Global Acceleration"+name_nice);
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'simplePendulum_Q_Plot_'+name+'.png');
end


%% Origin Plots
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,simulation.q(1,:));
plot(simulation.t,simulation.q(2,:));
plot(simulation.t,simulation.q(3,:));
title("Simple Pendulum: O Global Position"+name_nice);
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
ylim([-2.125, 2.125]);
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(simulation.t,simulation.q_dot(1,:));
plot(simulation.t,simulation.q_dot(2,:));
plot(simulation.t,simulation.q_dot(3,:));
title("Simple Pendulum: O Global Velocity"+name_nice);
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
ylim([-3.25, 3.25]);
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(simulation.t,simulation.q_ddot(1,:));
plot(simulation.t,simulation.q_ddot(2,:));
plot(simulation.t,simulation.q_ddot(3,:));
title("Simple Pendulum: O Global Acceleration"+name_nice);
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
ylim([-5.25, 5.25]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'simplePendulum_O_Plot_'+name+'.png');
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
	title(['Simple Pendulum: \omega Body ' num2str(i)]+name_nice);
	xlabel("t (s)");
	ylabel("\omega (rad/s)");
	legend('x','y','z');
	ylim([-1.6, 1.6]);
	hold off;
	if SAVE_PLOTS
		saveas(gcf,'simplePendulum_omega_Plot_'+name+'.png');
	end
end

if SolveType == "InverseDynamics"
	%% Torque Plot
	for i = 1:simulation.N_Bodies
		for tt = 1:simulation.N_t
			% Driving GCon
			GCon_idx = 6;
			tau(:,tt) = simulation.tau_rxn{i}{GCon_idx,tt};
		end
		figure;
		hold on;
		plot(simulation.t, tau(1,:));
		plot(simulation.t, tau(2,:));
		plot(simulation.t, tau(3,:));
		title("Driving Torque"+name_nice);
		xlabel("t (s)");
		ylabel("torque (N-m)");
		legend('x','y','z');
		ylim([-250, 250]);
		hold off;
		if SAVE_PLOTS
			saveas(gcf,"Torque_Plot_"+num2str(i)+"_"+name+".png");
		end
	end
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
ylim([-2.125, 2.125]);
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
ylim([-3.25, 3.25]);
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
ylim([-5.25, 5.25]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'simplePendulum_O_Analytical_Plot.png');
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
ylim([-1.6, 1.6]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'simplePendulum_omega_Analytical_Plot.png');
end

%{
if case_id == 10
	figure;
	theta = atan2(simulation.q(3,:),simulation.q(2,:));
	for ii=1:length(theta)
		if theta <0
			theta=theta+2*pi;
		end
	end
	theta = theta+pi/2;
	plot(simulation.t,theta)

	figure;
	theta_anal = 1/2*10/104*simulation.t.^2+pi/4;
	plot(simulation.t,theta_anal)
end
%}

toc;
%profile viewer