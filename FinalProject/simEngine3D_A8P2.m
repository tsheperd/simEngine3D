%% simEngine3D_A8P2 Driver Function
clear; close all; clc;
% Profiler shows the timings, calls, etc.
%profile on

% Timer
tic;

% Booleans
SAVE_PLOTS = 0;


%% Initial state
% Body 1
disp('Body 1');
theta_ini = pi/2;
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


% Body 2
disp('Body 2');
theta_ini = 0;
A_ini = [0,				0,					1;...
		sin(theta_ini), cos(theta_ini),		0;...
		-cos(theta_ini), sin(theta_ini),	0;];
e0_ini = ((trace(A_ini)+1)/4)^(1/2)
e1_ini = ((2*A_ini(1,1)-trace(A_ini)+1)/4)^(1/2)
e2_ini = ((2*A_ini(2,2)-trace(A_ini)+1)/4)^(1/2)
e3_ini = ((2*A_ini(3,3)-trace(A_ini)+1)/4)^(1/2)
p_ini = [e0_ini, e1_ini, e2_ini, e3_ini]
p_ini_norm = norm(p_ini)
r_ini = A_ini*[1,0,0]'+2*r_ini


L = 2;
rho = 7800;
a = 1*L;
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
%simulation.ReadInputDeck("revJoint_A8P2.mdl");
simulation.ReadInputDeck("revJoint_A8P2_2.mdl");

% Run the kinematic solver: (t_initial, dt, t_final, tolerance)
%simulation.InverseDynamicsSolver(0, 0.01, 10, 1e-6);
simulation.DynamicsSolver(0, 0.005, 10, 1e-3);




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
for i = 1:simulation.N_Bodies
	r_idx = 3*(i-1)+1+0:3*(i-1)+3;
	
	% O Position plot
	figure;
	subplot(3,1,1);
	hold on;
	plot(simulation.t,simulation.r(r_idx,:));
	title(['O Global Position for Body: ',num2str(i)']);
	xlabel("t (s)");
	ylabel("position (m)");
	legend('x','y','z');
	hold off;

	% O Velocity plot
	%figure;
	subplot(3,1,2);
	hold on;
	plot(simulation.t,simulation.r_dot(r_idx,:));
	title(['O Global Velocity for Body: ',num2str(i)']);
	xlabel("t (s)");
	ylabel("velocity (m/s)");
	legend('x','y','z');
	hold off;

	% O Acceleration plot
	%figure;
	subplot(3,1,3);
	hold on;
	plot(simulation.t,simulation.r_ddot(r_idx,:));
	title(['O Global Acceleration for Body: ',num2str(i)']);
	xlabel("t (s)");
	ylabel("acceleration (m/s^2)");
	legend('x','y','z');
	hold off;
	if SAVE_PLOTS
		saveas(gcf,['O_Body_',num2str(i),'_Plot.png']);
	end
end


%{
%% Velocity Violation Plot
i = 2;
figure;
hold on;
plot(simulation.t, simulation.v_violation{i});
title(['Velocity Violation of Joint ' num2str(i)]);
xlabel("t (s)");
ylabel("Velocity Violation (m/s)");
hold off;
if SAVE_PLOTS
	saveas(gcf,'Velocity_Violation_Plot.png');
end
%}


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
	title(['\omega Body ' num2str(i)]);
	xlabel("t (s)");
	ylabel("\omega (rad/s)");
	legend('x','y','z');
	hold off;
	if SAVE_PLOTS
		saveas(gcf,['omega_Body_',num2str(i),'_Plot.png']);
	end
end



toc;
%profile viewer