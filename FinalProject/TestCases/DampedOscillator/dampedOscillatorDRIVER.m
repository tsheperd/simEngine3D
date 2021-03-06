%% Damped Oscillator Driver Function
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
m = 4.1888;
k = 5;
x0 = z_0;
v0 = 0;
c = 1;


%% Add library of functions to path
addpath('./Functions');
addpath('./Functions/Constraints');


%% Solver
% Initialize class for the simEngine3D
simulation = simEngine3D;

% Read the input deck
simulation.ReadInputDeck("dampedOscillatorINPUT.mdl");


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
title("Damped Oscillator: O Global Position");
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
title("Damped Oscillator: O Global Velocity");
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
title("Damped Oscillator: O Global Acceleration");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'dampedOscillator_O_Plot.png');
end

% Quaternian Plot
figure;
hold on;
plot(simulation.t,simulation.p(1,:));
plot(simulation.t,simulation.p(2,:));
plot(simulation.t,simulation.p(3,:),'-.');
plot(simulation.t,simulation.p(4,:),'--');
title("Damped Oscillator: Quaternian Components");
xlabel("t (s)");
ylabel("Quaternian (-)");
legend('1','2','3','4');
ylim([-1,2]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'dampedOscillator_Quaternian_Plot.png');
end


%% Analytical Solution
ttt = simulation.t;
x_anal = g.*k.^(-1).*m+exp(1).^((-1/2).*c.*m.^(-1).*ttt).*(k.^(-1).*((-1).* ...
  g.*m+k.*x0).*cos((1/2).*m.^(-1).*((-1).*c.^2+4.*k.*m).^(1/2).*ttt)+( ...
  (-1).*c.^2+4.*k.*m).^(-1/2).*(2.*m.*v0+c.*k.^(-1).*((-1).*g.*m+k.* ...
  x0)).*sin((1/2).*m.^(-1).*((-1).*c.^2+4.*k.*m).^(1/2).*ttt));

x_dot_anal = exp(1).^((-1/2).*c.*m.^(-1).*ttt).*((1/2).*m.^(-1).*(2.*m.*v0+c.* ...
  k.^(-1).*((-1).*g.*m+k.*x0)).*cos((1/2).*m.^(-1).*((-1).*c.^2+4.* ...
  k.*m).^(1/2).*ttt)+(-1/2).*k.^(-1).*m.^(-1).*((-1).*c.^2+4.*k.*m).^( ...
  1/2).*((-1).*g.*m+k.*x0).*sin((1/2).*m.^(-1).*((-1).*c.^2+4.*k.*m) ...
  .^(1/2).*ttt))+(-1/2).*c.*exp(1).^((-1/2).*c.*m.^(-1).*ttt).*m.^(-1).* ...
  (k.^(-1).*((-1).*g.*m+k.*x0).*cos((1/2).*m.^(-1).*((-1).*c.^2+4.* ...
  k.*m).^(1/2).*ttt)+((-1).*c.^2+4.*k.*m).^(-1/2).*(2.*m.*v0+c.*k.^( ...
  -1).*((-1).*g.*m+k.*x0)).*sin((1/2).*m.^(-1).*((-1).*c.^2+4.*k.*m) ...
  .^(1/2).*ttt));

x_ddot_anal = (-1).*c.*exp(1).^((-1/2).*c.*m.^(-1).*ttt).*m.^(-1).*((1/2).*m.^(-1) ...
  .*(2.*m.*v0+c.*k.^(-1).*((-1).*g.*m+k.*x0)).*cos((1/2).*m.^(-1).*( ...
  (-1).*c.^2+4.*k.*m).^(1/2).*ttt)+(-1/2).*k.^(-1).*m.^(-1).*((-1).* ...
  c.^2+4.*k.*m).^(1/2).*((-1).*g.*m+k.*x0).*sin((1/2).*m.^(-1).*(( ...
  -1).*c.^2+4.*k.*m).^(1/2).*ttt))+(1/4).*c.^2.*exp(1).^((-1/2).*c.* ...
  m.^(-1).*ttt).*m.^(-2).*(k.^(-1).*((-1).*g.*m+k.*x0).*cos((1/2).* ...
  m.^(-1).*((-1).*c.^2+4.*k.*m).^(1/2).*ttt)+((-1).*c.^2+4.*k.*m).^( ...
  -1/2).*(2.*m.*v0+c.*k.^(-1).*((-1).*g.*m+k.*x0)).*sin((1/2).*m.^( ...
  -1).*((-1).*c.^2+4.*k.*m).^(1/2).*ttt))+exp(1).^((-1/2).*c.*m.^(-1) ...
  .*ttt).*((-1/4).*k.^(-1).*m.^(-2).*((-1).*c.^2+4.*k.*m).*((-1).*g.* ...
  m+k.*x0).*cos((1/2).*m.^(-1).*((-1).*c.^2+4.*k.*m).^(1/2).*ttt)+( ...
  -1/4).*m.^(-2).*((-1).*c.^2+4.*k.*m).^(1/2).*(2.*m.*v0+c.*k.^(-1) ...
  .*((-1).*g.*m+k.*x0)).*sin((1/2).*m.^(-1).*((-1).*c.^2+4.*k.*m).^( ...
  1/2).*ttt));

% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(simulation.t,0*simulation.t);
plot(simulation.t,0*simulation.t);
plot(simulation.t,x_anal);
title("Damped Oscillator: O Global Position (Analytical)");
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
plot(simulation.t,x_dot_anal);
title("Damped Oscillator: O Global Velocity (Analytical)");
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
plot(simulation.t,x_ddot_anal);
title("Damped Oscillator: O Global Acceleration (Analytical)");
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,'dampedOscillator_O_Plot_Analytical.png');
end

% Quaternian Plot
figure;
hold on;
plot(simulation.t,0*simulation.t+1);
plot(simulation.t,0*simulation.t);
plot(simulation.t,0*simulation.t,'-.');
plot(simulation.t,0*simulation.t,'--');
title("Damped Oscillator: Quaternian Components (Analytical)");
xlabel("t (s)");
ylabel("Quaternian (-)");
legend('1','2','3','4');
ylim([-1,2]);
hold off;
if SAVE_PLOTS
	saveas(gcf,'dampedOscillator_Quaternian_Plot_Analytical.png');
end

toc;
%profile viewer