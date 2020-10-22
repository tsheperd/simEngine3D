clear;close all; clc;
% Booleans
SAVE_PLOTS = 0;


%% Inputs
% Parameters
alpha = 1;
beta = 1;

% Initial conditions
x0 = 0;
y0 = 2;

% Time range
t0 = 0;
tf = 20;

% Stepsize
h = 0.001;

% Solve over specified times
t = t0:h:tf;

% Initialize solution space
x = t*0;
y = x;

% Initial conditions
x(1) = x0;
y(1) = y0;

%% Backwards Euler for each time step
for n = 2:length(t)
	% Newton-Raphson iterative approach
	delta = [1;1];
	tol = 1e-8;
	k = 0;
	while abs(norm(delta)) > tol
		% g, as given
		g = [x(n)*(1+h)+(4*h*x(n)*y(n))/(1 + x(n)^2) - x(n-1) - h*alpha;...
			-h*beta*x(n) + y(n) + (h*beta*x(n)*y(n))/(1 + x(n)^2) - y(n-1);];
		
		% J, as given
		J = [1 + h + 4*h*y(n)*(1 - x(n)^2)/(1 + x(n)^2)^2,	(4*h*x(n))/(1 + x(n)^2);...
			-h*beta + h*beta*(1 - x(n)^2)/(1 + x(n)^2)^2,	1 + beta*h*(x(n))/(1 + x(n)^2)];

		% Calculate the residual
		delta = J\g;

		% Update the guess
		x(n) = x(n) - delta(1);
		y(n) = y(n) - delta(2);

		% Breakout counter
		k = k+1;
		if (k > 10000)
			disp("BREAK");
			break;
		end
	end
end


%% Plots
figure;
hold on;
plot(t,x);
title("x");
xlabel("t");
ylabel("x");
hold off;
if SAVE_PLOTS
	saveas(gcf,'HW7_P2_x.png');
end


figure;
hold on;
plot(t,y);
title("y");
xlabel("t");
ylabel("y");
hold off;
if SAVE_PLOTS
	saveas(gcf,'HW7_P2_y.png');
end


figure;
hold on;
plot(t,x);
plot(t,y);
title("x and y");
xlabel("t");
ylabel("x, y");
legend('x','y');
hold off;
if SAVE_PLOTS
	saveas(gcf,'HW7_P2_xy.png');
end