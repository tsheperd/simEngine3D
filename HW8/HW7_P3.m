clear;close all; clc;
% Booleans
SAVE_PLOTS = 0;


%% Inputs
% Initial conditions
y0 = 1;

% Time range
t0 = 1;
tf = 10;


%% Scaling analysis by looking at different h sizes
h_v = logspace(-1,-3,25);%linspace(0.01,0.00001,1000);
for h_i = 1:length(h_v)

	% Stepsize
	h = h_v(h_i);
	
	% Solve over specified times
	t = t0:h:tf;

	% Initialize solution space
	y = t*0;

	% Initial conditions
	y(1) = y0;

	% Analytical solution for comparison
	y_analytical = 1./t + 1./t.^2.*tan(1./t+pi-1);
	
	
	%% BE
	for n = 2:length(t)
		% Newton-Raphson iterative approach
		delta = 1;
		tol = 1e-5;
		k = 0;
		while abs(norm(delta)) > tol
			% g
			g = y(n) - y(n-1) + h*y(n)^2 + h/t(n)^4;

			% J
			J = 1 + 2*h*y(n);

			% Calculate the residual
			delta = g/J;

			% Update the guess
			y(n) = y(n) - delta;

			% Breakout counter
			k = k+1;
			if (k > 10000)
				disp("BREAK");
				break;
			end
		end
	end
	y_BE = y;
	y = y*0; 
	y(1) = y0;


	%% BDF
	% Seed the first 4 y values with the analytical solution
	y(1:4) = y_analytical(1:4);

	for n = 5:length(t)
		% Newton-Raphson iterative approach
		delta = 1;
		tol = 1e-12;
		k = 0;
		while abs(norm(delta)) > tol
			% g
			g = y(n) - 48/25*y(n-1) + 36/25*y(n-2) - 16/25*y(n-3) + 3/25*y(n-4) - h*12/25*(-y(n)^2 - 1/t(n)^4);

			% J
			J = 1 + 24/25*h*y(n);

			% Calculate the residual
			delta = g/J;

			% Update the guess
			y(n) = y(n) - delta;

			% Breakout counter
			k = k+1;
			if (k > 10000)
				disp("BREAK");
				break;
			end
		end
	end
	
	y_BDF = y;
	y = y*0; 
	y(1) = y0;

	error_BE(h_i) = abs(y_analytical(end) - y_BE(end));
	error_BDF(h_i) = abs(y_analytical(end) - y_BDF(end));
end


%% Plot the finest h answer
figure;
hold on;
plot(t,y_BE);
plot(t,y_BDF);
plot(t,y_analytical);
title("y");
xlabel("t");
ylabel("y");
legend('y_{BE}','y_{BDF}','y_{analytical}');
hold off;
if SAVE_PLOTS
	saveas(gcf,'HW7_P3_y.png');
end


%% Plot the errors
figure;
hold on;
plot(h_v,error_BE);
plot(h_v,error_BDF);
title("Error");
xlabel("h");
ylabel("Error");
legend('BE','BDF');
hold off;
if SAVE_PLOTS
	saveas(gcf,'HW7_P3_Error.png');
end


%% LogLog plot the errors with line of best fit
C_BE = polyfit(log(h_v), log(error_BE), 1);
y_fit_BE = exp(C_BE(1)*log(h_v) + C_BE(2));

C_BDF = polyfit(log(h_v), log(error_BDF), 1);
y_fit_BDF = exp(C_BDF(1)*log(h_v) + C_BDF(2));

figure;
loglog(h_v,error_BE, h_v,error_BDF, h_v,y_fit_BE, '-.', h_v,y_fit_BDF, '-.');
hold on;
title("Error (Log-Log Scale)");
xlabel("h");
ylabel("error");
le = legend({'BE','BDF',['exp(',num2str(C_BE(1)),' ln(h) + ',num2str(C_BE(2)),')'],...
			['exp(',num2str(C_BDF(1)),' ln(h) + ',num2str(C_BDF(2)),')']}, 'Location', 'southeast');
hold off;
if SAVE_PLOTS
	saveas(gcf,'HW7_P3_Error_LogLog.png');
end
