%% Main class for simEngine3D
% Contains all simulation requirements
classdef simEngine3D < handle
	properties
		Args
		inputDeckFileName
		inputDeckFile
		input
		t_i
		t_f
		dt
		t
		N_t
		
		N_Bodies
		N_GCons
		N_EPCons
		N_Cons
		
		N_params
		q
		q_dot
		q_ddot
		
		f
		f_dot
		f_ddot
		
		Phi_G
		nu_G
		gamma_G
		Jacobian_G
		
		tol
	end
	methods
		
		
		%% Constructor function
		function obj = simEngine3D(val)
			if nargin > 0
				obj.Args = val;
			end
		end
		
		
		%% Input Deck Parsing Function
		% Function that takes a file name "inputDeck.mdl", opens the file,
		% and parses it. Input deck must be in JSON format
		function obj = ReadInputDeck(obj,val)
			% Read the input deck file
 			obj.inputDeckFileName = val;
			obj.inputDeckFile = fileread(val);
			
			% Remove all comments of the form /*...*/
			inputDeckFile_Mod = obj.inputDeckFile;
			while ~isempty(strfind(inputDeckFile_Mod,'/*'))
				i1 = strfind(inputDeckFile_Mod,"/*");
				i2 = strfind(inputDeckFile_Mod,"*/");
				inputDeckFile_Mod = strcat(inputDeckFile_Mod(1:i1(1)-1),inputDeckFile_Mod(i2(1)+2:end));
			end
			obj.inputDeckFile = inputDeckFile_Mod;
			
			% Parse the input as JSON
			obj.input = jsondecode(obj.inputDeckFile);
			
			% Need to convert the contraints to a cell because matlab
			% sometimes makes them a structure and sometimes a cell but
			% these are accessed differently
			if isstruct(obj.input.constraints)
				struct2cell(obj.input.constraints);
			end
			
			% Parse input times
			%obj.t_i = obj.input.time(1);
			%obj.dt = obj.input.time(2);
			%obj.t_f = obj.input.time(3);
			
			
		end
		
		
		%% Function to act as a kinematic solver
		function obj = KinematicSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp)
			% Parse input times
			obj.t_i = t_i_temp;
			obj.dt = dt_temp;
			obj.t_f = t_f_temp;
			
			% Create a time vector
			obj.t = (obj.t_i:obj.dt:obj.t_f);
			obj.N_t = length(obj.t);
			
			obj.tol = tol_temp;
			
			
			% Create bodies' generalized parameters, q vector, for all time
			% Find the number of bodies
			obj.N_Bodies = size(obj.input.bodies,1);
			% Initialize the q, q_dot, q_ddot vectors as a matrix for all time
			obj.N_params = 7*obj.N_Bodies;
			obj.q = zeros(obj.N_params,obj.N_t);
			obj.q_dot = zeros(obj.N_params,obj.N_t);
			obj.q_ddot = zeros(obj.N_params,obj.N_t);
			% Set the initial conditons of q from the input deck
			for i = 1:obj.N_Bodies
				% q for each body at time t_i
				q_i = [obj.input.bodies(i).r; obj.input.bodies(i).p];
				obj.q(7*(i-1)+1:7*(i-1)+7,1) = q_i;
				
				% q_dot for each body at time t_i
				q_i_dot = [obj.input.bodies(i).r_dot; obj.input.bodies(i).p_dot];
				obj.q_dot(7*(i-1)+1:7*(i-1)+7,1) = q_i_dot;
			end
			
			
			% Create forcing function (and derivatives) vectors
			% Find the number of GCons
			obj.N_GCons = size(obj.input.constraints,1);
			% Initialize the f, f_dot, f_ddot vectors as a matrix for all time
			obj.f = zeros(obj.N_GCons,obj.N_t);
			obj.f_dot = zeros(obj.N_GCons,obj.N_t);
			obj.f_ddot = zeros(obj.N_GCons,obj.N_t);
			for i = 1:obj.N_GCons
				% Calculate f for each time
				f_t = str2func(['@(t)', obj.input.constraints{i}.f]);
				for j = 1:obj.N_t
					obj.f(i,j) = f_t(obj.t(j));
				end
				% Calculate f_dot for each time
				f_t_dot = str2func(['@(t)', obj.input.constraints{i}.f_dot]);
				for j = 1:obj.N_t
					obj.f_dot(i,j) = f_t_dot(obj.t(j));
				end
				% Calculate f_ddot for each time
				f_t_ddot = str2func(['@(t)', obj.input.constraints{i}.f_ddot]);
				for j = 1:obj.N_t
					obj.f_ddot(i,j) = f_t_ddot(obj.t(j));
				end
			end
			
			
			% Number of euler parameter constraints
			obj.N_EPCons = obj.N_Bodies;
			% Number of constraints in total (Note: N_GCons contains both
			% kinematic and driving)
			obj.N_Cons = obj.N_GCons + obj.N_EPCons;

			% Initialize global Phi, nu, gamma, and Jacobian
			obj.Phi_G = zeros(obj.N_Cons, 1);
			obj.nu_G = zeros(obj.N_Cons, 1);
			obj.gamma_G = zeros(obj.N_Cons, 1);
			obj.Jacobian_G = zeros(obj.N_Cons, obj.N_params);
			
			
			% Compute Phi_G, nu_G, gamma_G, Jacobian_G
			Global_Phi_nu_gamma_Jacobian(obj, 1);
			
			%%{
			% Get Phi_G, nu_G, gamma_G, Jacobian_G
			for tt = 1:obj.N_t
				if tt > 1
					% For all but the first timestep the initial guess for
					% q and q_dot are obtained from the previous time step
					obj.q(:,tt) = obj.q(:,tt-1);
					obj.q_dot(:,tt) = obj.q_dot(:,tt-1);
				end
				
				% Newton-Raphson iterative approach to solving the
				% current q
				dq = 1;
				%tol = 1e-3;
				k = 0;
				while abs(norm(dq,2)) > obj.tol
					% Compute Phi_G, nu_G, gamma_G, Jacobian_G
					Global_Phi_nu_gamma_Jacobian(obj, tt);
					
					% Calculate the residual
					dq = obj.Jacobian_G\obj.Phi_G;
					
					% Update the guess for q
					obj.q(:,tt) = obj.q(:,tt) - dq;

					% Breakout counter
					k = k+1;
					if (k > 10000)
						break;
						disp("BREAK");
					end
				end
				%disp(k);
				%disp(abs(norm(dq)));
				%disp(obj.tol);
				% Compute Phi_G, nu_G, gamma_G, Jacobian_G
				Global_Phi_nu_gamma_Jacobian(obj, tt);
				
				% Solve for q_dot and q_ddot at this timestep
				obj.q_dot(:,tt) = obj.Jacobian_G\obj.nu_G;
				obj.q_ddot(:,tt) = obj.Jacobian_G\obj.gamma_G;
			end
			%%}
		end
		
		
		%% Compute Phi_G, nu_G, gamma_G, Jacobian_G
		function obj = Global_Phi_nu_gamma_Jacobian(obj, tt)
			% Calculate from each GCon
			for k = 1:obj.N_GCons
				i = obj.input.constraints{k}.i;
				j = obj.input.constraints{k}.j;

				% First body
				% If the first body is ground
				if i == 0
					q_i = [0, 0, 0, 1, 0, 0, 0]';
					q_i_dot = [0, 0, 0, 0, 0, 0, 0]';
				% Otherwise coordinates are from q
				else
					q_i = obj.q(7*(i-1)+1:7*(i-1)+7,tt);
					q_i_dot = obj.q_dot(7*(i-1)+1:7*(i-1)+7,tt);
				end

				% Second body
				% If the second body is ground
				if j == 0
					q_j = [0, 0, 0, 1, 0, 0, 0]';
					q_j_dot = [0, 0, 0, 0, 0, 0, 0]';
				% Otherwise coordinates are from q
				else
					q_j = obj.q(7*(j-1)+1:7*(j-1)+7,tt);
					q_j_dot = obj.q_dot(7*(j-1)+1:7*(j-1)+7,tt);
				end

				% Get a handle for the GCon type call
				Phi_handle = str2func(obj.input.constraints{k}.type);

				% Populate global Phi, nu, gamma, Jacobian
				obj.Phi_G(k,1) = Phi_handle(obj.input.constraints{k}, q_i, q_j,  q_i_dot, q_j_dot, obj.f(k,tt), obj.f_dot(k,tt), obj.f_ddot(k,tt), 'Phi');
				obj.nu_G(k,1) = Phi_handle(obj.input.constraints{k}, q_i, q_j,  q_i_dot, q_j_dot, obj.f(k,tt), obj.f_dot(k,tt), obj.f_ddot(k,tt), 'nu');
				obj.gamma_G(k,1) = Phi_handle(obj.input.constraints{k}, q_i, q_j,  q_i_dot, q_j_dot, obj.f(k,tt), obj.f_dot(k,tt), obj.f_ddot(k,tt), 'gamma');
				Jacobian_temp = Phi_handle(obj.input.constraints{k}, q_i, q_j,  q_i_dot, q_j_dot, obj.f(k,tt), obj.f_dot(k,tt), obj.f_ddot(k,tt), 'Jacobian');

				% Jacobian expanded to global size
				% If the first body is ground
				if i == 0
					obj.Jacobian_G(k,7*(j-1)+1:7*(j-1)+3) = Jacobian_temp(1:3);
					obj.Jacobian_G(k,7*(j-1)+4:7*(j-1)+7) = Jacobian_temp(4:7);
				% If the second body is ground
				elseif j == 0
					obj.Jacobian_G(k,7*(i-1)+1:7*(i-1)+3) = Jacobian_temp(1:3);
					obj.Jacobian_G(k,7*(i-1)+4:7*(i-1)+7) = Jacobian_temp(4:7);
				else
					obj.Jacobian_G(k,7*(i-1)+1:7*(i-1)+3) = Jacobian_temp(1:3);
					obj.Jacobian_G(k,7*(i-1)+4:7*(i-1)+7) = Jacobian_temp(4:7);

					obj.Jacobian_G(k,7*(j-1)+1:7*(j-1)+3) = Jacobian_temp(8:10);
					obj.Jacobian_G(k,7*(j-1)+4:7*(j-1)+7) = Jacobian_temp(11:14);
				end
			end

			% Calculate from each Euler parameter constraint
			for i = 1:obj.N_EPCons
				p_i = obj.q(7*(i-1)+1+3:7*(i-1)+7,tt);
				p_i_dot = obj.q_dot(7*(i-1)+1+3:7*(i-1)+7,tt);
				
				k = obj.N_GCons + i;
				obj.Phi_G(k, 1) = p_i'*p_i - 1;
				obj.nu_G(k, 1) = 0;
				obj.gamma_G(k, 1) = -2*p_i_dot'*p_i_dot;
				obj.Jacobian_G(k,7*(i-1)+1+3:7*(i-1)+7) = 2*p_i';
			end
		end
		
		
		
	end
	%methods(Static)
	%end
end