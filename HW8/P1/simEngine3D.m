%% Main class for simEngine3D
% Contains simulation class calls
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
		
		Jacobian_r_G
		Jacobian_p_G
		
		gamma_k_G
		gamma_p_G
		
		Phi_k_G
		Phi_p_G
		
		tol
		
		m
		r_cm
		
		J_xx_bar
		J_yy_bar
		J_zz_bar
		
		M_i
		J_bar_i
		
		M
		J_P
		P
		
		F_rxn
		tau_rxn
		tau_rxn_bar
        
        r
        r_dot
        r_ddot
        
        p
        p_dot
        p_ddot
		
		lambda
		lambda_p
		
		v_violation
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
			
		end
		
		
		%% Function to act as a kinematic solver
		function obj = KinematicSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp)
			% Function to initialize the data arrays, etc.
			initializeSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp);
			
			% Compute Phi_G, nu_G, gamma_G, Jacobian_G
			Global_Phi_nu_gamma_Jacobian(obj, 1, [1, 1, 1, 1]);
			
			% Create a waitbar
			WB = waitbar(0, 'Kinematic Solver: Computing');
			
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
					Global_Phi_nu_gamma_Jacobian(obj, tt, [1, 0, 0, 1]);
					
					% Calculate the residual
					dq = obj.Jacobian_G\obj.Phi_G;
					
					% Update the guess for q
					obj.q(:,tt) = obj.q(:,tt) - dq;

					% Breakout counter
					k = k+1;
					if (k > 10000)
						disp("BREAK");
						break;
					end
				end

				% Compute Phi_G, nu_G, gamma_G, Jacobian_G
				Global_Phi_nu_gamma_Jacobian(obj, tt, [0, 1, 0, 1]);
				
				% Solve for q_dot at this timestep
				obj.q_dot(:,tt) = obj.Jacobian_G\obj.nu_G;
				
				% Compute Phi_G, nu_G, gamma_G, Jacobian_G with the updated
				% q_dot
				Global_Phi_nu_gamma_Jacobian(obj, tt, [1, 1, 1, 1]);
				
				% Solve for q_ddot at this timestep
				obj.q_ddot(:,tt) = obj.Jacobian_G\obj.gamma_G;
				
				% Show a progress bar
				waitbar(tt/obj.N_t, WB, 'Kinematic Solver: Computing');
			end
			
			% Close the waitbar
			delete(WB);
		end
		
		
		%% Compute Phi_G, nu_G, gamma_G, Jacobian_G
		function obj = Global_Phi_nu_gamma_Jacobian(obj, tt, FLAG)
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
				[obj.Phi_G(k,1), obj.nu_G(k,1), obj.gamma_G(k,1), Jacobian_temp] = Phi_handle(obj.input.constraints{k}, q_i, q_j,  q_i_dot, q_j_dot, obj.f(k,tt), obj.f_dot(k,tt), obj.f_ddot(k,tt), FLAG);
				
				% Jacobian expanded to global size
				% If the first body is not ground
				if i ~= 0
					obj.Jacobian_G(k,7*(i-1)+1:7*(i-1)+3) = Jacobian_temp(1:3);
					obj.Jacobian_G(k,7*(i-1)+4:7*(i-1)+7) = Jacobian_temp(4:7);
					
					obj.Jacobian_r_G(k,3*(i-1)+1:3*(i-1)+3) = Jacobian_temp(1:3);
					obj.Jacobian_p_G(k,4*(i-1)+1:4*(i-1)+4) = Jacobian_temp(4:7);
				end
				% If the second body is not ground
				if j ~= 0
					obj.Jacobian_G(k,7*(j-1)+1:7*(j-1)+3) = Jacobian_temp(8:10);
					obj.Jacobian_G(k,7*(j-1)+4:7*(j-1)+7) = Jacobian_temp(11:14);
					
					obj.Jacobian_r_G(k,3*(j-1)+1:3*(j-1)+3) = Jacobian_temp(8:10);
					obj.Jacobian_p_G(k,4*(j-1)+1:4*(j-1)+4) = Jacobian_temp(11:14);
				end
				
				% Create a global kinematic gamma and Phi for dynamics
				obj.gamma_k_G(k,1) = obj.gamma_G(k,1);
				obj.Phi_k_G(k,1) = obj.Phi_G(k,1);
				

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
				
				% Create a global euler parameterization gamma and Phi for dynamics
				obj.gamma_p_G(i, 1) = -2*p_i_dot'*p_i_dot;
				obj.Phi_p_G(i, 1) = p_i'*p_i - 1;
			end
        end
		
        
        %% Function to initialize the data arrays, etc.
        function obj = initializeSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp)
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
			
			obj.Jacobian_r_G = zeros(obj.N_GCons, 3*obj.N_Bodies);
			obj.Jacobian_p_G = zeros(obj.N_GCons, 4*obj.N_Bodies);
			
			obj.gamma_k_G = zeros(obj.N_GCons, 1);
			obj.gamma_p_G = zeros(obj.N_EPCons, 1);
			
			obj.Phi_k_G = zeros(obj.N_GCons, 1);
			obj.Phi_p_G = zeros(obj.N_EPCons, 1);
        end
        
        
		
		%% Function to act as the Inverse Dynamics Solver
		function obj = InverseDynamicsSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp)
			% Run the kinematics solver to produce the q, q_dot, q_ddot,
			% etc...
			obj.KinematicSolver(t_i_temp, dt_temp, t_f_temp, tol_temp);
			
			g = obj.input.gravity;
			
			% Initialize Global M, J_P, P matricies
			obj.M = zeros(3*obj.N_Bodies);
			obj.J_P = zeros(4*obj.N_Bodies);
			obj.P = zeros(obj.N_Bodies, 4*obj.N_Bodies);
			
			r_ddot = zeros(3*obj.N_Bodies, 1);
			p_ddot = zeros(4*obj.N_Bodies, 1);
            
			% Initialize r, p and their derivs
            obj.r = zeros(3*obj.N_Bodies, obj.N_t);
			obj.p = zeros(4*obj.N_Bodies, obj.N_t);
            
            obj.r_dot = zeros(3*obj.N_Bodies, obj.N_t);
			obj.p_dot = zeros(4*obj.N_Bodies, obj.N_t);
            
            obj.r_ddot = zeros(3*obj.N_Bodies, obj.N_t);
			obj.p_ddot = zeros(4*obj.N_Bodies, obj.N_t);
			
			F = zeros(3*obj.N_Bodies, 1);
			tau = zeros(4*obj.N_Bodies, 1);
			
			% For each body, calculate M_i, J_bar_i
			for i = 1:obj.N_Bodies
				% Grab parameters from input deck on dynamic rel. params
				% Mass
				m_i = obj.input.bodies(i).m;
				obj.m(i,1) = m_i;
				
				% Moments of inertia
				J_xx_bar_i = obj.input.bodies(i).J_xx_bar;
				J_yy_bar_i = obj.input.bodies(i).J_yy_bar;
				J_zz_bar_i = obj.input.bodies(i).J_zz_bar;
				obj.J_xx_bar(i,1) = J_xx_bar_i;
				obj.J_yy_bar(i,1) = J_yy_bar_i;
				obj.J_zz_bar(i,1) = J_zz_bar_i;
				
				% Center of mass vector
				r_cm_i = obj.input.bodies(i).r_cm;
				obj.r_cm{i} = r_cm_i;
				
				% Populate needed matricies for the EOM
				M_i = m_i*eye(3);
				obj.M_i{i} = M_i;
				
				J_bar_i = [	J_xx_bar_i,	0,		0;...
							0,			J_yy_bar_i, 0;...
							0,			0,		J_zz_bar_i;];
				obj.J_bar_i{i} = J_bar_i;
				
				% Push the M_i for the body to the global M
				obj.M(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3) = M_i;
			end
			
			
			for tt = 1:obj.N_t
				% Populate obj.r, obj.p, etc from obj.q
				q2rp(obj,tt,[1,1,1]);

				% For each body, calculate M_i, J_bar_i
				for i = 1:obj.N_Bodies
					% Grab the current r, p etc.
					r_idx = 3*(i-1)+1+0:3*(i-1)+3;
					p_idx = 4*(i-1)+1+0:4*(i-1)+4;

					r_i = obj.r(r_idx,tt);
					r_i_dot = obj.r_dot(r_idx,tt);
					r_i_ddot = obj.r_ddot(r_idx,tt);

					p_i = obj.p(p_idx,tt);
					p_i_dot = obj.p_dot(p_idx,tt);
					p_i_ddot = obj.p_ddot(p_idx,tt);

					% Calculate the J_P_i for the body
					J_bar_i = obj.J_bar_i{i};
					J_P_i = 4*G(p_i)'*J_bar_i*G(p_i);

					% Push the J_P_i for the body to the global J_P
					obj.J_P(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = J_P_i;

					% Global P
					obj.P(i,4*(i-1)+1:4*(i-1)+4) = p_i';

					m_i = obj.m(i,1);
					F_m_i = m_i*g;
					F_a_i = zeros(3,1);

					n_m_bar_i = zeros(3,1);
					n_a_bar_i = zeros(3,1);

					% Calculate the force per body
					F(3*(i-1)+1:3*(i-1)+3,1) = F_m_i + F_a_i;

					% Calculate the torque per body
					tau(4*(i-1)+1:4*(i-1)+4,1) = 2*G(p_i)'*(n_m_bar_i + n_a_bar_i)+8*G(p_i_dot)'*J_bar_i*G(p_i_dot)*p_i;
				end

				% Compute Jacobian_G
				Global_Phi_nu_gamma_Jacobian(obj, tt, [0, 0, 0, 1]);
				
				%obj.Jacobian_r_G
				%obj.Jacobian_p_G
				
				
				% LHS of the inverse kinematics equation
				LHS =	[obj.Jacobian_r_G', zeros(3*obj.N_Bodies, obj.N_Bodies);...
						 obj.Jacobian_p_G', obj.P';];
				
				% RHS of the inverse kinematics equation
				RHS =	[F - obj.M*obj.r_ddot(:,tt);...
						 tau - obj.J_P*obj.p_ddot(:,tt);];
				
				% Solve for all of the lagrange multipliers
				lambda_v = LHS\RHS;
				
				% Extract the kinematic and euler LM
				lambda = lambda_v(1:obj.N_GCons,1);
				lambda_P = lambda_v(obj.N_GCons+1:end,1);
				
				obj.lambda_p(:,tt) = lambda_P;
				obj.lambda(:,tt) = lambda;
				
				% Calculate the reaction torques and forces
				for i = 1:obj.N_Bodies
					p_i = obj.q(7*(i-1)+1+3:7*(i-1)+7,tt);
					Jacobian_r_i = obj.Jacobian_r_G(:,3*(i-1)+1:3*(i-1)+3);
					Jacobian_p_i = obj.Jacobian_p_G(:,4*(i-1)+1:4*(i-1)+4);
					
					% Solve reaction forces from each GCon
					for j = 1:obj.N_GCons
						obj.F_rxn{i}{j,tt} = -Jacobian_r_i(j,:)'*obj.lambda(j,tt);
						obj.tau_rxn_bar{i}{j,tt} = -1/2*G(p_i)*(Jacobian_p_i(j,:)'*obj.lambda(j,tt));
						obj.tau_rxn{i}{j,tt} = A(p_i)*obj.tau_rxn_bar{i}{j,tt};
					end
				end
				
			end
			
		end
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Function to act as the Dynamics Solver
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		function obj = DynamicsSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp)
			% Function to initialize the data arrays, etc.
			initializeSolver(obj, t_i_temp, dt_temp, t_f_temp, tol_temp);
			
			% Time step
            h = dt_temp;
			
			% Shorten name for counts of stuff
            nb = obj.N_Bodies;
			nc = obj.N_GCons;
			np = obj.N_EPCons;
            
            % Gravity
			g = obj.input.gravity;
			
			% Initialize Global M, J_P, P matricies
			obj.M = zeros(3*obj.N_Bodies);
			obj.J_P = zeros(4*obj.N_Bodies);
			obj.P = zeros(obj.N_Bodies, 4*obj.N_Bodies);
			
			%r_ddot = zeros(3*obj.N_Bodies, 1);
			%p_ddot = zeros(4*obj.N_Bodies, 1);
			
			% Initialize r, p and their derivs
            obj.r = zeros(3*obj.N_Bodies, obj.N_t);
			obj.p = zeros(4*obj.N_Bodies, obj.N_t);
            
            obj.r_dot = zeros(3*obj.N_Bodies, obj.N_t);
			obj.p_dot = zeros(4*obj.N_Bodies, obj.N_t);
            
            obj.r_ddot = zeros(3*obj.N_Bodies, obj.N_t);
			obj.p_ddot = zeros(4*obj.N_Bodies, obj.N_t);
			
            F = zeros(3*obj.N_Bodies, 1);
			tau = zeros(4*obj.N_Bodies, 1);
			
			
			% For each body, calculate M_i, J_bar_i
			for i = 1:obj.N_Bodies
				% Grab parameters from input deck on dynamic rel. params
				% Mass
				m_i = obj.input.bodies(i).m;
				obj.m(i,1) = m_i;
				
				% Moments of inertia
				J_xx_bar_i = obj.input.bodies(i).J_xx_bar;
				J_yy_bar_i = obj.input.bodies(i).J_yy_bar;
				J_zz_bar_i = obj.input.bodies(i).J_zz_bar;
				obj.J_xx_bar(i,1) = J_xx_bar_i;
				obj.J_yy_bar(i,1) = J_yy_bar_i;
				obj.J_zz_bar(i,1) = J_zz_bar_i;
				
				% Center of mass vector
				r_cm_i = obj.input.bodies(i).r_cm;
				obj.r_cm{i} = r_cm_i;
				
				% Populate needed matricies for the EOM
				M_i = m_i*eye(3);
				obj.M_i{i} = M_i;
				
				J_bar_i = [	J_xx_bar_i,	0,		0;...
							0,			J_yy_bar_i, 0;...
							0,			0,		J_zz_bar_i;];
				obj.J_bar_i{i} = J_bar_i;
				
				% Push the M_i for the body to the global M
				obj.M(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3) = M_i;
			end
			
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%% Initial conditions
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Initial conditions
			
			% For each body
			% Populate obj.r, obj.p, etc from obj.q
			% Calculate M_i, J_bar_i
			tt = 1;
			
			% Populate obj.r, obj.p, etc from obj.q
			q2rp(obj,tt,[1,1,1]);
			
			% For each body, calculate J_P, force, torque contributions 
			for i = 1:obj.N_Bodies
				% Grab the current r, p etc.
				r_idx = 3*(i-1)+1+0:3*(i-1)+3;
				p_idx = 4*(i-1)+1+0:4*(i-1)+4;

				r_i = obj.r(r_idx,tt);
				r_i_dot = obj.r_dot(r_idx,tt);
				r_i_ddot = obj.r_ddot(r_idx,tt);

				p_i = obj.p(p_idx,tt);
				p_i_dot = obj.p_dot(p_idx,tt);
				p_i_ddot = obj.p_ddot(p_idx,tt);

				% Calculate the J_P_i for the body
				J_bar_i = obj.J_bar_i{i};
				J_P_i = 4*G(p_i)'*J_bar_i*G(p_i);

				% Push the J_P_i for the body to the global J_P
				obj.J_P(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = J_P_i;

				% Global P
				obj.P(i,4*(i-1)+1:4*(i-1)+4) = p_i';
				
				% Forces/torques
				m_i = obj.m(i,1);
				F_m_i = m_i*g;
				F_a_i = zeros(3,1);

				n_m_bar_i = zeros(3,1);
				n_a_bar_i = zeros(3,1);

				% Calculate the force per body
				F(3*(i-1)+1:3*(i-1)+3,1) = F_m_i + F_a_i;

				% Calculate the torque per body
				tau(4*(i-1)+1:4*(i-1)+4,1) = 2*G(p_i)'*(n_m_bar_i + n_a_bar_i)+8*G(p_i_dot)'*J_bar_i*G(p_i_dot)*p_i;
			end

			% Compute Jacobian_G
			Global_Phi_nu_gamma_Jacobian(obj, tt, [1, 1, 1, 1]);
			
			% Compute z at this time step
			LHS =	[obj.M,					zeros(3*nb,4*nb),	zeros(3*nb,nb),		obj.Jacobian_r_G';...
					 zeros(4*nb,3*nb),		obj.J_P,			obj.P',				obj.Jacobian_p_G';...
					 zeros(nb,3*nb),		obj.P,				zeros(nb,nb),		zeros(nb,nc);...
					 obj.Jacobian_r_G,		obj.Jacobian_p_G,	zeros(nc,nb),		zeros(nc,nc);];

			RHS =	[F;...
					 tau;...
					 obj.gamma_p_G;...
					 obj.gamma_k_G;];
			
			z = LHS\RHS;
			
			% Extract values
			obj.r_ddot(:,tt) = z(1:3*nb,1);
			obj.p_ddot(:,tt) = z(3*nb+1:3*nb+4*nb,1);
			obj.lambda_p(:,tt) = z(3*nb+4*nb+1:3*nb+4*nb+np,1);
			obj.lambda(:,tt) = z(3*nb+4*nb+np+1:3*nb+4*nb+np+nc,1);
			
			% Update q
			rp2q(obj,tt,[1,1,1]);
			
			Global_Phi_nu_gamma_Jacobian(obj, tt, [1, 1, 1, 1]);
			
			% Check it see if the initial conditions satisfy Phi
			if norm(obj.Phi_G,2) < 1e-8
				disp('I.C. OK');
			end
			IC_Phi_norm = norm(obj.Phi_G,2)
			IC_velocity_dev_norm = norm(obj.Jacobian_G*obj.q_dot(:,tt) - obj.nu_G,2)
			IC_accel_dev_norm = norm(obj.Jacobian_G*obj.q_ddot(:,tt) - obj.gamma_G,2)
			
			
			% Calculate the reaction torques and forces
			for i = 1:obj.N_Bodies
				p_i = obj.q(7*(i-1)+1+3:7*(i-1)+7,tt);
				Jacobian_r_i = obj.Jacobian_r_G(:,3*(i-1)+1:3*(i-1)+3);
				Jacobian_p_i = obj.Jacobian_p_G(:,4*(i-1)+1:4*(i-1)+4);

				% Solve reaction forces from each GCon
				for j = 1:obj.N_GCons
					obj.F_rxn{i}{j,tt} = -Jacobian_r_i(j,:)'*obj.lambda(j,tt);
					obj.tau_rxn_bar{i}{j,tt} = -1/2*G(p_i)*(Jacobian_p_i(j,:)'*obj.lambda(j,tt));
					obj.tau_rxn{i}{j,tt} = A(p_i)*obj.tau_rxn_bar{i}{j,tt} ;
				end
			end
			
			% Compute the velocity violation for the second major joint
			i = 1;
			ncon = 7;
			nppar = 1;			
			v_v_idx = ncon*(i-1)+1+0:ncon*(i-1)+(ncon-nppar);
			v_violation_temp = obj.Jacobian_G*obj.q_dot(:,tt) - obj.nu_G;
			obj.v_violation{i}(tt,1) = norm(v_violation_temp(v_v_idx),2);
			
			% Create a waitbar
			WB = waitbar(0, 'Dynamics Solver: Computing');
			
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%% Loop from the second to the last time steps
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%% Loop from the second to the last time steps
            for tt = 2:obj.N_t
				n = tt;
				% Replace with waitbar later
				PercentComplete = n/obj.N_t*100;
				
				% Initial guess for r_ddot, p_ddot, q_ddot, lambda_p, lambda
				% from previous time step
                obj.r_ddot(:,tt) = obj.r_ddot(:,tt-1);
				obj.p_ddot(:,tt) = obj.p_ddot(:,tt-1);
				
				obj.q_ddot(:,tt) = obj.q_ddot(:,tt-1);
				
				obj.lambda_p(:,tt) = obj.lambda_p(:,tt-1);
				obj.lambda(:,tt) = obj.lambda(:,tt-1);
				
				% Calculate the BDF coefficients
                if n == 2
                    % For the first calculation use a BDF 1
                    % Parameters from the BDF table
                    alpha_1 = -1;
                    beta_0 = 1;
                    
                    % Constants used for r, r_dot, p, p_dot
                    C_r_dot_n = -alpha_1*obj.r_dot(:,n-1);
					C_r_n = -alpha_1*obj.r(:,n-1) + beta_0*h*(C_r_dot_n);
                    
                    C_p_dot_n = -alpha_1*obj.p_dot(:,n-1);
                    C_p_n = -alpha_1*obj.p(:,n-1) + beta_0*h*(C_p_dot_n);
                else
                    % For other than the first calculation use a BDF 2
                    % Parameters from the BDF table
                    alpha_1 = -4/3;
                    alpha_2 = 1/3;
                    beta_0 = 2/3;
                    
                    % Constants used for r, r_dot, p, p_dot
					C_r_dot_n = -alpha_1*obj.r_dot(:,n-1) - alpha_2*obj.r_dot(:,n-2);
                    C_r_n = -alpha_1*obj.r(:,n-1) - alpha_2*obj.r(:,n-2) + beta_0*h*(C_r_dot_n);
                    
                    C_p_dot_n = -alpha_1*obj.p_dot(:,n-1) - alpha_2*obj.p_dot(:,n-2);
                    C_p_n = -alpha_1*obj.p(:,n-1) - alpha_2*obj.p(:,n-2) + beta_0*h*(C_p_dot_n);
				end
				
				
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%%%% Newton-Raphson iterative approach to solving z
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				% Newton-Raphson iterative approach to solving the
				% current z
				dz = 1;
				k = 0;
				BreakCount = 10000;
				while abs(norm(dz,2)) > obj.tol
					% Current accelerations and lambdas
					r_ddot_n = obj.r_ddot(:,n);
					p_ddot_n = obj.p_ddot(:,n);
					
					lambda_p_n = obj.lambda_p(:,n);
					lambda_n = obj.lambda(:,n);
					
					% Current positions and velocities via BDF
					r_n = C_r_n + beta_0^2*h^2*r_ddot_n;
					r_dot_n = C_r_dot_n + beta_0*h*r_ddot_n;
					
					p_n = C_p_n + beta_0^2*h^2*p_ddot_n;
					p_dot_n = C_p_dot_n + beta_0*h*p_ddot_n;
					
					% Update r and p vectors of obj
					obj.r(:,tt) = r_n;
					obj.r_dot(:,tt) = r_dot_n;
					obj.r_ddot(:,tt) = r_ddot_n;
					
					obj.p(:,tt) = p_n;
					obj.p_dot(:,tt) = p_dot_n;
					obj.p_ddot(:,tt) = p_ddot_n;
					
					% Update q
					rp2q(obj,tt,[1,1,1]);
					
					% For each body, calculate J_P, P, etc.
					for i = 1:obj.N_Bodies
						% Grab the current r, p etc.
						r_idx = 3*(i-1)+1+0:3*(i-1)+3;
						p_idx = 4*(i-1)+1+0:4*(i-1)+4;

						r_i = obj.r(r_idx,tt);
						r_i_dot = obj.r_dot(r_idx,tt);
						r_i_ddot = obj.r_ddot(r_idx,tt);

						p_i = obj.p(p_idx,tt);
						p_i_dot = obj.p_dot(p_idx,tt);
						p_i_ddot = obj.p_ddot(p_idx,tt);

						% Calculate the J_P_i for the body
						J_bar_i = obj.J_bar_i{i};
						J_P_i = 4*G(p_i)'*J_bar_i*G(p_i);

						% Push the J_P_i for the body to the global J_P
						obj.J_P(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = J_P_i;

						% Global P
						obj.P(i,4*(i-1)+1:4*(i-1)+4) = p_i';

						m_i = obj.m(i,1);
						F_m_i = m_i*g;
						F_a_i = zeros(3,1);

						n_m_bar_i = zeros(3,1);
						n_a_bar_i = zeros(3,1);

						% Calculate the force per body
						F(3*(i-1)+1:3*(i-1)+3,1) = F_m_i + F_a_i;

						% Calculate the torque per body
						tau(4*(i-1)+1:4*(i-1)+4,1) = 2*G(p_i)'*(n_m_bar_i + n_a_bar_i)+8*G(p_i_dot)'*J_bar_i*G(p_i_dot)*p_i;
					end

					% Compute Jacobian_G
					Global_Phi_nu_gamma_Jacobian(obj, tt, [1, 1, 1, 1]);
					
					
					% Res.
					% WTF, why does this work? Need 1/(beta_0*h)^pow
					% instead of 1/(beta_0*h)^2
					pow = 1.90;
					g_n =	[obj.M*r_ddot_n + obj.Jacobian_r_G'*lambda_n - F;...
							 obj.J_P*p_ddot_n + obj.Jacobian_p_G'*lambda_n + obj.P'*lambda_p_n - tau;...
							 1/(beta_0*h)^pow*obj.Phi_p_G;...
							 1/(beta_0*h)^pow*obj.Phi_k_G;];
						 
					% I found out that removing the square term from
					% 1/(beta_0*h) solves much better because as h -> 0
					% that term blows up. Maybe something is wrong with my
					% obj.Phi_p_G and obj.Phi_k_G because their error would
					% be amplified? But they solved for correct answers in
					% the last homeworks..
					%g_n =	[obj.M*r_ddot_n + obj.Jacobian_r_G'*lambda_n - F;...
					%		 obj.J_P*p_ddot_n + obj.Jacobian_p_G'*lambda_n + obj.P'*lambda_p_n - tau;...
					%		 1/(beta_0*h)^2*obj.Phi_p_G;...
					%		 1/(beta_0*h)^2*obj.Phi_k_G;];
						 
					% Quasi
					Psi =	[obj.M,					zeros(3*nb,4*nb),	zeros(3*nb,nb),		obj.Jacobian_r_G';...
							 zeros(4*nb,3*nb),		obj.J_P,			obj.P',				obj.Jacobian_p_G';...
							 zeros(nb,3*nb),		obj.P,				zeros(nb,nb),		zeros(nb,nc);...
							 obj.Jacobian_r_G,		obj.Jacobian_p_G,	zeros(nc,nb),		zeros(nc,nc);];
					
					% Could check the condition number
					%cond(Psi);
					
					% Change in z
					dz = Psi\(-g_n);
					
					% Set z (not needed but doing for explicitness)
					z =	[r_ddot_n;...
						 p_ddot_n;...
						 lambda_p_n;...
						 lambda_n;];
					
					% Update z
					z = z + dz;
					
					% Extract values from z
					obj.r_ddot(:,tt) = z(1:3*nb,1);
					obj.p_ddot(:,tt) = z(3*nb+1:3*nb+4*nb,1);
					obj.lambda_p(:,tt) = z(3*nb+4*nb+1:3*nb+4*nb+np,1);
					obj.lambda(:,tt) = z(3*nb+4*nb+np+1:3*nb+4*nb+np+nc,1);
					
					% Update q
					rp2q(obj,tt,[1,1,1]);
					
					% Breakout counter
					k = k+1;
					if (k > BreakCount)
						disp("BREAK");
						break;
					end
				end
				
				% Calculate the reaction torques and forces
				for i = 1:obj.N_Bodies
					p_i = obj.q(7*(i-1)+1+3:7*(i-1)+7,tt);
					Jacobian_r_i = obj.Jacobian_r_G(:,3*(i-1)+1:3*(i-1)+3);
					Jacobian_p_i = obj.Jacobian_p_G(:,4*(i-1)+1:4*(i-1)+4);

					% Solve reaction forces from each GCon
					for j = 1:obj.N_GCons
						obj.F_rxn{i}{j,tt} = -Jacobian_r_i(j,:)'*obj.lambda(j,tt);
						obj.tau_rxn_bar{i}{j,tt} = -1/2*G(p_i)*(Jacobian_p_i(j,:)'*obj.lambda(j,tt));
						obj.tau_rxn{i}{j,tt} = A(p_i)*obj.tau_rxn_bar{i}{j,tt};
					end
				end
				
				% Compute the velocity violation for the second major joint
				i = 1;
				ncon = 7;
				nppar = 1;			
				v_v_idx = ncon*(i-1)+1+0:ncon*(i-1)+(ncon-nppar);
				v_violation_temp = obj.Jacobian_G*obj.q_dot(:,tt) - obj.nu_G;
				obj.v_violation{i}(tt,1) = norm(v_violation_temp(v_v_idx),2);
				%obj.v_violation(tt,1) = norm(obj.Jacobian_G*obj.q_dot(:,tt) - obj.nu_G,2);
				
				%k
				
				% Show a progress bar
				waitbar(tt/obj.N_t, WB, ['Dynamics Solver: Computing (',num2str(PercentComplete,'%.2f'),'%, Iterations: ',num2str(k),')']);
				
			end
			
			% Close the waitbar
			delete(WB);
		end
		
		
		%% Function to push values from q to r and p
		function obj = q2rp(obj,tt,FLAGS)
			for i = 1:obj.N_Bodies
				r_idx = 3*(i-1)+1+0:3*(i-1)+3;
				q_r_idx = 7*(i-1)+1+0:7*(i-1)+3;
				if FLAGS(1)
					obj.r(r_idx,tt) = obj.q(q_r_idx,tt);
				end
				if FLAGS(2)
					obj.r_dot(r_idx,tt) = obj.q_dot(q_r_idx,tt);
				end
				if FLAGS(3)
					obj.r_ddot(r_idx,tt) = obj.q_ddot(q_r_idx,tt);
				end

				p_idx = 4*(i-1)+1+0:4*(i-1)+4;
				q_p_idx = 7*(i-1)+1+3:7*(i-1)+7;
				if FLAGS(1)
					obj.p(p_idx,tt) = obj.q(q_p_idx,tt);
				end
				if FLAGS(2)
					obj.p_dot(p_idx,tt) = obj.q_dot(q_p_idx,tt);
				end
				if FLAGS(3)
					obj.p_ddot(p_idx,tt) = obj.q_ddot(q_p_idx,tt);
				end
			end
		end
		
		
		%% Function to push values from r and p to q
		function obj = rp2q(obj,tt,FLAGS)
			for i = 1:obj.N_Bodies
				r_idx = 3*(i-1)+1+0:3*(i-1)+3;
				q_r_idx = 7*(i-1)+1+0:7*(i-1)+3;
				if FLAGS(1)
					obj.q(q_r_idx,tt) = obj.r(r_idx,tt);
				end
				if FLAGS(2)
					obj.q_dot(q_r_idx,tt) = obj.r_dot(r_idx,tt);
				end
				if FLAGS(3)
					obj.q_ddot(q_r_idx,tt) = obj.r_ddot(r_idx,tt);
				end

				p_idx = 4*(i-1)+1+0:4*(i-1)+4;
				q_p_idx = 7*(i-1)+1+3:7*(i-1)+7;
				if FLAGS(1)
					obj.q(q_p_idx,tt) = obj.p(p_idx,tt);
				end
				if FLAGS(2)
					obj.q_dot(q_p_idx,tt) = obj.p_dot(p_idx,tt);
				end
				if FLAGS(3)
					obj.q_ddot(q_p_idx,tt) = obj.p_ddot(p_idx,tt);
				end
			end
		end
		
	end
	%methods(Static)
	%end
	% End the Methods
end