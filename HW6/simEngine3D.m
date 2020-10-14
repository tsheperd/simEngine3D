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
			% Parse the input as JSON
			obj.input = jsondecode(obj.inputDeckFile);
			
			% Need to convert the contraints to a cell because matlab
			% sometimes makes them a structure and sometimes a cell but
			% these are accessed differently
			if isstruct(obj.input.constraints)
				struct2cell(obj.input.constraints);
			end
			
			% Parse input times
			obj.t_i = obj.input.time(1);
			obj.dt = obj.input.time(2);
			obj.t_f = obj.input.time(3);
			
			
		end
		
		
		%% Function to act as a kinematic solver
		function obj = KinematicSolver(obj)
			% Create a time vector
			obj.t = (obj.t_i:obj.dt:obj.t_f);
			obj.N_t = length(obj.t);
			
			
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
		
		end
		
		
		
	end
	%methods(Static)
	%end
end