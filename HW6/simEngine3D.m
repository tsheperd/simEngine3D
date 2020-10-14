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
	end
	%methods(Static)
	%end
end