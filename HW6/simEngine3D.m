%% Main class for simEngine3D
% Contains all simulation requirements
classdef simEngine3D < handle
	properties

	end
	methods
		%% Constructor function
		function obj = simEngine3D(val)
			if nargin > 0
				obj.Args = val;
			end
		end
	end
	%methods(Static)
	%end
end