%% Main class for simEngine3D
% Contains all simulation requirements
classdef simEngine3D < handle
	properties
		Args
		inputDeckFileName
		inputDeckFile
		input
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
 			obj.inputDeckFileName = val;
			obj.inputDeckFile = fileread(val);
			obj.input = jsondecode(obj.inputDeckFile);
			Normalize([1 2])
		end
		
		
		
	end
	methods(Static)
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%% Function to calculate A, the rotation matrix
		function out = A(p)
			e0 = p(1);
			e1 = p(2);
			e2 = p(3);
			e3 = p(4);

			e = [e1,e2,e3]';

			out_1 = (2*e0^2 - 1)*eye(3) + 2*(e*e' + e0*ToTilde(e));

			out_2 = 2*[e0^2+e1^2-1/2, e1*e2-e0*e3, e1*e3+e0*e2;...
					 e1*e2+e0*e3, e0^2+e2^2-1/2, e2*e3-e0*e1;...
					 e1*e3-e0*e2, e2*e3+e0*e1, e0^2+e3^2-1/2;];

			out = out_1;
		end

		function out = B(p, a_bar)
			e0 = p(1);
			e1 = p(2);
			e2 = p(3);
			e3 = p(4);

			e = [e1,e2,e3]';    

			out = 2*[(e0*eye(3) + ToTilde(e))*a_bar, e*a_bar'-(e0*eye(3) + ToTilde(e))*ToTilde(a_bar)];
		end


		function out = ToTilde(a)
			out = [ 0,    -a(3),  a(2);...
					a(3),  0,    -a(1);...
				   -a(2),  a(1),  0;];
		end

		function out = FromTilde(a)
			out = [a(3,2);...
				   a(1,3);...
				   a(2,1);];
		end

		function out = Normalize(a)
			out = a./norm(a);
		end
	end
end