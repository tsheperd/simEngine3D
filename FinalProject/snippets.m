% Number of compound constraints (CCons)
N_CCons = size(simulation.input.constraints_Compound,1)

% If there are compound constraints extract their GCons
if N_CCons > 0
	% Initialize variable to store the number of GCons per CCon
	N_CCons_GCons = zeros(N_CCons,1);
	
	% For each compound constraint extract its GCons
	for CC = 1:N_CCons
		% Number of GCons per CCon
		N_CCons_GCons(CC,1) = size(simulation.input.constraints_Compound(CC).GCons,1);
				
		% For each GCon per CCon
		for CC_CG = 1:N_CCons_GCons(CC,1)
			% JSONDecode is annoying with structs and cells
			if isstruct(simulation.input.constraints_Compound(CC).GCons)
				simulation.input.constraints{end+1,1} = simulation.input.constraints_Compound(CC).GCons(CC_CG);
			else
				simulation.input.constraints{end+1,1} = simulation.input.constraints_Compound(CC).GCons{CC_CG};
			end
		end
	end

	% Total number of GCons from compound constraints
	N_CCons_GCons_tot = sum(N_CCons_GCons);
	
end
