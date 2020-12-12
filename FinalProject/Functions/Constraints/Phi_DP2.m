function [Phi, nu, gamma, Jacobian] = Phi_DP2(constraint, q_i, q_j,  q_i_dot, q_j_dot,...
						f,f_dot,f_ddot,FLAG)
	
	i = constraint.i;
	j = constraint.j;
	
	a_i_bar = constraint.a_i_bar;
	%a_j_bar = constraint.a_j_bar;
	
	s_i_bar = constraint.s_i_bar;
	s_j_bar = constraint.s_j_bar;
	
	r_i = q_i(1:3,1);
	r_j = q_j(1:3,1);

	r_i_dot = q_i_dot(1:3,1);
	r_j_dot = q_j_dot(1:3,1);

	p_i = q_i(4:7,1);
	p_j = q_j(4:7,1);
	
	p_i_dot = q_i_dot(4:7,1);
	p_j_dot = q_j_dot(4:7,1);

	
    % Calculate tranformation matrices
    A_i = A(p_i);
    A_j = A(p_j);
    
    % Calculate global a vectors
    a_i = A_i*a_i_bar;
    %a_j = A_j*a_j_bar;
    
    % Calculate global a vector time derivatives
    a_i_dot = B(p_i, a_i_bar)*p_i_dot;
    %a_j_dot = B(p_j, a_j_bar)*p_j_dot;
    
    % Calculate distance vector and its time derivative
    d_ij = r_j + A_j*s_j_bar - r_i - A_i*s_i_bar;
    d_ij_dot = r_j_dot + B(p_j, s_j_bar)*p_j_dot - r_i_dot - B(p_i, s_i_bar)*p_i_dot;

    % Initialize return variables
	Phi = 0;
	nu = 0;
	gamma = 0;
	Jacobian = [zeros(1,3), zeros(1,4), zeros(1,3), zeros(1,4)];
	
	% Calculate either Phi, nu, gamma, or the jacobian based on a flag
    if FLAG(1)
        Phi = a_i_bar'*A_i'*(r_j + A_j*s_j_bar - r_i - A_i*s_i_bar) - f;
    end
    if FLAG(2)
        nu = f_dot;
    end
    if FLAG(3)
        gamma = -a_i'*B(p_j_dot, s_j_bar)*p_j_dot + a_i'*B(p_i_dot, s_i_bar)*p_i_dot - d_ij'*B(p_i_dot, a_i_bar)*p_i_dot - 2*a_i_dot'*d_ij_dot + f_ddot;
    end
    if FLAG(4)
        % For i and j not ground
        if j ~= 0 && i ~= 0
            Jacobian = [-a_i', d_ij'*B(p_i, a_i_bar)-a_i'*B(p_i,s_i_bar), a_i', a_i'*B(p_j, s_j_bar)];
        % Drop ground jacobian terms
        elseif j == 0
            Jacobian = [-a_i', d_ij'*B(p_i, a_i_bar)-a_i'*B(p_i,s_i_bar)];
        elseif i == 0
            Jacobian = [a_i', a_i'*B(p_j, s_j_bar)];
        end
    end

end