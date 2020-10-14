function out = Phi_DP1(i,a_i_bar,r_i,p_i,p_i_dot,...
                      j,a_j_bar,r_j,p_j,p_j_dot,...
                      f,f_dot,f_ddot,FLAG)
    % Calculate tranformation matrices
    A_i = A(p_i);
    A_j = A(p_j);
    
    % Calculate global a vectors
    a_i = A_i*a_i_bar;
    a_j = A_j*a_j_bar;
    
    % Calculate global a vector time derivatives
    a_i_dot = B(p_i, a_i_bar)*p_i_dot;
    a_j_dot = B(p_j, a_j_bar)*p_j_dot;
    
    % Calculate either Phi, nu, gamma, or the jacobian based on a flag
    if FLAG == "Phi"
        out = a_i_bar'*A_i'*A_j*a_j_bar - f;
    elseif FLAG == "nu"
        out = f_dot;
    elseif FLAG == "gamma"
        out = -a_i'*B(p_j_dot, a_j_bar)*p_j_dot - a_j'*B(p_i_dot, a_i_bar)*p_i_dot - 2*a_i_dot'*a_j_dot + f_ddot;
    elseif FLAG == "jacobian"
        % For i and j not ground
        if j ~= 0 && i ~= 0
            out = [zeros(1,3), a_j'*B(p_i,a_i_bar), zeros(1,3), a_i'*B(p_j, a_j_bar)];
        % Drop ground jacobian terms
        elseif j == 0
            out = [zeros(1,3), a_j'*B(p_i,a_i_bar)];
        elseif i == 0
            out = [zeros(1,3), a_i'*B(p_j, a_j_bar)];
        end
    end

end