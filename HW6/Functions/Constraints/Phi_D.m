function out = Phi_D(i,a_i_bar,s_i_bar,r_i_dot,r_i,p_i,p_i_dot,...
                      j,a_j_bar,s_j_bar,r_j_dot,r_j,p_j,p_j_dot,...
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
    
    % Calculate distance vector and its time derivative
    d_ij = r_j + A_j*s_j_bar - r_i - A_i*s_i_bar;
    d_ij_dot = r_j_dot + B(p_j, s_j_bar)*p_j_dot - r_i_dot - B(p_i, s_i_bar)*p_i_dot;

    % Calculate either Phi, nu, gamma, or the jacobian based on a flag
    if FLAG == "Phi"
        out = d_ij'*d_ij - f;
    elseif FLAG == "nu"
        out = f_dot;
    elseif FLAG == "gamma"
        out = -2*d_ij'*B(p_j_dot, s_j_bar)*p_j_dot + 2*d_ij'*B(p_i_dot, s_i_bar)*p_i_dot - 2*d_ij_dot'*d_ij_dot + f_ddot;
    elseif FLAG == "jacobian"
        % For i and j not ground
        if j ~= 0 && i ~= 0
            out = [-2*d_ij', -2*d_ij'*B(p_i, s_i_bar), 2*d_ij', 2*d_ij'*B(p_j, s_j_bar)];
        % Drop ground jacobian terms
        elseif j == 0
            out = [-2*d_ij', -2*d_ij'*B(p_i, s_i_bar)];
        elseif i == 0
            out = [2*d_ij', 2*d_ij'*B(p_j, s_j_bar)];
        end
    end

end