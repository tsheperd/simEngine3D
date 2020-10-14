function out = Phi_CD(c,...
                      i,s_i_bar,r_i,p_i,p_i_dot,...
                      j,s_j_bar,r_j,p_j,p_j_dot,...
                      f,f_dot,f_ddot,FLAG)
    % Calculate tranformation matrices
    A_i = A(p_i);
    A_j = A(p_j);
    
    % Calculate either Phi, nu, gamma, or the jacobian based on a flag
    if FLAG == "Phi"
        out = c'*(r_j + A_j*s_j_bar - r_i - A_i*s_i_bar) - f;
    elseif FLAG == "nu"
        out = f_dot;
    elseif FLAG == "gamma"
        out = c'*B(p_i_dot, s_i_bar)*p_i_dot - c'*B(p_j_dot, s_j_bar)*p_j_dot + f_ddot;
    elseif FLAG == "jacobian"
        % For i and j not ground
        if j ~= 0 && i ~= 0
            out = [-c', -c'*B(p_i, s_i_bar), c', c'*B(p_j, s_j_bar)];
        % Drop ground jacobian terms
        elseif j == 0
            out = [-c', -c'*B(p_i, s_i_bar)];
        elseif i == 0
            out = [c', c'*B(p_j, s_j_bar)];
        end
    end

end