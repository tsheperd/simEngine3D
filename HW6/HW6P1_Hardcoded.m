clear;close all;clc;
   

%% Given
r_i = [8, 6, -3]';
p_i = Normalize([4, 3, -5, 1]');
p_i_dot = [-0.2, 1.3, 3.4, 0]';
%p_i_dot(4) = -p_i_dot'*p_i/p_i(4);
p_i_dot(4) = -dot(p_i_dot,p_i)/p_i(4);
p_i_dot = Normalize(p_i_dot);


a_i_bar = [-1.2, 1 ,0.3]';
s_i_P_bar = [0.1, -0.3, 6.0]';

r_j = [-0.5, 1.6, -6.3]';
p_j = Normalize([3.3, -4, 5.1, 6]');

p_j_dot = [0.6, -3.7, 5.1, 0]';
%p_j_dot(4) = -p_j_dot'*p_j/p_j(4)
p_j_dot(4) = -dot(p_j_dot,p_j)/p_j(4);
p_j_dot = Normalize(p_j_dot);

a_j_bar = [1.2, 4.5, 3.1]';
s_j_Q_bar = [0.2, -1.0, 1.5]';

c = [0.3, 0.4, -6]';

f = 1.2;
f_dot = 2.5;
f_ddot = 0.2;

% Added these
r_i_dot = [8, 6, -3]';
r_j_dot = [8, 6, -3]';



%% GCon Calculations
% Flag used to define the output type of the GCON
% Can be "Phi", "nu", "gamma", or "jacobian"
FLAG = "Phi"

disp("DP1");
% Calculate the desired DP1 constrain value
Phi_DP1(1,a_i_bar,r_i,p_i,p_i_dot,...
        2,a_j_bar,r_j,p_j,p_j_dot,...
        f,f_dot,f_ddot,FLAG)

disp("CD");
% Calculate the desired CD constrain value
Phi_CD(c,...
       1,s_i_P_bar,r_i,p_i,p_i_dot,...
       2,s_j_Q_bar,r_j,p_j,p_j_dot,...
       f,f_dot,f_ddot,FLAG)

disp("DP2");
% Calculate the desired DP2 constrain value
Phi_DP2(1,a_i_bar,s_i_P_bar,r_i,r_i_dot,p_i,p_i_dot,...
        2,a_j_bar,s_j_Q_bar,r_j,r_j_dot,p_j,p_j_dot,...
        f,f_dot,f_ddot,FLAG)

disp("D");
% Calculate the desired D constrain value
Phi_D(1,a_i_bar,s_i_P_bar,r_i,r_i_dot,p_i,p_i_dot,...
    2,a_j_bar,s_j_Q_bar,r_j,r_j_dot,p_j,p_j_dot,...
    f,f_dot,f_ddot,FLAG)



%% GCons
% DP1
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


% DP2
function out = Phi_DP2(i,a_i_bar,s_i_bar,r_i_dot,r_i,p_i,p_i_dot,...
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
        out = a_i_bar'*A_i'*(r_j + A_j*s_j_bar - r_i - A_i*s_i_bar) - f;
    elseif FLAG == "nu"
        out = f_dot;
    elseif FLAG == "gamma"
        out = -a_i'*B(p_j_dot, s_j_bar)*p_j_dot + a_i'*B(p_i_dot, s_i_bar)*p_i_dot - d_ij'*B(p_i_dot, a_i_bar)*p_i_dot - 2*a_i_dot'*d_ij_dot + f_ddot;
    elseif FLAG == "jacobian"
        % For i and j not ground
        if j ~= 0 && i ~= 0
            out = [-a_i', d_ij'*B(p_i, a_i_bar)-a_i'*B(p_i,s_i_bar), a_i', a_i'*B(p_j, s_j_bar)];
        % Drop ground jacobian terms
        elseif j == 0
            out = [-a_i', d_ij'*B(p_i, a_i_bar)-a_i'*B(p_i,s_i_bar)];
        elseif i == 0
            out = [a_i', a_i'*B(p_j, s_j_bar)];
        end
    end

end


% D
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

% CD
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


%% Other Functions
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

