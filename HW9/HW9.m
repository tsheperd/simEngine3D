clear;close all;clc;
% S_xi(-1, 0, 0, 2, 1, 1);
% S_xi_xi(-1, 0, 0, 2, 1, 1);
% S_xi_eta(-1, 0, 0, 2, 1, 1);
% S_xi_zeta(-1, 0, 0, 2, 1, 1);

%% Given
rho = 7700;
L = 0.5;
H = 0.003;
W = 0.003;

e_0 =[	0;...
		0;...
		0;...
		1;...
		0;...
		0;...
		0;...
		1;...
		0;...
		0;...
		0;...
		1;...
		L;...
		0;...
		0;...
		1;...
		0;...
		0;...
		0;...
		1;...
		0;...
		0;...
		0;...
		1;];


%% (b)
GQ = [6, 2, 2];

M = Mass(rho, L, W, H, e_0, GQ);

disp("M Diagonals = ");
for ii = 1:size(M)
    disp([num2str(M(ii,ii),'%.2e')]);
end

%% (c)
F_g = [0; 0; -9.81];
M_tot = rho*L*W*H;
F_g_unitMass = F_g/M_tot;


Q_g = GravF(rho, F_g_unitMass, L, W, H, e_0, GQ)


%% (d)
E = 2.0e11;
nu = 0.3;

k_1 = 1;
k_2 = 10*((1+nu)/(12+11*nu));
k_3 = k_2;

% D_test = D_mtx(E, nu, k_1, k_2, k_3)
% D_test_swapped = D_mtx(E, nu, k_3, k_1, k_2)
% S_xiF_test = S_xiF_func(1, 0, 0, L, W, H, e_0)

GQ2 = [5,3,3];

e = [0;...
	0;...
	0;...
	1;...
	0;...
	0;...
	0;...
	1;...
	0;...
	0;...
	0;...
	1;...
	0.5;...
	0;...
	0;...
	0;...
	0;...
	-1;...
	0;...
	1;...
	0;...
	1;...
	0;...
	0];

Q_inter = InternalF(E, nu, k_1, k_2, k_3, L, W, H, e_0, e, GQ2)


%% (e)
eta = 0;
zeta = 0;
ii = 1;
for xi = -1:0.1:1
	r_i(:,ii) = nodal2Global(xi, eta, zeta, L, W, H, e_0);
	r(:,ii) = nodal2Global(xi, eta, zeta, L, W, H, e);
	ii = ii + 1;
end
figure;
hold on;
plot(r_i(1,:), r_i(3,:));
plot(r(1,:), r(3,:));
axis equal;
title("Global Beam Deflection");
xlabel("x (m)");
ylabel("z (m)");
legend('Undeformed','Deformed');
hold off;
hold off;


%% (f)
F_app =	[	10*cosd(45);...
			0;...
			10*sind(45)];

xi_app = 1;
eta_app = 0;
zeta_app = 0;

Q_point = PointF(xi_app, eta_app, zeta_app, L, W, H, F_app)


%% Function to calculate the mass matrix
function M_out = Mass(rho, L, W, H, e_0, GQ)
	M_temp = 0;
	
	% Gaussian Quadrature Weights
	W_xi = Gauss_Quadrature_Weights(GQ(1));
	W_eta = Gauss_Quadrature_Weights(GQ(2));
	W_zeta = Gauss_Quadrature_Weights(GQ(3));
	
	% Gaussian Quadrature loop
	i = 1;
	for xi = Gauss_Quadrature_Locations(GQ(1))
		j = 1;
		for eta = Gauss_Quadrature_Locations(GQ(2))
			k = 1;
			for zeta = Gauss_Quadrature_Locations(GQ(3))
				% Jacobian
				J_xi_0 = [S_xi_xi(xi,eta,zeta,L,W,H)*e_0, S_xi_eta(xi,eta,zeta,L,W,H)*e_0, S_xi_zeta(xi,eta,zeta,L,W,H)*e_0];
				J_xi_0_det = det(J_xi_0);
				
				% Sum mass matrix
				M_temp = M_temp + rho*S_xi(xi,eta,zeta,L,W,H)'*S_xi(xi,eta,zeta,L,W,H)*J_xi_0_det*W_xi(i)*W_eta(j)*W_zeta(k);
				
				k = k+1;
			end
			j = j+1;
		end
		i = i+1;
	end
	
	% Output
	M_out = M_temp;
end


%% Function to calculate the generalized gravity force vector
function G_out = GravF(rho, F_g, L, W, H, e_0, GQ)
	G_temp = 0;
	
	% Gaussian Quadrature Weights
	W_xi = Gauss_Quadrature_Weights(GQ(1));
	W_eta = Gauss_Quadrature_Weights(GQ(2));
	W_zeta = Gauss_Quadrature_Weights(GQ(3));
	
	% Gaussian Quadrature loop
	i = 1;
	for xi = Gauss_Quadrature_Locations(GQ(1))
		j = 1;
		for eta = Gauss_Quadrature_Locations(GQ(2))
			k = 1;
			for zeta = Gauss_Quadrature_Locations(GQ(3))
				% Jacobian
				J_xi_0 = [S_xi_xi(xi,eta,zeta,L,W,H)*e_0, S_xi_eta(xi,eta,zeta,L,W,H)*e_0, S_xi_zeta(xi,eta,zeta,L,W,H)*e_0];
				J_xi_0_det = det(J_xi_0);
				
				% Sum mass matrix
				G_temp = G_temp + rho*S_xi(xi,eta,zeta,L,W,H)'*F_g*J_xi_0_det*W_xi(i)*W_eta(j)*W_zeta(k);
				
				k = k+1;
			end
			j = j+1;
		end
		i = i+1;
	end
	
	% Output
	G_out = G_temp;
end


%% Function to calculate the internal force vector
% Linear Hookian
function Q_int_out = InternalF(E, nu, k_1, k_2, k_3, L, W, H, e_0, e, GQ)
	Q_int_temp = 0;
	Q_integrand = 0;
	
	% Gaussian Quadrature Weights
	W_xi = Gauss_Quadrature_Weights(GQ(1));
	W_eta = Gauss_Quadrature_Weights(GQ(2));
	W_zeta = Gauss_Quadrature_Weights(GQ(3));
	
	% Get constitutive matrix
	D = D_mtx(E, nu, k_1, k_2, k_3);
	
	% Gaussian Quadrature loop
	l = 1;
	for xi = Gauss_Quadrature_Locations(GQ(1))
		m = 1;
		for eta = Gauss_Quadrature_Locations(GQ(2))
			n = 1;
			for zeta = Gauss_Quadrature_Locations(GQ(3))
				% Jacobian
				J_xi_0 = [S_xi_xi(xi,eta,zeta,L,W,H)*e_0, S_xi_eta(xi,eta,zeta,L,W,H)*e_0, S_xi_zeta(xi,eta,zeta,L,W,H)*e_0];
				J_xi_0_det = det(J_xi_0);
				
				% Shape function derivatives accounting for deformed shape
				S_xiF = S_xiF_func(xi, eta, zeta, L, W, H, e_0);
				
				% Find the integrand (at each GQ point)
				Q_integrand = 0;
				for i = 1:3
					for j = 1:3
						Q_integrand = Q_integrand + (S_xiF{i}'*S_xiF{i}*e)*(D(i,j)/2*(e'*S_xiF{j}'*S_xiF{j}*e - 1));
					end
				end
				
				Q_integrand = Q_integrand + ((S_xiF{2}'*S_xiF{3} + S_xiF{3}'*S_xiF{2})*e)*(D(4,4)*(e'*S_xiF{2}'*S_xiF{3}*e));
				Q_integrand = Q_integrand + ((S_xiF{1}'*S_xiF{3} + S_xiF{3}'*S_xiF{1})*e)*(D(5,5)*(e'*S_xiF{1}'*S_xiF{3}*e));
				Q_integrand = Q_integrand + ((S_xiF{1}'*S_xiF{2} + S_xiF{2}'*S_xiF{1})*e)*(D(6,6)*(e'*S_xiF{1}'*S_xiF{2}*e));
				
				Q_integrand = -Q_integrand*J_xi_0_det;
				
				% Sum mass matrix
				Q_int_temp = Q_int_temp + Q_integrand*W_xi(l)*W_eta(m)*W_zeta(n);
				
				n = n+1;
			end
			m = m+1;
		end
		l = l+1;
	end
	
	% Output
	Q_int_out = Q_int_temp;
end


%% Function to convert nodal coordinates to global
function r_global_out = nodal2Global(xi, eta, zeta, L, W, H, e)
	r_global_out = S_xi(xi,eta,zeta,L,W,H)*e;
end


%% Function to get force vector from point load
function Q_point_out = PointF(xi, eta, zeta, L, W, H, F_app)
	Q_point_out = S_xi(xi,eta,zeta,L,W,H)'*F_app;
end


%% Function giving shape function derivatives accounting for deformed shape
function S_xiF_out = S_xiF_func(xi, eta, zeta, L, W, H, e_0)
	% Jacobian and its inverse
	J_xi_0 = [S_xi_xi(xi,eta,zeta,L,W,H)*e_0, S_xi_eta(xi,eta,zeta,L,W,H)*e_0, S_xi_zeta(xi,eta,zeta,L,W,H)*e_0];
	J_xi_0_inv = inv(J_xi_0);
	
	% Get the shape function derivatives for the point
	S_xi_xi_temp = S_xi_xi(xi,eta,zeta,L,W,H);
	S_xi_eta_temp = S_xi_eta(xi,eta,zeta,L,W,H);
	S_xi_zeta_temp = S_xi_zeta(xi,eta,zeta,L,W,H);
	
	% Calculate shape function derivatives while accounting for deformation
	S_xiF_temp{1} = J_xi_0_inv(1,1)*S_xi_xi_temp + J_xi_0_inv(2,1)*S_xi_eta_temp + J_xi_0_inv(3,1)*S_xi_zeta_temp;
	S_xiF_temp{2} = J_xi_0_inv(1,2)*S_xi_xi_temp + J_xi_0_inv(2,2)*S_xi_eta_temp + J_xi_0_inv(3,2)*S_xi_zeta_temp;
	S_xiF_temp{3} = J_xi_0_inv(1,3)*S_xi_xi_temp + J_xi_0_inv(2,3)*S_xi_eta_temp + J_xi_0_inv(3,3)*S_xi_zeta_temp;
	
	% Output
	S_xiF_out = S_xiF_temp;
end


%% Function giving Constitutive Law
function D_out = D_mtx(E, nu, k_1, k_2, k_3)
	% Define some constants for brevity
	C1 = (1-nu)/nu;
	C2 = (1-2*nu)/(2*nu);
	C3 = (E*nu)/((1+nu)*(1-2*nu));

	% Constituative matrix
	D_out = C3*[C1,		1,		1,		0,		0,			0;...
				1,		C1,		1,		0,		0,			0;...
				1,		1,		C1,		0,		0,			0;...
				0,		0,		0,		C2*k_1,	0,			0;...
				0,		0,		0,		0,		C2*k_2,		0;...
				0,		0,		0,		0,		0,			C2*k_3;];
end


%% Function giving locations of Gauss Quadrature locations
function out = Gauss_Quadrature_Locations(n)
	if n == 1
		out = [	0.0000000000000000];
	elseif n == 2
		out = [	-0.5773502691896257,...
				0.5773502691896257];
	elseif n == 3
		out = [	-0.7745966692414834,...
				0.0000000000000000,...
				0.7745966692414834];
	elseif n == 4
		out = [	-0.8611363115940526,...
				-0.3399810435848563,...
				0.3399810435848563,...
				0.8611363115940526];
	elseif n == 5
		out = [	-0.9061798459386640,...
				-0.5384693101056831,...
				0.0000000000000000,...
				0.5384693101056831,...
				0.9061798459386640];
	elseif n == 6
		out = [	-0.9324695142031521,...
				-0.6612093864662645,...
				-0.2386191860831969,...
				0.2386191860831969,...
				0.6612093864662645,...
				0.9324695142031521];
	end
end


%% Function giving locations of Gauss Quadrature locations
function out = Gauss_Quadrature_Weights(n)
	if n == 1
		out = [	2.0000000000000000];
	elseif n == 2
		out = [	1.0000000000000000,...
				1.0000000000000000];
	elseif n == 3
		out = [	0.5555555555555556,...
				0.8888888888888888,...
				0.5555555555555556];
	elseif n == 4
		out = [	0.3478548451374538,...
				0.6521451548625461,...
				0.6521451548625461,...
				0.3478548451374538];
	elseif n == 5
		out = [	0.2369268850561891,...
				0.4786286704993665,...
				0.5688888888888889,...
				0.4786286704993665,...
				0.2369268850561891];
	elseif n == 6
		out = [	0.1713244923791704,...
				0.3607615730481386,...
				0.4679139345726910,...
				0.4679139345726910,...
				0.3607615730481386,...
				0.1713244923791704];
	end
end


%% Function giving normalized nodal coordinates
function S_xi_out = S_xi(xi, eta, zeta, L, W, H)
	S = [(1/4).*(2+(-3).*xi+xi.^3);...
		(1/8).*L.*((-1)+xi).^2.*(1+xi);...
		(-1/4).*eta.*W.*((-1)+xi);...
		(-1/4).*H.*((-1)+xi).*zeta;...
		(1/4).*(2+3.*xi+(-1).*xi.^3);...
		(1/8).*L.*((-1)+xi).*(1+xi).^2;...
		(1/4).*eta.*W.*(1+xi);...
		(1/4).*H.*(1+xi).*zeta];

	S_xi_out = [S(1)*eye(3), S(2)*eye(3), S(3)*eye(3), S(4)*eye(3), S(5)*eye(3), S(6)*eye(3), S(7)*eye(3), S(8)*eye(3)];
end


%% Function giving normalized nodal coordinates derivative wrt xi
function S_xi_xi_out = S_xi_xi(xi, eta, zeta, L, W, H)
	S = [(3/4).*((-1)+xi.^2);...
		(1/8).*L.*((-1)+(-2).*xi+3.*xi.^2);...
		(-1/4).*eta.*W;...
		(-1/4).*H.*zeta;...
		(-3/4).*((-1)+xi.^2);...
		(1/8).*L.*((-1)+2.*xi+3.*xi.^2);...
		(1/4).*eta.*W;...
		(1/4).*H.*zeta];

	S_xi_xi_out = [S(1)*eye(3), S(2)*eye(3), S(3)*eye(3), S(4)*eye(3), S(5)*eye(3), S(6)*eye(3), S(7)*eye(3), S(8)*eye(3)];
end


%% Function giving normalized nodal coordinates derivative wrt eta
function S_xi_eta_out = S_xi_eta(xi, eta, zeta, L, W, H)
	S = [0;...
		0;...
		(-1/4).*W.*((-1)+xi);...
		0;...
		0;...
		0;...
		(1/4).*W.*(1+xi);...
		0];

	S_xi_eta_out = [S(1)*eye(3), S(2)*eye(3), S(3)*eye(3), S(4)*eye(3), S(5)*eye(3), S(6)*eye(3), S(7)*eye(3), S(8)*eye(3)];
end


%% Function giving normalized nodal coordinates derivative wrt zeta
function S_xi_zeta_out = S_xi_zeta(xi, eta, zeta, L, W, H)
	S = [0;...
		0;...
		0;...
		(-1/4).*H.*((-1)+xi);...
		0;...
		0;...
		0;...
		(1/4).*H.*(1+xi)];

	S_xi_zeta_out = [S(1)*eye(3), S(2)*eye(3), S(3)*eye(3), S(4)*eye(3), S(5)*eye(3), S(6)*eye(3), S(7)*eye(3), S(8)*eye(3)];
end

