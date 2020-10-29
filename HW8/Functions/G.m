function out = G(p)
	e0 = p(1);
	e1 = p(2);
	e2 = p(3);
	e3 = p(4);

	e = [e1,e2,e3]';

	out_1 = [-e, -ToTilde(e) + e0*eye(3)];

	out_2 = [-e1, e0, e3, -e2;...
			 -e2, -e3, e0, e1;...
			 -e3, e2, -e1, e0;];

	out = out_1;
end