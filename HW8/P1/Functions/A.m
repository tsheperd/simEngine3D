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