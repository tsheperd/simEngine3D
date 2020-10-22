function out = B(p, a_bar)
	e0 = p(1);
	e1 = p(2);
	e2 = p(3);
	e3 = p(4);

	e = [e1,e2,e3]';    

	out = 2*[(e0*eye(3) + ToTilde(e))*a_bar, e*a_bar'-(e0*eye(3) + ToTilde(e))*ToTilde(a_bar)];
end