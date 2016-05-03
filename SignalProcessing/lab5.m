
data = load('lab5exp1.txt');

f = data(:,1);		
w = 2*pi.*f./1000;	
V_in = 5*ones(size(f,1),1);
V_out = data(:,3);	

f_L = f(23,1);
f_0 = f(32,1);
f_H = f(43,1);
fx = log10(f./f_H);
vee = V_out./V_in;
V = 20*log10(vee);
w = 2*pi*f;

f = fit(w, vee, '(x^2)/(x^2 + (a/b)*x + 1/(b*c))', 'StartPoint', [0, 0, 0])
% a = R, b = L, c = C
