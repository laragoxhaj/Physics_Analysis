% EXP 1

data = load('lab4exp1.txt');

f = data(:,1)*1000;		% Hz
V_c = data(:,2)./1000;	% V
R_s = 390;			% ohms
L = .050;			% H
C = .005*10^(-6);		% F
V_i = 2;			% V

% Q1

V = V_c ./ V_i;

f_L = f(2,1);
f_0 = f(11,1);
f_H = f(31,1);
V_max = V(11,1);
V_2 = V(2,1);

df = f_H-f_L;
Q = f_0 / df;

w_0 = 2*pi*f_0;
r = w_0 * L / Q;

figure(1)
hold on
plot(f(1:32)./1000,V(1:32),'o');
set(gca, 'XTick', [f_L, f_0, f_H]./1000);
set(gca, 'YTick', [V(2,1), V(11,1)]);
plot(f_L*ones(2,1)./1000,V(1:2),':');
plot(f_0*ones(11,1)./1000,V(1:11),':');
plot(f_H*ones(2,1)./1000,V(1:2),':');
plot(f(1:11,:)./1000,V_max*ones(11,1),':');
plot(f(1:31,:)./1000,V_2*ones(31,1),':');
title('Resonance Diagram of Exp. 1 Series RLC Circuit')
xlabel('frequency f, kHz')
ylabel('Voltage ratio Vc/Vi')
hold off


% EXP 2


% Q3
f02 = [10.250000, 10.280000]*1000;	% Hz
fL2 = [10.091000, 9.7500000]*1000;	% Hz
fH2 = [10.368000, 10.710000]*1000;	% Hz

df2 = fH2-fL2;
Q2 = f02 ./ df2;

Rp = 47000;		% ohms
Q2bt = Rp / (2*pi*f02(1,1)*L);



