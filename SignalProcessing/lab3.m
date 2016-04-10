data = load('q2.txt');
f = data(:,1);
Vi = data(:,2);
Vr = data(:,3);
y = 20*log10(Vr./Vi)

% Figure 1 plotting + formatting
figure(1)
hold on
%loglog(f,y,'mo','MarkerSize',8)
fi = fit(f,y,'exp1')
plot(fi,f,y,'mo')
h = legend('Bode points','Fitted curve')
set(h,'Box','Off')
set(h,'Location','SouthEast')
title('CR Circuit - Bode Magnitude Plot')
xlabel('Frequency, Hz')
ylabel('Magnitude, dB')
hold off

z_r = 10000;
z_c = 1/(2*pi*1592*(.01e-6));
z_i = sqrt(z_r^2 + z_c^2)
phi = 180-atan(1/(2*pi*1592*10000*(.01e-6)))*180

data = load('q4.txt');
f = data(:,1);
Vi = data(:,2);
Vr = data(:,3);
y = 20*log10(Vr./Vi)

% Figure 2 plotting + formatting
figure(2)
hold on
%semilogx(f,y,'mo','MarkerSize',8)
fi = fit(f,y,'exp1')
plot(fi,f,y,'mo')
h = legend('Bode points','Fitted curve')
set(h,'Box','Off')
set(h,'Location','SouthWest')
title('LR Circuit - Bode Magnitude Plot')
xlabel('Frequency, Hz')
ylabel('Magnitude, dB')
hold off