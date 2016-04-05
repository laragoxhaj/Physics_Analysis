% CLEMENT-DESORMES
% theoretical gamma values
yta = 1.403;
ytc = 1.304;
ytr = 1.668;
y_err = .001;

data = load('data.txt');
air = data(:,1:4);
co2 = data(:,5:8);
arg = data(:,9:12);
h_err = .01*ones(11,1);
h_Err = .01*ones(10,1);

% Figure 1 plotting + formatting
figure(1)
hold on
errorbar(air(1:10,3),air(1:10,1),h_Err,'bo');
errorbar(co2(1:10,3),co2(1:10,1),h_Err,'m*');
errorbar(arg(1:10,3),arg(1:10,1),h_Err,'rx');
a = fit(air(1:11,3),air(1:11,1),'poly1');
c = fit(co2(1:11,3),co2(1:11,1),'poly1');
r = fit(arg(1:11,3),arg(1:11,1),'poly1');
ac = coeffvalues(a);		% coefficients of fit
ae = confint(a);
aerr = (ae(:,2)-ae(:,1))/2;		% errors on coefficients
cc = coeffvalues(c);
ce = confint(c);
cerr = (ce(:,2)-ce(:,1))/2;
rc = coeffvalues(r);
re = confint(r);
rerr = (re(:,2)-re(:,1))/2;
afit = ac(1)*air(1:11,3) + ac(2);
cfit = cc(1)*co2(1:11,3) + cc(2);
rfit = rc(1)*arg(1:11,3) + cc(2);
plot(air(1:11,3),afit,'b');
plot(co2(1:11,3),cfit,'m');
plot(arg(1:11,3),rfit,'r');
title('Manometer readings h0 vs. h0-h1 for the Determination of Respective Adiabatic Coefficients')
xlabel('h0-h1, cmHg')
ylabel('h0, cmHg')
legend('Air','CO_2','Argon');
m = sprintf('m_{Air} = %0.1f(h0-h1)+%0.3f',ac(1),ac(2));
text(15,10,m);
m = sprintf('m_{CO_2} = %0.1f(h0-h1)+%0.1f',cc(1),cc(2));
text(15,7,m);
m = sprintf('m_{Argon} = %0.1f(h0-h1)+%0.2f',rc(1),rc(2));
text(15,4,m);
axis([0,25,0,30]);
hold off

% Residuals plotting + formatting
figure(2)
hold on
ares = afit(1:10,:)-air(1:10,1);	% residuals 
plot(air(1:10,1),ares,'bo','MarkerSize',12);
cres = cfit(1:10,:)-co2(1:10,1);	% residuals 
plot(co2(1:10,3),cres,'m*','MarkerSize',12);
rres = rfit(1:10,:)-arg(1:10,1);	% residuals 
plot(arg(1:10,3),rres,'rx','MarkerSize',12);
legend('Air','CO_2','Argon');
title('Plot of Air, CO_2, and Argon Fit Residuals');
xlabel('h0-h1, cmHg');
ylabel('Residual');
axis([0,20,-1.5,1.5]);
hold off

ay = ac(1);     % adiabatic coefficient of air
ay_err = aerr(1);
cy = cc(1);     % adiabatic coefficient of co2
cy_err = cerr(1);
ry = rc(1);     % adiabatic coefficient of argon
ry_err = rerr(1);

ay_sig = abs(ay-yta)/sqrt(ay_err^2 + y_err^2);
cy_sig = abs(cy-ytc)/sqrt(cy_err^2 + y_err^2);
ry_sig = abs(ry-ytr)/sqrt(ry_err^2 + y_err^2);

% RUCHHARDT
T_dat = [5.90,5.94,5.97,6.15,6.02,5.94,6.02,5.94,6.08,5.96];
T = mean(T_dat);
T_err = std(T_dat);
m = .005620;		% mass of bob, kg
m_err = .000001;
V = 9.8*.001;   % volume of jar, m^3
V_err = .1*.001;
d = .0166;		% diameter of tube, m
d_err = .0001;
A = pi*(d/2)^2;		% area of tube, m
A_err = d_err*(pi*d/2);
P0 = 753.06*(133.322365);   % atmospheric pressure on day of experiment, Pa
P0_err = sqrt((.01*(133.322365))^2 + (.000001*753.06)^2);
g = 9.81;		% acceleration due to gravity, m/s^2
g_err = .01;
P = P0 + (m*g/A); % pressure in jar, Pa
P=P/100;
P_err = sqrt((P0_err/100)^2 + (m_err*g/A/100)^2 + (g_err*m/A/100)^2 + (A_err*(-m*g/(A^2)/100))^2);
ay2 = (4*(pi^2)*m*V)/((A^2)*P*(T^2));		% adiabatic coefficient
dydm = (4*(pi^2)*V)/((A^2)*P*(T^2));
dydV = (4*(pi^2)*m)/((A^2)*P*(T^2));
dydA = -2*(4*(pi^2)*m*V)/((A^3)*P*(T^2));
dydP = -(4*(pi^2)*m*V)/((A^2)*(P^2)*(T^2));
dydT = -2*(4*(pi^2)*m*V)/((A^2)*P*(T^3));
ay2_err = sqrt((m_err*dydm)^2 + (V_err*dydV)^2 + (A_err*dydA)^2 + (P_err*dydP)^2 + (T_err*dydT)^2);

y12_sig = abs(ay2-ay)/sqrt(ay2_err^2 + ay_err^2);
yt2_sig = abs(ay2-yta)/sqrt(ay2_err^2 + y_err^2);