% Part 1 ------------------------------------------------------------------
T_room = round(23.0 + 273.15,1);    % room temperatre
w_calorimeter = 208.3;              % weight, g
w_water = 202.0;                    % weight, g
w_err = 0.1;                        % error on weight, g
c_w = 4.18 * w_water;               % water heat capacity, J/K
dc_ww = 4.18;
c_err = w_err * dc_ww;              % error on heat capacity
c_c = 40;                           % calorimeter heat capcity, J/K
c_T = 0;                            % thermometer heat capcity, J/K - we neglect, very small
c_v = c_w + c_c + c_T;              % system heat capacity, J/K
Ti = 286;                           % initial temp
Tf = 306;                           % final temp
T_err = 1;                          % error in temp
dT = Tf - Ti;                       % temp change
dQ = dT * c_v;                      % system heat change
dQdT = c_v;
dQdc_v = dT;
dQ_err = sqrt((c_err*dQdc_v)^2 + (T_err*dQdT)^2);    % error on Q

I = 2.92;                           % current, A
I_err = 0.01;                       % error on A
V = 6.00;                           % voltage, V
V_err = .01;                        % error on V;
t = 17*60 + 49.08;                  % time from Ti->Tf, s
t_err = .01;                        % error on time
W = I * V * t;                      % electrical work on system
dWdI = V * t;
dWdV = t;
dWdt = V;
W_err = sqrt((I_err*dWdI)^2 + (V_err*dWdV)^2 + (t_err*dWdt)^2);  % error on W

% see if W and dQ overlap within their errors
% then for % difference, there is no theoretical value per se, but we
% choose to designate dQ as the 'theoretical' value, as only 1/2 factors in
% calculating dQ are experimentally measured, whereas 2/2 are measured in
% calculating W
per_diff = abs(W-dQ)/dQ;            % % difference

% Part 2 ------------------------------------------------------------------
dat = load('part2.txt');

ni = dat(1,1);                      % initial crank count
nf = dat(size(dat,1),1);            % final crank count
N_tot = nf - ni;                    % total cranks
N = dat(:,1) - ni;                  % adjusted cranks
T = dat(:,2) + 273.15;              % adjust temp, K
Ti = T(1,1);                        % initial temp, K 
Tf = T(size(T,1),1);                % final temp, K
dT = Tf - Ti;                       % change in temp, K
T_err = .1;                         % error in temp
d_cal = 47.03 / 1000;               % calorimeter diameter, m
d_err = .01 / 1000;                 % error in calorimeter diameter, m
c_v = 264;                          % heat capacity of calorimter, J/K
dQ = dT * c_v;                      % change of heat
dQ = T_err;                         % error in chagne of heat
m_sw = 5;                           % mass of suspended weight, kg
alpha = m_sw * 9.81 * d_cal * pi;   % alpha factor
g_err = .01;                        % error on g (9.81)
W = alpha * N_tot;                  % work
dWdg = m_sw * d_cal * pi * N_tot;
dWdd = m_sw * 9.81 * pi * N_tot;
W_err = sqrt((g_err * dWdg)^2 + (d_err * dWdd)^2);  % error on Work


% Figure 1 plotting + formatting
figure(1)
hold on
errorbar(N,T,T_err*ones(size(T,1),1),'.','MarkerSize',7);
cf = fit(N,T,'poly3');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
x = linspace(N(1,1),N(size(N,1),1),20);
y = polyval(cf_coeff,x);
plot(x,y,'Color','m');
plot([N(1,1) N(size(N,1),1)],[T_room T_room],'r--','Linewidth',1)
axis([N(1,1),N(size(N,1),1),Ti-1,Tf+1])
rt = sprintf('Room temp = %0.1f%s%0.1f K', T_room, char(177), T_err);
text(550,295,rt);
m = sprintf('Quadratic fit:\n%0.2fx^2 + %0.2fx + %0.2f', cf_coeff(1), cf_coeff(2), cf_coeff(3));
text(250,290,m);
t = sprintf('Temperature of Calorimeter Per Total Number of Cranks\n');
title(t);
yl = sprintf('\nTemperature, K');
ylabel(yl)
xl = sprintf('Number of rotations\n');
xlabel(xl);
hold off

% Figure 2 plotting + formatting
figure(2)
hold on
res = T-reshape(y,size(y,2),size(y,1));
plot(x,res,'+','MarkerSize',7);
axis([1,1000,-0.6,0.6])
title('Plot of the Residuals');
xl = sprintf('\nNumber of rotations\n');
xlabel(xl);
yl = sprintf('Residual\n');
ylabel(yl);
hold off

delta_TN = cf_coeff(2);             % ignore cf_coeff(1), ~= 0
delta_T = delta_TN * N_tot;         % change in temp
delta_T_err = (cf_confint(2,2)-cf_confint(1,2))/2*N_tot;  % error on delta_T
dQ = delta_T * c_v;                 % change in heat
dQ_err = delta_T_err * c_v;         % error in change in heat

% see if W and dQ overlap within their errors
% then for % difference, there is no theoretical value per se, but we
% choose to designate W as the 'theoretical' value, as only 1/5 factors in
% calculating dQ is experimentally measured (though 2/5 have non-infinite
% precisions), whereas 2/2 are measured (though 1/2 have non-infinite
% precisions) in calculating dQ
per_diff = abs(dQ-W)/W*100;         %  % difference