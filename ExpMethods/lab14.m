data = load('data.txt');

% raw data
I = data(:,6);
I_err = .00001;		% precision of measurements
V = data(:,7);
V_err = .00001;		% precision of measurements

% intermediate calculations: resistance, temperature
RE = V./I
RE_err = sqrt((V_err*1./I).^2 + (I_err*V./(I.^2)).^2);
R0 = .271*ones(size(RE,1),1);
R0_err = .001*ones(size(RE,1),1);
T0 = 22.6*ones(size(RE,1),1);
T0_err = .1*ones(size(RE,1),1);
a = .0045;
a_err = .0001;
T_f = (RE - R0)./(R0*a) + T0;
T_f_err = sqrt((RE_err.*(1./(R0*a))).^2 + (R0_err.*(RE./(a*R0.^2))).^2 + (a_err*(RE-R0)./(R0*a^2)).^2 + T0_err.^2);

% final calculations: power, temparture difference
W = I.*V;
W_err = sqrt((I_err*V).^2 + (V_err*I).^2);
T_room = 22.2*ones(size(T_f,1),1);
T_room_err = .1*ones(size(T_f,1),1);
T = T_f - T_room;
T_err = sqrt(T_f_err.^2 + T_room_err.^2);
Bf = fit(T(1:7,:),W(1:7,:),'poly1');
B = coeffvalues(Bf);
Be = confint(Bf);
B_err =  (Be(2,:)-Be(1,:))/2;
WNR = B(1,1)*T + B(1,2)*ones(size(T0,1),1);
WNR_err = sqrt((B_err(1,1)*T).^2 + (B(1,2)*ones(size(T0,1),1)).^2);
WR = W-WNR;
WR_err = sqrt(W_err.^2 + WNR_err.^2);
sigma = 5.67*10^(-12);	% W/(cm^2 * K^4)
sigma_err = .01*10^(-12);
C = 5E-2;	% convenient constant
WB = C*sigma*(T_f.^4);
WB_err = sqrt((T_f_err.*(C*sigma*4*T_f.^3)).^2 + (sigma_err*C*T_f.^4).^2);

f = fit(T(19:end,:),WR(19:end,:),'poly1');
m = coeffvalues(f);
me = confint(f);
me = (me(2,1)-me(1,1))/2;

% Figure 1 plotting + formatting
figure(1)
hold on
loglog(T,W,'bo')
loglog(T,WNR,'m-.')			% full
loglog(T_f(19:end,:),WR(19:end,:),'r.')
WR_line = m(1,1)*T_f + m(1,2);
loglog(T,WR_line,'r-')
loglog(T(1:7,:),WNR(1:7,:),'m')		% non-extrapolated only
axis([0,3500,-5,35])
title('Power at Various Temperature Difference Points for Stefan-Boltzmann Lamp')
xlabel('T-Ts, C')
ylabel('W, Watts')
legend('T-Ts vs. W','T-Ts vs W_{NR}','T-Ts vs W_R','Location','northwest');
errorbar(T,W,W_err,'bo')
errorbar(T,WNR,WNR_err,'m-.')
errorbar(T_f(19:end,:),WR(19:end,:),WR_err(19:end,:),'r.')
mt = sprintf('m = %0.3f %c %0.3f', m(1), char(177), me(1));
text(2600,13,mt);
hold off

W_ratio = WR./WB
W_ratio_err = sqrt((WR_err./WB).^2 + (WB_err.*WR_err./(WB_err.^2)).^2)

% Figure 2 plotting + formatting - WR/WB log-log plot
figure(2)
hold on
loglog(T_f(2:end,:)+273.15*ones(size(T_f,1)-1,1),W_ratio(2:end,:),'b*')
errorbar(T_f(19:end,:)+273.15*ones(size(T_f,1)-18,1),W_ratio(19:end,:),W_ratio_err(19:end,:),'b*')
loglog(T_f(19:end,:)+273.15*ones(size(T_f,1)-18,1),ones(size(T,1)-18,1),'c-')
legend('W_R/W_B','Expected Ratio','Location','Southeast')
title('Log-Log Power Ratio at Various Temperature Points for Stefan-Boltzmann Lamp')
xlabel('T, K')
ylabel('W_{R}/W_{B}, Watts')
hold off

% Figure 3 plotting + formatting - WR/WB linear plot
figure(3)
hold on
plot(T_f(2:end,:)+273.15*ones(size(T_f,1)-1,1),W_ratio(2:end,:),'b*')
errorbar(T_f(19:end,:)+273.15*ones(size(T_f,1)-18,1),W_ratio(19:end,:),W_ratio_err(19:end,:),'b*')
plot(T_f(19:end,:)+273.15*ones(size(T_f,1)-18,1),ones(size(T,1)-18,1),'c-')
legend('W_R/W_B','Expected Ratio','Location','Southeast')
title('Linear Power Ratio at Various Temperature Points for Stefan-Boltzmann Lamp')
xlabel('T, K')
ylabel('W_{R}/W_{B}, Watts')
hold off
