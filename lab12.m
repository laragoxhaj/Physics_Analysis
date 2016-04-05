P_atm = 76.77;   % atmospheric pressure, cmHg
P_atm_err = .01;

data = load('data.txt');    % loading + formatting data

T = data(:,1);      % K
T_err = 0.1;
Pr = data(:,2);     % cmHg
Pr_err = .01;
Pl = data(:,3);     % cmHg
Pl_err = Pr_err;
diff = Pr-Pl;       % cmHg
diff_err = Pr_err;
data(:,4) = diff;
Pabs = P_atm - diff;    % cmHg
Pabs_err = Pr_err*ones(13,1);
Ti = 1./T;          % 1/K
Ti_err = T_err*(1/(T_err^2));

R = 8.314;      % molar gas constant, J/(K*mol)
R_err = .001;

% Figure 1 plotting + formatting
figure(1)
hold on
errorbar(Ti,exp(2.85)*log(Pabs),Pabs_err,'o');
[p1,S1] = fit(Ti,log(Pabs),'poly1');
c1 = coeffvalues(p1);           % coefficients of fit
e1 = confint(p1);               % 95% certainty bounds on coefficients
err1 = (e1(2,:) - e1(1,:))/2;   % errors on coefficients
pfit = c1(1)*Ti+c1(2);          % predicted values of amplitude according to fit
plot(Ti,exp(2.85)*pfit,'r');
title('Absolute Pressure vs. 1/T')
xlabel('Temperature^{-1}, K^{-1}')
ylabel('Absolute Pressure, cmHg')
m = sprintf('slope = (-5.3%c.2)x10^3 K^{-1}',char(177));
text(3.2E-3,60,m);
hold off
c1(2)

% Figure 2 plotting + formatting
figure(2)
hold on
res = pfit - log(Pabs);          % residuals 
plot(Ti,res,'ro','MarkerSize',12);
title('Plot of Linear Fit Residuals');
xlabel('Temperature^{-1}, K^{-1}');
ylabel('Residual');
hold off

% Latent heat of vaporization calculations
m = c1(1,1);        % slope, K
m_err = err1(1,1);
kgmol = .018;       % molar mass, kg/mol
L = -R * m * 1/kgmol / 1000;    % Latent heat, J/g
L_err = sqrt((R_err*m*(1/kgmol)*(1/1000))^2 + (m_err*R*(1/kgmol)*(1/1000))^2);
L_acc = 2258;   % kJ/kg
Lacc_err = 1;
L_sig = abs(L-L_acc)/sqrt(L_err^2 + Lacc_err^2);    % L significance

% Triple point pressure calculations
T0 = 1/273.16;      % 1 / triple point temp, 1/K
T0_err = .01*(1/(273.16^2));
tpp_exp = exp(m*(1/273.16) + c1(2));    % triple point pressure (experimental), cmHg
tppexp_err = sqrt((m_err*T0*tpp_exp)^2 + (T0_err*m*tpp_exp)^2 + (err1(2)*tpp_exp)^2);
tpp0 = .485;    % triple point pressure (accepted), cmHg
tpp0_err = .001;
tpp_sig = abs(tpp_exp-tpp0)/sqrt(tpp0_err^2 + tppexp_err^2);    % TPP significance
c1
L_err
L
