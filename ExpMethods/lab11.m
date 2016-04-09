a0 = 45.10;      % equilibrium position on ruler
a0_err = 0.01;

% data
data = sortrows(load('data11.txt'),5);
T5 = data(:,1);
a = data(:,2);
a_err = 0.05;
t = data(:,3);
t_err = 0.001;
T = data(:,4);
T_err = 0.01;
f = data(:,5);
f_err = 0.01;
a_net = data(:,6);      % = |a - a0|
a_net_err = a_err;
theta = data(:,7);
theta_err = 0.01;       % since theta was provided to us, we determine its error
                        % by the precision of the data provided, rather than by its formula 
w = 2*pi*f;
w_err = 2*pi*f_err;

% get parameters from custom fit - gamma, p, w_0
afittype = fittype('((p*w0^2)/sqrt(((w0)^2-w^2)^2+(gamma*w)^2))','dependent',{'a'},...
                   'independent',{'w'},'coefficients',{'p','gamma','w0'});
a_fit = fit(w,a_net,afittype);
a_coeff = coeffvalues(a_fit)';
a_confint = confint(a_fit);
a_err = ((a_confint(2,:) - a_confint(1,:))/2)';
p = a_coeff(1);
p_err = a_err(1);
gamma_a = abs(a_coeff(2));  % use abs val as fit may produce gamma < 0 (legitimately)
gamma_a_err = a_err(2);
w0_a = a_coeff(3);          % use abs val as fit may producem w0 < 0 (legitimately)
w0_a_err = a_err(3);

% get parameteres from custom fit - gamma, w_0
tfittype = fittype('atan2((-gamma*w),((w_0)^2-w^2))',...
                   'dependent',{'theta'},'independent',{'w'},...
                   'coefficients',{'gamma','w_0'});
theta_fit = fit(w,theta,tfittype,'StartPoint',[gamma_a,w0_a]);
t_coeff = abs(coeffvalues(theta_fit))'; % use abs val as fit may produce gamma, w0 < 0 (legitimately)
t_confint = confint(theta_fit);
t_err = ((t_confint(2,:) - t_confint(1,:))/2)';
gamma_t = t_coeff(1);
gamma_t_err = t_err(1);
w0_t = t_coeff(2);
w0_t_err = t_err(2);

% use parameters from fit of amplitude graph
gamma = mean([gamma_a, gamma_t]);
gamma_err = sqrt(gamma_a_err^2 + gamma_t_err^2);
w0 = mean([w0_a, w0_t]);
w0_err = sqrt(w0_a_err^2 + w0_t_err^2);

% plot A vs w, A vs theta (same graph)
figure(1)
hold on
grid on
plot(w,a_net,'bo');
plot(a_fit,'c');
plot(w,theta,'m*');
plot(theta_fit,'r');
title('Amplitude and Phase Values vs. Driving (Angular) Frequency of a Driven, Damped Oscillator')
xlabel('Angular Frequency w, s^{-1}')
ylabel('')
ampfit = sprintf('a = (%0.2f*%0.2f^2)/sqrt((%0.2f^2-w^2)^2+(%0.2fw)^2))',p,w0_a,w0_a,gamma_a);
tfit = sprintf('theta = arctan((-%0.2fw)/(%0.2f^2-w^2))',gamma_t,w0_t);
legend('Amplitude, a (cm)',ampfit,'Phase, theta (rad)',tfit);
hold off

% Figure 2 plotting + formatting
figure(2)
hold on
res = ((p*w0_a^2)./sqrt(((w0_a)^2-w.^2).^2+(gamma_a*w).^2)) - a_net;    % residuals 
plot(w,res,'o','MarkerSize',12);
title('Plot of Amplitude Fit Residuals');
xlabel('Amplitude, cm');
ylabel('Residual');
hold off

% Figure 3 plotting + formatting
figure(3)
hold on
res = abs(atan2((-gamma_t*w),((w0_t)^2-w.^2))) - theta;	% residuals 
plot(w,res,'o','MarkerSize',12);
title('Plot of Theta Fit Residuals');
xlabel('Theta, s^{-1}');
ylabel('Residual');
hold off

% three estimates of Q
ar = ((p*w0_a)/gamma);      % resonance amplitude
dardp = w0_a/gamma;
dardw0 = p/gamma;
dardg = -p*w0_a/(gamma^2);
ar_err = sqrt((p_err*dardp)^2 + (w0_err*dardw0)^2 + (gamma_err*dardg)^2);
Q1 = ar/p;              % Q estimated from resonance and natural amplitudes
dQ1dar = 1/p;
dQ1dp = -ar/(p^2);
Q1_err = sqrt((ar_err*dQ1dar)^2 + (p_err*dQ1dp)^2);

aa = 1.75;
b = (gamma_a^2 - 2*w0_a^2);
dbdg = 2*gamma_a;
dbdw = -4*w0_a;
b_err = sqrt((gamma_a_err*dbdg)^2 + (w0_a_err*dbdw)^2);
c = (gamma_a^4 - 4*(gamma_a*w0_a)^2);
dcdg = 4*gamma_a^3 - 8*gamma_a*w0_a^2;
dcdw = -8*w0_a*gamma_a^2;
c_err = sqrt((gamma_a_err*dcdg)^2 + (w0_a_err*dcdw)^2);
w2 = sqrt((-b+sqrt(b^2 -4*aa*c))/(2*aa));   % quadratic dunction
dw2da = -sqrt((-b+sqrt(b^2 -4*aa*c))/(2*aa^2))-c/(aa*sqrt(b^2-4*aa*c));
dw2db = ((b/sqrt(b^2-4*aa*c))-1)/(2*aa);
dw2dc = -1/sqrt(b^2-4*aa*c);
w2_err = sqrt((b_err*dw2db)^2 + (c_err*dw2dc)^2);
dw = 2*abs(w0_a-w2);
dw_err = sqrt((2*w0_a_err)^2 + (-2*w2_err)^2);
Q2 = sqrt(3) * w0_a / dw;    % Q estimated from FWHM
dQ2dw0 = sqrt(3)/dw;
dQ2ddw = -sqrt(3) * w0_a / dw^2;
Q2_err = sqrt((w0_a_err*dQ2dw0)^2 + (dw_err*dQ2ddw)^2);

Qtau = w0/gamma;        % Q given decay time 1/gamma
dQtdw0 = 1/gamma;
dQtdg = -w0/(gamma^2);
Qtau_err = sqrt((w0_err*dQtdw0)^2 + (gamma_err*dQtdg)^2);

% compare w/ Q from measurement of damping time in part 2
Q1_sig = abs(Q1 - Qtau) / sqrt( (Q1_err)^2 + (Qtau_err)^2);  % Q factor
Q2_sig = abs(Q2 - Qtau) / sqrt( (Q2_err)^2 + (Qtau_err)^2);  % Q factor

[Q2_err]
