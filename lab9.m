data = load('data.txt');
gamma = data(:,1);      % wavelength, nm
f = data(:,2)*10^(14);  % frequency, Hz
V = data(:,3:end);      % stopping potential, Volts
v_avg = mean(V,2);      % average stopping potential, volts
v_std = std(V,0,2)     % standard deviation fo stopping potential results

md = fit(f,v_avg,'poly1');          % linear fit
cf = coeffvalues(md);
m = cf(1,1);                        % slope = h/e
ci = confint(md);
ms(1,1) = (ci(2,1)-ci(1,1))/4;      % sigma h/e
ms(1,2) = (ci(2,2)-ci(1,2))/4;      % sigma work function W_o 
% as per the National Institute of Standards & Technology (NIST)
e = 1.6021766208*10^(-19);          % charge of electron, C
e_err = .0000000098*10^(-19);       % error on electron charge, e
h_exp = m*e;                        % experimental h (planck's constant), J*s
h_exp_err = sqrt((ms(1,1)*e)^2 + (e_err*m)^2);      % h error in quadrature

% as per the NIST
h_o = 6.626070040 * 10^(-34);                       % accepted h, J*s
h_o_err = .000000081 * 10^(-34);                    % error on accepted h

sig = abs(h_exp - h_o) / abs(h_exp_err - h_o_err);  % significant (num standard deviations)

W_o = cf(1,2) * e;                                  % work function, Volts
W_o_err = sqrt((ms(1,2)*e)^2 + (e_err*cf(1,2))^2);  % error work function

% Figure 1 plotting + formatting
figure(1)
hold on
errorbar(f,v_avg,v_std,'.','MarkerSize',7,'Color','b');
y = polyval(cf,f);
plot(f,y,'Color','m');
m5 = m*10^(15);
ft = sprintf('fit line = (%0.2fe-15)x - %0.2f', m5, -cf(1,2));
text(6.5*10^14,1,ft);
title('Stopping Potential vs Lightsource Frequency');
yl = sprintf('Potential, V\n');
ylabel(yl);
xl = sprintf('\nFrequency, Hz');
xlabel(xl);
hold off

% Figure 2 plotting + formatting
figure(2)
hold on
res = v_avg-y;
plot(f,res,'o','MarkerSize',12);
%axis([f(1,1),f(end,1),-1,1]);
title('Plot of Stopping Potential Residuals');
xl = sprintf('\nFrequency, Hz');
xlabel(xl);
yl = sprintf('Residual\n');
ylabel(yl);
hold off
