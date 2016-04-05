% Part 1 - Period over 10 measurements of free oscillator
T10 = [8.15, 8.43, 7.49];       % time over 10 periods
T10mean = mean(T10);            % mean
T10stdev = std(T10);            % standard deviation
T1 = T10mean/10;                % period
T1_err = T10stdev/10;           % error on period

% Part 2 - Damped oscillator
eq = 46.7;                      % equilibrium position
% position of oscillator
top = [51.9,51.3,50.5,49.7,49.25,48.95,48.65];  % top of oscillation - n/2 T
bot = [41.7,42.6,43.4,44,44.3,44.65,45];        % bottom of oscillation - n T
% normalize position data to equilibrium and combine position data
half = top-eq;
full = eq-bot;
A_T(1:2:14)=half;       
A_T(2:2:14)=full;
T = 0.5:0.5:7;
A_Terr = 0.01*ones(1,14);       % error on position data

% Figure 1 plotting + formatting
figure(1)
hold on
errorbar(T,log(A_T),A_Terr,'.');
[p1,S1] = fit(T',log(A_T'),'poly1');
c1 = coeffvalues(p1);           % coefficients of fit
e1 = confint(p1);               % 95% certainty bounds on coefficients
err1 = (e1(2,:) - e1(1,:))/2;   % errors on coefficients
Afit = c1(1)*T+c1(2);           % predicted values of amplitude according to fit
plot(T,Afit,'r-.');
title('Amplitude of Oscillator at Half-Periods, Linear Fit')
xlabel('Oscillation Period, T')
ylabel('ln(Amplitude)')
m = sprintf('A_T = (%0.2ft + %0.2f) cm', c1(1), c1(2));
text(1,0.5,m);
hold off

% Figure 2 plotting + formatting
figure(2)
hold on
res = Afit - log(A_T);          % residuals 
plot(T,res,'o','MarkerSize',12);
title('Plot of Linear Fit Residuals');
xlabel('Oscillation Period, T');
ylabel('Residual');
hold off

gamma_1 = -2*c1(1);             % gamma estimated using part 1 period data, linear fit
dgdp = -2;
gamma_1_err = sqrt((err1(1)*dgdp)^2); % error on estimated gamma

% Figure 3 plotting + formatting
figure(3)
hold on
errorbar(T,A_T,A_Terr,'x')
[p2,S2] = fit(T',A_T','exp1');
c2 = coeffvalues(p2);           % coefficients of fit
e2 = confint(p2);               % 95% certainty bounds on coefficients
err2 = (e2(2,:) - e2(1,:))/2;   % errors on coefficients
efit = c2(1)*exp(c2(2)*T);      % predicted values of amplitude according to fit
plot(T,efit,'r-.');
title('Amplitude of Oscillator at Half-Periods, Exponential Fit')
xlabel('Oscillation Period, T')
ylabel('Amplitude, cm')
m = sprintf('A_T = (%0.2fe^{(%0.2fs)}) cm', c2(1), c2(2));
text(0.5,2,m);
hold off

% Figure 4 plotting + formatting
figure(4)
hold on
res = efit - A_T;               % residuals
plot(T,res,'o','MarkerSize',12);
title('Plot of Exponential Fit Residuals');
xlabel('Oscillation Period, T');
ylabel('Residual');
hold off

gamma_2 = -2*c2(2);          % gamma estimated using part 1 period data, exponential fit
dgdp = -2;
gamma_2_err = sqrt((err2(2)*dgdp)^2);       % error on estimated gamma

tau = 1 / gamma_2;                          % decay time from exponential fit
dt2dg = -1 / (gamma_2)^2;
tau_err = sqrt((gamma_2_err * dt2dg)^2);    % error on decay time

w_0 = 1/T1;                     % undamped frequency
w_0_err = T1_err*1/T1^2;        % undamped frequency error

Q_0 = w_0 / gamma_2;            % Q factor estimated from exponenital fit, natural frequency
dQdw0 = 1 / gamma_2;
dQdg0 = -w_0 / (gamma_2)^2;
Q_0_err = sqrt( (w_0_err*dQdw0)^2 + (gamma_2_err*dQdg0)^2 );    % error on Q factor

w_1 = sqrt((w_0)^2 - ((gamma_2)/2)^2);          % damped frequency
dw1dw0 = w_0 / sqrt((w_0)^2 - (gamma_2^2)/4);
dw1dg = -gamma_2 / (2*sqrt(4*(w_0)^2 - gamma_2^2));
w_1_err = sqrt( (w_0_err*dw1dw0)^2 + (gamma_2_err*dw1dg)^2 );   % error on damped frequency

Q_1 = w_1 / gamma_2;    % Q factor estimated frrom exponential fit, damped frequency
dQdw1 = 1 / gamma_2;
dQdg1 = -w_1 / (gamma_2)^2;
Q_1_err = sqrt( (w_1_err*dQdw1)^2 + (gamma_2_err*dQdg1)^2 );    % error on Q factor

% relative significance values of
Q_sig = abs(Q_1 - Q_0) / sqrt( (Q_1_err)^2 + (Q_0_err)^2);  % Q factor
w_sig = abs(w_1 - w_0) / sqrt( (w_1_err)^2 + (w_0_err)^2);  % omega
fit_sig = abs(c1(1)-c2(2))/sqrt((err2(2))^2 + (err1(1))^2); % linear & exponential fits
gamma_sig = abs(gamma_1 - gamma_2) / sqrt(gamma_1_err^2 + gamma_2_err^2);   % gamma

w_0
w_0_err