% Part 1 - Coriolis Force ----------------------------

% raw data
r_circ = 238 / 1000;                % radius of curvature of the track, m
r_circ_err = 0.5 / 1000;
theta = 142.0;                      % angle bw first and last dots on track, degrees
theta_err = 0.5;

n = 13;                             % total num intervals

d_t = 20 / 1000;                    % time between dots, s
d_t_err = 0.5 / 1000;

s = load('part1.txt') / 100;        % distance traveled along track, m
for i=2:n+1,
    s(i,1) = s(i,1)+s(i-1,1);
end
s_err = 0.05 / 100;
t = reshape(linspace(0, d_t*n, n+1),n+1,1);     % time vector

% calculations
w = pi*(180-theta)/(180*n*d_t);                 % angular velocity, 1/s
dw_dtheta = -pi/(180*n*d_t);
dw_ddelta_t = -w/d_t;
w_err = sqrt((theta_err*dw_dtheta)^2 + (d_t_err*dw_ddelta_t)^2);

cf = fit(t,s,'poly1');
cf_coeff = coeffvalues(cf)
cf_confint = confint(cf)
v_rot = cf_coeff(1);                             % rotational velocity, m/s
v_rot_err = (cf_confint(2,1) - cf_confint(1,1))/2;

% Figure 1 plotting + formatting
figure(1)
% make graph of distance traveled along track versus elapsed time
% slope of graph at center is v_rot
hold on
errorbar(t,s,s_err*ones(size(s,1),1),'.','MarkerSize',7);
cf = fit(t,s,'poly1');
y = polyval(cf_coeff,t);
plot(t,y,'Color','m');
axis([t(1,1),t(size(t,1),1),s(1,1),s(size(s,1),1)])
m = sprintf('v_{rot} = %0.2f %c %0.2f m/s', v_rot, char(177), v_rot_err);
text(0.17,0.17,m);
ti = sprintf('Arc Length Traveled Over Time in Motion Dominated by Coriolis Force\n');
title(ti);
yl = sprintf('\nArc length traveled, m');
ylabel(yl);
xl = sprintf('Time, s\n');
xlabel(xl);
hold off

% Figure 2 plotting + formatting
figure(2)
hold on
res = s-y;
plot(t,res,'+','MarkerSize',7);
title('Plot of Coriolis Residuals');
xl = sprintf('\nTime, s\n');
xlabel(xl);
yl = sprintf('Residual\n');
ylabel(yl);
hold off

% verify expected relation a_rot = a_circ
a_rot = 2 * w * v_rot;
dardw = 2 * v_rot;
dardvr = 2 * w;
a_rot_err = sqrt((w_err*dardw)^2 + (v_rot_err*dardvr)^2);

a_circ = (v_rot)^2 / r_circ;
dacdvr = 2*v_rot / r_circ;
dacdrc = -a_circ/r_circ;
a_circ_err = sqrt((v_rot_err*dacdvr)^2 + (r_circ_err*dacdrc)^2);

% Part 2 - Centrifugal Acceleration -------------------

% measured values
s = load('part2.txt') / 100;            % arc length traveled, m
n = size(s,1);
for i=2:n,
    s(i,1) = s(i,1)+s(i-1,1);
end

theta = 150.0;                          % angle bw first and last dots on track, degrees
d_t = 30 / 1000;                        % time, s
t = reshape(linspace(0, d_t*n, n),n,1);  % time vector

% calculated values

w = (pi/180)*(180 + theta)/((n-1)*d_t); % angular velocity, 1/s
dw_dtheta = -pi/(180*n*d_t);
dw_ddelta_t = -w/d_t;
w_err = sqrt((theta_err*dw_dtheta)^2 + (d_t_err*dw_ddelta_t)^2);

r_p = load('part2rp.txt') / 100;        % radial component perpendicular to t, m
rp_err = 0.05 / 100;

cf = fit(t(2:end,1),r_p,'poly4');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a_err = (cf_confint(2,1) - cf_confint(1,1))/2;
b_err = (cf_confint(2,2) - cf_confint(1,2))/2;
c_err = (cf_confint(2,3) - cf_confint(1,3))/2;
d_err = (cf_confint(2,4) - cf_confint(1,4))/2;
e_err = (cf_confint(2,5) - cf_confint(1,5))/2;

% Figure 3 plotting + formatting
figure(3)
hold on
errorbar(t(2:end,1),r_p,rp_err*ones(size(r_p,1),1),'.','MarkerSize',7);
y = polyval(cf_coeff,t(2:end));
plot(t(2:end),y,'Color','m');
%axis([t(2,1),t(n,1),r_p(1,1),r_p(size(r_p,1),1)])
m = sprintf('fit = %0.2fx^4 + %0.2fx^3 + %0.2fx^2 + %0.2fx + %0.2f m', cf_coeff);
text(0.6,-0.04,m);
ti = sprintf('Angular Radius Over Time in Motion Dominated by Centrifugal Force\n');
title(ti);
yl = sprintf('\nAngular radius, m');
ylabel(yl);
xl = sprintf('Time, s\n');
xlabel(xl);
hold off

% Figure 4 plotting + formatting
figure(4)
hold on
res = r_p-y;
plot(t(2:end),res,'+','MarkerSize',7);
%axis([1,1000,-0.6,0.6])
title('Plot of Centrifugal Residuals for Angular Radius');
xl = sprintf('\nTime, s\n');
xlabel(xl);
yl = sprintf('Residual\n');
ylabel(yl);
hold off

fun = @(x)(cf_coeff(1,1)*x.^4 + cf_coeff(1,2)*x.^3 + cf_coeff(1,3)*x.^2 + cf_coeff(1,4)*x + cf_coeff(1,5));
t0 = solve(fun==0,x);
Fvovf = integral(fun,t0,t(end,1));     % integral, 0->f
Fvovs = integral(fun,t(1,1),t0);     % integral, s->0
dt = [t(end,1)-t0 t0-t(1,1)];
dFdt = cf_coeff(1,1)*dt.^4 + cf_coeff(1,2)*dt.^3 + cf_coeff(1,3)*dt.^2 + cf_coeff(1,4)*dt + cf_coeff(1,5);
F_err = sqrt((a_err/5*dt.^5).^2 + (b_err/4*dt.^4).^2 + (c_err/3*dt.^3).^2 + (d_err/2*dt.^2).^2 + (d_t_err*dFdt).^2);
vovf = w^2 * Fvovf;
dvovfdw = 2*w*Fvovf;
dvdint = w^2;
vovf_err = sqrt((w_err*dvovfdw)^2 + (F_err(1,1)*dvdint)^2);
vovs = w^2 * Fvovs;
dvovsdw = 2*w*Fvovs;
vovs_err = sqrt((w_err*dvovsdw)^2 + (F_err(1,2)*dvdint)^2);

cf = fit(t,s,'poly3');
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a_err = (cf_confint(2,1) - cf_confint(1,1))/2;
b_err = (cf_confint(2,2) - cf_confint(1,2))/2;
c_err = (cf_confint(2,3) - cf_confint(1,3))/2;
d_err = (cf_confint(2,4) - cf_confint(1,4))/2;

syms x
f = cf_coeff(1,1)*x^3 + cf_coeff(1,2)*x^2 + cf_coeff(1,3)*x + cf_coeff(1,4);
y = diff(f,x);
dtv = [t(4,1) t0 t(34,1)];
dsds = vpa(subs(y,x,dtv(1,1)));
dsdo = vpa(subs(y,x,dtv(1,2)));
dsdf = vpa(subs(y,x,dtv(1,3)));
dfda = 3*dtv.^2;
dfdb = 2*dtv;
dfdc = 1;
dfdt = 6*cf_coeff(1,1)*dtv + 2*cf_coeff(1,2);
y_err = sqrt((a_err*dfda).^2 + (b_err*dfdb).^2 + (c_err*dfdc).^2 + (d_t_err*dfdt).^2);

vs = dsds*t(1:9,1)+0.01;
vo = dsdo*t(15:22,1)+0.09;      % v_rotational
vf = dsdf*t(30:38,1)-0.17;
v_per_diff = (dsds-dsdf)/dsds*100

% Figure 5 plotting + formatting
figure(5)
% make graph of distance traveled along track versus elapsed time
% slope of graph at center is v_rot
hold on
errorbar(t,s,s_err*ones(size(s,1),1),'.','MarkerSize',7);
plot(t(1:9,1),vs,'Color','b');
plot(t(15:22,1),vo,'Color','b');
plot(t(30:38,1),vf,'Color','b');
axis([t(1,1),t(size(t,1),1),s(1,1),s(size(s,1),1)])
m = sprintf('v_s = %0.2f %c %0.2f m/s', 0.60, char(177), y_err(1,1));
text(0.2,0.08,m);
m = sprintf('v_{rot} = %0.2f %c %0.2f m/s', double(dsdo), char(177), y_err(1,2));
text(0.5,0.2,m);
m = sprintf('v_f = %0.2f %c %0.2f m/s', double(dsdf), char(177), y_err(1,3));
text(1,0.35,m);
ti = sprintf('Arc Length Traveled Over Time in Motion Dominated by Centrifugal Force\n');
title(ti);
yl = sprintf('\nArc length traveled, m');
ylabel(yl);
xl = sprintf('Time, s\n');
xlabel(xl);
hold off
