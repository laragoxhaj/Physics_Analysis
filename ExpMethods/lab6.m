% Part 1 ------------------------------------------------------------------
part1 = load('part1.txt');
N = part1(:,1);                                     % number fringes
d = (part1(:,3) - part1(:,2)) * 10^-4;              % displacement, cm
d_err = 0.1 * 10^-4;                                % error on displacement
lambda = 2 * d ./ N;                                % wavelength, cm
lambda_accepted = 632.8*(10^-7)*ones(size(lambda)); % cm
dldd = 2 ./ N;
dldN = - 2 * d ./ (N.^2);
lambda_err = abs(d_err*dldd);                       % error in wavelength
lambda_mean = mean(lambda);                         % mean wavelength
lambda_mean_err = lambda_err(1,1) / sqrt(size(lambda,1));     % mean wavelength error
per_diff = (lambda-lambda_accepted)./lambda_accepted * 100;   % % difference

lam = lambda*10^7
lam_err = lambda_err*10^7
lam_accepted = lambda_accepted*10^7
lam_mean = lambda_mean*10^7


% Figure 1 plotting + formatting
for i=1:size(part1,1),
    trials(1,i) = i;
end
figure(1)
errorbar(trials,lam,lam_err,'o');
box on
hold on
plot(get(gca,'xlim'), [lam_accepted lam_accepted]);   % plot accepted wavelength line
gtext = sprintf('true wavelength = %0.0f nm\n', lam_accepted(1,1));
text(1.5,633,gtext);
plot(get(gca,'xlim'), [lam_mean lam_mean]);           % plot mean wavelength
mtext = sprintf('mean experimental wavelength = %0.0f nm\n', lam_mean);
text(1.5, 657, mtext);
legend('experimental wavelength')
title('Part 1 Experimental vs. Accepted Wavelength');
xlabel('Trial')
ylabel('Wavelength, nm')
hold off


% Part 2 ------------------------------------------------------------------

part2 =  load('part2.txt');
fringes = part2(:,1);
for i=1:(size(part2,2)-2),
    dP(:,i) = (part2(:,i+2) - part2(:,2));          % pressure differential, cmHg
end
dP_err = 0.1*ones(size(dP,1),1);                    % error on pressure differential, cmHg

for i=1:(size(dP,1))
    dP_avg(i,1) = mean(dP(i,:));                    % average pressure diff, cmHg
end
dP_avg_err = dP_err./ (size(dP,2)-2);               % cmHg


% linear fit: slope = #fringes / change in pressure 
x = zeros(size(dP_avg,1)+1,1);                      % aliasing
x(2:end,:)=dP_avg;
y = zeros(size(fringes,1)+1,1);
y(2:end,:)= fringes;
N = size(fringes,1);
delta = N*sum(x.^2)-sum(x)^2;
m = (N*sum(x.*y)-sum(x)*sum(y))/(delta);            % slope
c = (sum(y)*sum(x.^2)-sum(x)*sum(x.*y))/(delta);    % intercept
a_cu = sqrt(sum((y-m*x-c).^2)/(N-2));               % common uncertainty
m_err = a_cu * sqrt(N/delta);                       % error in the gradient

% Figure 2 plotting + formatting
figure(2)
hold on
plot(dP_avg,fringes)
title('dP vs N_{fringes}');
slope = sprintf('m = %0.2f%s%0.2f fringes/cmHg', m,char(177),m_err);
legend(slope)
errorbar(dP_avg,fringes,dP_avg_err,'bo');
xlabel('Pressure differential, cmHg')
ylabel('Distance differential, cm')
hold off


% calculate n based on experimental wavelength values
l = 3.0;                                            % vacuum cell length, cm
l_err = 0.1;                                        % error on vacuum length cell
m_0 = m * lambda_mean / (2*l);                      % m prime
dm0dm = lambda_mean / (2*l);
dm0dlam = m / (2*l);
dm0dl = m * lambda_mean / (-2*l^2);
m_0_err = sqrt( (m_err*dm0dm).^2 + (lambda_mean_err*dm0dlam).^2 + (lambda_mean_err*dm0dlam).^2);    % error on m prime
alpha = 1;                                          % index of refraction of vacuum
P = 76.0;                                           % pressure - cmHg
P_err = 0.1;                                        % error on pressure
n_exp = m_0 * P + alpha;                            % experimental index of refraction of air
dnexpdm0 = P;
n_exp_err = abs(m_0_err*dnexpdm0);

n_air = 1.000263;
n_exp_per_diff = abs(n_air - n_exp) / n_air * 100;