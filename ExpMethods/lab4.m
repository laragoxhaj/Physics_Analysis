% sodium calibration
sod = load('sodium.txt');                   % load sodium beams dataset
center_angle = 85.000;                      % calibration angle
error_theta = (0.5 / 60) / 180 * pi;        % convert to radians
lambda = sod(1);                            % angstrom
theta1 = sod(2) + sod(3)/60 + sod(4)/60;    % sum deg,min,sec
theta2 = sod(5) + sod(6)/60 + sod(7)/60;
offset1 = abs(theta1 - center_angle);       % angle after offset
offset2 = abs(theta2 - center_angle);
mean_off = (offset1 + offset2) / 2;         % mean angle after offset
n = 1;
s = n * lambda / sind(mean_off);            % angstroms/slit
slits = n/(s*10e-8);                        % slits/mm
% partial derivative w/ respect to theta
error_s = sqrt((error_theta * n * lambda * (cosd(mean_off) / (sind(mean_off) ^ 2))) .^ 2)
error_s = roundn(error_s, 1)                % sig figs

% mercury beam angle dataset
dat = load('angles.txt');       % load mercury beam angles dataset
theta = zeros(7,4);             % initialization
lambda = zeros(7,4);
error_lambda = zeros(7,4);
n = repmat([1,1,2,2],7,1);      % for vectorization
for i = 1:4
    theta(:,i) = abs(dat(:,3*i-2)) + dat(:,3*i-1)/60 + dat(:,3*i)/60/60;    % sum deg,min,sec
    lambda(:,i) = s * sind(theta(:,i)) ./ n(:,i);                           % calculate lambda
    % calcualte error on lambda, which is function of theta, s
    ds = error_s * (sind(theta(:,i)) ./ n(:,i));
    dtheta = error_theta * s * cosd(theta(:,i)) ./ (n(:,i));
    error_lambda(:,i) = sqrt((ds .^ 2) + (dtheta .^ 2));
end
error_lambda

std_lambda = zeros(7,1);
mean_lambda = zeros(7,1);
for i = 1:7         % calculate mean and error on lambda
    std_lambda(i) = sqrt(1 / sum(1 ./ (error_lambda(i,:) .^ 2)));
    mean_lambda(i) = (std_lambda(i)^2) * sum(lambda(i,:) ./ (error_lambda(i,:) .^ 2));
end
std_lambda

% percent error of accepted values
true_lambda = [4047; 4358; 4916; 5461; 5770; 5791; 6234];
per_error = (true_lambda - mean_lambda) ./ true_lambda * 100;

% Figure 1 plotting + formatting
figure(1)
hold on
s = scatter(1:7,reshape(true_lambda(1:7,:),1,7),'x');
e = errorbar(1:7,mean_lambda(1:7,:),std_lambda(1:7,:),'o','MarkerSize',7);
text(1.2,true_lambda(1)+50,'Violet');
text(2.2,true_lambda(2)+50,'Blue');
text(3.2,true_lambda(3)+50,'Teal');
text(4.2,true_lambda(4)+50,'Green');
text(5.2,true_lambda(5)+50,'Yellow');
text(6.2,true_lambda(6)+50,'Yellow');
text(7.2,true_lambda(7)+50,'Red');
legend('True wavelengths', 'Experimental wavelengths','Location','southeast');
t = sprintf('Wavelengths of Mercury light beams seen through diffraction grating\n');
title(t);
x = sprintf('\nColor');
xlabel(x)
y = sprintf('Wavelength, %s\n', char(197));
ylabel(y);
hold off
