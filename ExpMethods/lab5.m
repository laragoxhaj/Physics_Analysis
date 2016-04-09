% ALL UNITS IN SI


% Part 1 a ----------------------------------------------------------------

% angular velocity w = v * x / r
x = load('part1a.txt');
r = (126.56/ 1000) / 2;             % upper disk radius, m
w = .002 * x(2:end,2:end) ./ r;     % rad/s

% error calculation
x_err = 1;                          % instrumental error
r_err = .01 / 1000 / 2;
w_err = zeros(size(x,1)-1,size(x,2)-1);
dwdx = .002 / r;
dwdr = -.002 * x(2:end,2:end) ./ (r^2);
w_err = sqrt((x_err^2)*(dwdx^2) + (r_err^2)*(dwdr .^ 2));

alpha = zeros(size(w));
alpha_err = zeros(size(w));

% error calculation
t = repmat(x(2:end,1),1,size(w,2));     % vectorized times
dadw = 1 ./ (t(2:end,:) - t(1:end-1,:));
alpha(2:end,:) = (w(2:end,:) - w(1:end-1,:)) .* dadw;
alpha_err(2:end,:) = sqrt((w_err(2:end,:) .^ 2) .* (dadw .^ 2));

% angular acceleration and averaging
alpha_avg = zeros(1,size(w,2));
alpha_avg_err = zeros(1,size(w,2));
a_end = [7,5,4,11];
for i = 1:size(w,2)
    % omit first and last measures (innacurate - low instrumental precision)
    alpha_avg(1,i) = mean(alpha(3:a_end(i),i));             % mean
    alpha_avg_err(1,i) = round(std(alpha(3:a_end(i),i)),2); % standard dev
end

Torque = zeros(2,4);
m_hang = [24.43 49.87 24.43 24.43] / 1000;      % hanging mass - kg
r_p = [24.99 24.99 50.00 24.99] / 1000 / 2;     % pulley radius
Torque(1,:) = m_hang * 9.81 .* r_p;             % F*r
m_err = 0.01 / 1000;
g_err = .01;
dTdm = 9.81 * r_p;
dTdg = m_hang .* r_p;
dTdr = m_hang * 9.81;
Fr_err = sqrt((m_err*dTdm).^2 + (g_err*dTdg).^2 + (r_err*dTdr).^2);     % F*r error

m_d = [1359.00 1359.00 1359.00 (1359.00+1345.17)]/1000;   % kg
m = m_d(1,1);
r_d = [126.56 126.56 126.56 ((126.56+126.32)/2)]/1000/2;
I_d = ((m_d) .* (r_d.^2)) ./ 2;                      % inertia of solid disk
I = I_d(1,1);
dIdm = (r_d.^2) ./ 2;
dIdr = r_d .* m_d;
Id_err = sqrt((m_err .* dIdm).^2 + (r_err .* dIdr).^2);
I_err = Id_err(1,1);
Torque(2,:) = I_d .* alpha_avg;                 % I*alpha
dTdI = alpha_avg;
Ia_err = sqrt((Id_err .* dTdI).^2 + (I_d .* alpha_avg_err).^2);    % I*alpha error

T_err = (Torque(1,:) - Torque(2,:)) ./ Torque(1,:) * 100;      % % difference

% Figure 1 plotting + formatting
sys = [1 2 3 4];
figure(1)
hold on
errorbar(sys,Torque(2,:),Ia_err(1,:),'ro')  % Torque = I * alpha
errorbar(sys,Torque(1,:),Fr_err(1,:),'bo')  % Torque = F * r
legend('I*\alpha','F*r')
title('Comparison of Torque Relation Values of Steel Cylindrical Systems');
xlabel('System')
ylabel('Torque, kg * m^2/s^2')
hold off


% Part 1 b ----------------------------------------------------------------
xvst = load('part1b.txt');
t = xvst(2:end,:);
for i=1:size(t,2)
    t_avg(1,i) = mean(t(:,i));
    t_err(1,i) = std(t(:,i));
end

h = xvst(1,:) / 100;
h_err = 1 / 1000;
v = 2 * h ./ t_avg;
v(1,1) = 0;
dvdt = -2 * h ./ (t_avg.^2);
dvdh = 2 ./ t_avg;
v_err = sqrt((t_err.^2).*(dvdt.^2) + (h_err*dvdh).^2);

r_sp = r_p(1,1);
w = v ./ r_sp;
dwdv = 1 / r_sp;
dwdr = v ./ (r_sp^2);
w_err = sqrt((dwdv*v_err).^2 + (r_err*dwdr).^2);
w_err(1,1) = 0;

m_sp = 24.43 / 1000;
for i=1:size(h,2)
    PE(1,i) = m_sp * (9.81) * -h(1,i);
    dPdm = 9.81 * -h(1,i);
    dPdg = m_sp * -h(1,i);
    dPdh = repmat((m_sp*9.81),size(h(1,i)));
    PE_err(1,i) = sqrt((m_err*dPdm).^2 + (g_err*dPdg).^2 + (h_err*dPdh).^2);
end

KE_mass = m_sp * (v.^2) / 2;
dKmdm = (v.^2)/2;
dKmdv = m_sp * v;
KE_mass_err = sqrt((m_err*dKmdm).^2 + (v_err.^2).*(dKmdv.^2));

KE_disc = I * (w.^2) / 2;
dKddI = (w.^2) / 2;
dKddw = I * w;
KE_disc_err = sqrt((I_err*dKddI).^2 + (w_err.^2).*(dKddw.^2));

E = KE_mass + KE_disc;
E_err = sqrt(KE_mass_err.^2 + KE_disc_err.^2);

E_PE_diff = (PE + E) ./ PE * 100;
E_PE_diff(1,1) = 0;
EvsPE = polyfit(PE,E,1);    % slope = rate of energy decrease

% Figure 2 plotting + formatting
figure(2)
hold on
plot(PE,E)
title('Comparison of Total and Kinetic Energies of System 1');
slope = sprintf('slope = %f', EvsPE(1,1));
legend(slope)
errorbar(PE,E,PE_err,'bo');
xlabel('Potental Energy, kg*m^2/s^2')
ylabel('Total Energy, kg*m^2/s^2')
hold off


% Part 2 ------------------------------------------------------------------

freq = load('part2.txt');

m_d = [1359.00 1345.17 470.27] / 1000;  % disc masses
r_d = [126.56 126.32 126.25] / 2 / 1000;% disc radii
I_d = (m_d .* (r_d .^ 2)) / 2;          % disc moments of iner
dIdm = (r_d .^ 2) / 2;
dIdr = r_d .* m_d;
I_err = sqrt((m_err * dIdm).^2 + (r_err * dIdr).^2);

% pre-collision alphas (apply for all cases)
alpha_d = 2 * 9.81 ./ r_d;
dadg = 2 ./ r_d;
dadr = -2 * 9.81 ./ (r_d .^ 2);
alpha_err = sqrt((g_err * dadg).^2 + (r_err * dadr).^2);

w = zeros(size(freq));
w_err = zeros(size(w));
L = zeros(size(freq,1)/3*2,6);
L_err = zeros(size(L));
L_loss = zeros(size(freq,1)/3,size(L,2));

% j=0: same direction; j=1: opposing directions
for j=0:1,    
    i=j*6+1;
    k=j*4+1;
    l=j*2+1;
    
    % two steel discs & drop pin
    w(i,:) = .002 * freq(i,:) ./ r_d(1,1);              % pre-collision - top steel
    dwdx = repmat((.002 / r_d(1,1)),1,6);
    dwdr = -.002 * freq(i,:) / (r_d(1,1)^2);
    w_err(i,:) = sqrt((x_err^2)*(dwdx.^2) + (r_err^2)*(dwdr.^2));
    
    w(i+1,:) = .002 * freq(i+1,:) ./ r_d(1,2);          % pre-collision - bottom steel
    dwdx = repmat((.002 / r_d(1,2)),1,6);
    dwdr = -.002 * freq(i+1,:) / (r_d(1,2)^2);
    w_err(i+1,:) = sqrt((x_err^2)*(dwdx.^2) + (r_err^2)*(dwdr.^2));
    
    L(k,:) = I_d(1,1) * w(i,:) + I_d(1,2) * w(i+1,:);   % angular momentum - pre-collision
    L_err(k,:) = sqrt((I_err(1,1)*w(i,:)).^2 + (I_err(1,2)*w(i+1,:)).^2 + (w_err(i,:).*(I_d(1,1))).^2  + (w_err(i+1,:).*(I_d(1,2))).^2);
    
    w(i+2,:) = .002 * freq(i+2,:) ./ mean(r_d(1,1:2));  % post-collision - system
    dwdx = repmat((.002 / mean(r_d(1,1:2))),1,6);
    dwdr = -.002 * freq(i+2,:) / (mean(r_d(1,1:2))^2);
    w_err(i+2,:) = sqrt((x_err^2)*(dwdx.^2) + (r_err^2)*(dwdr.^2));
    
    L(k+1,:) = sum(I_d(1,1:2)) * w(i+2,:);                % angular momentum - post-collision
    L_err(k+1,:) = sqrt((I_err(1,1)*w(i+2,:)).^2 + (I_err(1,2)*w(i+2,:)).^2 + (w_err(i+2,:).*sum(I_d(1,1:2))).^2);
    
    L_loss(l,:) = abs(L(k+1,:) - L(k,:)) ./ abs(L(k,:)) * 100;
    
    % aluminum + steel discs & drop pin
    w(i+3,:) = .002 * freq(i+3,:) ./ r_d(1,3);          % pre-collision - top aluminum
    dwdx = repmat((.002 / r_d(1,1)),1,6);
    dwdr = -.002 * freq(i+3,:) / (r_d(1,1)^2);
    w_err(i+3,:) = sqrt((x_err^2)*(dwdx.^2) + (r_err^2)*(dwdr.^2));
    
    w(i+4,:) = .002 * freq(i+4,:) ./ r_d(1,2);          % pre-collision - bottom steel
    dwdx = repmat((.002 / r_d(1,2)),1,6);
    dwdr = -.002 * freq(i+4,:) / (r_d(1,2)^2);
    w_err(i+4,:) = sqrt((x_err^2)*(dwdx.^2) + (r_err^2)*(dwdr.^2));
    
    L(k+2,:) = I_d(1,3) * w(i+3,:) + I_d(1,2) * w(i+4,:);   % angular momentum - pre-collision
    L_err(k+2,:) = sqrt((I_err(1,3)*w(i+3,:)).^2 + (I_err(1,2)*w(i+4,:)).^2 + (w_err(i+3,:).*(I_d(1,3))).^2  + (w_err(i+4,:).*(I_d(1,2))).^2);
   
    w(i+5,:) = .002 * freq(i+5,:) ./ mean(r_d(1,2:3));  % post-collision - system
    dwdx = repmat((.002 / mean(r_d(1,1:2))),1,6);
    dwdr = -.002 * freq(i+5,:) / (mean(r_d(1,1:2))^2);
    w_err(i+5,:) = sqrt((x_err^2)*(dwdx.^2) + (r_err^2)*(dwdr.^2));
    
    L(k+3,:) = sum(I_d(1,2:3)) * w(i+5,:);                % angular momentum - post-collision
    L_err(k+3,:) = sqrt((I_err(1,3)*w(i+3,:)).^2 + (I_err(1,2)*w(i+4,:)).^2 + (w_err(i+5,:).*sum(I_d(1,2:3))).^2);
    
    L_loss(l+1,:) = abs(L(k+3,:) - L(k+2,:)) ./ abs(L(k+2,:)) * 100;
end