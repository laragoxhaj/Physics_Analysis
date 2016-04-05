% load dataset of periods per distance
Tvsd = load('T_vs_d.txt')
dist = Tvsd(:,1);

% actual correct data
T1a = Tvsd(:,2) ./ 4;
T2a = Tvsd(:,3) ./ 4;
% fit linear polynomial
p1a = polyfit(dist,T1a,1);
p2a = polyfit(dist,T2a,1);
% calculate intersection
dist_ia = fzero(@(dist) polyval(p1a-p2a,dist),3)
T_ia = polyval(p2a,dist_ia);

% what we used - prematurely rounded
T1 = round(T1a,3)
T2 = round(T2a,3)
% fit linear polynomial
p1 = polyfit(dist,T1,1);
p2 = polyfit(dist,T2,1);
% calculate intersection
dist_i = fzero(@(dist) polyval(p1-p2,dist),3)
T_i = polyval(p2,dist_i);

% Figure 1 plotting + formatting
figure(1)
hold on
plot(dist,T1,'ro',dist,polyval(p1,dist),'--r')
plot(dist,T2,'b*',dist,polyval(p2,dist),'--b')
plot(dist_i,T_i,'-o','MarkerSize',10)
i_text = sprintf('Intersection:\ndist = %0.3f cm\nT = %0.3f s',dist_i,T_i);
text(dist_i + 0.05, T_i + 0.001, i_text)
legend('T1','T2','T1 linear fit','T2 linear fit','Location','northeast')
title('Period of Pendulum from Opposite Positionings, 1 and 2, for Chosen Distances')
xlabel('Distance, cm')
ylabel('Period, s')
hold off

% Period data for estimated centre of mass positioning
Tdata = load('Tdata.txt');
T1 = Tdata(:,1) ./ 4;
T2 = Tdata(:,2) ./ 4;
m1 = round(mean(T1),5);       % mean periods
m2 = round(mean(T2),5);
s1 = std(T1);	% stddevs
s2 = std(T2);
syse = 0.001 / 4;
binwidth1 = sqrt(s1^2 + syse^2)/2;
binwidth2 = sqrt(s2^2 + syse^2)/2;
nbins1 = round((max(T1) - min(T1)) / binwidth1);   % binwidth = 1/2 stddev
nbins2 = round((max(T2) - min(T2)) / binwidth2);

% Figure 2 plotting + formatting
figure(2)
h1 = histfit(T1,nbins1);    % T1
set(h1(1),'facecolor','r','edgecolor','w','facealpha',0.5);
set(h1(2),'color','r');
hold on
plot([m1,m1],ylim,'r--','Linewidth',2)  % plot mean line
h2 = histfit(T2,nbins2);    % T2
set(h2(1),'facecolor','b','edgecolor','w','facealpha',0.5);
set(h2(2),'color','b');
pm = setstr(177);
text1 = sprintf('mean = %0.5f s\nstdev = %s%0.5f s\nsys error = %s%0.5f s',m1,pm,round(s1,5),pm,syse);
text2 = sprintf('mean = %0.5f s\nstdev = %s%0.5f s\nsys error = %s%0.5f s',m2,pm,round(s2,5),pm,syse);
text(m1+0.00008,4.25,text1);
text(m2+0.00013,3.2,text2);
plot([m2,m2],ylim,'b--','Linewidth',2)  % plot mean line
legend('T1','T1 normal','T1 mean','T2','T2 normal','T2 mean','Location','northeast')
title('Normal Distribution of Pendulum Periods for Opposite Positionings, 1 and 2');
xlabel('Period, s')
ylabel('Frequency')
hold off
