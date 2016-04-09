S = load('10cm.mn.txt');    % load manually generated time dataset at 10 cm
set1 = S(:,1)       % remove distance labels
m = round(mean(set1),2)      % mean time
s = round(std( set1 ),2)       % stddev on times

set1sim = load('lab2sim.mn.txt'); % load simulated time datasets at 10 cm
full = cat(2,set1,set1sim);     % concatenate manual & simulated datasets

% set binwidth on histogram to 1/2 of the standard deviation
binwidth = s/2;
nbins = round((max(max(full)) - min(min(full))) / binwidth);

% Figure 1 plotting + formatting
figure(1)
h = histfit(reshape(full,[1,5010]),nbins);
set(h(1),'facecolor','b','edgecolor','w','facealpha',0.5);
set(h(2),'color','b');
hold on
plot([m,m],ylim,'b--','Linewidth',2)    % plot mean line
m2 = sprintf('mean = %0.2f s', m);
s2 = sprintf('standard deviation = %0.2f s', s);
text(m + 0.06, 700.5, m2);
text(m + 0.06, 600.5, s2);
title('Normal Time Distribution of Cylinder Rolling 10 cm');
xlabel('Time, s')
ylabel('Frequency')
hold off

fullmean = mean(full);      % mean of each of 501 sets of 10 trials
fmean = round(mean(fullmean),2);     % mean of the means
fstd = round(std(fullmean),2);       % stddev on the means

% set binwidth on histogram to 1/2 of the standard deviation
binwidth = fstd/2;
n2bins = round((max(fullmean) - min(fullmean)) / binwidth);

% Figure 2 plotting + formatting
m2 = sprintf('mean = %0.2f s', fmean);
s2 = sprintf('standard deviation = %0.2f s', fstd);
figure(2)
g = histfit(fullmean, n2bins);
set(g(1),'facecolor','b','edgecolor','w','facealpha',0.5);
set(g(2),'color','b');
hold on
plot([fmean,fmean],ylim,'b--','Linewidth',2)    % plot mean line
text(fmean + 0.015, 70.5, m2);
text(fmean + 0.015, 60.5, s2);
title('Normal Mean Time Distribution of Cylinder Rolling 10 cm');
xlabel('Time, s')
ylabel('Frequency')
hold off

flatfull = full(:);
fmean = round(mean(flatfull),2);     % mean on all 5010 data points
fstd = round(std(flatfull),2);       % stddev on all 5010 data points
for i= 1:7,                 % set bin edges at stddev markers, up to 3 stddevs
    edges(i) = fmean + fstd*(i-4)
end
[N,edges] = histcounts(flatfull,edges)  % get # bin edges on full dataset

sizefull = size(flatfull,1);    % = num rows
for i = 1:3,            % calculate and store % of values within 1,2,3 stddevs
    edges(4-i);
    edges(4+i);
    count = sum(flatfull >= edges(4-i) & flatfull <= edges(4+i));
    fracinstd(i) = count / sizefull;
end

% dsiplay % of values within 1,2,3 stddevs
fracinstd

% Figure 3 plotting + formatting
figure(3)
h = histogram(flatfull,edges);
set(h(1),'facecolor','b','edgecolor','w','facealpha',0.5);
m2 = sprintf('mean = %0.2f s', fmean);
s2 = sprintf('standard deviation = %0.2f s', fstd);
hold on
plot([fmean,fmean],ylim,'b--','Linewidth',2)    % plot mean line
plot(flatfull,normpdf(flatfull,0,fstd*1800));
text(fmean + 0.05, 1000, m2);
text(fmean + 0.05, 900, s2);
title('Normal Time Distribution of Cylinder Rolling 10 cm within 3 Standard Deviations');
xlabel('Time, s')
ylabel('Frequency')
hold off

% calculate accelerations
d = 10.0 * [1:5];       % distances, formatted for vectorization
times(:,1) = set1;      % put 10 cm time points in time vector
for i=2:5,              % put 20-50 cm time points in time vector
    file = sprintf('%d0cm.mn.txt',i);
    S = load(file);
    times(:,i) = S(:,1);
end
t_avg = mean(times)                         % mean times
a_avg = 2 * d ./ (t_avg .^ 2)               % average accelerations at mean times
t_std = std(times);                         % standard deviations on time per distance
dadt = -4 * d ./ (t_avg .^ 3);              % derivative of acceleration at t = t_avg
a_std = sqrt((t_std .^ 2) .* (dadt .^ 2))   % derived errors on accelerations

t_true = [2.0307, 2.8718, 3.5173, 4.0614, 4.5408];  % true time values, given
a_true = 2 * d ./ (t_true .^ 2);            % true acceleration values
a = round(mean(a_true),2);                  % mean of true acceleration values

% Figure 4 plotting + formatting
figure(4)
errorbar(d,a_avg,a_std,'o');
box on
hold on
plot(get(gca,'xlim'), [a a]);               % plot 'true acceleration' line
atext = sprintf('true acceleration = %0.2f cm/s^2', a);
text(40.5, 4.87, atext);
title('Average Acceleration of Rolling Cylinder Per Distance Travelled');
xlabel('Distance, cm')
ylabel('Acceleration, cm/s^2')
hold off

a_std(4) = 0;
a_avg(4) = 0;