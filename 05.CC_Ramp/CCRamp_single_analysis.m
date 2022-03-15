CC_Ramp = table2array(PMPulse1621Vmon1);

name = 'SingleNeuron';

tStart = tic;

CC_Ramp = CC_Ramp*1000;
% V_Raw = detrend(V_Raw);
% V_Raw = V_Raw - 60;

% Global definitions
fontsize = 8;

% Intrinsic properties
sample_rate_kHz = 40;
T = 1/sample_rate_kHz;
size_V = size(CC_Ramp,1);

% Time axes
t = 0:T:(size_V-1) * T;

% Index = t/T + 1
% t = (Index - 1)*T = Index*T - 1

N = length(CC_Ramp);

% Filter parameters
fc = 0.006;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)

% Positivity bias (peaks are positive)
r = 6;          % r : asymmetry parameter

% Regularization parameters
amp = 0.05;
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

% run BEADS function
tBEADS = tic;
[x1, f1, cost] = beads(CC_Ramp, d, fc, r, lam0, lam1, lam2);
fprintf(1, 'The neuron denoising and baseline correction takes %ss\n', num2str(toc(tBEADS)));


% constrain the coordinate range
tmax = max(t);
Raw_v_max = max(CC_Ramp);
Raw_v_min = min(CC_Ramp);
vlim1 = [Raw_v_min-10 Raw_v_max+10];
tlim1 = [0 tmax];

% plot seperate components
f = figure('WindowState','minimized');
pause(1);
clf

subplot(5, 1, [1 1])
plot(t,CC_Ramp)
title('Spontaneous Firing','FontSize',fontsize * 1.2)
xlim(tlim1)
ylim(vlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

subplot(5, 1, [2 2])
plot(CC_Ramp,'color', [1 1 1]*0.7)
line(1:N, f1, 'LineWidth', 1)
legend('Data', 'Baseline')
legend boxoff
title(['Baseline, as estimated by BEADS', ' (r = ', num2str(r), ', fc = ', num2str(fc), ', d = ', num2str(d),')'],'FontSize',fontsize * 1.2)
ylim(vlim1)
set(gca,'ytick', vlim1)
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

subplot(5, 1, [3 3])
plot(t,x1)
title('Peaks data','FontSize',fontsize * 1.2)
xlim(tlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

subplot(5, 1, [4 4])
plot(t,CC_Ramp - x1 - f1)
title('Noise','FontSize',fontsize * 1.2)
xlim(tlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

% define denoised curve as "V_Clean"
Ramp_Clean = x1+f1;

% plot "CleanV", marks the Peaks,Valleys,and Threshold in curve
subplot(5, 1, [5 5])

plot(t, Ramp_Clean)
hold on
title('Clean CC_Ramp','FontSize',fontsize * 1.2)
hold off
xlim(tlim1)
ylim(vlim1)
legend('Data')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

% -------------------------------------------------------------------------

fprintf(1, 'The whole process takes %ss\n', num2str(toc(tStart)));
