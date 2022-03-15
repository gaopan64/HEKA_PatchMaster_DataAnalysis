firing = table2array(DA2201171152StimSpike);

name = 'DA2201171152StimSpike';

%% set default parameters--------------------------------------------------

disp('Set default parameters')
% start a timer
tStart = tic;

% change unit to mV
V_Raw = firing*1000;
% V_Raw = firing; % For axograph export

% Global definitions
fontsize = 8;

% Intrinsic properties
sample_rate_input = "Please input the sampling rate(Only input interger and the unit is kHz): \n";
sample_rate_kHz = input(sample_rate_input); % set the sample rate, usually 40kHz for current clamp and 20kHz for voltage clamp
T = 1/sample_rate_kHz;
size_V = size(V_Raw,1);

% Time axes
t = 0:T:(size_V-1) * T;

% Index = int64(t/T + 1)
% t = (Index - 1)*T = Index*T - T

%% run BEADS to denoise----------------------------------------------------
disp('run BEADS to denoise, takes long time, be patient....')
% Filter parameters
fc = 0.006;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)

% Positivity bias (peaks are positive)
r = 6;          % r : asymmetry parameter

% Regularization parameters
amp = 0.08; % default value is 0.04,more noise give less amp
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

% run BEADS function
tBEADS = tic;
[x1, f1, cost] = beads(V_Raw, d, fc, r, lam0, lam1, lam2);
fprintf(1, 'The neuron denoising and baseline correction took %ss\n', num2str(toc(tBEADS)));

% define denoised curve as "V_Clean"
V_Clean = x1+f1;
disp('Denoise and correction finished.')

% Stim spike area
t_spike_start = 100;
Index_spike_start = int64(t_spike_start/T + 1);
t_spike_end = 200;
Index_spike_end =  int64(t_spike_end/T + 1);
t_AHP_end = 400;
Index_AHP_end =  int64(t_AHP_end/T + 1);

%% count spikes------------------------------------------------------------
disp('Start to count spikes.')
% find peaks in "V_Clean"
% "loc" means t value, "i or N" means index
[CleanV_pos_pks,CleanV_pos_loc,half_widths,peak_proms] = findpeaks(V_Clean(Index_spike_start:Index_spike_end),t(Index_spike_start:Index_spike_end),'MinPeakHeight',5,'MinPeakProminence',15);
% define AHP time after peak
t_peak = CleanV_pos_loc;
Index_peak = int64(t_peak/T + 1);
% calculate AHP
CleanV_neg_pks = min(V_Clean(Index_peak:Index_AHP_end));
CleanV_neg_pks_index = find(V_Clean == CleanV_neg_pks);
CleanV_neg_loc = CleanV_neg_pks_index*T - T;

%% Spike Parameters -------------------------------------------------------

disp('Analyze spikes.')
% peak amplitude see CleanV_pos_pks above

% half width see half_widths above

% Afterhyperpolarization see CleanV_neg_pks above

% define threshold 
dCleanV = diff(V_Clean)/T;
dCleanV = smooth(dCleanV,3);
ddCleanV = diff(dCleanV)/T;
ddCleanV = smooth(ddCleanV,3);
dddCleanV = diff(ddCleanV)/T;
dddCleanV = smooth(dddCleanV,3);
t_diff3 = t(1:end-3) + 2*T;
ref_point = CleanV_pos_loc - 1.3; % time
st_pt = int64(ref_point/T); % index
end_pt = int64(CleanV_pos_loc/T); % index
dddCleanV_sub = dddCleanV(st_pt:end_pt);
[dddCleanV_spike_p,dddCleanV_spike_i] = findpeaks(dddCleanV_sub,'MinPeakHeight',1000,'MinPeakProminence',400,'Npeaks',1);
dddCleanV_1st_spike_p = dddCleanV_spike_p;
dddCleanV_spike_i = dddCleanV_spike_i + st_pt - 1; % index
dddCleanV_spike_loc = t_diff3(dddCleanV_spike_i); % peak location on "t_diff3" value(ms)
spike_t = dddCleanV_spike_i + 2; % switch the index of the same time point from t_aavg_diff3 to t_avg axis
CleanV_threshold = V_Clean(spike_t);

%% Plot Waves -------------------------------------------------------------

% constrain the coordinate range
tmax = max(t);
Raw_v_max = max(V_Raw);
Raw_v_min = min(V_Raw);
vlim1 = [Raw_v_min-10 Raw_v_max+10];
tlim1 = [0 tmax];

% plot seperate components
f = figure('WindowState','minimized');
pause(1);
clf

subplot(6, 1, [1 1])
plot(t,V_Raw)
title('Spontaneous Firing','FontSize',fontsize * 1.2)
xlim(tlim1)
ylim(vlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off

subplot(6, 1, [2 2])
plot(V_Raw,'color', [1 1 1]*0.7)
line(1:size_V, f1, 'LineWidth', 1)
legend('Data', 'Baseline')
legend boxoff
title(['Baseline, as estimated by BEADS', ' (r = ', num2str(r), ', fc = ', num2str(fc), ', d = ', num2str(d),')'],'FontSize',fontsize * 1.2)
ylim(vlim1)
set(gca,'ytick', vlim1)
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off

subplot(6, 1, [3 3])
plot(t,x1)
title('Peaks data','FontSize',fontsize * 1.2)
xlim(tlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off

subplot(6, 1, [4 4])
plot(t,V_Raw - x1 - f1)
title('Noise','FontSize',fontsize * 1.2)
xlim(tlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off

% plot "1st_spike" in dddCleanV
subplot(6, 1, [5 5])

plot(t_diff3, dddCleanV)
hold on
plot(dddCleanV_spike_loc, dddCleanV_1st_spike_p, '*c')
title('d^3V/dt^3','FontSize',fontsize * 1.2)
dddCleanV_min = min(dddCleanV);
dddCleanV_max = max(dddCleanV);
dddlim1 = [dddCleanV_min-500 dddCleanV_max+500];
ylim(dddlim1)
hold off
legend('Data', '1st spike')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('d^3V/dt^3','FontSize',fontsize * 1.2);
box off

% plot "CleanV", marks the Peaks,Valleys,and Threshold in curve
subplot(6, 1, [6 6])

plot(t, V_Clean)
hold on
plot(CleanV_pos_loc, CleanV_pos_pks, '^r')
plot(CleanV_neg_loc, CleanV_neg_pks, 'og')
plot(dddCleanV_spike_loc,CleanV_threshold,'xk')
title('Clean Spontaneous Firing','FontSize',fontsize * 1.2)
hold off
xlim(tlim1)
ylim(vlim1)
legend('Data','Peaks','AHPs','Threshold')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off