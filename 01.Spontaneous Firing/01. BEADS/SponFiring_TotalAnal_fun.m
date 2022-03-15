function [half_widths_mean,half_widths_CV,firing_rate,spikes_inter_mean,spikes_inter_CV,peak_amp_mean,peak_amp_CV,AHP_mean,AHP_CV,threshold_mean,threshold_CV,avg_half_widths,avg_spike_pks,avg_spike_neg_pks,avg_spike_neg_pks_3ms,avg_spike_neg_pks_20ms,avg_spike_neg_pks_150ms,avg_spike_threshold,t_avg,avg_spike,avg_spike_dt,d_avg_spike_dt,slope_th,range] = SponFiring_TotalAnal_fun(V_Raw,AHP_length,PrePeak_length,AHP_extend,PrePeak_extend,T,name)

%% set default parameters--------------------------------------------------
fprintf(1, 'start to analysis: %s\n',name)
disp('Set default parameters')
V_Raw = V_Raw*1000;
V_Raw = detrend(V_Raw);
V_Raw = V_Raw - 60;

% Global definitions
fontsize = 8;

size_V = size(V_Raw,1);

% Time axes
t = 0:T:(size_V-1) * T;

% Index = t/T + 1
% t = (Index - 1)*T = Index*T - 1

%% run BEADS to denoise----------------------------------------------------
disp('run BEADS to denoise, takes long time, be patient....')

% Filter parameters
fc = 0.006;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)

% Positivity bias (peaks are positive)
r = 6;          % r : asymmetry parameter

% Regularization parameters
amp = 0.04; % default value is 0.04,more noise give less amp
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

% run BEADS function
tBEADS = tic;
[x1, f1, ~] = beads(V_Raw, d, fc, r, lam0, lam1, lam2);
fprintf(1, 'The neuron denoising and baseline correction takes %ss\n', num2str(toc(tBEADS)));

% define denoised curve as "V_Clean"
V_Clean = x1+f1;
disp('Denoise and correction finished.')

%% count spikes------------------------------------------------------------
disp('Start to count spikes.')
% find peaks in "V_Clean"
% "loc" means t value, "i" means index
[CleanV_pos_pks,CleanV_pos_loc,half_widths,~] = findpeaks(V_Clean,t,'MinPeakHeight',5,'MinPeakProminence',15,'MinPeakDistance',5);
% [CleanV_neg_pks_N,CleanV_neg_loc] = findpeaks(-V_Clean,t,'MinPeakHeight',60,'MinPeakProminence',10);
% CleanV_neg_pks = -CleanV_neg_pks_N; 
CleanV_neg_pks_count = numel(CleanV_pos_pks)-1;% AHP is always between in two peaks,so number of AHPs is one less than number of peaks
CleanV_neg_pks = zeros(CleanV_neg_pks_count,1);% preallocate the matrix
CleanV_neg_loc = zeros(1,CleanV_neg_pks_count);
for i = 1:CleanV_neg_pks_count 
    t_i = CleanV_pos_loc(i:i+1);
    N_i = int64(t_i/T + 1);
    CleanV_peaks_interval_i = V_Clean(N_i(1,1):N_i(1,2));
    CleanV_neg_pks(i,1) = min(CleanV_peaks_interval_i);
    CleanV_neg_loc(1,i) = find(CleanV_peaks_interval_i == CleanV_neg_pks(i,1),1,'last') + N_i(1,1) - 1; 
    % Without BEADS denoise, there could be more than one value in the minimum, we coose the last on here
end
CleanV_neg_loc = CleanV_neg_loc*T - T; % turn index into timepoint

%% Spike Parameters -------------------------------------------------------

% interspike interval 
if numel(CleanV_pos_loc) <= 1
    fprintf(1, 'Check if the firing has only one spike or less.\n');
    spikes_inter_mean = NaN;
    spikes_inter_CV = NaN;
else
    spikes_inter(1,:) = diff(CleanV_pos_loc);
    spikes_inter_mean = mean(spikes_inter,2);
    spikes_inter_CV = std(spikes_inter)/spikes_inter_mean;
end

% firing rate
if spikes_inter_CV > 1 
    fprintf(1, 'Firing is not stable.\n');
    firing_rate = NaN;
elseif isnan(spikes_inter_CV)
    firing_rate = 0;
else
    firing_rate = numel(CleanV_pos_pks)/(max(t)/1000);
end

% peak amplitude
peak_amp_mean = mean(CleanV_pos_pks,1);
peak_amp_CV = std(CleanV_pos_pks)/peak_amp_mean;

% half width
half_widths_mean = mean(half_widths,2);
half_widths_CV = std(half_widths)/half_widths_mean;

% Afterhyperpolarization
AHP_mean = mean(CleanV_neg_pks,1);
AHP_CV = std(CleanV_neg_pks)/abs(AHP_mean);

% define spike threshold 
[dddCleanV,t_diff3,dddCleanV_spike_loc,dddCleanV_1st_spike_p,CleanV_threshold] = SponFiring_thresh_fun(V_Clean,CleanV_pos_loc,t,T);
% threshold
threshold_mean = mean(CleanV_threshold,1);
threshold_CV = std(CleanV_threshold)/abs(threshold_mean);

%% Plot Waves -------------------------------------------------------------

% constrain the coordinate range
tmax = max(t);
Raw_v_max = max(V_Raw);
Raw_v_min = min(V_Raw);
vlim1 = [Raw_v_min-10 Raw_v_max+10];
tlim1 = [0 tmax];

% plot seperate components
figure('Name',sprintf('Firing of %s ', name),'WindowState','minimized');
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

%% average every single firing spike --------------------------------------
[t_avg,avg_spike,avg_spike_dt,d_avg_spike_dt,avg_spike_pks,~,avg_spike_neg_pks,avg_spike_neg_pks_3ms,avg_spike_neg_pks_20ms,avg_spike_neg_pks_150ms,~,avg_half_widths,~,avg_spike_threshold,slope_th,range] = Average_spikes_fun(V_Clean,CleanV_pos_loc,AHP_length,PrePeak_length,AHP_extend,PrePeak_extend,t,T,name);

end