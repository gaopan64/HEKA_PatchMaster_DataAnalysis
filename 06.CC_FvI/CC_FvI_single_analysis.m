% warning! open the parallel proccess
parpool(12) 

V_RAW = table2array(DA2101040604CCFvIV);

name = 'DA210107023502';

tStart = tic;

V_Raw = V_RAW*1000;
% V_Raw = detrend(V_Raw);
% V_Raw = V_Raw - 60;

% Global definitions
fontsize = 8;

% Intrinsic properties
sample_rate_kHz = 40;
T = 1/sample_rate_kHz;
size_V = size(V_Raw,1);

% Time axes
t = 0:T:(size_V-1) * T;

% Index = t/T + 1
% t = (Index - 1)*T = Index*T - T

N = length(V_Raw);

% Filter parameters
fc = 0.006;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)

% Positivity bias (peaks are positive)
r = 6;          % r : asymmetry parameter

% Regularization parameters
amp = 0.08;
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

% run BEADS function
tBEADS = tic;

% preallocate matrix
x1_column = size(V_Raw,2);
x1_row = size(V_Raw,1);
x1 = zeros(x1_row,x1_column); % x stands for peak data
f1 = zeros(x1_row,x1_column); % f stands for baseline, same size with xl_sweeps

parfor i = 1:size(V_Raw,2)
    [xi, fi, cost] = beads(V_Raw(:,i), d, fc, r, lam0, lam1, lam2);
    x1(:,i) = xi(:,1);
    f1(:,i) = fi(:,1);
end
fprintf(1, 'The neuron denoising and baseline correction takes %ss\n', num2str(toc(tBEADS)));

% define denoised curve as "V_Clean"
V_Clean = x1 + f1;

tDrawing = tic;

% constrain the coordinate range
tmax = max(t);
Raw_v_max = max(V_Raw,[],'all'); % Find the minimum value of all elements of A. This syntax is applicable to MATLAB® R2018b and later.
Raw_v_min = min(V_Raw,[],'all'); % Find the minimum value of all elements of A. This syntax is applicable to MATLAB® R2018b and later.
vlim1 = [Raw_v_min-10 Raw_v_max+10];
tlim1 = [0 tmax];

% plot seperate components
f = figure('WindowState','minimized');
pause(1);
clf

subplot(5, 1, [1 1])
plot(t,V_Raw)
title('CC_FvI','FontSize',fontsize * 1.2,'Interpreter', 'none')
xlim(tlim1)
ylim(vlim1)
legend('0pA','+20pA','+40pA','+60pA','+80pA','+100pA','+120pA','+140pA','+160pA','+180pA')
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

subplot(5, 1, [2 2])
plot(V_Raw,'color', [1 1 1]*0.7)
line(1:N, f1, 'LineWidth', 1)
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
plot(t,V_Raw - x1 - f1)
title('Noise','FontSize',fontsize * 1.2)
xlim(tlim1)
set(gca,'ytick', vlim1)
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

% plot "CleanV", marks the Peaks,Valleys,and Threshold in curve
subplot(5, 1, [5 5])

plot(t, V_Clean)
hold on
title('Clean CC_FvI','FontSize',fontsize * 1.2,'Interpreter', 'none')
hold off
xlim(tlim1)
ylim(vlim1)
legend('0pA','+20pA','+40pA','+60pA','+80pA','+100pA','+120pA','+140pA','+160pA','+180pA')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

fprintf(1, 'The drawing takes %ss\n', num2str(toc(tDrawing)));
% -------------------------------------------------------------------------
% measure the AHP following the spikes
valley_window_t = (3000:T:4000);
valley_window_i = round(valley_window_t/T + 1);
valley_V_Clean = V_Clean(valley_window_i,:);
valley_V = min(valley_V_Clean);
% -------------------------------------------------------------------------
% Find peaks for each sweep
V_Clean_Sub_t = (1000:T:3000);
V_Clean_Sub_i = round(V_Clean_Sub_t/T + 1);
V_Clean_Sub = V_Clean(V_Clean_Sub_i,:);

parfor i = 1:size(V_Clean_Sub,2)
    [CleanV_pos_pks_i,CleanV_pos_loc_i] = findpeaks(V_Clean_Sub(:,i),V_Clean_Sub_t,'MinPeakHeight',0); % "loc" means t value, "pks" means voltage value
    if isempty(CleanV_pos_pks_i)
        CleanV_pos_pks(:,i) = 0;
        CleanV_pos_loc(i,:) = 0;
        continue
    else
        CleanV_pos_pks(:,i) = CleanV_pos_pks_i(:,1);
        CleanV_pos_loc(i,:) = CleanV_pos_loc_i(1,:);
    end
end
% Number of peaks
parfor i = 1:size(CleanV_pos_loc,1)
    pks_n(i,1) = nnz(CleanV_pos_loc(i,:)); % numbers of nonzeros
end
% -------------------------------------------------------------------------
[VSC_row,VCS_col] = size(V_Clean);
V_Clean_S = zeros(VSC_row,VCS_col);
parfor i = 1:size(V_Clean,2)
    V_Clean_S(:,i) = smooth(V_Clean(:,i));
end
% -------------------------------------------------------------------------
delete(gcp('nocreate'));