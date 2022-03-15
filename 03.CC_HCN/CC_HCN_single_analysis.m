V_RAW = table2array(DA2111241624CCHCN);

name = 'SingleNeuron';

tStart = tic;

V_Raw = V_RAW*1000;
V_Raw = detrend(V_Raw);
V_Raw = V_Raw-70;

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
    fprintf(1, 'Start to denoise the sweep number: %s, please wait...\n', num2str(i));
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

% plot "CleanV", marks the Peaks,Valleys,and Threshold in curve
subplot(1, 1, [1 1])

plot(t, V_Clean)
hold on
title('Clean CC_HCN','FontSize',fontsize * 1.2,'Interpreter', 'none')
hold off
xlim(tlim1)
ylim(vlim1)
legend('0pA','-10pA','-20pA','-30pA','-40pA','-50pA','-60pA','-70pA','-80pA','-90pA','-100pA')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off
axis off

fprintf(1, 'The drawing takes %ss\n', num2str(toc(tDrawing)));
% -------------------------------------------------------------------------

% bridge balance (To be developed)

% calculate the sag ratio
% current clamped to make the voltage at -70mV from the begining

% measure the steady state decrease in the voltage
steady_window_t = (4500:T:5000); % because the CC_HCN stim stop at 5s, we set the average value during the time before 5s as the voltage of steady state
steady_window_i = int64(steady_window_t/T + 1); % Index = t/T + 1
steady_V_Clean = V_Clean(steady_window_i,:); 
steady_V = mean(steady_V_Clean);
steady_decrease = steady_V - (-70); % here we set the base line at -70mV

% measure the largest decrease in voltage following a hyperpolarizing current step
valley_window_t = (1000:T:2000);
valley_window_i = int64(valley_window_t/T + 1);
valley_V_Clean = V_Clean(valley_window_i,:);
valley_V = min(valley_V_Clean);
Valley_decrease = valley_V - (-70); % here we set the base line at -70mV

% calculate sag ratio
sag_ratio = steady_decrease./Valley_decrease;
% replace the fist value of sag ratio to 1
sag_ratio(:,1) = 1;

% -------------------------------------------------------------------------
% calculate rebound delay
% rebound_window_t = (5000:T:10000);
% rebound_window_i = round(rebound_window_t/T + 1);
% rebound_V_Clean = V_Clean(rebound_window_i,:); 
% [rebound_spike_first,rebound_spike_first_i] = findpeaks(rebound_V_Clean,'MinPeakHeight',0,'Npeaks',1);