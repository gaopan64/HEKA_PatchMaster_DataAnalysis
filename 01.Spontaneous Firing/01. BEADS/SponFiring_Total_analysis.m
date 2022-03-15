%% Specify the folder where are the files ---------------------------------
myFolder = '‪C:\Users\gaopa\Documents\Ongoing';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

%% Extract designated files -----------------------------------------------

tStart = tic;
disp('Extract designated files')
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*_SponFiring.txt'); % Change to whatever pattern you need.
% filePattern = fullfile(myFolder, '*_1.5HzFiring.txt'); % For fixed 1.5Hz
theFiles = dir(filePattern);
% If want to get a list of all files in the folder, and its subfolders, with the desired file name pattern.
% filePattern = fullfile(myFolder, '**/*.txt'); % Change to whatever pattern you need.

%% set default parameters -------------------------------------------------

% Intrinsic properties
sample_rate_input = "Please input the sampling rate(Only input interger and the unit is kHz): \n";
sample_rate_kHz = input(sample_rate_input); % set the sample rate, usually 40kHz for current clamp and 20kHz for voltage clamp
T = 1/sample_rate_kHz;

%% Collect neurons --------------------------------------------------------

% Preallocate matrix

neurons = length(theFiles);
% read one neuron to see how big a matrix need to be allocated
baseFileName_ = theFiles(1).name; 
fullFileName_ = fullfile(theFiles(1).folder, baseFileName_);
% sgl_table_ = readtable(fullFileName_,'PreserveVariableNames',true);
sgl_table_ = readtable(fullFileName_);
Height_table = height(sgl_table_);
% preallocate the matrix for collected neurons
FiringCollect = zeros(Height_table,neurons);

% preallocate the matrix for averaged spike
disp('In order to preallocate the matrix for averaged spike,')
AHP_length_input = "please input how long is the AHP that you want to average and plot (Please input the interger and the unit is ms): \n";
AHP_length = input(AHP_length_input);
PrePeak_length = 20; % 20ms is the time duration before the peak
row_spike = int64((PrePeak_length + AHP_length)/T + 1); 
avg_spike_Collect = zeros(row_spike,neurons); 

% preallocate the matrix for 'dV/dt vs Voltage' plot (averaged spike)
disp('Preallocate the matrix for dV/dt plot,')
AHP_extend_input = "please input how long you want to extend the AHP to plot dV/dt (Please input the interger and the unit is ms,suggest 300): \n";
AHP_extend = input(AHP_extend_input); 
PrePeak_extend_input = "please input how long you want to extend the length before the peak to plot dV/dt (Please input the interger and the unit is ms,suggest 300): \n";
PrePeak_extend = input(PrePeak_extend_input);
avg_spike_dt_row = int64((PrePeak_extend + AHP_extend)/T + 1);
avg_spike_dt_Collect = zeros(avg_spike_dt_row,neurons);
d_avg_spike_dt_Collect = zeros(avg_spike_dt_row,neurons);

for neurons_i = 1 : neurons
    baseFileName = theFiles(neurons_i).name;
    fullFileName = fullfile(theFiles(neurons_i).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an array with readtable and table2array
    % sgl_table = readtable(fullFileName,'PreserveVariableNames',true);
    sgl_table = readtable(fullFileName);
    FiringCollect(:,neurons_i) = table2array(sgl_table);
end

% trim preallocated matrix
fprintf(1, 'Collecting and analysing data from single neuron, please wait.....\n');


%% Analyzing single neuron ------------------------------------------------

% loop the analysis for different cells
for neurons_i = 1 : neurons
    baseFileName = theFiles(neurons_i).name;
    fprintf(1, 'Now analysing neuron of %s\n', baseFileName);
    [half_widths_mean,half_widths_CV,firing_rate,spikes_inter_mean,spikes_inter_CV,peak_amp_mean,peak_amp_CV,AHP_mean,AHP_CV,threshold_mean,threshold_CV,avg_half_widths,avg_spike_pks,avg_spike_neg_pks,avg_spike_neg_pks_3ms,avg_spike_neg_pks_20ms,avg_spike_neg_pks_150ms,avg_spike_threshold,t_avg,avg_spike,avg_spike_dt,d_avg_spike_dt,slope_th,range] = SponFiring_TotalAnal_fun(FiringCollect(:,neurons_i),AHP_length,PrePeak_length,AHP_extend,PrePeak_extend,T,baseFileName);
    tLoop = tic; % set a timer for loop
    % Analyze different features
    theFiles(neurons_i).half_widths = half_widths_mean;
    theFiles(neurons_i).hw_CV = half_widths_CV;
    theFiles(neurons_i).firing_rate = firing_rate;
    theFiles(neurons_i).spikes_interval = spikes_inter_mean;
    theFiles(neurons_i).si_CV = spikes_inter_CV;
    theFiles(neurons_i).peak_amp = peak_amp_mean;
    theFiles(neurons_i).pa_CV = peak_amp_CV;
    theFiles(neurons_i).AHP = AHP_mean;
    theFiles(neurons_i).AHP_CV = AHP_CV;
    theFiles(neurons_i).threshold = threshold_mean;
    theFiles(neurons_i).thres_CV = threshold_CV;
    theFiles(neurons_i).avg_half_widths = avg_half_widths;
    theFiles(neurons_i).avg_peak_amp = avg_spike_pks;
    theFiles(neurons_i).avg_AHP = avg_spike_neg_pks;
    theFiles(neurons_i).avg_threshold = avg_spike_threshold;
    theFiles(neurons_i).slope_th = slope_th;
    theFiles(neurons_i).range = range;
    theFiles(neurons_i).avg_fAHP = avg_spike_neg_pks_3ms;
    theFiles(neurons_i).avg_mAHP = avg_spike_neg_pks_20ms;
    theFiles(neurons_i).avg_sAHP = avg_spike_neg_pks_150ms;
    avg_spike_Collect(:,neurons_i) = avg_spike;
    avg_spike_dt_Collect(:,neurons_i) = avg_spike_dt;
    d_avg_spike_dt_Collect(:,neurons_i) = d_avg_spike_dt;
    fprintf(1, 'Analysing the neuron takes %ss\n', num2str(toc(tLoop))); 
end

fprintf(1, 'Analysis finished. Pleas check in "theFile"\n');


%% average all neurons

disp('Try average all neurons')
aavg_spike = mean(avg_spike_Collect,2); % "1" means average columns, "2" means average rows
aavg_spike_dt = mean(avg_spike_dt_Collect,2);
fprintf(1, 'Average all the data. Please wait...\n');

% find peaks in "aavg_spike"
[aavg_spike_pks,aavg_spike_loc,aavg_half_widths] = findpeaks(aavg_spike,t_avg,'MinPeakHeight',10); % Only one spike after average, no need to set many parameters
[aavg_spike_neg_pks,aavg_spike_neg_i] = min(aavg_spike);
aavg_spike_neg_loc = aavg_spike_neg_i*T - T;

% Voltage derivatives
d_aavg_spike = diff(aavg_spike)/T;
d_aavg_spike = smooth(d_aavg_spike,3);

dd_aavg_spike = diff(d_aavg_spike)/T;
dd_aavg_spike = smooth(dd_aavg_spike,3);

ddd_aavg_spike = diff(dd_aavg_spike)/T;
ddd_aavg_spike = smooth(ddd_aavg_spike,3);

% Time axes
t_aavg_diff  = t_avg(1:end-1) + T/2;
t_aavg_diff2 = t_avg(1:end-2) + T;
t_aavg_diff3 = t_avg(1:end-3) + 2*T;

% define "refer point" for the peak
ref_point = 15; % After average, the curve is smooth, don't need to strictly limit the window of onset

% Find the first peak dddV value to define the onset
% This code constrains the dddV maxumum search to time after the onset

ddd_aavg_spike_sub = ddd_aavg_spike(int64((ref_point+0.2)/T):end);

[~,first_ddd_aavg_spike_sub_i] = findpeaks(ddd_aavg_spike_sub,'MinPeakHeight',500,'MinPeakProminence',400,'Npeaks',1);
% 'MinPeakHeight'=1000,'MinPeakProminence'=200 used for averaged curve
first_ddd_aavg_spike_i = first_ddd_aavg_spike_sub_i + int64((ref_point+0.2)/T) - 1; % '-1': because when slicing ddd_aavg_spike, the sub contains the index '(ref_point+0.2)/T'

time_spike = t_aavg_diff3(1,first_ddd_aavg_spike_i); % time point of the first spike in t_aavg_diff3 axis 
i_spike = first_ddd_aavg_spike_i + 2; % switch the index of the same time point from t_aavg_diff3 to t_avg axis

% Calculate threshold
aavg_spike_threshold = aavg_spike(i_spike);

% Determine 20%-80% zone between the onset of depolarization and Voltage
% Threshold
range = time_spike - ref_point;
t20   = range * 0.2 + ref_point;
t80   = range * 0.8 + ref_point;


% Determine Slope threshold
integerIndex20 = max([1 round(t20/T)]);
integerIndex80 = max([1 round(t80/T)]);
slope_th = mean(d_aavg_spike(integerIndex20:integerIndex80));

%% Plots figures ----------------------------------------------------------

fprintf(1, 'Starting to plot data\n');


% Global definitions
fontsize = 10;

FigHandle = figure('Name', 'Total average spike','WindowState','minimized');
pause(1);
set(FigHandle, 'Position', [0, 0, 800, 800]);
  
% Plot time vs Voltage
hold on
subplot(4,1,[1 1]);
title('Total Average Spike')
set(gca,'XLim',[0 max(t_avg)])
set(gca,'YLim',[-80,max(aavg_spike)+10])
xlabel('Time (ms)','FontSize',fontsize*1.2);
ylabel('Voltage','FontSize',fontsize*1.2);

line([time_spike time_spike],[-80 20],...
    'color',[0.25 0.25 0.25],...
    'linestyle','--');

line([ref_point ref_point],[-80 20],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([t20  t20],[-80 20],...
    'color','k',...
    'linestyle','-');

line([t80  t80],[-80 20],...
    'color','k',...
    'linestyle','-');

line([time_spike - 0.08 time_spike + 0.08],[aavg_spike_threshold aavg_spike_threshold],...
    'color',[1 0 0],...
    'linestyle','-',...
    'linewidth',6);

%Place the yint value on the plot
text((t_avg(end)*0.48),aavg_spike_threshold,...
    ['V_t_h_r_e_s_h = ',num2str(aavg_spike_threshold,3),' mV'],...
    'VerticalAlignment','cap',...
    'FontSize',fontsize);
text(.1, -20,'^V^o^l^t^a^g^e/_t_i_m_e',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);
text((t20+((t80-t20)/2)),-20,...
    '\leftarrow  20%-80%\rightarrow',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize);

hold on
p_aavg_spike = plot(t_avg,aavg_spike);
set(p_aavg_spike,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t_avg,aavg_spike,'bo','markersize',3);


hold off

% Plot(tdiff, dV/dt);
hold on
subplot(4,1,[2 2]);


set(gca,'XLim',[0 max(t_avg)])
set(gca,'YLim',[min(d_aavg_spike)*1.3,max(d_aavg_spike)*1.3])
xlabel('Time (ms)','FontSize',fontsize*1.2);
ylabel('dV/dt','FontSize',fontsize*1.2);

line([time_spike time_spike],[-800 800],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([ref_point ref_point],[-800 800],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([t20  t20],[-800 800],...
    'color','k',...
    'linestyle','-');

line([t80 t80],[-800 800],...
    'color','k',...
    'linestyle','-');

line([0 100],[0 0],...
    'color',[.5 .5 .5],...
    'linestyle','-');

line([integerIndex20*T  integerIndex80*T],[slope_th slope_th],...
    'color','c',...
    'linestyle','-',...
    'linewidth',6);

hold on

p_d_aavg_spike = plot(t_aavg_diff,d_aavg_spike);
set(p_d_aavg_spike,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t_aavg_diff,d_aavg_spike,'bo','markersize',3);


%Place text on the plot
text((t_aavg_diff(end)*0.38),max(d_aavg_spike),...
    ['Slope_t_h_r_e_s_h_o_l_d = ',num2str(slope_th,3),'  ^m^V/_m_s'],...
    'VerticalAlignment','cap',...
    'FontSize',fontsize)
text(.1, max(d_aavg_spike),'dV/dt',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold off


% Plot (tdiff3, d^3/dt^3)

hold on
subplot(4,1,[3 3]);

%set(gca,'XLim',[0 max(t)])
set(gca,'XLim',[0 max(t_avg)])
set(gca,'YLim',[min(ddd_aavg_spike) * 1.3, max(ddd_aavg_spike) * 1.3])
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('d^3V/dt^3','FontSize',fontsize * 1.2);

line([(time_spike) (time_spike)],[-50000 50000],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([ref_point ref_point],[-50000 50000],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([t20  t20],[-50000 50000],...
    'color','k',...
    'linestyle','-');

line([t80 t80],[-50000 50000],...
    'color','k',...
    'linestyle','-');

line([0 100],[0 0],...
    'color',[.5 .5 .5],...
    'linestyle','-');

%Place text on the plot
text((time_spike + .05),max(ddd_aavg_spike),...
    ['\leftarrow Time of spike = ',num2str(time_spike),' ms'],...
    'VerticalAlignment','middle',...
    'FontSize',fontsize)
text((t20+((t80-t20)/2)),-1000,...
    ['Integ. period = ',num2str(range),' ms'],...
    'VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize)
text(.1, max(ddd_aavg_spike),'d^3V/dt^3',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold on

p_ddd_aavg_spike = plot(t_aavg_diff3,ddd_aavg_spike);
set(p_ddd_aavg_spike,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t_aavg_diff3,ddd_aavg_spike,'bo','markersize',3);

hold off

% plot total average spike, marks the Peaks,Valleys,Onset and Threshold in curve

hold on 
subplot(4,1,[4 4]);

plot(t_avg,aavg_spike)
hold on
plot(aavg_spike_loc, aavg_spike_pks, '^r')
plot(aavg_spike_neg_loc, aavg_spike_neg_pks, 'og')
plot(time_spike,aavg_spike_threshold,'xk')
vlim1 = [aavg_spike_neg_pks-10 aavg_spike_pks+10];
ylim(vlim1)
tliml = [0 max(t_avg)];
xlim(tliml)
legend('Data','Peak','AHP','Threshold')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);

text(.1, -20,'^V^o^l^t^a^g^e/_t_i_m_e',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold off


% Voltage derivatives
d_aavg_spike_dt = diff(aavg_spike_dt)/T;
d_aavg_spike_dt = smooth(d_aavg_spike_dt,3);
d_aavg_spike_dt(end + 1) = 0;% make the d_avg_spike_dt length is as the same as avg_spike_dt


% Plot dv/dt vs Voltage
FigHandle = figure('Name', sprintf('dV/dt versus voltage of total'),'WindowState','minimized');
pause(1);
set(FigHandle, 'Position', [0, 0, 800, 800]);
  
% Plot time vs Voltage
title(sprintf('dV/dt versus voltage'))

hold on

set(gca,'XLim',[min(aavg_spike_dt)-10,max(aavg_spike_dt)*2])
set(gca,'YLim',[min(d_aavg_spike_dt)*2,max(d_aavg_spike_dt)*2])
xlabel('Membrane Voltage (mV)','FontSize',fontsize*1.2);
ylabel('dV/dt (mV/ms)','FontSize',fontsize*1.2);

p_d_aavg_spike_dt = plot(aavg_spike_dt,d_aavg_spike_dt);
set(p_d_aavg_spike_dt,'Color',[1 0 0],'LineWidth',1);
p_d_avg_spike_dt_Collect = plot(avg_spike_dt_Collect,d_avg_spike_dt_Collect);
set(p_d_avg_spike_dt_Collect,'Color',[.5 .5 .5],'LineWidth',1);


hold on
plot(aavg_spike_dt,d_aavg_spike_dt,'bo','markersize',3)
hold off

box off
hold off

% d_aavg_spike_dt(end) = [];% remove the 0 value added above

fprintf(1, 'The whole process takes %ss\n', num2str(toc(tStart)));

writetable(struct2table(theFiles), 'theFiles.csv')
