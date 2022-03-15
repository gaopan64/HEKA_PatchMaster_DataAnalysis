function [t_avg,avg_spike,avg_spike_dt,d_avg_spike_dt,avg_spike_pks,avg_spike_loc,avg_spike_neg_pks,avg_spike_neg_pks_3ms,avg_spike_neg_pks_20ms,avg_spike_neg_pks_150ms,avg_spike_neg_loc,avg_half_widths,time_spike,avg_spike_threshold,slope_th,range] = Average_spikes_fun(V_Clean,CleanV_pos_loc,AHP_length,PrePeak_length,AHP_extend,PrePeak_extend,t,T,name)

%% split and average every single spike

% preallocate a matrix
sgl_spike_column = length(CleanV_pos_loc);
sgl_spike_row = int64((PrePeak_length + AHP_length)/T + 1); 
sgl_spike = zeros(sgl_spike_row,sgl_spike_column); 

% split every single spike values
sgl_spike_index = 1;
for i_CleanV_pos_loc = CleanV_pos_loc(1,:)
    if i_CleanV_pos_loc - PrePeak_length < 0 || i_CleanV_pos_loc + AHP_length > max(t)
        disp('Check if any spike crosses the limit.')
        continue
    end
    st_pt = int64((i_CleanV_pos_loc - PrePeak_length)/T) + 1;
    end_pt = int64((i_CleanV_pos_loc + AHP_length)/T) + 1; 

    sgl_spike(:,sgl_spike_index) = V_Clean(st_pt:end_pt,1);
    sgl_spike_index = sgl_spike_index + 1;
end

% trim preallocated matrix 
sgl_spike = sgl_spike(:,1:(sgl_spike_index-1)); 

% average every single firing spike
avg_spike = mean(sgl_spike,2); % "1" means average columns, "2" means average rows

% set time scale
size_avg_spike = size(avg_spike,1);
t_avg = 0:T:(size_avg_spike-1) * T;

% find peaks in "avg_spike"
[avg_spike_pks,avg_spike_loc,avg_half_widths,~] = findpeaks(avg_spike,t_avg,'MinPeakHeight',5); % Only one spike after average, no need to set many parameters
[avg_spike_neg_pks,avg_spike_neg_i] = min(avg_spike);
avg_spike_neg_loc = avg_spike_neg_i*T - T;

%% split and average every single spike and calculate avg_fAHP(3ms)
% preallocate a matrix
sgl_spike_column_3ms = length(CleanV_pos_loc);
sgl_spike_row_3ms = int64(3/T + 1); % AHP_length = 3,start from peak,so no PrePeak value
sgl_spike_3ms = zeros(sgl_spike_row_3ms,sgl_spike_column_3ms); 

% split every single spike values
sgl_spike_3ms_index = 1;
for i_CleanV_pos_loc = CleanV_pos_loc(1,:)
    if i_CleanV_pos_loc < 0 || i_CleanV_pos_loc + 3 > max(t)
        disp('Check if any spike crosses the limit.')
        continue
    end
    st_pt_3ms = int64(i_CleanV_pos_loc/T) + 1;
    end_pt_3ms = int64((i_CleanV_pos_loc + 3)/T) + 1; 

    sgl_spike_3ms(:,sgl_spike_3ms_index) = V_Clean(st_pt_3ms:end_pt_3ms,1);
    sgl_spike_3ms_index = sgl_spike_3ms_index + 1;
end

% trim preallocated matrix 
sgl_spike_3ms = sgl_spike_3ms(:,1:(sgl_spike_3ms_index-1)); 

% average every single firing spike
avg_spike_3ms = mean(sgl_spike_3ms,2); % "1" means average columns, "2" means average rows

% find fAHP in "avg_spike_3ms"
[avg_spike_neg_pks_3ms,~] = min(avg_spike_3ms);

%% split and average every single spike and calculate avg_mAHP(20ms)
% preallocate a matrix
sgl_spike_column_20ms = length(CleanV_pos_loc);
sgl_spike_row_20ms = int64(20/T + 1); % AHP_length = 20,start from peak,so no PrePeak value
sgl_spike_20ms = zeros(sgl_spike_row_20ms,sgl_spike_column_20ms); 

% split every single spike values
sgl_spike_20ms_index = 1;
for i_CleanV_pos_loc = CleanV_pos_loc(1,:)
    if i_CleanV_pos_loc < 0 || i_CleanV_pos_loc + 20 > max(t)
        disp('Check if any spike crosses the limit.')
        continue
    end
    st_pt_20ms = int64(i_CleanV_pos_loc/T) + 1;
    end_pt_20ms = int64((i_CleanV_pos_loc + 20)/T) + 1; 

    sgl_spike_20ms(:,sgl_spike_20ms_index) = V_Clean(st_pt_20ms:end_pt_20ms,1);
    sgl_spike_20ms_index = sgl_spike_20ms_index + 1;
end

% trim preallocated matrix 
sgl_spike_20ms = sgl_spike_20ms(:,1:(sgl_spike_20ms_index-1)); 

% average every single firing spike
avg_spike_20ms = mean(sgl_spike_20ms,2); % "1" means average columns, "2" means average rows

% find fAHP in "avg_spike_20ms"
[avg_spike_neg_pks_20ms,~] = min(avg_spike_20ms);

%% split and average every single spike and calculate avg_sAHP(150ms)
% preallocate a matrix
sgl_spike_column_150ms = length(CleanV_pos_loc);
sgl_spike_row_150ms = int64(150/T + 1); % AHP_length = 150,start from peak,so no PrePeak value
sgl_spike_150ms = zeros(sgl_spike_row_150ms,sgl_spike_column_150ms); 

% split every single spike values
sgl_spike_150ms_index = 1;
for i_CleanV_pos_loc = CleanV_pos_loc(1,:)
    if i_CleanV_pos_loc < 0 || i_CleanV_pos_loc + 150 > max(t)
        disp('Check if any spike crosses the limit.')
        continue
    end
    st_pt_150ms = int64(i_CleanV_pos_loc/T) + 1;
    end_pt_150ms = int64((i_CleanV_pos_loc + 150)/T) + 1; 

    sgl_spike_150ms(:,sgl_spike_150ms_index) = V_Clean(st_pt_150ms:end_pt_150ms,1);
    sgl_spike_150ms_index = sgl_spike_150ms_index + 1;
end

% trim preallocated matrix 
sgl_spike_150ms = sgl_spike_150ms(:,1:(sgl_spike_150ms_index-1)); 

% average every single firing spike
avg_spike_150ms = mean(sgl_spike_150ms,2); % "1" means average columns, "2" means average rows

% find fAHP in "avg_spike_150ms"
[avg_spike_neg_pks_150ms,~] = min(avg_spike_150ms);

%% Calculate single spike threshold

% Voltage derivatives
d_avg_spike = diff(avg_spike)/T;
d_avg_spike = smooth(d_avg_spike,3);

dd_avg_spike = diff(d_avg_spike)/T;
dd_avg_spike = smooth(dd_avg_spike,3);

ddd_avg_spike = diff(dd_avg_spike)/T;
ddd_avg_spike = smooth(ddd_avg_spike,3);

% Time axes
t_avg_diff  = t_avg(1:end-1) + T/2;
% t_avg_diff2 = t_avg(1:end-2) + T;
t_avg_diff3 = t_avg(1:end-3) + 2*T;

% define "refer point" for the peak
ref_point = 15; % After average, the curve is smooth, don't need to strictly limit the window of onset

% Find the first spike of ddd_avg_spike value to define the onset
% This code constrains the ddd_avg_spike maxumum search to time after onset
ddd_avg_spike_sub = ddd_avg_spike(int64((ref_point+0.2)/T):end);

[~,first_ddd_avg_spike_sub_i] = findpeaks(ddd_avg_spike_sub,'MinPeakHeight',500,'MinPeakProminence',400,'Npeaks',1); 
% 'MinPeakHeight'=1000,'MinPeakProminence'=200 used for averaged curve
first_ddd_avg_spike_i = first_ddd_avg_spike_sub_i + int64((ref_point+0.2)/T) - 1;

time_spike = t_avg_diff3(1,first_ddd_avg_spike_i); % time point of the first spike in t_aavg_diff3 axis 
i_spike = first_ddd_avg_spike_i + 2; % switch the index of the same time point from t_aavg_diff3 to t_avg axis

% Calculate threshold voltage as the voltage at the Time of spike
avg_spike_threshold = avg_spike(i_spike);

% Determine 20%-80% zone between the onset of depolarization and Voltage
% threshold
range = time_spike - ref_point;
t20   = range * 0.2 + ref_point;
t80   = range * 0.8 + ref_point;


% Determine Slope threshold
integerIndex20 = max([1 round(t20/T)]);
integerIndex80 = max([1 round(t80/T)]);
slope_th = mean(d_avg_spike(integerIndex20:integerIndex80));

%% Prepare the data for plot "dv/dt vs Voltage" 

% preallocate a matrix
sgl_spike_column_dt = length(CleanV_pos_loc);
sgl_spike_row_dt = int64((PrePeak_extend + AHP_extend)/T + 1); 
sgl_spike_dt = zeros(sgl_spike_row_dt,sgl_spike_column_dt); 

% split every single spike values
sgl_spike_index_dt = 1;
for i_CleanV_pos_loc = CleanV_pos_loc(1,:)
    if i_CleanV_pos_loc - PrePeak_extend < 0 || i_CleanV_pos_loc + AHP_extend > max(t) 
        disp('Check if any spike crosses the limit.')
        continue
    end
    st_pt = int64((i_CleanV_pos_loc - PrePeak_extend)/T) + 1;
    end_pt = int64((i_CleanV_pos_loc + AHP_extend)/T) + 1; 
    sgl_spike_dt(:,sgl_spike_index_dt) = V_Clean(st_pt:end_pt,1);
    sgl_spike_index_dt = sgl_spike_index_dt + 1;
end

% trim preallocated matrix 
sgl_spike_dt = sgl_spike_dt(:,1:(sgl_spike_index_dt-1)); 

% average every single firing spike
avg_spike_dt = mean(sgl_spike_dt,2); % "1" means average columns, "2" means average rows

% Voltage derivatives
d_avg_spike_dt = diff(avg_spike_dt)/T;
d_avg_spike_dt = smooth(d_avg_spike_dt,3);

%% Plot figures

% Global definitions
fontsize = 10;

% Plots____________________________________________________________________
FigHandle = figure('Name', sprintf('Average spike of %s ', name),'WindowState','minimized');
pause(1);
set(FigHandle, 'Position', [0, 0, 800, 800]);
  
% Plot time vs Voltage
title(sprintf('Average Spike of %s', name))

hold on
subplot(4,1,[1 1]);

set(gca,'XLim',[0 max(t_avg)])
set(gca,'YLim',[min(avg_spike)-20,max(avg_spike)+20])
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

line([time_spike - 0.08 time_spike + 0.08],[avg_spike_threshold avg_spike_threshold],...
    'color',[1 0 0],...
    'linestyle','-',...
    'linewidth',6);

% Place the yint value on the plot
text((t_avg(end)*0.48),avg_spike_threshold,...
    ['V_t_h_r_e_s_h = ',num2str(avg_spike_threshold,3),' mV'],...
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
p_avg_spike = plot(t_avg,avg_spike);
set(p_avg_spike,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t_avg,avg_spike,'bo','markersize',3);


hold off

% Plot(tdiff, dV/dt);
hold on
subplot(4,1,[2 2]);


set(gca,'XLim',[0 max(t_avg)])
set(gca,'YLim',[min(d_avg_spike)*1.3,max(d_avg_spike)*1.3])
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
    'linewidth',10);

hold on

p_d_avg_spike = plot(t_avg_diff,d_avg_spike);
set(p_d_avg_spike,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t_avg_diff,d_avg_spike,'bo','markersize',3);


% Place text on the plot
text((t_avg_diff(end)*0.38),max(d_avg_spike),...
    ['Slope_t_h_r_e_s_h_o_l_d = ',num2str(slope_th,3),'  ^m^V/_m_s'],...
    'VerticalAlignment','cap',...
    'FontSize',fontsize)
text(.1, max(d_avg_spike),'dV/dt',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold off


% Plot (tdiff3, dddV/dt3)

hold on
subplot(4,1,[3 3]);

% set(gca,'XLim',[0 max(t)])
set(gca,'XLim',[0 max(t_avg)])
set(gca,'YLim',[min(ddd_avg_spike) * 1.3, max(ddd_avg_spike) * 1.3])
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

% Place text on the plot
text((time_spike + .05),max(ddd_avg_spike),...
    ['\leftarrow Time of spike = ',num2str(time_spike),' ms'],...
    'VerticalAlignment','middle',...
    'FontSize',fontsize)
text((t20+((t80-t20)/2)),-1000,...
    ['Integ. period = ',num2str(range),' ms'],...
    'VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize)
text(.1, max(ddd_avg_spike),'d^3V/dt^3',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold on

p_ddd_avg_spike = plot(t_avg_diff3,ddd_avg_spike);
set(p_ddd_avg_spike,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t_avg_diff3,ddd_avg_spike,'bo','markersize',3);

hold off

% plot average spike, marks the Peaks,Valleys,and Threshold in curve
hold on
subplot(4,1,[4 4]);

plot(t_avg,avg_spike)
hold on
plot(avg_spike_loc, avg_spike_pks, '^r')
plot(avg_spike_neg_loc, avg_spike_neg_pks, 'og')
plot(time_spike,avg_spike_threshold,'xk')
hold off
V_Clean_min = min(V_Clean);
V_Clean_max = max(V_Clean);
vlim = [V_Clean_min-10 V_Clean_max+10];
ylim(vlim)
tLim = [0 max(t_avg)];
xlim(tLim)
legend('Data','Peak','AHP','Threshold')
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('Voltage (mV)','FontSize',fontsize * 1.2);
box off

text(.1, -20,'^V^o^l^t^a^g^e/_t_i_m_e',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

%% Plots dv/dt vs Voltage

FigHandle = figure('Name', sprintf('dV/dt versus voltage of %s', name),'WindowState','minimized');
pause(1);
set(FigHandle, 'Position', [0, 0, 800, 800]);
  
% Plot time vs Voltage
title(sprintf('dV/dt versus voltage of %s', name),'Interpreter', 'none') % set the interpreter to 'none' to keep the underscore format

hold on

set(gca,'XLim',[min(avg_spike_dt)-10,max(avg_spike_dt)+10])
set(gca,'YLim',[min(d_avg_spike_dt)*1.5,max(d_avg_spike_dt)*1.5])
xlabel('Membrane Voltage (mV)','FontSize',fontsize*1.2);
ylabel('dV/dt (mV/ms)','FontSize',fontsize*1.2);

d_avg_spike_dt(end + 1) = 0;% make the d_avg_spike length is as the same as avg_spike

p_d_avg_spike_dt = plot(avg_spike_dt,d_avg_spike_dt);
set(p_d_avg_spike_dt,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(avg_spike_dt,d_avg_spike_dt,'bo','markersize',3)
hold off

box off
hold off

end