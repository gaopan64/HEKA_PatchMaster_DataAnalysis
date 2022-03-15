%% set variables
aavg_spike_dt = mean(avg_spike_dt_Collect,2);

% Voltage derivatives
d_aavg_spike_dt = diff(aavg_spike_dt)/T;
d_aavg_spike_dt = smooth(d_aavg_spike_dt,3);
d_aavg_spike_dt(end + 1) = 0;% make the d_avg_spike_dt length is as the same as avg_spike_dt

%% Plot dv/dt vs Voltage
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

% p_c_aavg_spike_dt = plot(wt_aavg_spike_dt,wt_d_aavg_spike_dt);
% set(p_c_aavg_spike_dt,'Color',[.5 .5 .5],'LineWidth',1);

% p_c_aavg_spike_dt = plot(het_aavg_spike_dt,het_d_aavg_spike_dt);
% set(p_c_aavg_spike_dt,'Color',[0 0 .9],'LineWidth',1);

% p_k_aavg_spike_dt = plot(ko_aavg_spike_dt,ko_d_aavg_spike_dt);
% set(p_k_aavg_spike_dt,'Color',[.9 0 0],'LineWidth',1);


hold on
plot(aavg_spike_dt,d_aavg_spike_dt,'bo','markersize',3)
% plot(wt_aavg_spike_dt,wt_d_aavg_spike_dt,'ko','markersize',1)
% plot(het_aavg_spike_dt,het_d_aavg_spike_dt,'bs','markersize',1)
% plot(ko_aavg_spike_dt,ko_d_aavg_spike_dt,'r^','markersize',1)
hold off

box off
hold off
