% Global definitions
fontsize = 10;

% Plots____________________________________________________________________
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

fprintf(1, 'Plotting finished.\n');