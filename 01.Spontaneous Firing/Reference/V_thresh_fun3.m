function [ slope_th,range, V_th] = V_thresh_fun3( sample_rate_kHz, V, onset )
%V_THRESH_FUN Summary of this function goes here
%   Detailed explanation goes here

% Global definitions
fontsize = 10;


% Intrinsic properties
T = 1/sample_rate_kHz;
size_V = size(V,1);


% Time axes
t = 0:T:(size_V-1) * T;
tdiff  = t(1:end-1) + T/2;
tdiff2 = t(1:end-2) + T;


% Voltage derivatives
dV  = diff(V)  * sample_rate_kHz;
dV = smooth(dV,3);

ddV = diff(dV) * sample_rate_kHz;
ddV = smooth(ddV,3);


% Find the maximum ddV value to define the Time of spike
% This code constrains the ddV maxumum search to time after the onset

ddV_sub = ddV((onset+.2)/T:end);

[max_ddV_sub,max_ddV_sub_i] = max(ddV_sub);
max_ddV = max_ddV_sub;
max_ddV_i = max_ddV_sub_i + (onset+.2)/T;


time_spike = max_ddV_i * T - T;
i_spike = max_ddV_i;

% Calculate Voltage threshold as the voltage at the Time of spike
V_th = V(i_spike);


% Determine 20%-80% zone between the onset of depolarization and Voltage
% threshold
range = time_spike - onset;
t20   = range * .2 + onset;
t80   = range * .8 + onset;


% Determine Slope threshold
integerIndex20 = max([1 round(t20 * sample_rate_kHz)]);
integerIndex80 = max([1 round(t80 * sample_rate_kHz)]);
slope_th = mean(dV(integerIndex20:integerIndex80));






% Plots____________________________________________________________________
FigHandle = figure;
  set(FigHandle, 'Position', [0, 0, 800, 800]);
  
% Plot time vs Voltage
hold on
subplot(3,1,[1 1]);

set(gca,'XLim',[0 time_spike*1.3])
set(gca,'YLim',[-75,max(V)+5])
xlabel('Time (ms)','FontSize',fontsize*1.2);
ylabel('Voltage','FontSize',fontsize*1.2);

line([time_spike time_spike],[-80 20],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([onset onset],[-80 20],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([t20  t20],[-80 20],...
    'color','k',...
    'linestyle','-');

line([t80  t80],[-80 20],...
    'color','k',...
    'linestyle','-');

line([time_spike - .08 time_spike + .08],[V_th V_th],...
    'color',[1 0 0],...
    'linestyle','-',...
    'linewidth',3);

%Place the yint value on the plot
text((t(end)*0.7),10,...
    ['V_t_h_r_e_s_h_o_l_d = ',num2str(V_th,3),' mV'],...
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
p_V = plot(t,V);
set(p_V,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(t,V,'bo','markersize',3);


hold off






% Plot(tdiff, dV/dt);
hold on
subplot(3,1,[2 2]);


set(gca,'XLim',[0 time_spike*1.3])
set(gca,'YLim',[min(dV)*1.3,max(dV)*1.3])
xlabel('Time (ms)','FontSize',fontsize*1.2);
ylabel('dV/dt','FontSize',fontsize*1.2);

line([time_spike time_spike],[-800 800],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([onset onset],[-800 800],...
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

p_dV = plot(tdiff,dV);
set(p_dV,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(tdiff,dV,'bo','markersize',3);


%Place text on the plot
text((tdiff(end)*.7),max(dV),...
    ['Slope_t_h_r_e_s_h_o_l_d = ',num2str(slope_th,3),'  ^m^V/_m_s'],...
    'VerticalAlignment','cap',...
    'FontSize',fontsize)
text(.1, max(dV),'^d^V/_d_t',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold off





% Plot (tdiff2, ddV/dt2)

hold on
subplot(3,1,[3 3]);

%set(gca,'XLim',[0 max(t)])
set(gca,'XLim',[0 time_spike*1.3])
set(gca,'YLim',[min(ddV) * 1.3, max(ddV) * 1.3])
xlabel('Time (ms)','FontSize',fontsize * 1.2);
ylabel('ddV/dt^2','FontSize',fontsize * 1.2);

line([(time_spike) (time_spike)],[-5000 5000],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([onset onset],[-5000 5000],...
    'color',[.25 .25 .25],...
    'linestyle','--');

line([t20  t20],[-5000 5000],...
    'color','k',...
    'linestyle','-');

line([t80 t80],[-5000 5000],...
    'color','k',...
    'linestyle','-');

line([0 100],[0 0],...
    'color',[.5 .5 .5],...
    'linestyle','-');

%Place text on the plot
text((time_spike + .05),max(ddV),...
    ['\leftarrow Time of spike = ',num2str(time_spike),' ms'],...
    'VerticalAlignment','middle',...
    'FontSize',fontsize)
text((t20+((t80-t20)/2)),-1000,...
    ['Integ. period = ',num2str(range),' ms'],...
    'VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize)
text(.1, max(ddV),'^d^d^V/_d_t2',...
    'VerticalAlignment','cap',...
    'FontSize',fontsize*1.5,...
    'BackgroundColor',[.7 .9 .7],...
	'Margin',3);

hold on

p_ddV = plot(tdiff2,ddV);
set(p_ddV,'Color',[.5 .5 .5],'LineWidth',1);

hold on
plot(tdiff2,ddV,'bo','markersize',3);

hold off



end
