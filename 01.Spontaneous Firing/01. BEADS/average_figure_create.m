plot(t_avg,aavg_spike,'color','blue')
hold on
plot(aavg_spike_loc, aavg_spike_pks, '^c')
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
