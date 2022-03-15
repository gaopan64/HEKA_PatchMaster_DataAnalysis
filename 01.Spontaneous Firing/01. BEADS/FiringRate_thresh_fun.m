function [dddCleanV,t_diff3,dddCleanV_spike_loc,dddCleanV_1st_spike_p,CleanV_threshold] = FiringRate_thresh_fun(V_Clean,CleanV_pos_loc,t,T)

% V_THRESH_FUN Summary of this function goes here
% Detailed explanation goes here
% Voltage derivatives
dCleanV = diff(V_Clean)/T;
dCleanV = smooth(dCleanV,3);
ddCleanV = diff(dCleanV)/T;
ddCleanV = smooth(ddCleanV,3);
dddCleanV = diff(ddCleanV)/T;
dddCleanV = smooth(dddCleanV,3);

% Time axes
% t_diff  = t(1:end-1) + T/2;
% t_diff2 = t(1:end-2) + T;
t_diff3 = t(1:end-3) + 2*T;

% define "ref_point" for each peak
ref_point = CleanV_pos_loc - 1.6; % define ref_point at 1.6ms before each peak

% preallocate a matrix
dddCleanV_sub_column = length(CleanV_pos_loc); % number of spikes
dddCleanV_sub_row = int64(1/T + 1); % number of values, 1ms/T+1
dddCleanV_sub = zeros(dddCleanV_sub_row,dddCleanV_sub_column);% For 40kHz sampling 

% limit the dddCleanV maxumum search to the time window between the onset(ref_point: usually 1ms before peak) and the peak
ref_index = 1;
for i_ref = ref_point(1,:)
    if i_ref < 0
        ref_point(1) = [];
        disp('ref_point value before the start point, please check if any spike goes too early')
        continue
    elseif i_ref+1 > max(t)
        ref_point(end) = [];
        disp('ref_point value exceeds the end point, please check if any spike goes too late')
        continue
    end
    st_pt = int64(i_ref/T);
    end_pt = int64((i_ref+1)/T);
    dddCleanV_sub(:,ref_index) = dddCleanV(st_pt:end_pt,1); % dddCleanV_sub stores dddCleanV values for each peak in different columns
    ref_index = ref_index + 1;
end

% trim preallocated matrix
dddCleanV_sub = dddCleanV_sub(:,1:(ref_index-1)); % matrix need to be trimed because there may be spikes cross the limit

% define the peak point 
[dddCleanV_spike_p,dddCleanV_spike_i] = arrayfun(@(col)findpeaks(dddCleanV_sub(:,col),'MinPeakHeight',500,'MinPeakProminence',400,'Npeaks',1),1:size(dddCleanV_sub,2),'UniformOutput',false); 
% 'MinPeakProminence' = 400 and 'MinPeakHeight'=500 used for only the coarse curve of firingrate
% For SponFiring and 1.5Hz Firing,apply 'MinPeakHeight'=800

% dddCleanV_spike_p, dddCleanV_spike_i would be cell arrays. Index will correspond to column 1 and so on
dddCleanV_1st_spike_p = cellfun(@(v)v(1),dddCleanV_spike_p); % select first peak as threshold
dddCleanV_1st_spike_sub_i = cellfun(@(v)v(1),dddCleanV_spike_i); % location on "ddCleanV_sub" index 
dddCleanV_1st_spike_i = dddCleanV_1st_spike_sub_i + round(ref_point/T) - 1; % switch "ddCleanV_sub" index into "t" or "V_Clean" index; '-1': because when slicing dddCleanV_sub, the sub contains the index 'ref_point/T'
dddCleanV_spike_loc = t_diff3(1,dddCleanV_1st_spike_i); % peak location on "t_diff3" value(ms)
first_spike = dddCleanV_1st_spike_i + 2; % switch the index of the same time point from t_aavg_diff3 to t_avg axis

% preallocate a matrix
CleanV_thres_row = length(CleanV_pos_loc);
CleanV_threshold = zeros(CleanV_thres_row,1);

% Calculate threshold voltage as the voltage at the Time of first spike
ref_index = 1;
for i_spike = first_spike(1,:)
    CleanV_threshold(ref_index,1) = V_Clean(i_spike,1);
    ref_index = ref_index +1;
end

% trim preallocated matrix 
CleanV_threshold = CleanV_threshold(1:(ref_index-1),1);

end