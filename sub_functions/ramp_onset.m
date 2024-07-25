function [time_of_onset_d1r,time_of_onset_non] = ramp_onset(x,y,x_time,y_time,start_time,end_time)
number_of_change_points = 1;

clear ipt ipt_1_d1r time_of_onset_d1r
for i = 1:size(x,1)
    
    if isempty(findchangepts(x(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==1
        ipt_1_d1r(i) = NaN;
        time_of_onset_d1r(i)= NaN;
        
    elseif isempty(findchangepts(x(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==0
        [ipt{i}] = findchangepts(x(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear');
        time_of_onset_d1r(i)= x_time(i,ipt{1,i}(1)+(start_time-1));
    end
end
%
clear ipt ipt_1_non time_of_onset_non
for i = 1:size(y,1)
    
    
    if isempty(findchangepts(y(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==1
        ipt_1_non(i) = NaN;
        time_of_onset_non(i)= NaN;
        
    elseif isempty(findchangepts(y(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==0
        [ipt{i}] = findchangepts(y(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear');
        ipt_1_non(i) = ipt{1,i}(1);
        time_of_onset_non(i)= y_time(i,ipt{1,i}(1)+(start_time-1));
    end
    
end
time_of_onset_non=(time_of_onset_non(~isnan(time_of_onset_non)));
%
time_of_onset_d1r=(time_of_onset_d1r(~isnan(time_of_onset_d1r)));

end