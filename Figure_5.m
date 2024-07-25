%% Optoinhibition analyses

clear 
clc
% close all
load('opto_delay_summary.mat')
load('opto_sample_summary.mat')

%
close all
clearvars stim_ipsi_perf stim_contra_perf nostim_ipsi_perf nostim_contra_perf p_values_combined_hemispheres
stim_ipsi_perf.delay = mean([delay_performance(1).stim_right_average_per_mouse; delay_performance(2).stim_left_average_per_mouse]);

stim_contra_perf.delay = mean([delay_performance(1).stim_left_average_per_mouse; delay_performance(2).stim_right_average_per_mouse]);

nostim_ipsi_perf.delay = mean([delay_performance(1).no_stim_right_average_per_mouse; delay_performance(2).no_stim_left_average_per_mouse]);

nostim_contra_perf.delay = mean([delay_performance(1).no_stim_left_average_per_mouse; delay_performance(2).no_stim_right_average_per_mouse]);

[~,p_values_combined_hemispheres.delay.ipsi] = ttest(stim_ipsi_perf.delay,nostim_ipsi_perf.delay);
[~,p_values_combined_hemispheres.delay.contra] = ttest(stim_contra_perf.delay,nostim_contra_perf.delay);
[~,p_values_combined_hemispheres.delay.control_trials] = ttest(nostim_ipsi_perf.delay,nostim_contra_perf.delay);

stim_ipsi_perf.sample = mean([sample_performance(1).stim_right_average_per_mouse; sample_performance(2).stim_left_average_per_mouse]);

stim_contra_perf.sample = mean([sample_performance(1).stim_left_average_per_mouse; sample_performance(2).stim_right_average_per_mouse]);

nostim_ipsi_perf.sample = mean([sample_performance(1).no_stim_right_average_per_mouse; sample_performance(2).no_stim_left_average_per_mouse]);

nostim_contra_perf.sample = mean([sample_performance(1).no_stim_left_average_per_mouse; sample_performance(2).no_stim_right_average_per_mouse]);
[~,p_values_combined_hemispheres.sample.ipsi] = ttest(stim_ipsi_perf.sample,nostim_ipsi_perf.sample);
[~,p_values_combined_hemispheres.sample.contra] = ttest(stim_contra_perf.sample,nostim_contra_perf.sample);

figure(30);clf;
% subplot(1,2,1)
marker_size = 25;
plot(1,mean(nostim_ipsi_perf.sample),'ro','MarkerSize',marker_size-17); hold on;
std(nostim_ipsi_perf.sample)
plot(1,mean(nostim_contra_perf.sample),'bo','MarkerSize',marker_size-17);
std(nostim_contra_perf.sample)

plot(2,mean(stim_ipsi_perf.sample),'r.','MarkerSize',marker_size);
std(stim_ipsi_perf.sample)

plot(2,mean(stim_contra_perf.sample),'b.','MarkerSize',marker_size);
std(stim_contra_perf.sample)

line([1 2],[mean(nostim_ipsi_perf.sample) mean(stim_ipsi_perf.sample)],'Color','r','LineWidth',2);
line([1 2],[mean(nostim_contra_perf.sample) mean(stim_contra_perf.sample)],'Color','b','LineWidth',2)
for i = 1:size(stim_ipsi_perf.sample,2)
    lh = line([1 2],[nostim_contra_perf.sample(i) stim_contra_perf.sample(i)],'Color','b','LineWidth',2);
    lh.Color=[0,0,1,0.3];
    lh = line([1 2],[nostim_ipsi_perf.sample(i) stim_ipsi_perf.sample(i)],'Color','r','LineWidth',2);
    lh.Color=[1,0,0,0.3];
end
plot(3,mean(nostim_ipsi_perf.delay),'ro','MarkerSize',marker_size-17); hold on;
std(nostim_ipsi_perf.delay)
plot(3,mean(nostim_contra_perf.delay),'bo','MarkerSize',marker_size-17);
std(nostim_contra_perf.delay)

plot(4,mean(stim_ipsi_perf.delay),'r.','MarkerSize',marker_size);
std(stim_ipsi_perf.delay)

plot(4,mean(stim_contra_perf.delay),'b.','MarkerSize',marker_size);
std(stim_contra_perf.delay)

line([3 4],[mean(nostim_ipsi_perf.delay) mean(stim_ipsi_perf.delay)],'Color','r','LineWidth',2);
line([3 4],[mean(nostim_contra_perf.delay) mean(stim_contra_perf.delay)],'Color','b','LineWidth',2)
for i = 1:size(stim_ipsi_perf.delay,2)
    lh = line([3 4],[nostim_contra_perf.delay(i) stim_contra_perf.delay(i)],'Color','b','LineWidth',2);
    lh.Color=[0,0,1,0.3];
    lh = line([3 4],[nostim_ipsi_perf.delay(i) stim_ipsi_perf.delay(i)],'Color','r','LineWidth',2);
    lh.Color=[1,0,0,0.3];
end
line([0 5],[0.5 0.5],'Color','k','LineStyle','--');
xlim([0 5]);ylim([0.4 1]); set(gca,'Xtick',0:1:4)

ylabel('Performance (%)');xticklabels({'','Control', 'Sample','Control','Delay'});
text(2.5,0.45,['delay ipsi: ' num2str(p_values_combined_hemispheres.delay.ipsi)])
text(2.5,0.49,['delay contra: ' num2str(p_values_combined_hemispheres.delay.contra)])
text(0.7,0.45,['sample ipsi: ' num2str(p_values_combined_hemispheres.sample.ipsi)])
text(0.7,0.49,['sample contra: ' num2str(p_values_combined_hemispheres.sample.contra)])
if p_values_combined_hemispheres.delay.contra < 0.001
    text(4.1,mean(stim_contra_perf.delay),'***','FontSize',14)
elseif p_values_combined_hemispheres.delay.ipsi > 0.05
    text(4.1,mean(stim_contra_perf.delay),'N.S.','FontSize',14)

end
if p_values_combined_hemispheres.sample.ipsi > 0.05
    text(2.1,mean(stim_ipsi_perf.sample),'N.S.','FontSize',14)
end
box off;set(gca,'LineWidth',1);set(gca,'TickDir','out'); % The only other option is 'in'


% subplot(1,2,2)
% plot(1,0,'r.','MarkerSize',marker_size); hold on;
% plot(1,0,'b.','MarkerSize',marker_size);

[~,p_values_combined_hemispheres.changeinperf_between_periods.contra] = ttest2((stim_contra_perf.delay - nostim_contra_perf.delay),(stim_contra_perf.sample - stim_contra_perf.sample));
if p_values_combined_hemispheres.changeinperf_between_periods.contra < 0.001
    text(2.4,mean((stim_contra_perf.delay - nostim_contra_perf.delay))-0.06,'***','FontSize',14)
    line([2 3],[mean((stim_contra_perf.delay - nostim_contra_perf.delay)) mean((stim_contra_perf.delay - nostim_contra_perf.delay))]-0.04,'Color','k')
end
% Optostimulation analyses
clear
clc
% close all
load('opto_delay_summary_stimulation_mice.mat')
load('opto_sample_summary_stimulation_mice.mat')
%
% close all
clc;
clearvars stim_ipsi_perf stim_contra_perf nostim_ipsi_perf nostim_contra_perf p_values_combined_hemispheres
stim_ipsi_perf.delay = mean([delay_performance(1).stim_right_average_per_mouse; delay_performance(2).stim_left_average_per_mouse]);
std(stim_ipsi_perf.delay)

stim_contra_perf.delay = mean([delay_performance(1).stim_left_average_per_mouse; delay_performance(2).stim_right_average_per_mouse]);
std(stim_contra_perf.delay)

nostim_ipsi_perf.delay = mean([delay_performance(1).no_stim_right_average_per_mouse; delay_performance(2).no_stim_left_average_per_mouse]);
std(nostim_ipsi_perf.delay)

nostim_contra_perf.delay = mean([delay_performance(1).no_stim_left_average_per_mouse; delay_performance(2).no_stim_right_average_per_mouse]);
std(nostim_contra_perf.delay)

[~,p_values_combined_hemispheres.delay.ipsi] = ttest(stim_ipsi_perf.delay,nostim_ipsi_perf.delay);
[~,p_values_combined_hemispheres.delay.contra] = ttest(stim_contra_perf.delay,nostim_contra_perf.delay);
[~,p_values_combined_hemispheres.delay.control_trials] = ttest(nostim_ipsi_perf.delay,nostim_contra_perf.delay);

stim_ipsi_perf.sample = mean([sample_performance(1).stim_right_average_per_mouse; sample_performance(2).stim_left_average_per_mouse]);
% std(stim_ipsi_perf.sample)
stim_contra_perf.sample = mean([sample_performance(1).stim_left_average_per_mouse; sample_performance(2).stim_right_average_per_mouse]);
% std(stim_contra_perf.sample)

nostim_ipsi_perf.sample = mean([sample_performance(1).no_stim_right_average_per_mouse; sample_performance(2).no_stim_left_average_per_mouse]);
% std(nostim_ipsi_perf.sample)

nostim_contra_perf.sample = mean([sample_performance(1).no_stim_left_average_per_mouse; sample_performance(2).no_stim_right_average_per_mouse]);
% std(nostim_contra_perf.sample)

[~,p_values_combined_hemispheres.sample.ipsi] = ttest(stim_ipsi_perf.sample,nostim_ipsi_perf.sample);
[~,p_values_combined_hemispheres.sample.contra] = ttest(stim_contra_perf.sample,nostim_contra_perf.sample);

figure(31);clf;
% subplot(1,2,1)
marker_size = 25;
plot(1,mean(nostim_ipsi_perf.sample),'ro','MarkerSize',marker_size-17); hold on;
plot(1,mean(nostim_contra_perf.sample),'bo','MarkerSize',marker_size-17);
plot(2,mean(stim_ipsi_perf.sample),'r.','MarkerSize',marker_size);
plot(2,mean(stim_contra_perf.sample),'b.','MarkerSize',marker_size);
line([1 2],[mean(nostim_ipsi_perf.sample) mean(stim_ipsi_perf.sample)],'Color','r','LineWidth',2);
line([1 2],[mean(nostim_contra_perf.sample) mean(stim_contra_perf.sample)],'Color','b','LineWidth',2)
for i = 1:size(stim_ipsi_perf.sample,2)
    lh = line([1 2],[nostim_contra_perf.sample(i) stim_contra_perf.sample(i)],'Color','b','LineWidth',2);
    lh.Color=[0,0,1,0.3];
    lh = line([1 2],[nostim_ipsi_perf.sample(i) stim_ipsi_perf.sample(i)],'Color','r','LineWidth',2);
    lh.Color=[1,0,0,0.3];
end
plot(3,mean(nostim_ipsi_perf.delay),'ro','MarkerSize',marker_size-17); hold on;
plot(3,mean(nostim_contra_perf.delay),'bo','MarkerSize',marker_size-17);
plot(4,mean(stim_ipsi_perf.delay),'r.','MarkerSize',marker_size);
plot(4,mean(stim_contra_perf.delay),'b.','MarkerSize',marker_size);
line([3 4],[mean(nostim_ipsi_perf.delay) mean(stim_ipsi_perf.delay)],'Color','r','LineWidth',2);
line([3 4],[mean(nostim_contra_perf.delay) mean(stim_contra_perf.delay)],'Color','b','LineWidth',2)
for i = 1:size(stim_ipsi_perf.delay,2)
    lh = line([3 4],[nostim_contra_perf.delay(i) stim_contra_perf.delay(i)],'Color','b','LineWidth',2);
    lh.Color=[0,0,1,0.3];
    lh = line([3 4],[nostim_ipsi_perf.delay(i) stim_ipsi_perf.delay(i)],'Color','r','LineWidth',2);
    lh.Color=[1,0,0,0.3];
end
line([0 5],[0.5 0.5],'Color','k','LineStyle','--');
xlim([0 5]);ylim([0.4 1]); set(gca,'Xtick',0:1:4)

ylabel('Performance (%)');xticklabels({'','Control', 'Sample','Control','Delay'});
text(2.5,0.45,['delay ipsi: ' num2str(p_values_combined_hemispheres.delay.ipsi)])
text(2.5,0.49,['delay contra: ' num2str(p_values_combined_hemispheres.delay.contra)])
text(0.7,0.45,['sample ipsi: ' num2str(p_values_combined_hemispheres.sample.ipsi)])
text(0.7,0.49,['sample contra: ' num2str(p_values_combined_hemispheres.sample.contra)])
if p_values_combined_hemispheres.delay.ipsi < 0.01
    text(4.1,mean(stim_ipsi_perf.delay),'**','FontSize',14)
elseif p_values_combined_hemispheres.delay.ipsi > 0.05
    text(4.1,mean(stim_ipsi_perf.delay),'N.S.','FontSize',14)

end
if p_values_combined_hemispheres.sample.ipsi > 0.05
    text(2.1,mean(stim_ipsi_perf.sample),'N.S.','FontSize',14)
end
box off;set(gca,'LineWidth',1);set(gca,'TickDir','out'); % The only other option is 'in'


