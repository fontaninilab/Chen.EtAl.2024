 %% Figure 1: Behavioral data analyses
clear
 
load('all_behavioral_sessions.mat')

%%
close all; clc;
suc_trials = cell(1,size(trial_1_valve,2));
for i=1:length(trial_1_valve)
    if isfield(trial_1_valve{1,i},'TasteID')
        for j = 1:size(trial_1_valve{1,i},2)
            if strcmp(trial_1_valve{1,i}(j).TasteID,'Sucrose') && trial_1_valve{1,i}(j).correct_choice == 1
                suc_perf{j,i} = trial_1_valve{1,i}(j).correct_choice;
                suc_dur{j,i} = trial_1_valve{1,i}(j).centSp(6:end);
            elseif strcmp(trial_1_valve{1,i}(j).TasteID,'Sucrose') && trial_1_valve{1,i}(j).correct_choice == 0
                suc_perf{j,i} = nan;
                suc_dur_error{j,i} = trial_1_valve{1,i}(j).centSp(6:end);
                left_lick_error{j,i} = trial_1_valve{1,i}(j).LeftSp;
            end
            
        end
    else
        for j = 1:size(trial_1_valve{1,i},2)
            if strcmp(trial_1_valve{1,i}(j).tasteID,'Suc') && trial_1_valve{1,i}(j).choice_correct == 1
                suc_perf{j,i} = trial_1_valve{1,i}(j).choice_correct;
                suc_dur{j,i} = trial_1_valve{1,i}(j).all_TDT_central_licks(7:end);
            elseif strcmp(trial_1_valve{1,i}(j).tasteID,'Suc') && trial_1_valve{1,i}(j).choice_correct == 0
                suc_perf{j,i} = nan;
                suc_dur_error{j,i} = trial_1_valve{1,i}(j).all_TDT_central_licks(7:end);
            end
        end
    end
end
%
for i=1:length(trial_1_valve)
    if isfield(trial_1_valve{1,i},'TasteID')
        for j = 1:size(trial_1_valve{1,i},2)
            if strcmp(trial_1_valve{1,i}(j).TasteID,'NaCl') && trial_1_valve{1,i}(j).correct_choice == 1
                salt_perf{j,i} = trial_1_valve{1,i}(j).correct_choice;
                salt_dur{j,i} = trial_1_valve{1,i}(j).centSp(6:end);
            elseif strcmp(trial_1_valve{1,i}(j).TasteID,'NaCl') && trial_1_valve{1,i}(j).correct_choice == 0
                salt_perf{j,i} = NaN;
                salt_dur_error{j,i} = trial_1_valve{1,i}(j).centSp(6:end);
            end
            
        end
    else
        for j = 1:size(trial_1_valve{1,i},2)
            if strcmp(trial_1_valve{1,i}(j).tasteID,'NaCl') && trial_1_valve{1,i}(j).choice_correct == 1
                salt_perf{j,i} = trial_1_valve{1,i}(j).choice_correct;
                salt_dur{j,i} = trial_1_valve{1,i}(j).all_TDT_central_licks(7:end);
            elseif strcmp(trial_1_valve{1,i}(j).tasteID,'NaCl') && trial_1_valve{1,i}(j).choice_correct == 0
                salt_perf{j,i} = NaN;
                salt_dur_error{j,i} = trial_1_valve{1,i}(j).all_TDT_central_licks(7:end);
            end
        end
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear left_lick_correct left_lick_error right_lick_correct right_lick_error right_dura2 left_dura2 left_correct_ra left_ra2
for i=1:length(trial_1_valve)
    if isfield(trial_1_valve{1,i},'TasteID')
        for j = 1:size(trial_1_valve{1,i},2)
            if trial_1_valve{1,i}(j).L_R_trial == 1 && trial_1_valve{1,i}(j).correct_choice == 1
                left_lick_correct{j,i} = trial_1_valve{1,i}(j).LeftSp;
                left_correct_ra{j,i} = trial_1_valve{1,i}(j).LeftSp(1) -  trial_1_valve{1,i}(j).Up;
                
            elseif  trial_1_valve{1,i}(j).L_R_trial == 1 && trial_1_valve{1,i}(j).correct_choice == 0
                left_lick_error{j,i} = trial_1_valve{1,i}(j).LeftSp;
                
            elseif  trial_1_valve{1,i}(j).L_R_trial == 2 && trial_1_valve{1,i}(j).correct_choice == 1
                right_lick_correct{j,i} = trial_1_valve{1,i}(j).RightSp;
                right_correct_ra{j,i} = trial_1_valve{1,i}(j).RightSp(1) - trial_1_valve{1,i}(j).Up;
                
            elseif  trial_1_valve{1,i}(j).L_R_trial == 2 && trial_1_valve{1,i}(j).correct_choice == 0
                right_lick_error{j,i} = trial_1_valve{1,i}(j).RightSp;
            end
            
        end
    else
        for j = 1:size(trial_1_valve{1,i},2)
            if trial_1_valve{1,i}(j).choice_correct == 1 && trial_1_valve{1,i}(j).choice == 1
                left_lick_correct{j,i} = trial_1_valve{1,i}(j).all_TDT_left_licks;
                left_correct_ra{j,i} = (trial_1_valve{1,i}(j).left_lick(1) -  trial_1_valve{1,i}(j).lateral_motor(1))*0.0051;
                
            elseif trial_1_valve{1,i}(j).choice_correct == 0 && trial_1_valve{1,i}(j).choice == 1
                left_lick_error{j,i} = trial_1_valve{1,i}(j).all_TDT_left_licks;
                
            elseif trial_1_valve{1,i}(j).choice_correct == 1 && trial_1_valve{1,i}(j).choice == 2
                right_lick_correct{j,i} = trial_1_valve{1,i}(j).all_TDT_right_licks;
                right_correct_ra{j,i} = (trial_1_valve{1,i}(j).right_lick(1) -  trial_1_valve{1,i}(j).lateral_motor(1))*0.0051;
                
            elseif trial_1_valve{1,i}(j).choice_correct == 0 && trial_1_valve{1,i}(j).choice == 2
                right_lick_error{j,i} = trial_1_valve{1,i}(j).all_TDT_right_licks;
            end
        end
    end
end
%
i = [];
j = [];
for i=1:length(trial_1_valve)
    a = cell2mat(suc_perf(:,i));
    a(:,2) = ~isnan(a);
    suc_corr(i) = (sum(a(:,2)))/length(a(:,2));
    b = cell2mat(salt_perf(:,i));
    b(:,2) = ~isnan(b);
    salt_corr(i) = (sum(b(:,2)))/length(b(:,2));
    
    clear c
    c = suc_dur(:,i);
    c = c(~cellfun('isempty',c));
    
    for k = 1:size(c,1)
        
        if length(c{k,1}(1:end)) == 1
            diff1(k) = length(c{k,1}(1:end))/1;
        else
            diff1(k) = length(c{k,1}(1:end))/(max(c{k,1}(1:end))-min(c{k,1}(1:end)));
        end
    end
    suc_dura2(1,i) = mean(diff1,'omitnan');
    
    clear d
    d = salt_dur(:,i);
    d = d(~cellfun('isempty',d));
    
    for q = 1:size(d,1)
        
        if length(d{q,1}(1:end)) == 1
            diff2(q) = length(d{q,1}(1:end))/1;
        else
            diff2(q) = length(d{q,1}(1:end))/(max(d{q,1}(1:end))-min(d{q,1}(1:end)));
        end
    end
    
    salt_dura2(1,i) = mean(diff2,'omitnan');
    
    clear e
    e = left_lick_correct(:,i);
    e = e(~cellfun('isempty',e));
    clear q diff2
    for q = 1:size(e,1)
        
        if length(e{q,1}(1:end)) == 1
%        
            diff2(q) = nan;

        else
            diff2(q) = length(e{q,1}(1:end))/(max(e{q,1}(1:end))-min(e{q,1}(1:end)));
%        
        end
    end
    if ~isempty(e)
        left_dura2(1,i) = mean(diff2,'omitnan');
    else
            left_dura2(1,i) = nan;
    end

    clear f diff2
    f = right_lick_correct(:,i);
    f = f(~cellfun('isempty',f));
    clear q
    for q = 1:size(f,1)
        
        if length(f{q,1}(1:end)) == 1
            diff2(q) = nan;
        else
            diff2(q) = length(f{q,1}(1:end))/(max(f{q,1}(1:end))-min(f{q,1}(1:end)));
%             
        end
    end
    if ~isempty(f)
        right_dura2(1,i) = mean(diff2,'omitnan');
    else
            right_dura2(1,i) = nan;
    end
    
    clear g diff2
    g = left_correct_ra(:,i);
    g = g(~cellfun('isempty',g));

    left_react2(1,i) = mean(cell2mat(g),'omitnan');
    
     clear h diff2
    h = right_correct_ra(:,i);
    h = h(~cellfun('isempty',h));

    right_react2(1,i) = mean(cell2mat(h),'omitnan');
end
%
disp('Plotting summary figures')
clear suc_corr1 salt_corr1


suc_corr1 = [mean(suc_corr(1:3)) mean(suc_corr(4:6)) mean(suc_corr(7:9)) mean(suc_corr(10:12))...
    mean(suc_corr(13:15)) mean(suc_corr(16:19)) mean(suc_corr(20:23)) mean(suc_corr(24:26))...
    mean(suc_corr(27:29)) mean(suc_corr(30:35)) mean(suc_corr(36:39)) mean(suc_corr(40:43))...
    mean(suc_corr(44:48)) mean(suc_corr(49:53)) mean(suc_corr(54:57))];
salt_corr1 = [mean(salt_corr(1:3)) mean(salt_corr(4:6)) mean(salt_corr(7:9)) mean(salt_corr(10:12))...
    mean(salt_corr(13:15)) mean(salt_corr(16:19)) mean(salt_corr(20:23)) mean(salt_corr(24:26))...
    mean(salt_corr(27:29)) mean(salt_corr(30:35)) mean(salt_corr(36:39)) mean(salt_corr(40:43))...
    mean(salt_corr(44:48)) mean(salt_corr(49:53)) mean(salt_corr(54:57))];

figure(1);clf
 subplot(1,5,1) %% Performance
std_suc = std(suc_corr1); std_salt = std(salt_corr1);
avg_suc = mean(suc_corr1);
avg_salt = mean(salt_corr1);
p = bar(1:2,[avg_suc,avg_salt]);hold on
p.FaceColor = 'flat';
p.CData(1,:) = [0 0 1];
p.CData(2,:) = [1 0 0];
er = errorbar([1,2], [avg_suc,avg_salt],zeros(1,2),[std_suc/sqrt(size(suc_corr1,2)),std_salt/sqrt(size(salt_corr1,2))],'k','LineWidth',1.2);hold on
ylabel('Performance');
scatter(ones(size(suc_corr1)),suc_corr1,200,'k.'); hold on;
scatter(ones(size(salt_corr1))+1,salt_corr1,200,'k.'); hold on;
for i = 1:size(salt_corr1,2)
    line([1,2],[suc_corr1(1,i) salt_corr1(1,i)],'Color','k')
end

xlim([0 3]);xticklabels({'Suc', 'NaCl'});
% xtickangle(45);
ylim([0.5 1]); yticks([-0.5:0.1:1]);
line([0 4],[0.7 0.7],'Color','k','LineStyle','--');
[~,pvalues.performance,~,stats_table.performance] = ttest(suc_corr1,salt_corr1);
set(er,'linestyle','none');
box off
% axis tight
disp({'suc perf: ',num2str(avg_suc)}); disp({'nacl perf: ',num2str(avg_salt)})
disp({'suc std: ',num2str(std_suc/sqrt(size(suc_corr1,2)))}); disp({'nacl std: ',num2str(std_salt/sqrt(size(salt_corr1,2)))})

%
clear suc_dura3 salt_dura3


suc_dura3 = [mean(suc_dura2(1:3)) mean(suc_dura2(4:6)) mean(suc_dura2(7:9)) mean(suc_dura2(10:12))...
    mean(suc_dura2(13:15)) mean(suc_dura2(16:19)) mean(suc_dura2(20:23)) mean(suc_dura2(24:26))...
    mean(suc_dura2(27:29)) mean(suc_dura2(30:35)) mean(suc_dura2(36:39)) mean(suc_dura2(40:43))...
    mean(suc_dura2(44:48)) mean(suc_dura2(49:53)) mean(suc_dura2(54:57))];
salt_dura3 = [mean(salt_dura2(1:3)) mean(salt_dura2(4:6)) mean(salt_dura2(7:9)) mean(salt_dura2(10:12))...
    mean(salt_dura2(13:15)) mean(salt_dura2(16:19)) mean(salt_dura2(20:23)) mean(salt_dura2(24:26))...
    mean(salt_dura2(27:29)) mean(salt_dura2(30:35)) mean(salt_dura2(36:39)) mean(salt_dura2(40:43))...
    mean(salt_dura2(44:48)) mean(salt_dura2(49:53)) mean(salt_dura2(54:57))];

 subplot(1,5,2)     %central lick frequency
std_suc = std(suc_dura3); std_salt = std(salt_dura3);
avg_suc = mean(suc_dura3);
avg_salt = mean(salt_dura3);
p = bar(1:2,[avg_suc,avg_salt]);hold on
p.FaceColor = 'flat';
p.CData(1,:) = [0 0 1];
p.CData(2,:) = [1 0 0];
ylabel('Central lick frequency (Hz)');xticklabels({'Suc', 'NaCl'});
er = errorbar([1,2], [avg_suc,avg_salt],zeros(1,2),[std_suc/sqrt(size(suc_dura3,2)),std_salt/sqrt(size(salt_dura3,2))],'k','LineWidth',1.2);hold on
scatter(ones(size(suc_dura3)),suc_dura3,200,'k.'); hold on;
scatter(ones(size(salt_dura3))+1,salt_dura3,200,'k.'); hold on;
set(er,'linestyle','none');
for i = 1:size(suc_dura3,2)
    line([1,2],[suc_dura3(1,i) salt_dura3(1,i)],'Color','k')
end
box off
xlim([0 3]);

[~,pvalues.central_lick_freq,~,stats_table.central_lick_freq] = ttest(suc_dura3,salt_dura3);
disp({'suc perf: ',num2str(avg_suc)}); disp({'nacl perf: ',num2str(avg_salt)})
disp({'suc std: ',num2str(std_suc/sqrt(size(suc_corr1,2)))}); disp({'nacl std: ',num2str(std_salt/sqrt(size(salt_corr1,2)))})



clear left_dura3 right_dura3 std_left std_right avg_left avg_right

left_dura3 = [mean(left_dura2(1:3)) mean(left_dura2(4:6)) mean(left_dura2(7:9)) mean(left_dura2(10:12))...
    mean(left_dura2(13:15)) mean(left_dura2(16:19)) mean(left_dura2(20:23)) mean(left_dura2(24:26))...
    mean(left_dura2(27:29)) mean(left_dura2(30:35)) mean(left_dura2(36:39)) mean(left_dura2(40:43))...
    mean(left_dura2(44:48)) mean(left_dura2(49:53)) mean(left_dura2(54:57))];
right_dura3 = [mean(right_dura2(1:3)) mean(right_dura2(4:6)) mean(right_dura2(7:9)) mean(right_dura2(10:12))...
    mean(right_dura2(13:15)) mean(right_dura2(16:19)) mean(right_dura2(20:23)) mean(right_dura2(24:26))...
    mean(right_dura2(27:29)) mean(right_dura2(30:35)) mean(right_dura2(36:39)) mean(right_dura2(40:43))...
    mean(right_dura2(44:48)) mean(right_dura2(49:53)) mean(right_dura2(54:57))];
std_left = std(left_dura3); std_right = std(right_dura3);
avg_left = mean(left_dura3);
avg_right = mean(right_dura3);
 subplot(1,5,3) % Lateral lick frequency
p = bar(1:2,[avg_right,avg_left]);hold on
p.FaceColor = 'flat';
p.CData(1,:) = [0 0 1];
p.CData(2,:) = [1 0 0];
er = errorbar([1,2], [avg_right,avg_left],zeros(1,2),[std_right/sqrt(size(right_dura3,2)),std_left/sqrt(size(left_dura3,2))],'k','LineWidth',1.2);hold on
scatter(ones(size(right_dura3)),right_dura3,200,'k.'); hold on;
scatter(ones(size(left_dura3))+1,left_dura3,200,'k.'); hold on;
for i = 1:size(left_dura3,2)
    line([1,2],[right_dura3(1,i) left_dura3(1,i)],'Color','k')
end
ylabel('Lateral lick frequency (Hz)');xticklabels({'Right', 'Left'});

set(er,'linestyle','none');
box off
xlim([0 3]);
[~,pvalues.lateral_lick_freq,~,stats_table.lateral_lick_freq] = ttest(right_dura3,left_dura3);
disp({'left lick rate: ',num2str(avg_left)}); disp({'right lick rate',num2str(avg_right)})
disp({'leftlick rate std: ',num2str(std_right/sqrt(size(right_dura3,2)))}); disp({'rightlick rate std: ',num2str(std_left/sqrt(size(left_dura3,2)))})

clearvars x y

xlim([0 3]);

clear left_react3 right_react3


left_react3 = [mean(left_react2(1:3)) mean(left_react2(4:6)) mean(left_react2(7:9)) mean(left_react2(10:12))...
    mean(left_react2(13:15)) mean(left_react2(16:19)) mean(left_react2(20:23)) mean(left_react2(24:26))...
    mean(left_react2(27:29)) mean(left_react2(30:35)) mean(left_react2(36:39)) mean(left_react2(40:43))...
    mean(left_react2(44:48)) mean(left_react2(49:53)) mean(left_react2(54:57))];
right_react3 = [mean(right_react2(1:3)) mean(right_react2(4:6)) mean(right_react2(7:9)) mean(right_react2(10:12))...
    mean(right_react2(13:15)) mean(right_react2(16:19)) mean(right_react2(20:23)) mean(right_react2(24:26))...
    mean(right_react2(27:29)) mean(right_react2(30:35)) mean(right_react2(36:39)) mean(right_react2(40:43))...
    mean(right_react2(44:48)) mean(right_react2(49:53)) mean(right_react2(54:57))];

avg_left = mean(left_react3);
avg_right = mean(right_react3);
std_left = std(left_react3); std_right = std(right_react3);

 subplot(1,5,4)
p = bar(1:2,[avg_right,avg_left]);hold on
p.FaceColor = 'flat';
p.CData(1,:) = [0 0 1];
p.CData(2,:) = [1 0 0];
er = errorbar([1,2], [avg_right,avg_left],zeros(1,2),[std_right/sqrt(size(right_react3,2)),std_left/sqrt(size(left_react3,2))],'k','LineWidth',1.2);hold on
scatter(ones(size(right_react3)),right_react3,200,'k.'); hold on;
scatter(ones(size(left_react3))+1,left_react3,200,'k.'); hold on;
for i = 1:size(suc_dura3,2)
    line([1,2],[right_react3(1,i) left_react3(1,i)],'Color','k')
end
ylabel('Reaction time (s)');xticklabels({'Right', 'Left'});
ylim([0 1.5]);
box off
xlim([0 3]);
set(er,'linestyle','none');
[~,pvalues.reaction_time,~,stats_table.reaction_time] = ttest(left_react3,right_react3);
disp({'left_react: ',num2str(avg_left)}); disp({'right_react',num2str(avg_right)})
disp({'left_react std: ',num2str(std_right/sqrt(size(right_react3,2)))}); disp({'right_react std: ',num2str(std_left/sqrt(size(left_react3,2)))})


 subplot(1,5,5)
 num_trials = 107;
trial_FC = trial_1_valve{1,19}; 
for n =1:num_trials
    if trial_FC(n).L_R_trial == 1
        h1 =  plot([-1 -1],[n-0.3 n+0.3],'r','LineWidth', 1.5);hold on;
        h1a = plot([trial_FC(n).T_1 trial_FC(n).T_1],[n-0.3 n+0.3],'r','LineWidth', 1.5);hold on;
    elseif trial_FC(n).L_R_trial == 2
        h2 =  plot([-1 -1],[n-0.3 n+0.3],'b','LineWidth', 1.5);hold on;
        h2a = plot([trial_FC(n).T_9 trial_FC(n).T_9],[n-0.3 n+0.3],'b','LineWidth', 1.5);hold on;
        
    end
end

for n = 1:num_trials
    j = ones(1,length(trial_FC(n).centSp(1:5)))*n;
    
    h3=scatter(trial_FC(n).centSp(1:5),j,4,'filled','k'); %central licks are red
    
    hold on
end
hold on
for n = 1:num_trials
    j = ones(1,length(trial_FC(n).centSp(6:end)))*n;
    
    h3=scatter(trial_FC(n).centSp(6:end),j,4,'filled','c'); %central licks are red
    
    hold on
end
hold on


for n = 1:num_trials
    j = ones(1,length(trial_FC(n).RightSp))*n;
    h4=scatter(trial_FC(n).RightSp,j,4,'filled','b'); %central licks are red
    
    hold on
end
hold on

for n = 1:num_trials
    j = ones(1,length(trial_FC(n).LeftSp))*n;
    h5=scatter(trial_FC(n).LeftSp,j,4,'filled','r'); %central licks are red
    
    hold on
end
hold on

% Plot Green/Magenta to display incorrect or correct trial
for n = 1:num_trials
    
    if trial_FC(n).correct_choice == 1 && trial_FC(n).L_R_trial == 1
        h6 =  plot([7 7.5],[n n],'g', 'LineWidth', 1.5);hold on;
    elseif trial_FC(n).correct_choice == 0 && trial_FC(n).L_R_trial == 1
        h6 = plot([7 7.5],[n n],'m', 'LineWidth', 1.5);hold on;
    elseif trial_FC(n).correct_choice == 1 && trial_FC(n).L_R_trial == 2
        h6 = plot([7.5 8],[n n],'g', 'LineWidth', 1.5);hold on;
    elseif trial_FC(n).correct_choice == 0 && trial_FC(n).L_R_trial == 2
        h6 = plot([7.5 8],[n n],'m', 'LineWidth', 1.5);hold on;
    end
end
box off
set(gcf, 'Position',  [-1479 59 620 929])
% legend([h1,h2,h3,h4,h5],'Sucrose','NaCl','Central', 'Right', 'Left','Location','northwest')
title('Example Session')
ylabel('Trial #');
xlabel('Time (sec)')
xlim([-2 9]);
ylim([0 num_trials+5]);
set(gcf, 'Position',  [-1630 549 1406 183])

 
 %% Figure 1: Fiber Photometry data analyses

clear

load('all_photometry_sessions.mat')

load('photometry_summary_control_mice.mat') %%Load in GRAB_DAmut data


% addpath('D:\JohnData\Matlab\Fiber_photometry_matlab')


addpath(genpath('D:\JohnData\Matlab\bounded_line'));

close all;
subplot_dim = [4,4];
clims = [-7 7];xlimits = [-8 2]; ylimits = [-2 7];
map1 = [255,255,255
        255,245,240
        254,224,210
        252,187,161
        252,146,114
        251,106,74
        239,59,44]./255;
    
clear all_resp
all_resp = all_resp_1(2, 2);


all_trials = [all_resp.summarydata.LL.correct_contra_zdff; all_resp.summarydata.LL.correct_ipsi_zdff*-1];
figure
subplot(subplot_dim(1),subplot_dim(2),[4,8])
imagesc(all_resp.summarydata.LL.x_time,1:size(all_trials,1),all_trials,clims);
colormap(flipud(redblue))
colorbar
hold on;
line([0 0],[-1 size(all_trials,1)],'Color','k','LineWidth',1.5,'LineStyle','--'); hold on;

rectangle('Pos',[-8 0 1 size(all_resp.summarydata.LL.correct_contra_zdff,1)],'FaceColor','b'); hold on;
rectangle('Pos',[-8 size(all_resp.summarydata.LL.correct_contra_zdff,1) 1 size(all_resp.summarydata.LL.correct_ipsi_zdff,1)],'FaceColor','r'); hold on;
rectangle('Pos',[-8 size(all_resp.summarydata.LL.correct_contra_zdff,1)+size(all_resp.summarydata.LL.correct_ipsi_zdff,1) 1 size(all_resp.summarydata.LL.incorrect_contra_zdff,1)],'FaceColor',[0 0 1],'LineStyle','--'); hold on;
rectangle('Pos',[-8 size(all_resp.summarydata.LL.correct_contra_zdff,1)+size(all_resp.summarydata.LL.correct_ipsi_zdff,1)+size(all_resp.summarydata.LL.incorrect_contra_zdff,1) 1 size(all_resp.summarydata.LL.incorrect_ipsi_zdff,1)],'FaceColor',[1 0 0],'LineStyle','--'); hold on;

title('Example Mouse #1'); ylabel('Trial #');xlabel('Time (s)');
xlim(xlimits);
box off;
yticks(gca,[1 size(all_trials,1)])
yticklabels(gca,{'1',num2str(size(all_trials,1))})
subplot(subplot_dim(1),subplot_dim(2),12)
boundedline(all_resp.summarydata.LL.x_time,mean(all_resp.summarydata.LL.correct_contra_zdff,1),std(all_resp.summarydata.LL.correct_contra_zdff)/sqrt(size(all_resp.summarydata.LL.correct_contra_zdff,1)),'alpha','b','LineWidth',2); hold on;
hold on;
boundedline(all_resp.summarydata.LL.x_time,mean(all_resp.summarydata.LL.correct_ipsi_zdff,1),std(all_resp.summarydata.LL.correct_ipsi_zdff)/sqrt(size(all_resp.summarydata.LL.correct_ipsi_zdff,1)),'alpha','r','LineWidth',2); hold on;
hold on;
line([0 0],[-3 10],'Color','k','LineWidth',1.5,'LineStyle','--'); hold on;
line([-10 10],[0 0],'Color','k','LineWidth',1.5,'LineStyle','--')

xlim(xlimits);ylim(ylimits);
xlabel('Time (s)'); ylabel('DA z-\DeltaF/F')
set(gca,'LineWidth',1); % The only other option is 'in'
set(gca,'TickDir','out'); % The only other option is 'in'


N = 200;
% 
% 
GFP_ds  = arrayfun(@(i) mean(deltaFoverF(i:i+N-1)), 1:N:length(deltaFoverF)-N+1);

time_ds  = arrayfun(@(i) mean(time(i:i+N-1)), 1:N:length(time)-N+1);
%
clear x
y = [trial_by_trial.TDT_central_licks];
figure(1)
h1 = subplot(subplot_dim(1),subplot_dim(2),[1,2,3]);
plot(time_ds,GFP_ds,'k'); hold on;


for i = 1:size(trial_by_trial,2)
    x = trial_by_trial(i).all_TDT_central_licks(1:6);
    line([x x],[0.6 0.8],'Color','k','LineWidth',0.5); hold on
     x1 = trial_by_trial(i).all_TDT_central_licks(7:end);
    line([x1 x1],[0.6 0.8],'Color','c','LineWidth',0.5); hold on
    rectangle('Pos',[x(1) -0.7 (x(end)-x(1)) 1.7],'FaceColor', [0.7, 0.7, 0.7, 0.5],'LineStyle','none');
    rectangle('Pos',[x1(1) -0.7 (x1(end)-x1(1)) 1.7],'FaceColor', [0, 1, 1, 0.5],'LineStyle','none');
    if ~isempty(trial_by_trial(i).all_TDT_left_licks)
        y = trial_by_trial(i).all_TDT_left_licks;
        line([y y],[0.6 0.8],'Color','r','LineWidth',0.5); hold on
        rectangle('Pos',[y(1) -0.7 (y(end)-y(1)) 1.7],'FaceColor', [1, 0, 0, 0.5],'LineStyle','none');
    end
    
    if ~isempty(trial_by_trial(i).all_TDT_right_licks)
        z = trial_by_trial(i).all_TDT_right_licks;
        line([z z],[0.6 0.8],'Color','b','LineWidth',0.5); hold on
        rectangle('Pos',[z(1) -0.7 (z(end)-z(1)) 1.7],'FaceColor', [0, 0, 1, 0.5],'LineStyle','none');
    end
    
    hold on;
end
line([950 955],[-0.6 -0.6],'Color','k','LineWidth',1.5); hold on
line([955 955],[-0.6 -0.1],'Color','k','LineWidth',1.5)

xlim([840 955]); ylim([-0.7 1])
box off
ylabel('\Delta F/F'); xlabel('Time (s)');
set(gca,'TickDir','out'); % The only other option is 'in'
set(gca,'LineWidth',1); % The only other option is 'in'

%
clear all_resp
% 



%
clear output1 output2 output3 output4 output5 output6 output7 output8 mean_by_animal AUC_by_epoch
output1 = vertcat(mouse_average_1.stimulus_group1);%%sucrose
output2 = vertcat(mouse_average_1.stimulus_group2);%%nacl
output3 = vertcat(mouse_average_1.direction_group1);
output4 = vertcat(mouse_average_1.direction_group2);
output5 = vertcat(mouse_average_1.error_group1);
output6 = vertcat(mouse_average_1.error_group2);
output7 = vertcat(mouse_average_1.stimulus_group1_FL);
output8 = vertcat(mouse_average_1.stimulus_group2_FL);

output_time_FL = vertcat(mouse_average_1.time_FL);
output_time_LL = vertcat(mouse_average_1.time_LL);
time_FL = mean(output_time_FL);
time_LL = mean(output_time_LL);
figure(1)
clims = [0 3];
    map1 = [255,255,255
        255,245,240
        254,224,210
        252,187,161
        252,146,114
        251,106,74
        239,59,44]./255;
    
    map2 = [255,255,255
        247,251,255
        222,235,247
        198,219,239
        158,202,225
        107,174,214
        66,146,198
        33,113,181]./255;
    
subplot(subplot_dim(1),subplot_dim(2),5)
rectangle('Pos',[0 min(mean(output1,1))-0.5 2 max(mean(output1,1))+1.5],'FaceColor', [0.7, 0.7, 0.7, 0.4],'LineStyle','none')
b1 = boundedline(mean(output_time_FL,1),mean(output1,1),std(output1)/sqrt(size(output1,1)),'alpha','b','LineWidth',1); hold on;
b2 = boundedline(mean(output_time_FL,1),mean(output2,1),std(output2)/sqrt(size(output2,1)),'alpha','r','LineWidth',1); hold on;
line([0 0],[-1 7],'Color','k','LineWidth',1.5,'LineStyle','--');hold on;
line([-8 8],[0 0],'Color','k','LineStyle','--','LineWidth',1.5);
text(0,5.5,'First dry lick')
xlim([-3 3]); ylim([min(mean(output1,1))-0.5 6]);
xlabel('Time (s)'); ylabel('DA z-\DeltaF/F'); legend([b1;b2],{'NaCl','Suc'},'Location','northwest'); legend boxoff  
    set(gca,'TickDir','out'); set(gca,'linewidth',1)
yticks([0:1:6])

subplot(subplot_dim(1),subplot_dim(2),6)
rectangle('Pos',[-1.5 min(mean(output3,1))-0.5 3 max(mean(output3,1))+1.5],'FaceColor', [0.7, 0.7, 0.7, 0.4],'LineStyle','none')
b1 = boundedline(mean(output_time_LL,1),mean(output3,1),std(output3)/sqrt(size(output3,1)),'alpha','b','LineWidth',1); hold on;
b2 = boundedline(mean(output_time_LL,1),mean(output4,1),std(output4)/sqrt(size(output4,1)),'alpha','r','LineWidth',1); hold on;
line([0 0],[-1 7],'Color','k','LineWidth',1.5,'LineStyle','--');hold on;
line([-8 8],[0 0],'Color','k','LineStyle','--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('DA z-\DeltaF/F')
xlim([-3 2]); ylim([min(mean(output1,1))-0.5 6]);title('Averaged dopamine response (n = 7, sessions = 31)');
text(-1,5.5,'First lateral lick');

legend([b1;b2],{'Contra','Ipsi'},'Location','northwest'); legend boxoff
    set(gca,'TickDir','out'); set(gca,'linewidth',1);yticks([0:1:6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT CONTROL MICE DATA ON SAME PLOT AS EXPERIMENTAL

subplot(subplot_dim(1),subplot_dim(2),5)
clear x y
x = control_mice_photometry.stim1;
y = control_mice_photometry.stim2;
boundedline(mean(control_mice_photometry.output_time_FL,1),mean([x;y],1),std([x;y])/sqrt(size([x;y],1)),'alpha','k','LineWidth',1); hold on;

subplot(subplot_dim(1),subplot_dim(2),6)
clear x y
x = control_mice_photometry.direc1;
y = control_mice_photometry.direc2;
boundedline(mean(control_mice_photometry.output_time_LL,1),mean([x;y],1),std([x;y])/sqrt(size([x;y],1)),'alpha','k','LineWidth',1); hold on;
% boundedline(mean(control_mice_photometry.output_time_LL,1),mean(y,1),std(y)/sqrt(size([x;y],1)),'alpha','c','LineWidth',1); hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
time_FL = mean(output_time_FL);
time_LL = mean(output_time_LL);
mean_by_animal.stimulus.NaCl_pre = mean(output1(:,find(time_FL(1,:)>-1 & time_FL(1,:)<0)),2);
mean_by_animal.stimulus.Suc_pre = mean(output2(:,find(time_FL(1,:)>-1 & time_FL(1,:)<0)),2);

mean_by_animal.stimulus.NaCl = max(output1(:,find(time_FL(1,:)>0 & time_FL(1,:)<1.5)),[],2);
mean_by_animal.stimulus.Suc = max(output2(:,find(time_FL(1,:)>0 & time_FL(1,:)<1.5)),[],2);

mean_by_animal.direction.contra = max(output3(:,find(time_LL(1,:)>-1.5 & time_LL(1,:)<0)),[],2);
mean_by_animal.direction.ipsi = max(output4(:,find(time_LL(1,:)>-1.5 & time_LL(1,:)<0)),[],2);

mean_by_animal.direction.contra_choice = max(output3(:,find(time_LL(1,: )>0 & time_LL(1,:)<1.5)),[],2);
mean_by_animal.direction.ipsi_choice= max(output4(:,find(time_LL(1,:)>0 & time_LL(1,:)<1.5)),[],2);

mean_by_animal.outcome.correct = max(output5(:,find(time_LL(1,:)>=2 & time_LL(1,:)<3)),[],2);
mean_by_animal.outcome.error = max(output6(:,find(time_LL(1,:)>=2 & time_LL(1,:)<3)),[],2);
for i = 1:size(output3,1)
     AUC_by_epoch.baseline(i,1) = trapz(output1(i,time_FL(1,:) > -4 & time_FL(1,:) <-2.5));
    AUC_by_epoch.baseline(i,2) = trapz(output2(i,time_FL(1,:) > -4 & time_FL(1,:) <-2.5));
    AUC_by_epoch.sample(i,1) = trapz(output1(i,time_FL(1,:) > 0 & time_FL(1,:) <1.5));
    AUC_by_epoch.sample(i,2) = trapz(output2(i,time_FL(1,:) > 0 & time_FL(1,:) <1.5));
    AUC_by_epoch.delay(i,1) = trapz(output3(i,time_LL(1,:) > -1.5 & time_LL(1,:) <0));
    AUC_by_epoch.delay(i,2) = trapz(output4(i,time_LL(1,:) > -1.5 & time_LL(1,:) <0));
    AUC_by_epoch.choice(i,1) = trapz(output3(i,time_LL(1,:) > 0 & time_LL(1,:) <1.5));
    AUC_by_epoch.choice(i,2) = trapz(output4(i,time_LL(1,:) > 0 & time_LL(1,:) <1.5));
end


figure(1)

%
stim_selec = (output1-output2)/2;
direct_selec = (output3-output4)/2;
outcome_selec = (output5-output6)/2;
time_FL = mean(output_time_FL);
time_LL = mean(output_time_LL);

subplot(subplot_dim(1),subplot_dim(2),10)
rectangle('Pos',[-1.5 -0.2 3 max(mean(direct_selec,1))+2],'FaceColor', [0.7, 0.7, 0.7, 0.4],'LineStyle','none')
hold on;
boundedline(mean(output_time_LL),mean(direct_selec),std(direct_selec)/sqrt(size(direct_selec,1)),'alpha','k','LineWidth',1); hold on;

line([0 0],[-1 7],'Color','k','LineWidth',1.5,'LineStyle','--');hold on;
line([-8 8],[0 0],'Color','k','LineStyle','--','LineWidth',1.5);
xlim([-3 2]);ylim([-0.2 0.4]);yticks([-0.2:0.2:0.4])
title('Direction Selectivity');xlabel('Time (s)'); ylabel('Selectivity (DA z-\DeltaF/F)')
set(gca,'LineWidth',1); % The only other option is 'in'
set(gca,'TickDir','out'); % The only other option is 'in'

subplot(subplot_dim(1),subplot_dim(2),9)
rectangle('Pos',[0 min(mean(output1,1))-0.5 2 max(mean(output1,1))+1.5],'FaceColor', [0.7, 0.7, 0.7, 0.4],'LineStyle','none')

hold on;
boundedline(mean(output_time_FL),mean(stim_selec),std(stim_selec)/sqrt(size(stim_selec,1)),'alpha','k','LineWidth',1); hold on;

line([0 0],[-1 7],'Color','k','LineWidth',1.5,'LineStyle','--');hold on;
line([-8 8],[0 0],'Color','k','LineStyle','--','LineWidth',1.5);
xlim([-3 3]);ylim([-0.2 0.4]);yticks([-0.2:0.2:0.4])
title('Stimulus Selectivity');xlabel('Time (s)'); ylabel('Selectivity (DA z-\DeltaF/F)')
set(gca,'LineWidth',1); % The only other option is 'in'
set(gca,'TickDir','out'); % The only other option is 'in'



selectivity.sensory.baseline = mean(stim_selec(:,find(time_FL(1,:)>-5.5 & time_FL(1,:)<-4)),2);

selectivity.sensory.sample = mean(stim_selec(:,find(time_FL(1,:)>0 & time_FL(1,:)<1.5)),2);


selectivity.direction.baseline = mean(direct_selec(:,find(time_LL(1,:)>-10.5 & time_LL(1,:)<-9)),2);

selectivity.direction.sample = mean(direct_selec(:,find(time_LL(1,:)>-5.6& time_LL(1,:)<-4.1)),2);

selectivity.direction.delay = mean(direct_selec(:,find(time_LL(1,:)>-1.5 & time_LL(1,:)<0)),2);

selectivity.direction.choice = mean(direct_selec(:,find(time_LL(1,:)>0 & time_LL(1,:)<1.5)),2);

clear p
[h,p.sample] = ttest(selectivity.sensory.baseline,selectivity.sensory.sample,'tail','both')

[h,p.delay] = ttest(selectivity.direction.baseline,selectivity.direction.delay,'tail','both')

[h,p.choice] = ttest(selectivity.direction.baseline,selectivity.direction.choice,'tail','both')

[h,p.sample] = ttest(selectivity.sensory.baseline,selectivity.sensory.sample,'tail','both')

[h,p.delay_sample] = ttest(selectivity.sensory.sample,selectivity.direction.delay,'tail','both')

[h,p.delay_choice] = ttest(selectivity.direction.delay,selectivity.direction.choice,'tail','both')

[h,p.sample_choice] = ttest(selectivity.sensory.sample,selectivity.direction.choice,'tail','both')




subplot(subplot_dim(1),subplot_dim(2),14)
b = bar([1:3],[mean(selectivity.direction.sample),mean(selectivity.direction.delay),mean(selectivity.direction.choice)],'FaceColor','none','EdgeColor','k','LineWidth',1.5); 

hold on;
e = errorbar([1:3],[mean(selectivity.direction.sample),mean(selectivity.direction.delay),mean(selectivity.direction.choice)],[std(selectivity.direction.sample),std(selectivity.direction.delay),std(selectivity.direction.choice)]./sqrt(size(selectivity.sensory.sample,1)),'.','Color', 'k'); 
scatter(ones(size(selectivity.direction.sample)),selectivity.sensory.sample,70,'k.', 'MarkerFaceAlpha',.1);
hold on;
scatter(ones(size(selectivity.direction.delay))+1,selectivity.direction.delay,70,'k.');
hold on;
scatter(ones(size(selectivity.direction.choice))+2,selectivity.direction.choice,70,'k.');

for i = 1:size(selectivity.sensory.baseline,1)
    plot([1,2],[selectivity.direction.sample(i),selectivity.direction.delay(i)],'Color',[0.7 0.7 0.7 0.5]);hold on;
end
for i = 1:size(selectivity.sensory.baseline,1)
    plot([2,3],[selectivity.direction.delay(i),selectivity.direction.choice(i)],'Color',[0.7 0.7 0.7 0.5]);hold on;
end
ylabel('Selectivity (DA z-\DeltaF/F)');
ylim([-1 1]);xlim([0.5 3.5]);
xticks([1 2 3 4])
xticklabels({'Sample', 'Delay', 'Choice'});xtickangle(0)
yticks([-1:0.5:1])
set(gca,'LineWidth',1); % The only other option is 'in'
set(gca,'TickDir','out'); % The only other option is 'in'

box off

if p.sample < 0.05
    text(0.9,0.95,'*','FontSize',12)
end
if p.delay < 0.05
    text(1.9,0.95,'*','FontSize',20)
end
if p.choice < 0.05
    text(2.9,0.95,'*','FontSize',20)
end
if p.delay_choice < 0.05
    text(2.4,1.3,'+','FontSize',12);
    line([2 3],[1.2 1.2],'Color','k')
end
if p.sample_choice < 0.05
    text(2,1.55,'+','FontSize',12);
    line([1 3],[1.45 1.45],'Color','k')
end


marker_size = 25;
subplot(subplot_dim(1),subplot_dim(2),13)
clearvars x y 

x = AUC_by_epoch.sample(:,1);
y = AUC_by_epoch.sample(:,2);

b = bar([1,2],[mean(x),mean(y)],'LineWidth',1); 
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
hold on;
scatter(ones(size(mean_by_animal.stimulus.NaCl,1))*1,x,marker_size,'k.');
hold on;
scatter(ones(size(mean_by_animal.stimulus.Suc,1))*2,y,marker_size,'k.');
for k = 1:size(mean_by_animal.direction.contra_choice,1)
    line([1,2],[x(k) y(k)],'Color','k','LineWidth',1);
    
end
error_1 = [std(x)/sqrt(size(x,1));...
    std(y)/sqrt(size(y,1))];
e = errorbar([1],[mean(x)],error_1(1)','k');
hold on;
e = errorbar([2],[mean(y)],error_1(2)','k');
%
e.LineStyle = 'none';
[~,figures_stats.subplot13.pval(1),~,figures_stats.subplot13.table(1)] = ttest(x,y,'tail','both');
if figures_stats.subplot13.pval(1)>0.05
   line([1 2],[9500 9500],'Color','k','LineWidth',1);
   hold on;
   text(1.3,9800,'ns')
end
subplot(subplot_dim(1),subplot_dim(2),13)
clearvars x y 

x = AUC_by_epoch.delay(:,1);
y = AUC_by_epoch.delay(:,2);

b = bar([4,5],[mean(x),mean(y)],'LineWidth',1); 
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
hold on;
scatter(ones(size(mean_by_animal.direction.contra,1))*4,x,marker_size,'k.');
hold on;
scatter(ones(size(mean_by_animal.direction.ipsi,1))*5,y,marker_size,'k.');
for k = 1:size(mean_by_animal.direction.contra_choice,1)
    line([4,5],[x(k) y(k)],'Color','k','LineWidth',1);
    
end
error_1 = [std(x)/sqrt(size(x,1));...
    std(y)/sqrt(size(y,1))];
e = errorbar([4],[mean(x)],error_1(1)','k');
hold on;
e = errorbar([5],[mean(y)],error_1(2)','k');
e.LineStyle = 'none';
title('Delay Epoch');ylabel('Mean DA z-\DeltaF/F');
xticklabels({'Correct Trials','Error Trials'})
box off
[~,figures_stats.subplot13.pval(2),~,figures_stats.subplot13.table(2)] = ttest(x,y,'tail','both');
if figures_stats.subplot13.pval(2)<0.01
   line([4 5],[9500 9500],'Color','k','LineWidth',1);
   hold on;
   text(4.3,9600,'**')
end
subplot(subplot_dim(1),subplot_dim(2),13)
clearvars x y 

x = AUC_by_epoch.choice(:,1);
y = AUC_by_epoch.choice(:,2);

b = bar([7,8],[mean(x),mean(y)],'LineWidth',1); 
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
hold on;
scatter(ones(size(mean_by_animal.direction.contra_choice,1))*7,x,marker_size,'k.');
hold on;
scatter(ones(size(mean_by_animal.direction.ipsi_choice,1))*8,y,marker_size,'k.');
hold on;
for k = 1:size(mean_by_animal.direction.contra_choice,1)
    line([7,8],[x(k) y(k)],'Color','k','LineWidth',1);
    
end
error_1 = [std(x)/sqrt(size(x,1));...
    std(y)/sqrt(size(y,1))];
e = errorbar([7],[mean(x)],error_1(1)','k');
hold on;
e = errorbar([8],[mean(y)],error_1(2)','k');

e.LineStyle = 'none';
[~,figures_stats.subplot13.pval(3),~,figures_stats.subplot13.table(3)] = ttest(x,y,'tail','both');
if figures_stats.subplot13.pval(3)<0.001
   line([7 8],[9500 9500],'Color','k','LineWidth',1);
   hold on;
   text(7.2,9600,'***')
end

 xlim([0 9]);
 set(gca,'LineWidth',1); % The only other option is 'in'
set(gca,'TickDir','out'); % The only other option is 'in'

title('Average DA response by epoch');ylabel('AUC (A.U.)');
yticks([0:5000:10000]);

xticks([1.5 4.5 7.5]);
xticklabels({'Sample','Delay','Choice'})
box off
%
