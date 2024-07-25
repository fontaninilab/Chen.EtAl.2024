function [d1r_mean_activity, non_mean_activity] = plot_mean_response_for_all_neurons(group1a,group2a,group1, group2,number_of_error_trials,randomized);
figure_startup
start_index = 1;
% clims = [0 1];
% index_rank = 23;
% threshold = 0;
% x_limits = [-3 2];
% y_limits = [-0.5 8];
% clear A
% A = 1:size(group1,2);
clear correct_trials_L incorrect_trials_L correct_trials_R incorrect_trials_R direction_e_right direction_e_left 
for j =1:size(group1a,2)
   
    correct_trials_L = find(group1(group1a(j)).Left.correct==1);
    incorrect_trials_L = find(group1(group1a(j)).Left.incorrect==1);
    correct_trials_R = find(group1(group1a(j)).Right.correct==1);
    incorrect_trials_R = find(group1(group1a(j)).Right.incorrect==1);
    
    correct_trials_suc = find(group1(group1a(j)).Sucrose.correct==1);
    incorrect_trials_suc = find(group1(group1a(j)).Sucrose.incorrect==1);
    correct_trials_nacl = find(group1(group1a(j)).NaCl.correct==1);
    incorrect_trials_nacl = find(group1(group1a(j)).NaCl.incorrect==1);
    
    if randomized == 1
        correct_trials_L_random = randsample(size(correct_trials_L,1),round(size(correct_trials_L,1)/2),true);
        
        correct_trials_R_random = randsample(size(correct_trials_R,1),round(size(correct_trials_R,1)/2),true);
        
        correct_trials_suc_random = randsample(size(correct_trials_suc,1),round(size(correct_trials_suc,1)/2),true);
        
        correct_trials_nacl_random = randsample(size(correct_trials_nacl,1),round(size(correct_trials_nacl,1)/2),true);
    else
        correct_trials_L_random = correct_trials_L;
        
        correct_trials_R_random = correct_trials_R;
        
        correct_trials_suc_random = correct_trials_suc;
        
        correct_trials_nacl_random = correct_trials_nacl;
    end
    
    directionA11(j,:) = mean(group1(group1a(j)).Left.dff_aligned(correct_trials_L_random,start_index:end));
    time_directionA1(j,:) = mean(group1(group1a(j)).Left.frames_aligned(correct_trials_L_random,start_index:end));
    directionB11(j,:) = mean(group1(group1a(j)).Right.dff_aligned(correct_trials_R_random,start_index:end));
    time_directionB1(j,:) = mean(group1(group1a(j)).Right.frames_aligned(correct_trials_R_random,start_index:end));
    
    directionA11_not(j,:) = mean(group1(group1a(j)).Left.dff_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
    time_directionA1_not(j,:) = mean(group1(group1a(j)).Left.frames_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
    
    directionB11_not(j,:) = mean(group1(group1a(j)).Right.dff_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));
    time_directionB1_not(j,:) = mean(group1(group1a(j)).Right.frames_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));
    

    
    if size(incorrect_trials_L,1) >= number_of_error_trials && size(incorrect_trials_R,1) >= number_of_error_trials
        direction_e_left(j,:) = mean(group1(group1a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end),1);
        direction_e_right(j,:) = mean(group1(group1a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end),1);
    else
        direction_e_right(j,1:82) = NaN;
        direction_e_left(j,1:82) = NaN;
    end
    taste_correct_trials_suc(j,:) = mean(group1(group1a(j)).Sucrose.dff_aligned(correct_trials_suc_random,start_index:end));
    taste_correct_trials_nacl(j,:) = mean(group1(group1a(j)).NaCl.dff_aligned(correct_trials_nacl_random,start_index:end));
    taste_incorrect_trials_suc(j,:) = mean(group1(group1a(j)).Sucrose.dff_aligned(incorrect_trials_suc,start_index:end),1);
    taste_incorrect_trials_nacl(j,:) = mean(group1(group1a(j)).NaCl.dff_aligned(incorrect_trials_nacl,start_index:end),1);
    
    taste_correct_trials_suc_time(j,:) = mean(group1(group1a(j)).Sucrose.frames_aligned(correct_trials_suc_random,start_index:end));
    taste_correct_trials_nacl_time(j,:) = mean(group1(group1a(j)).NaCl.frames_aligned(correct_trials_nacl_random,start_index:end));
    taste_incorrect_trials_suc_time(j,:) = mean(group1(group1a(j)).Sucrose.frames_aligned(incorrect_trials_suc,start_index:end),1);
    taste_incorrect_trials_nacl_time(j,:) = mean(group1(group1a(j)).NaCl.frames_aligned(incorrect_trials_nacl,start_index:end),1);

%     directionA11(j,:) = mean([group1(group1a(j)).Left.dff_aligned(correct_trials_L,start_index:end);group1(group1a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end)]);
%     if abs(max(directionA11(j,:))) > abs(min(directionA11(j,:)))
        directionA1(j,:) = directionA11(j,:)./max(directionA11(j,:));
%     elseif abs(max(directionA11(j,:))) < abs(min(directionA11(j,:)))
%         directionA1(j,:) = directionA11(j,:)./min(directionA11(j,:));
%     end

%     time_directionA1(j,:) = mean([group1(group1a(j)).Left.frames_aligned(correct_trials_L,start_index:end);group1(group1a(j)).Right.frames_aligned(incorrect_trials_R,start_index:end)]);
%     central_lick_left_trial(j) = mean(group1(group1a(j)).Left.central_licks(:,1));


%         directionB11(j,:) = mean([group1(group1a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end);group1(group1a(j)).Right.dff_aligned(correct_trials_R,start_index:end)]);

%     if abs(max(directionB11(j,:))) > abs(min(directionB11(j,:)))
        directionB1(j,:) = directionB11(j,:)./max(directionB11(j,:));
%     elseif abs(max(directionB11(j,:))) < abs(min(directionB11(j,:)))
%         directionB1(j,:) = directionB11(j,:)./min(directionB11(j,:));
%     end
%     time_directionB1(j,:) = mean([group1(group1a(j)).Left.frames_aligned(incorrect_trials_L,start_index:end);group1(group1a(j)).Right.frames_aligned(correct_trials_R,start_index:end)]);
%     central_lick_right_trial(j) = mean(group1(group1a(j)).Right.central_licks(:,1));
    all_first_central_lick_time1(j) = mean([group1(group1a(j)).Left.central_licks(:,1)]);
    all_first_central_lick_time2(j) = mean([group1(group1a(j)).Right.central_licks(:,1)]);
    all_first_central_lick_time_1(j) = mean([group1(group1a(j)).Left.central_licks(:,6)]);
    all_first_central_lick_time_2(j) = mean([group1(group1a(j)).Right.central_licks(:,6)]);
    all_first_wetcentral_lick_time1(j) = mean(nonzeros([group1(group1a(j)).Left.central_licks(:,end-1:end)]));
    all_first_wetcentral_lick_time2(j) = mean(nonzeros([group1(group1a(j)).Right.central_licks(:,end-1:end)]));
%     d1r_all_trials(j,:) = mean([group1(group1a(j)).Left.dff_aligned(group1(group1a(j)).Left.correct ==1,start_index:end); group1(group1a(j)).Right.dff_aligned(group1(group1a(j)).Right.correct ==1,start_index:end)]);
%     directiond1r_all(j,:) = mean([group1(group1a(j)).Left.dff_aligned(group1(group1a(j)).Left.correct==1,start_index:end);...
%         group1(group1a(j)).Left.dff_aligned(group1(group1a(j)).Left.incorrect==1,start_index:end);...
%         group1(group1a(j)).Right.dff_aligned(group1(group1a(j)).Right.correct==1,start_index:end);...
%         group1(group1a(j)).Right.dff_aligned(group1(group1a(j)).Right.incorrect==1,start_index:end)]);
    
    directiond1r_all(j,:) = mean([mean(group1(group1a(j)).Left.dff_aligned(group1(group1a(j)).Left.correct==1,start_index:end),1);...
        mean(group1(group1a(j)).Right.dff_aligned(group1(group1a(j)).Right.correct==1,start_index:end),1)]);
    
    directiond1r_all_time(j,:) = mean([mean(group1(group1a(j)).Left.frames_aligned(group1(group1a(j)).Left.correct==1,start_index:end),1);...
        mean(group1(group1a(j)).Right.frames_aligned(group1(group1a(j)).Right.correct==1,start_index:end),1)]);
%     directiond1r_all_taste_aligned(j,:) = mean([mean(group1(group1a(j)).Sucrose.dff_aligned(group1(group1a(j)).Sucrose.correct==1,start_index:end));...
%         mean(group1(group1a(j)).NaCl.dff_aligned(group1(group1a(j)).NaCl.correct==1,start_index:end))]);
%     
    directiond1r_all_taste_aligned(j,:) = mean([mean(group1(group1a(j)).Sucrose.dff_aligned(group1(group1a(j)).Sucrose.correct==1,start_index:end),1);...
        mean(group1(group1a(j)).NaCl.dff_aligned(group1(group1a(j)).NaCl.correct==1,start_index:end),1)]);
    directiond1r_all_taste_aligned_time(j,:) = mean([mean(group1(group1a(j)).Sucrose.frames_aligned(group1(group1a(j)).Sucrose.correct==1,start_index:end),1);...
        mean(group1(group1a(j)).NaCl.frames_aligned(group1(group1a(j)).NaCl.correct==1,start_index:end),1)]);
%     d1r_all_trials(j,:) = mean([directionA1(j,:) ; directionB1(j,:)]);

      d1r_mean_activity.all_trials_drylick_aligned_activity(j,:) = group1(group1a(j)).all_trials.lick_aligned_traces_avg;
      d1r_mean_activity.all_trials_drylick_aligned_time(j,:) = mean(group1(group1a(j)).all_trials.lick_aligned_frames_matrix,1);
%        d1r_mean_activity.laterallicks.sucrose(j,:) = mean(group1(group1a(j)).Sucrose.right_licks(:,1)');
%     d1r_mean_activity.laterallicks.nacl(j,:) = mean(group1(group1a(j)).NaCl.left_licks(:,1)');
    clear correct_trials_L incorrect_trials_L correct_trials_R incorrect_trials_R correct_trials_suc incorrect_trials_suc correct_trials_nacl incorrect_trials_nacl
end

d1r_mean_activity.all_trials_lateral_aligned =  directiond1r_all;
d1r_mean_activity.all_trials_lateral_aligned_time =  directiond1r_all_time;


d1r_mean_activity.all_trials_taste_aligned =  directiond1r_all_taste_aligned;
d1r_mean_activity.all_trials_taste_aligned_time =  directiond1r_all_taste_aligned_time;

d1r_mean_activity.left_mean_resp = directionA11;
d1r_mean_activity.right_mean_resp =  directionB11;
d1r_mean_activity.left_mean_time = time_directionA1;
d1r_mean_activity.right_mean_time =  time_directionB1;
d1r_mean_activity.direction_e_right = direction_e_right;
d1r_mean_activity.direction_e_left = direction_e_left;


d1r_mean_activity.left_mean_resp_not = directionA11_not;
d1r_mean_activity.right_mean_resp_not =  directionB11_not;
d1r_mean_activity.left_mean_time_not = time_directionA1_not;
d1r_mean_activity.right_mean_time_not =  time_directionB1_not;

d1r_mean_activity.correct_trials_suc = taste_correct_trials_suc;
d1r_mean_activity.correct_trials_nacl =  taste_correct_trials_nacl;
d1r_mean_activity.incorrect_trials_suc = taste_incorrect_trials_suc;
d1r_mean_activity.incorrect_trials_nacl = taste_incorrect_trials_nacl;

d1r_mean_activity.correct_trials_suc_time = taste_correct_trials_suc_time;
d1r_mean_activity.correct_trials_nacl_time =  taste_correct_trials_nacl_time;
d1r_mean_activity.incorrect_trials_suc_time = taste_incorrect_trials_suc_time;
d1r_mean_activity.incorrect_trials_nacl_time = taste_incorrect_trials_nacl_time;


%%
clearvars -except group1a group2a group1 group2 y_limits d1r_mean_activity number_of_error_trials correct_trials_R_random correct_trials_R_random randomized

start_index = 1;
clims = [0 1];
index_rank = 23;
threshold = 0;
x_limits = [-3 2];
% unresp1 = setdiff(A,group2a);
for j =1:size(group2a,2)
    correct_trials_L = find(group2(group2a(j)).Left.correct==1);
    incorrect_trials_L = find(group2(group2a(j)).Left.incorrect==1);
    correct_trials_R = find(group2(group2a(j)).Right.correct==1);
    incorrect_trials_R = find(group2(group2a(j)).Right.incorrect==1);
    correct_trials_suc = find(group2(group2a(j)).Sucrose.correct==1);
    incorrect_trials_suc = find(group2(group2a(j)).Sucrose.incorrect==1);
    correct_trials_nacl = find(group2(group2a(j)).NaCl.correct==1);
    incorrect_trials_nacl = find(group2(group2a(j)).NaCl.incorrect==1);
    
    if randomized == 1
        correct_trials_L_random = randsample(size(correct_trials_L,1),round(size(correct_trials_L,1)/2),true);
        
        correct_trials_R_random = randsample(size(correct_trials_R,1),round(size(correct_trials_R,1)/2),true);
        
        correct_trials_suc_random = randsample(size(correct_trials_suc,1),round(size(correct_trials_suc,1)/2),true);
        
        correct_trials_nacl_random = randsample(size(correct_trials_nacl,1),round(size(correct_trials_nacl,1)/2),true);
    else
        correct_trials_L_random = correct_trials_L;
        
        correct_trials_R_random = correct_trials_R;
        
        correct_trials_suc_random = correct_trials_suc;
        
        correct_trials_nacl_random = correct_trials_nacl;
    end
   
     
    
   
    directionA21(j,:) = mean(group2(group2a(j)).Left.dff_aligned(correct_trials_L_random,start_index:end));
    time_directionA2(j,:) = mean(group2(group2a(j)).Left.frames_aligned(correct_trials_L_random,start_index:end));
    directionB21(j,:) = mean(group2(group2a(j)).Right.dff_aligned(correct_trials_R_random,start_index:end));
    time_directionB2(j,:) = mean(group2(group2a(j)).Right.frames_aligned(correct_trials_R_random,start_index:end));
    
    directionA21_not(j,:) = mean(group2(group2a(j)).Left.dff_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
    time_directionA2_not(j,:) = mean(group2(group2a(j)).Left.frames_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
    directionB21_not(j,:) = mean(group2(group2a(j)).Right.dff_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));
    time_directionB2_not(j,:) = mean(group2(group2a(j)).Right.frames_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));

    if size(incorrect_trials_L,1) >= number_of_error_trials && size(incorrect_trials_R,1) >= number_of_error_trials
        direction_e_left(j,:) = mean(group2(group2a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end),1);
        direction_e_right(j,:) = mean(group2(group2a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end),1); %%actual right licks
    else
        direction_e_right(j,1:82) = NaN;
        direction_e_left(j,1:82) = NaN;
    end
%     directionA21(j,:) = mean([group2(group2a(j)).Left.dff_aligned(correct_trials_L,start_index:end);group2(group2a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end)]);
    directionA2(j,:) = directionA21(j,:)./max(directionA21(j,:));
    
%     time_directionA2(j,:) = mean([group2(group2a(j)).Left.frames_aligned(correct_trials_L,start_index:end);group2(group2a(j)).Right.frames_aligned(incorrect_trials_R,start_index:end)]);
    

%     directionB21(j,:) = mean([group2(group2a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end);group2(group2a(j)).Right.dff_aligned(correct_trials_R,start_index:end)]);
    directionB2(j,:) = directionB21(j,:)./max(directionB21(j,:));
%         time_directionB2(j,:) = mean([group2(group2a(j)).Right.frames_aligned(sampled_trials,start_index:end);group2(group2a(j)).Left.frames_aligned(incorrect_trials,start_index:end)]);
%     directionnon_all(j,:) = mean([group2(group2a(j)).Left.dff_aligned(group2(group2a(j)).Left.correct==1,start_index:end);...
%         group2(group2a(j)).Left.dff_aligned(group2(group2a(j)).Left.incorrect==1,start_index:end);...
%         group2(group2a(j)).Right.dff_aligned(group2(group2a(j)).Right.correct==1,start_index:end);...
%         group2(group2a(j)).Right.dff_aligned(group2(group2a(j)).Right.incorrect==1,start_index:end)]);
    
    directionnon_all(j,:) = mean([mean(group2(group2a(j)).Left.dff_aligned(group2(group2a(j)).Left.correct==1,start_index:end),1);...
        mean(group2(group2a(j)).Right.dff_aligned(group2(group2a(j)).Right.correct==1,start_index:end),1)]);
    
     directionnon_all_time(j,:) = mean([mean(group2(group2a(j)).Left.frames_aligned(group2(group2a(j)).Left.correct==1,start_index:end),1);...
        mean(group2(group2a(j)).Right.frames_aligned(group2(group2a(j)).Right.correct==1,start_index:end),1)]);
        
    directionnon_all_taste_aligned(j,:) = mean([mean(group2(group2a(j)).Sucrose.dff_aligned(group2(group2a(j)).Sucrose.correct==1,start_index:end),1);...
        mean(group2(group2a(j)).NaCl.dff_aligned(group2(group2a(j)).NaCl.correct==1,start_index:end),1)]);
    
    directionnon_all_taste_aligned_time(j,:) = mean([mean(group2(group2a(j)).Sucrose.frames_aligned(group2(group2a(j)).Sucrose.correct==1,start_index:end),1);...
        mean(group2(group2a(j)).NaCl.frames_aligned(group2(group2a(j)).NaCl.correct==1,start_index:end),1)]);
%     non_all_trials(j,:) = mean([group2(group2a(j)).Left.dff_aligned(group2(group2a(j)).Left.correct ==1,start_index:end);group2(group2a(j)).Right.dff_aligned(group2(group2a(j)).Right.correct ==1,start_index:end)]);
%     non_all_trials(j,:) = mean([directionA2(j,:); directionB2(j,:)]);
    all_first_central_lick_time(j) = mean([group2(group2a(j)).Left.central_licks(:,1); group2(group2a(j)).Right.central_licks(:,1)]);
    all_first_wetcentral_lick_time(j) = mean(nonzeros([group2(group2a(j)).Left.central_licks(:,end-1:end); group2(group2a(j)).Right.central_licks(:,end-1:end)]));
    
    taste_correct_trials_suc(j,:) = mean(group2(group2a(j)).Sucrose.dff_aligned(correct_trials_suc_random,start_index:end));
    taste_correct_trials_nacl(j,:) = mean(group2(group2a(j)).NaCl.dff_aligned(correct_trials_nacl_random,start_index:end));
    taste_incorrect_trials_suc(j,:) = mean(group2(group2a(j)).Sucrose.dff_aligned(incorrect_trials_suc,start_index:end));
    taste_incorrect_trials_nacl(j,:) = mean(group2(group2a(j)).NaCl.dff_aligned(incorrect_trials_nacl,start_index:end));
    
    taste_correct_trials_suc_time(j,:) = mean(group2(group2a(j)).Sucrose.frames_aligned(correct_trials_suc_random,start_index:end));
    taste_correct_trials_nacl_time(j,:) = mean(group2(group2a(j)).NaCl.frames_aligned(correct_trials_nacl_random,start_index:end));
    taste_incorrect_trials_suc_time(j,:) = mean(group2(group2a(j)).Sucrose.frames_aligned(incorrect_trials_suc,start_index:end));
    taste_incorrect_trials_nacl_time(j,:) = mean(group2(group2a(j)).NaCl.frames_aligned(incorrect_trials_nacl,start_index:end));

    non_mean_activity.all_trials_drylick_aligned_activity(j,:) = group2(group2a(j)).all_trials.lick_aligned_traces_avg;
      non_mean_activity.all_trials_drylick_aligned_time(j,:) = mean(group2(group2a(j)).all_trials.lick_aligned_frames_matrix,1);
%     non_mean_activity.laterallicks.sucrose(j,:) = mean(group2(group2a(j)).Sucrose.right_licks(:,1)');
%     non_mean_activity.laterallicks.nacl(j,:) = mean(group2(group2a(j)).NaCl.left_licks(:,1)');
    
    clear correct_trials_L incorrect_trials_L correct_trials_R incorrect_trials_R correct_trials_suc incorrect_trials_suc correct_trials_nacl incorrect_trials_nacl
end 
% clear r I t
% for t = 1:size(non_all_trials,1)
%     r(t) = length(find(non_all_trials(t,index_rank:end)>threshold));
% end
% [~,I]=sort(r);
% directionA2 = directionA2(flipud(I'),:);
% directionB2 = directionB2(flipud(I'),:);
% 
% ax(1) = subplot(2,3,4);
% imagesc(mean(time_directionA2),1:size(directionA2,1),directionA2,clims)
% line([0 0],[0 size(time_directionA2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% line([mean(all_first_central_lick_time) mean(all_first_central_lick_time)],[0 size(directionA2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% line([mean(all_first_wetcentral_lick_time) mean(all_first_wetcentral_lick_time)],[0 size(directionA2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% colorbar('eastoutside'); 
% xlim(x_limits)
% ylabel('D1R- neurons');
% yticks([1 size(directionA2,1)])
% box off
% % axis square
% colormap(ax(1),J);
% ax(1) = subplot(2,3,5);
% % imagesc(mean(time_directionB2),1:size(directionB2,1),directionB2,clims)
% line([0 0],[0 size(directionB2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% line([mean(all_first_central_lick_time) mean(all_first_central_lick_time)],[0 size(directionB2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% line([mean(all_first_wetcentral_lick_time) mean(all_first_wetcentral_lick_time)],[0 size(directionB2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% xlim(x_limits)
% colorbar('eastoutside'); colormap(ax(1),J);
% yticks([1 size(directionB2,1)])
% xlabel('Time from first lateral lick (s)');
% box off
% % axis square
% 
% clear sem1 sem2
% ax(1) = subplot(2,3,6);
% sem1 = std(directionB21)/sqrt(size(directionB21,1));
% sem2 = std(directionA21)/sqrt(size(directionA21,1));
% line([0 0],[-0.5 7],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% line([mean(all_first_central_lick_time) mean(all_first_central_lick_time)],[-0.5 7],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% line([mean(all_first_wetcentral_lick_time) mean(all_first_wetcentral_lick_time)],[-0.5 7],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% 
% boundedline(mean(time_directionA2),mean(directionA21),sem2,'alpha','r','LineWidth',2.5); hold on;
% boundedline(mean(time_directionB2),mean(directionB21),sem1,'alpha','b','LineWidth',2.5); 
% title('D1R- population PSTH')
% ylabel('z \DeltaF/F');xlabel('Time from first lateral lick (s)');
% ylim(y_limits);xlim(x_limits)
% axis square
% box off
% axis square
% set(gcf,'Position',[544 319 1103 566]);
non_mean_activity.all_trials_lateral_aligned =  directionnon_all;
non_mean_activity.all_trials_lateral_aligned_time =  directionnon_all_time;

non_mean_activity.all_trials_taste_aligned =  directionnon_all_taste_aligned;
non_mean_activity.all_trials_taste_aligned_time =  directionnon_all_taste_aligned_time;

non_mean_activity.left_mean_resp = directionA21;
non_mean_activity.right_mean_resp =  directionB21;
non_mean_activity.left_mean_time = time_directionA2;
non_mean_activity.right_mean_time =  time_directionB2;
non_mean_activity.direction_e_right = direction_e_right;
non_mean_activity.direction_e_left = direction_e_left;


non_mean_activity.left_mean_resp_not = directionA21_not;
non_mean_activity.right_mean_resp_not =  directionB21_not;
non_mean_activity.left_mean_time_not = time_directionA2_not;
non_mean_activity.right_mean_time_not =  time_directionB2_not;

non_mean_activity.correct_trials_suc = taste_correct_trials_suc;
non_mean_activity.correct_trials_nacl =  taste_correct_trials_nacl;
non_mean_activity.incorrect_trials_suc = taste_incorrect_trials_suc;
non_mean_activity.incorrect_trials_nacl = taste_incorrect_trials_nacl;


non_mean_activity.correct_trials_suc_time = taste_correct_trials_suc_time;
non_mean_activity.correct_trials_nacl_time =  taste_correct_trials_nacl_time;
non_mean_activity.incorrect_trials_suc_time = taste_incorrect_trials_suc_time;
non_mean_activity.incorrect_trials_nacl_time = taste_incorrect_trials_nacl_time;


%%
% clims = [-2 2];

% axis square

% clear directionA2 time_directionA2 directionB2 time_directionB2 non_all_trials
% for j =1:size(unresp1,2)
%     
%     directionA2(j,:) = mean(group2(unresp1(j)).Left.dff_aligned(group2(unresp1(j)).Left.correct ==1,start_index:end));
%     time_directionA2(j,:) = mean(group2(unresp1(j)).Left.frames_aligned(group2(unresp1(j)).Left.correct ==1,start_index:end));
%     directionB2(j,:) = mean(group2(unresp1(j)).Right.dff_aligned(group2(unresp1(j)).Right.correct ==1,start_index:end));
%     time_directionB2(j,:) = mean(group2(unresp1(j)).Right.frames_aligned(group2(unresp1(j)).Right.correct ==1,start_index:end));
%     non_all_trials(j,:) = mean([group2(unresp1(j)).Left.dff_aligned(group2(unresp1(j)).Left.correct ==1,start_index:end);group2(unresp1(j)).Right.dff_aligned(group2(unresp1(j)).Right.correct ==1,start_index:end)]);
% end
% clear r I t
% for t = 1:size(non_all_trials,1)
%     r(t) = length(find(non_all_trials(t,index_rank:end)>threshold));
% end
% [~,I]=sort(r);
% directionA2 = directionA2(flipud(I'),:);
% directionB2 = directionB2(flipud(I'),:);
% 
% ax(1) = subplot(2,4,7);
% imagesc(mean(time_directionA2),1:size(directionA2,1),directionA2,clims)
% line([0 0],[0 size(directionB2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% colorbar('eastoutside'); 
% xlim(x_limits)
% ylabel('D1R- neurons');
% yticks([1 size(time_directionA2,1)])
% 
% % J = customcolormap([0 0.2 1], [255,255,204;  253,141,60; 50,0,13]./255);
% colormap(ax(1),J);
% ax(1) = subplot(2,4,8);
% imagesc(mean(time_directionB2),1:size(directionB2,1),directionB2,clims)
% % J = customcolormap([0 0.2 1], [255,255,204;  253,141,60; 50,0,13]./255);
% colorbar('eastoutside'); colormap(ax(1),J);
% line([0 0],[0 size(directionB2,1)],'color',[0.7 0.7 0.7],'LineWidth',2); hold on
% xlim(x_limits)
% yticks([1 size(directionB2,1)])


end