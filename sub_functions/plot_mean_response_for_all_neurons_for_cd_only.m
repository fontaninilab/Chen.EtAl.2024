function [d1r_mean_activity, non_mean_activity] = plot_mean_response_for_all_neurons_for_cd_only(group1a,group2a,group1, group2,number_of_error_trials);
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
    
%     correct_trials_L_random = randsample(size(correct_trials_L,1),round(size(correct_trials_L,1)/2),true);
    
%     correct_trials_R_random = randsample(size(correct_trials_R,1),round(size(correct_trials_R,1)/2),true);
    
%     correct_trials_suc_random = randsample(size(correct_trials_suc,1),round(size(correct_trials_suc,1)/2),true);
    
%     correct_trials_nacl_random = randsample(size(correct_trials_nacl,1),round(size(correct_trials_nacl,1)/2),true);
    
    directionA11(j,:) = mean([group1(group1a(j)).Left.dff_aligned(correct_trials_L,start_index:end);group1(group1a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end)]);
    time_directionA1(j,:) = mean(group1(group1a(j)).Left.frames_aligned(correct_trials_L,start_index:end));
    
    directionB11(j,:) = mean([group1(group1a(j)).Right.dff_aligned(correct_trials_R,start_index:end);group1(group1a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end)]);
    time_directionB1(j,:) = mean(group1(group1a(j)).Right.frames_aligned(correct_trials_R,start_index:end));
    
%     directionA11_not(j,:) = mean(group1(group1a(j)).Left.dff_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
%     time_directionA1_not(j,:) = mean(group1(group1a(j)).Left.frames_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
    directionA11_not(j,:) = directionA11(j,:);
    time_directionA1_not(j,:) = time_directionA1(j,:);
    
%     directionB11_not(j,:) = mean(group1(group1a(j)).Right.dff_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));
%     time_directionB1_not(j,:) = mean(group1(group1a(j)).Right.frames_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));
    directionB11_not(j,:) = directionB11(j,:);
    time_directionB1_not(j,:) = time_directionB1(j,:);
    
    
%     if size(incorrect_trials_L,1) >= number_of_error_trials && size(incorrect_trials_R,1) >= number_of_error_trials
%         direction_e_right(j,:) = mean(group1(group1a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end),1);
%         direction_e_left(j,:) = mean(group1(group1a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end),1);
%     else
%         direction_e_right(j,1:82) = NaN;
%         direction_e_left(j,1:82) = NaN;
%     end
    taste_correct_trials_suc(j,:) = mean([group1(group1a(j)).Sucrose.dff_aligned(correct_trials_suc,start_index:end);group1(group1a(j)).Sucrose.dff_aligned(incorrect_trials_suc,start_index:end)]);
    taste_correct_trials_nacl(j,:) = mean([group1(group1a(j)).NaCl.dff_aligned(correct_trials_nacl,start_index:end);group1(group1a(j)).NaCl.dff_aligned(incorrect_trials_nacl,start_index:end)]);
    
%     taste_correct_trials_suc_not(j,:) = mean(group1(group1a(j)).Sucrose.dff_aligned(setdiff(correct_trials_suc,correct_trials_suc_random),start_index:end));
%     taste_correct_trials_nacl_not(j,:) = mean(group1(group1a(j)).NaCl.dff_aligned(setdiff(correct_trials_nacl,correct_trials_nacl_random),start_index:end));
    taste_correct_trials_suc_not(j,:) = taste_correct_trials_suc(j,:);
    taste_correct_trials_nacl_not(j,:) = taste_correct_trials_nacl(j,:);
    
%     taste_incorrect_trials_suc(j,:) = mean(group1(group1a(j)).Sucrose.dff_aligned(incorrect_trials_suc,start_index:end));
% %     taste_incorrect_trials_nacl(j,:) = mean(group1(group1a(j)).NaCl.dff_aligned(incorrect_trials_nacl,start_index:end));
    
    taste_correct_trials_suc_time(j,:) = mean(group1(group1a(j)).Sucrose.frames_aligned(correct_trials_suc,start_index:end));
    taste_correct_trials_nacl_time(j,:) = mean(group1(group1a(j)).NaCl.frames_aligned(correct_trials_nacl,start_index:end));
    
%     taste_correct_trials_suc_time_not(j,:) = mean(group1(group1a(j)).Sucrose.frames_aligned(setdiff(correct_trials_suc,correct_trials_suc_random),start_index:end));
%     taste_correct_trials_nacl_time_not(j,:) = mean(group1(group1a(j)).NaCl.frames_aligned(setdiff(correct_trials_nacl,correct_trials_nacl_random),start_index:end));
    
    taste_correct_trials_suc_time_not(j,:) = taste_correct_trials_suc_time(j,:);
    taste_correct_trials_nacl_time_not(j,:) = taste_correct_trials_nacl_time(j,:);
    
    
%     taste_incorrect_trials_suc_time(j,:) = mean(group1(group1a(j)).Sucrose.frames_aligned(incorrect_trials_suc,start_index:end));
%     taste_incorrect_trials_nacl_time(j,:) = mean(group1(group1a(j)).NaCl.frames_aligned(incorrect_trials_nacl,start_index:end));



    clear correct_trials_L incorrect_trials_L correct_trials_R incorrect_trials_R correct_trials_suc incorrect_trials_suc correct_trials_nacl incorrect_trials_nacl
end

% d1r_mean_activity.all_trials =  directiond1r_all;
d1r_mean_activity.left_mean_resp = directionA11;
d1r_mean_activity.right_mean_resp =  directionB11;
d1r_mean_activity.left_mean_time = time_directionA1;
d1r_mean_activity.right_mean_time =  time_directionB1;
% d1r_mean_activity.direction_e_right = direction_e_right;
% d1r_mean_activity.direction_e_left = direction_e_left;


d1r_mean_activity.left_mean_resp_not = directionA11_not;
d1r_mean_activity.right_mean_resp_not =  directionB11_not;
d1r_mean_activity.left_mean_time_not = time_directionA1_not;
d1r_mean_activity.right_mean_time_not =  time_directionB1_not;

d1r_mean_activity.correct_trials_suc = taste_correct_trials_suc;
d1r_mean_activity.correct_trials_nacl =  taste_correct_trials_nacl;
% d1r_mean_activity.incorrect_trials_suc = taste_incorrect_trials_suc;
% d1r_mean_activity.incorrect_trials_nacl = taste_incorrect_trials_nacl;

d1r_mean_activity.correct_trials_suc_time = taste_correct_trials_suc_time;
d1r_mean_activity.correct_trials_nacl_time =  taste_correct_trials_nacl_time;
% d1r_mean_activity.incorrect_trials_suc_time = taste_incorrect_trials_suc_time;
% d1r_mean_activity.incorrect_trials_nacl_time = taste_incorrect_trials_nacl_time;

d1r_mean_activity.correct_trials_suc_not = taste_correct_trials_suc_not;
d1r_mean_activity.correct_trials_nacl_not =  taste_correct_trials_nacl_not;
d1r_mean_activity.correct_trials_suc_time_not = taste_correct_trials_suc_time_not;
d1r_mean_activity.correct_trials_nacl_time_not =  taste_correct_trials_nacl_time_not;
%%
clearvars -except group1a group2a group1 group2 y_limits d1r_mean_activity number_of_error_trials correct_trials_R_random correct_trials_R_random

start_index = 1;
% clims = [0 1];
% index_rank = 23;
% threshold = 0;
% x_limits = [-3 2];
% unresp1 = setdiff(A,group2a);
for j =1:size(group2a,2)
    
    correct_trials_L = find(group2(group2a(j)).Left.correct==1);
    incorrect_trials_L = find(group2(group2a(j)).Left.incorrect==1);
    
    correct_trials_R = find(group2(group2a(j)).Right.correct==1);
    incorrect_trials_R = find(group2(group2a(j)).Right.incorrect==1);
    
%     correct_trials_L_random = randsample(size(correct_trials_L,1),round(size(correct_trials_L,1)/2),true);
%     correct_trials_R_random = randsample(size(correct_trials_R,1),round(size(correct_trials_R,1)/2),true);
%     
%     correct_trials_suc_random = randsample(size(correct_trials_suc,1),round(size(correct_trials_suc,1)/2),true); 
%     correct_trials_nacl_random = randsample(size(correct_trials_nacl,1),round(size(correct_trials_nacl,1)/2),true);
  
    correct_trials_suc = find(group2(group2a(j)).Sucrose.correct==1);
    incorrect_trials_suc = find(group2(group2a(j)).Sucrose.incorrect==1);
    correct_trials_nacl = find(group2(group2a(j)).NaCl.correct==1);
    incorrect_trials_nacl = find(group2(group2a(j)).NaCl.incorrect==1);
    
   
    directionA21(j,:) = mean([group2(group2a(j)).Left.dff_aligned(correct_trials_L,start_index:end);group2(group2a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end)]);
    time_directionA2(j,:) = mean(group2(group2a(j)).Left.frames_aligned(correct_trials_L,start_index:end));
    directionB21(j,:) = mean([group2(group2a(j)).Right.dff_aligned(correct_trials_R,start_index:end);group2(group2a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end)]);
    time_directionB2(j,:) = mean(group2(group2a(j)).Right.frames_aligned(correct_trials_R,start_index:end));
    
%     directionA21_not(j,:) = mean(group2(group2a(j)).Left.dff_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
%     time_directionA2_not(j,:) = mean(group2(group2a(j)).Left.frames_aligned(setdiff(correct_trials_L,correct_trials_L_random),start_index:end));
%     directionB21_not(j,:) = mean(group2(group2a(j)).Right.dff_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));
%     time_directionB2_not(j,:) = mean(group2(group2a(j)).Right.frames_aligned(setdiff(correct_trials_R,correct_trials_R_random),start_index:end));

    directionA21_not(j,:) = directionA21(j,:);
    time_directionA2_not(j,:) = time_directionA2(j,:);
    directionB21_not(j,:) = directionB21(j,:);
    time_directionB2_not(j,:) = time_directionB2(j,:);


%     if size(incorrect_trials_L,1) >= number_of_error_trials && size(incorrect_trials_R,1) >= number_of_error_trials
%         direction_e_right(j,:) = mean(group2(group2a(j)).Right.dff_aligned(incorrect_trials_R,start_index:end),1);
%         direction_e_left(j,:) = mean(group2(group2a(j)).Left.dff_aligned(incorrect_trials_L,start_index:end),1);
%     else
%         direction_e_right(j,1:82) = NaN;
%         direction_e_left(j,1:82) = NaN;
%     end
%  
    
    taste_correct_trials_suc(j,:) = mean([group2(group2a(j)).Sucrose.dff_aligned(correct_trials_suc,start_index:end);group2(group2a(j)).Sucrose.dff_aligned(incorrect_trials_suc,start_index:end)]);
    taste_correct_trials_nacl(j,:) = mean([group2(group2a(j)).NaCl.dff_aligned(correct_trials_nacl,start_index:end);group2(group2a(j)).NaCl.dff_aligned(incorrect_trials_nacl,start_index:end)]);
%     taste_incorrect_trials_suc(j,:) = mean(group2(group2a(j)).Sucrose.dff_aligned(incorrect_trials_suc,start_index:end));
%     taste_incorrect_trials_nacl(j,:) = mean(group2(group2a(j)).NaCl.dff_aligned(incorrect_trials_nacl,start_index:end));
    
    taste_correct_trials_suc_time(j,:) = mean(group2(group2a(j)).Sucrose.frames_aligned(correct_trials_suc,start_index:end));
    taste_correct_trials_nacl_time(j,:) = mean(group2(group2a(j)).NaCl.frames_aligned(correct_trials_nacl,start_index:end));
%     taste_incorrect_trials_suc_time(j,:) = mean(group2(group2a(j)).Sucrose.frames_aligned(incorrect_trials_suc,start_index:end));
%     taste_incorrect_trials_nacl_time(j,:) = mean(group2(group2a(j)).NaCl.frames_aligned(incorrect_trials_nacl,start_index:end));

%     taste_correct_trials_suc_not(j,:) = mean(group2(group2a(j)).Sucrose.dff_aligned(setdiff(correct_trials_suc,correct_trials_suc_random),start_index:end));
%     taste_correct_trials_nacl_not(j,:) = mean(group2(group2a(j)).NaCl.dff_aligned(setdiff(correct_trials_nacl,correct_trials_nacl_random),start_index:end));
%     
%     taste_correct_trials_suc_time_not(j,:) = mean(group2(group2a(j)).Sucrose.frames_aligned(setdiff(correct_trials_suc,correct_trials_suc_random),start_index:end));
%     taste_correct_trials_nacl_time_not(j,:) = mean(group2(group2a(j)).NaCl.frames_aligned(setdiff(correct_trials_nacl,correct_trials_nacl_random),start_index:end));

    taste_correct_trials_suc_not(j,:) = taste_correct_trials_suc(j,:);
    taste_correct_trials_nacl_not(j,:) = taste_correct_trials_nacl(j,:);
    
    taste_correct_trials_suc_time_not(j,:) = taste_correct_trials_suc_time(j,:);
    taste_correct_trials_nacl_time_not(j,:) = taste_correct_trials_nacl_time(j,:);
    clear correct_trials_L incorrect_trials_L correct_trials_R incorrect_trials_R correct_trials_suc incorrect_trials_suc correct_trials_nacl incorrect_trials_nacl
end 
% non_mean_activity.all_trials =  directionnon_all;
non_mean_activity.left_mean_resp = directionA21;
non_mean_activity.right_mean_resp =  directionB21;
non_mean_activity.left_mean_time = time_directionA2;
non_mean_activity.right_mean_time =  time_directionB2;
% non_mean_activity.direction_e_right = direction_e_right;
% non_mean_activity.direction_e_left = direction_e_left;


non_mean_activity.left_mean_resp_not = directionA21_not;
non_mean_activity.right_mean_resp_not =  directionB21_not;
non_mean_activity.left_mean_time_not = time_directionA2_not;
non_mean_activity.right_mean_time_not =  time_directionB2_not;

non_mean_activity.correct_trials_suc = taste_correct_trials_suc;
non_mean_activity.correct_trials_nacl =  taste_correct_trials_nacl;
% non_mean_activity.incorrect_trials_suc = taste_incorrect_trials_suc;
% non_mean_activity.incorrect_trials_nacl = taste_incorrect_trials_nacl;


non_mean_activity.correct_trials_suc_time = taste_correct_trials_suc_time;
non_mean_activity.correct_trials_nacl_time =  taste_correct_trials_nacl_time;
% non_mean_activity.incorrect_trials_suc_time = taste_incorrect_trials_suc_time;
% non_mean_activity.incorrect_trials_nacl_time = taste_incorrect_trials_nacl_time;

non_mean_activity.correct_trials_suc_not = taste_correct_trials_suc_not;
non_mean_activity.correct_trials_nacl_not =  taste_correct_trials_nacl_not;
non_mean_activity.correct_trials_suc_time_not = taste_correct_trials_suc_time_not;
non_mean_activity.correct_trials_nacl_time_not =  taste_correct_trials_nacl_time_not;

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