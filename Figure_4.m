%% Figure 4 plots -  Population coding
%%
clearvars -except d1r_neurons non_neurons d1r_resp non_resp
include_error_trials = 1;
figure_startup
cd 'D:\John_data_figure_materials\Fig4'
% close all;
clc;
disp('Running, please wait')
group1 = d1r_neurons;
group2 = non_neurons;
group1a = 1:length(d1r_neurons);
group2a = 1:length(non_neurons);
y_limits = [-2 6];
number_of_error_trials = 3;  
random = 0;
figure;clf

[d1r_mean_activity, non_mean_activity] = plot_mean_response_for_all_neurons(group1a,group2a,group1, group2,number_of_error_trials,random);
d1r_mean_activity_copy = d1r_mean_activity;
non_mean_activity_copy = non_mean_activity;


if include_error_trials == 1
    d1r_mean_activity_copy.direction_e_right(isnan(d1r_mean_activity_copy.direction_e_right))=0;
    non_mean_activity_copy.direction_e_right(isnan(non_mean_activity_copy.direction_e_right))=0;
    d1r_mean_activity_copy.direction_e_left(isnan(d1r_mean_activity_copy.direction_e_left))=0;
    non_mean_activity_copy.direction_e_left(isnan(non_mean_activity_copy.direction_e_left))=0;
    
    clearvars r r1
    
    r = find(any(~d1r_mean_activity_copy.direction_e_right,2));
    r1 = find(any(~d1r_mean_activity_copy.direction_e_left,2));
    no_error_d1r_neurons = unique([r;r1],'sorted');
    
    clearvars r r1
    
    r = find(any(~non_mean_activity_copy.direction_e_right,2));
    r1 = find(any(~non_mean_activity_copy.direction_e_left,2));
    no_error_non_neurons = unique([r;r1],'sorted');
    f1 = fieldnames(non_mean_activity_copy);
    
    for i = 1:size(f1,1)
        non_mean_activity_copy.(f1{i})(no_error_non_neurons,:) = [];
        d1r_mean_activity_copy.(f1{i})(no_error_d1r_neurons,:) = [];
    end
end
    close all
    disp('Finished averaging responses')

%%
% clearvars -except d1r_neurons non_neurons d1r_resp non_resp
disp('Running, please wait')

clearvars d1r_mean_activity non_mean_activity
[d1r_mean_activity, non_mean_activity] = plot_mean_response_for_all_neurons_for_cd_only(group1a,group2a,group1, group2,number_of_error_trials);
disp('Finished averaging responses')
    %%
    clearvars d1r_selec non_selec p_value_time perm_selectivity_diff true_selectivity_diff
figure(1000);clf;
clc;
disp('Running')
delay_window = 48:80; 


rng('default')
clearvars diff_selec_perm diff_selec_true 
BW = 0.01;
 
    d1r_selec = d1r_mean_activity.right_mean_resp - d1r_mean_activity.left_mean_resp;
    non_selec = non_mean_activity.right_mean_resp - non_mean_activity.left_mean_resp;


num_of_permutations = 1000;
for i = 1:num_of_permutations
  

 
all_selectivity = [d1r_selec;non_selec];
        y = randperm(length(all_selectivity));

    diff_selec_perm(i,:) = mean(mean(all_selectivity(y(1:size(d1r_selec,1)),delay_window),2)) - mean(mean(all_selectivity(y(size(d1r_selec,1)+1:end),delay_window),2));
    
    diff_selec_true(i,:) = mean(mean(all_selectivity(randsample((1:size(d1r_selec,1)),size(d1r_selec,1),true),delay_window),2)) - mean(mean(all_selectivity(randsample((size(d1r_selec,1)+1:size(non_selec,1)),size(non_selec,1),true),delay_window),2));

end

%
p_observed = (size(find(diff_selec_perm>=(mean(mean(d1r_selec(:,delay_window),2))-mean(mean(non_selec(:,delay_window),2)))),1)+1)/(num_of_permutations+1);

histogram(diff_selec_true,'normalization','probability','BinWidth',BW,'EdgeColor','none','FaceColor',[0.2 0.2 0.2]);
hold on;

histogram(diff_selec_perm,'normalization','probability','BinWidth',BW,'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);

set(gca,'linewidth',1)

hold on;
line([mean(mean(d1r_selec(:,delay_window),2))-mean(mean(non_selec(:,delay_window),2)) mean(mean(d1r_selec(:,delay_window),2))-mean(mean(non_selec(:,delay_window),2))],[0 1000],'Color',[0.2 0.2 0.2],'LineWidth',1);
line([mean(diff_selec_perm) mean(diff_selec_perm)],[0 1000],'Color',[0.7 0.7 0.7],'LineWidth',1);
clearvars y
y =diff_selec_perm;
yMean = mean(y);                                    
ySEM = std(y)/sqrt(1);                              
CI95 = tinv([0.05 0.95], size(diff_selec_perm,1)-1);             
yCI95 = bsxfun(@times, ySEM, CI95(:));             
% 
line([yMean + yCI95(1) yMean + yCI95(1)],[0 100],'Color',[0.7 0.7 0.7],'LineStyle','--');
line([yMean + yCI95(2) yMean + yCI95(2)],[0 100],'Color',[0.7 0.7 0.7],'LineStyle','--');
legend('True','Shuffled','Location','northwest')

text(0.15,0.09,strcat('p =  ', num2str(p_observed)))
box off
ylim([0 0.1]); 
ylabel('Probability')

set(gcf, 'Renderer', 'painters');
set(gcf,'Position', [-1086 652 606 257]);


%%
number_of_error_trials = 3;  
rng('default')
clc;
group1 = d1r_neurons;
group2 = non_neurons;
group1a = 1:length(d1r_neurons);
group2a = 1:length(non_neurons);
clearvars d1r_mean_activity non_mean_activity
[d1r_mean_activity, non_mean_activity] = plot_mean_response_for_all_neurons_for_cd_only(group1a,group2a,group1, group2,number_of_error_trials);

%
CloseVariableWindows

number_of_resamples = 1000;

ramp_window_central_base = 10:15;
ramp_window_central = 15:35;

sample_window_base = 15:35;
delay_window_base = 35:42;

sample_window = 38:47;
delay_window = 57:68;
outcome_window_base =57:68;
choice_window = 69:80;
outcome_window = 69:82;

ramp_window_base = 15:35;
ramp_window = 65:68;

response_window = 70:82;
projection_multiplier = 1;
clearvars stimulus_projections direction_projections outcome_projections all_left_times_not all_right_times_not
clearvars ramp_mode_central response_mode_central stimulus_mode direction_mode outcome_mode ramp_mode_lateral response_mode_lateral
clearvars selectivity_across_runs


    cell_type = {'d1r','non'};
disp('Calculating modes and projecting activity modes')

for k = 1:number_of_resamples
%     tic
    
    ran_num_d1r = randsample(size(d1r_mean_activity.left_mean_resp,1),round(size(d1r_mean_activity.left_mean_resp,1)/1),true);
    
    ran_num_non = randsample(size(non_mean_activity.left_mean_resp,1),round(size(d1r_mean_activity.left_mean_resp,1)/1),true);
    
    ran_num_d1r_outcome = randsample(size(d1r_mean_activity_copy.left_mean_resp,1),round(size(d1r_mean_activity_copy.left_mean_resp,1)/1),true);
    
    ran_num_non_outcome = randsample(size(non_mean_activity_copy.left_mean_resp,1),round(size(d1r_mean_activity_copy.left_mean_resp,1)/1),true);
    
    %
    clearvars  ramp_mode_central
    
    
    ramp_mode_central.d1r.combined_trials = (d1r_mean_activity.correct_trials_suc(ran_num_d1r,ramp_window_central)+d1r_mean_activity.correct_trials_nacl(ran_num_d1r,ramp_window_central))/2;
    ramp_mode_central.d1r.combined_trials_base = (d1r_mean_activity.correct_trials_suc(ran_num_d1r,ramp_window_central_base)+d1r_mean_activity.correct_trials_nacl(ran_num_d1r,ramp_window_central_base))/2;
    ramp_mode_central.d1r.selectivity= (mean(ramp_mode_central.d1r.combined_trials,2) - mean(ramp_mode_central.d1r.combined_trials_base,2))/2;
    
    ramp_mode_central.non.combined_trials = (non_mean_activity.correct_trials_suc(ran_num_non,ramp_window_central)+non_mean_activity.correct_trials_nacl(ran_num_non,ramp_window_central))/2;
    ramp_mode_central.non.combined_trials_base = (non_mean_activity.correct_trials_suc(ran_num_non,ramp_window_central_base)+non_mean_activity.correct_trials_nacl(ran_num_non,ramp_window_central_base))/2;
    ramp_mode_central.non.selectivity = (mean(ramp_mode_central.non.combined_trials,2) - mean(ramp_mode_central.non.combined_trials_base,2))/2;
    
    clearvars  response_mode_central
    
    
    response_mode_central.d1r.combined_trials = (d1r_mean_activity.correct_trials_suc(ran_num_d1r,sample_window)+d1r_mean_activity.correct_trials_nacl(ran_num_d1r,sample_window));
    response_mode_central.d1r.combined_trials_base = (d1r_mean_activity.correct_trials_suc(ran_num_d1r,sample_window_base)+d1r_mean_activity.correct_trials_nacl(ran_num_d1r,sample_window_base));
    response_mode_central.d1r.selectivity= (mean(response_mode_central.d1r.combined_trials,2) - mean(response_mode_central.d1r.combined_trials_base,2))/2;
    
    response_mode_central.non.combined_trials = (non_mean_activity.correct_trials_suc(ran_num_non,sample_window)+non_mean_activity.correct_trials_nacl(ran_num_non,sample_window));
    response_mode_central.non.combined_trials_base = (non_mean_activity.correct_trials_suc(ran_num_non,sample_window_base)+non_mean_activity.correct_trials_nacl(ran_num_non,sample_window_base));
    response_mode_central.non.selectivity = (mean(response_mode_central.non.combined_trials,2) - mean(response_mode_central.non.combined_trials_base,2))/2;
    
    
    
    clearvars  stimulus_mode
%     stimulus_mode.d1r.stim1 = (d1r_mean_activity_copy.correct_trials_suc + d1r_mean_activity_copy.incorrect_trials_suc);
%     stimulus_mode.d1r.stim2 = (d1r_mean_activity_copy.correct_trials_nacl + d1r_mean_activity_copy.incorrect_trials_nacl);
    stimulus_mode.d1r.stim1 = (d1r_mean_activity.correct_trials_suc(ran_num_d1r,:));
    stimulus_mode.d1r.stim2 = (d1r_mean_activity.correct_trials_nacl(ran_num_d1r,:));
    stimulus_mode.d1r.selectivity = (stimulus_mode.d1r.stim1 - stimulus_mode.d1r.stim2)/2;
    selectivity_across_runs.stimulus.d1r(k,:) = mean(stimulus_mode.d1r.selectivity);
    
    stimulus_mode.d1r.stim1_not = (d1r_mean_activity.correct_trials_suc_not(ran_num_d1r,:));
    stimulus_mode.d1r.stim2_not = (d1r_mean_activity.correct_trials_nacl_not(ran_num_d1r,:));
%     stimulus_mode.non.stim1 =  (non_mean_activity_copy.correct_trials_suc + non_mean_activity_copy.incorrect_trials_suc);
%     stimulus_mode.non.stim2 =  (non_mean_activity_copy.correct_trials_nacl + non_mean_activity_copy.incorrect_trials_nacl);
    stimulus_mode.non.stim1 = (non_mean_activity.correct_trials_suc(ran_num_non,:));
    stimulus_mode.non.stim2 = (non_mean_activity.correct_trials_nacl(ran_num_non,:));
    stimulus_mode.non.selectivity = (stimulus_mode.non.stim1 - stimulus_mode.non.stim2)/2;
    selectivity_across_runs.stimulus.non(k,:) = mean(stimulus_mode.non.selectivity);
    
    stimulus_mode.non.stim1_not = (non_mean_activity.correct_trials_suc_not(ran_num_non,:));
    stimulus_mode.non.stim2_not = (non_mean_activity.correct_trials_nacl_not(ran_num_non,:));
    
    
    clearvars  direction_mode
    
    
%     direction_mode.d1r.contra = (d1r_mean_activity_copy.right_mean_resp + d1r_mean_activity_copy.direction_e_left);
%     direction_mode.d1r.ipsi = (d1r_mean_activity_copy.left_mean_resp + d1r_mean_activity_copy.direction_e_right);
    direction_mode.d1r.contra = (d1r_mean_activity.right_mean_resp(ran_num_d1r,:));
    direction_mode.d1r.ipsi = (d1r_mean_activity.left_mean_resp(ran_num_d1r,:));
    direction_mode.d1r.selectivity = (direction_mode.d1r.contra - direction_mode.d1r.ipsi)/2;
    selectivity_across_runs.direction.d1r(k,:) = mean(direction_mode.d1r.selectivity);
    
    direction_mode.d1r.contra_not = (d1r_mean_activity.right_mean_resp_not(ran_num_d1r,:));
    direction_mode.d1r.ipsi_not = (d1r_mean_activity.left_mean_resp_not(ran_num_d1r,:));
    direction_mode.d1r.selectivity_not= (direction_mode.d1r.contra_not - direction_mode.d1r.ipsi_not)/2;

    selectivity_across_runs.direction.d1r_not(k,:) = mean(direction_mode.d1r.selectivity_not);
    
%     direction_mode.non.contra = (non_mean_activity_copy.right_mean_resp + non_mean_activity_copy.direction_e_left);
%     direction_mode.non.ipsi = (non_mean_activity_copy.left_mean_resp + non_mean_activity_copy.direction_e_right);
    direction_mode.non.contra = (non_mean_activity.right_mean_resp(ran_num_non,:));
    direction_mode.non.ipsi = (non_mean_activity.left_mean_resp(ran_num_non,:));
    direction_mode.non.selectivity= (direction_mode.non.contra - direction_mode.non.ipsi)/2;
    selectivity_across_runs.direction.non(k,:) = mean(direction_mode.non.selectivity);
    
    direction_mode.non.contra_not = (non_mean_activity.right_mean_resp_not(ran_num_non,:));
    direction_mode.non.ipsi_not = (non_mean_activity.left_mean_resp_not(ran_num_non,:));
    direction_mode.non.non_selectivity_not= (direction_mode.non.contra_not - direction_mode.non.ipsi_not)/2;
    selectivity_across_runs.direction.non_not(k,:) = mean(direction_mode.non.non_selectivity_not);
    clearvars  outcome_mode
    
    outcome_mode.d1r.correct = (d1r_mean_activity_copy.left_mean_resp(ran_num_d1r_outcome,:) + d1r_mean_activity_copy.right_mean_resp(ran_num_d1r_outcome,:));
    outcome_mode.d1r.error = (d1r_mean_activity_copy.direction_e_right(ran_num_d1r_outcome,:) + d1r_mean_activity_copy.direction_e_left(ran_num_d1r_outcome,:));
    outcome_mode.d1r.selectivity= (outcome_mode.d1r.correct - outcome_mode.d1r.error)/2;
    selectivity_across_runs.outcome.d1r(k,:) = mean(outcome_mode.d1r.selectivity);
    
    outcome_mode.non.correct = (non_mean_activity_copy.left_mean_resp(ran_num_non_outcome,:) + non_mean_activity_copy.right_mean_resp(ran_num_non_outcome,:));
    outcome_mode.non.error = (non_mean_activity_copy.direction_e_right(ran_num_non_outcome,:) + non_mean_activity_copy.direction_e_left(ran_num_non_outcome,:));
    outcome_mode.non.selectivity= (outcome_mode.non.correct - outcome_mode.non.error)/2;
    selectivity_across_runs.outcome.non(k,:) = mean(outcome_mode.non.selectivity);
    clearvars  ramp_mode_lateral
    
    
    ramp_mode_lateral.d1r.combined_trials = (d1r_mean_activity.right_mean_resp(ran_num_d1r,ramp_window)+d1r_mean_activity.left_mean_resp(ran_num_d1r,ramp_window));
    ramp_mode_lateral.d1r.combined_trials_base = (d1r_mean_activity.right_mean_resp(ran_num_d1r,ramp_window_base)+d1r_mean_activity.left_mean_resp(ran_num_d1r,ramp_window_base));
    ramp_mode_lateral.d1r.selectivity= (mean(ramp_mode_lateral.d1r.combined_trials,2) - mean(ramp_mode_lateral.d1r.combined_trials_base,2))/2;
    
    ramp_mode_lateral.non.combined_trials = (non_mean_activity.right_mean_resp(ran_num_non,ramp_window)+non_mean_activity.left_mean_resp(ran_num_non,ramp_window));
    ramp_mode_lateral.non.combined_trials_base = (non_mean_activity.right_mean_resp(ran_num_non,ramp_window_base)+non_mean_activity.left_mean_resp(ran_num_non,ramp_window_base));
    ramp_mode_lateral.non.selectivity = (mean(ramp_mode_lateral.non.combined_trials,2) - mean(ramp_mode_lateral.non.combined_trials_base,2))/2;
    
    clearvars  response_mode_lateral
    response_mode_lateral.d1r.combined_trials = (d1r_mean_activity.right_mean_resp(ran_num_d1r,response_window)+d1r_mean_activity.left_mean_resp(ran_num_d1r,response_window));
    response_mode_lateral.d1r.combined_trials_base = (d1r_mean_activity.right_mean_resp(ran_num_d1r,delay_window)+d1r_mean_activity.left_mean_resp(ran_num_d1r,delay_window));
    response_mode_lateral.d1r.selectivity= (mean(response_mode_lateral.d1r.combined_trials,2) - mean(response_mode_lateral.d1r.combined_trials_base,2))/2;
    
    response_mode_lateral.non.combined_trials = (non_mean_activity.right_mean_resp(ran_num_non,response_window)+non_mean_activity.left_mean_resp(ran_num_non,response_window));
    response_mode_lateral.non.combined_trials_base = (non_mean_activity.right_mean_resp(ran_num_non,delay_window)+non_mean_activity.left_mean_resp(ran_num_non,delay_window));
    response_mode_lateral.non.selectivity = (mean(response_mode_lateral.non.combined_trials,2) - mean(response_mode_lateral.non.combined_trials_base,2))/2;
    
    
    for i = 1:2
        
        stimulus_mode.(cell_type{i}).selectivity_norm_base = mean(stimulus_mode.(cell_type{i}).selectivity(:,sample_window_base),2)./norm(mean(stimulus_mode.(cell_type{i}).selectivity(:,sample_window_base),2),1);
        stimulus_mode.(cell_type{i}).selectivity_norm = mean(stimulus_mode.(cell_type{i}).selectivity(:,sample_window),2)./norm(mean(stimulus_mode.(cell_type{i}).selectivity(:,sample_window),2),1);
        
        direction_mode.(cell_type{i}).selectivity_norm_delay_base = mean(direction_mode.(cell_type{i}).selectivity(:,delay_window_base),2)./norm(mean(direction_mode.(cell_type{i}).selectivity(:,delay_window_base),2),1);
        direction_mode.(cell_type{i}).selectivity_norm_delay = mean(direction_mode.(cell_type{i}).selectivity(:,delay_window),2)./norm(mean(direction_mode.(cell_type{i}).selectivity(:,delay_window),2),1);
        direction_mode.(cell_type{i}).selectivity_norm_choice = mean(direction_mode.(cell_type{i}).selectivity(:,choice_window),2)./norm(mean(direction_mode.(cell_type{i}).selectivity(:,choice_window),2),1);
        
        outcome_mode.(cell_type{i}).selectivity_norm_baseline = mean(outcome_mode.(cell_type{i}).selectivity(:,outcome_window_base),2)./norm(mean(outcome_mode.(cell_type{i}).selectivity(:,outcome_window_base),2),1);
        outcome_mode.(cell_type{i}).selectivity_norm = mean(outcome_mode.(cell_type{i}).selectivity(:,outcome_window),2)./norm(mean(outcome_mode.(cell_type{i}).selectivity(:,outcome_window),2),1);
        
        ramp_mode_lateral.(cell_type{i}).selectivity_norm_base = mean(ramp_mode_lateral.(cell_type{i}).selectivity(:,1),2)./norm(mean(ramp_mode_lateral.(cell_type{i}).selectivity(:,1),2),1);
        ramp_mode_lateral.(cell_type{i}).selectivity_norm = mean(ramp_mode_lateral.(cell_type{i}).selectivity(:,1),2)./norm(mean(ramp_mode_lateral.(cell_type{i}).selectivity(:,1),2),1);
        
        ramp_mode_central.(cell_type{i}).selectivity_norm_base = mean(ramp_mode_central.(cell_type{i}).selectivity(:,1),2)./norm(mean(ramp_mode_central.(cell_type{i}).selectivity(:,1),2),1);
        ramp_mode_central.(cell_type{i}).selectivity_norm = mean(ramp_mode_central.(cell_type{i}).selectivity(:,1),2)./norm(mean(ramp_mode_central.(cell_type{i}).selectivity(:,1),2),1);
        
%         response_mode.(cell_type{i}).selectivity_norm_base = mean(response_mode.(cell_type{i}).selectivity(:,1),2)./norm(mean(response_mode.(cell_type{i}).selectivity(:,1),2),1);
        response_mode_lateral.(cell_type{i}).selectivity_norm = mean(response_mode_lateral.(cell_type{i}).selectivity(:,1),2)./norm(mean(response_mode_lateral.(cell_type{i}).selectivity(:,1),2),1);
        response_mode_central.(cell_type{i}).selectivity_norm = mean(response_mode_central.(cell_type{i}).selectivity(:,1),2)./norm(mean(response_mode_central.(cell_type{i}).selectivity(:,1),2),1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ramp_mode_central.(cell_type{i}).selectivity_norm_ortho = ramp_mode_central.(cell_type{i}).selectivity_norm;
%         ramp_mode_central.(cell_type{i}).selectivity_norm_ortho = ramp_mode_central.(cell_type{i}).selectivity_norm - dot(ramp_mode_central.(cell_type{i}).selectivity_norm,ramp_mode_central.(cell_type{i}).selectivity_norm_base)*ramp_mode_central.(cell_type{i}).selectivity_norm_base/(norm(ramp_mode_central.(cell_type{i}).selectivity_norm_base)^2);
        stimulus_mode.(cell_type{i}).selectivity_norm_ortho = stimulus_mode.(cell_type{i}).selectivity_norm - dot(stimulus_mode.(cell_type{i}).selectivity_norm,stimulus_mode.(cell_type{i}).selectivity_norm_base)*stimulus_mode.(cell_type{i}).selectivity_norm_base/(norm(stimulus_mode.(cell_type{i}).selectivity_norm_base)^2);
        direction_mode.(cell_type{i}).selectivity_norm_ortho_delay = direction_mode.(cell_type{i}).selectivity_norm_delay - dot(direction_mode.(cell_type{i}).selectivity_norm_delay,direction_mode.(cell_type{i}).selectivity_norm_delay_base)*direction_mode.(cell_type{i}).selectivity_norm_delay_base/(norm(direction_mode.(cell_type{i}).selectivity_norm_delay_base)^2);
        direction_mode.(cell_type{i}).selectivity_norm_ortho_choice = direction_mode.(cell_type{i}).selectivity_norm_choice - dot(direction_mode.(cell_type{i}).selectivity_norm_choice,direction_mode.(cell_type{i}).selectivity_norm_delay)*direction_mode.(cell_type{i}).selectivity_norm_delay/(norm(direction_mode.(cell_type{i}).selectivity_norm_delay)^2);   
        outcome_mode.(cell_type{i}).selectivity_norm_ortho = outcome_mode.(cell_type{i}).selectivity_norm - dot(outcome_mode.(cell_type{i}).selectivity_norm,outcome_mode.(cell_type{i}).selectivity_norm_baseline)*outcome_mode.(cell_type{i}).selectivity_norm_baseline/(norm(outcome_mode.(cell_type{i}).selectivity_norm_baseline)^2);
        ramp_mode_lateral.(cell_type{i}).selectivity_norm_ortho = ramp_mode_lateral.(cell_type{i}).selectivity_norm;
%         ramp_mode_lateral.(cell_type{i}).selectivity_norm_ortho = ramp_mode_lateral.(cell_type{i}).selectivity_norm - dot(ramp_mode_lateral.(cell_type{i}).selectivity_norm,ramp_mode_lateral.(cell_type{i}).selectivity_norm_base)*ramp_mode_lateral.(cell_type{i}).selectivity_norm_base/(norm(ramp_mode_lateral.(cell_type{i}).selectivity_norm_base)^2);
        response_mode_lateral.(cell_type{i}).selectivity_norm_ortho = response_mode_lateral.(cell_type{i}).selectivity_norm;
        response_mode_central.(cell_type{i}).selectivity_norm_ortho = response_mode_central.(cell_type{i}).selectivity_norm;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ramp_projections_central.(cell_type{i}).stim1(k,:) = mean(ramp_mode_central.(cell_type{i}).selectivity_norm.*stimulus_mode.(cell_type{i}).stim1_not,1)*projection_multiplier;
        ramp_projections_central.(cell_type{i}).stim2(k,:) = mean(ramp_mode_central.(cell_type{i}).selectivity_norm.*stimulus_mode.(cell_type{i}).stim2_not,1)*projection_multiplier;
        ramp_projections_central.(cell_type{i}).selectivity(k,:) = ramp_projections_central.(cell_type{i}).stim1(k,:) - ramp_projections_central.(cell_type{i}).stim2(k,:)*projection_multiplier;
        
        stimulus_projections.(cell_type{i}).stim1(k,:) = mean(stimulus_mode.(cell_type{i}).selectivity_norm_ortho.*stimulus_mode.(cell_type{i}).stim1_not,1)*projection_multiplier;
        stimulus_projections.(cell_type{i}).stim2(k,:) = mean(stimulus_mode.(cell_type{i}).selectivity_norm_ortho.*stimulus_mode.(cell_type{i}).stim2_not,1)*projection_multiplier;
        stimulus_projections.(cell_type{i}).selectivity(k,:) = stimulus_projections.(cell_type{i}).stim1(k,:) - stimulus_projections.(cell_type{i}).stim2(k,:)*projection_multiplier;
        
        direction_projections.(cell_type{i}).delay.contra(k,:) = mean(direction_mode.(cell_type{i}).selectivity_norm_ortho_delay.*direction_mode.(cell_type{i}).contra_not,1)*projection_multiplier;
        direction_projections.(cell_type{i}).delay.ipsi(k,:) = mean(direction_mode.(cell_type{i}).selectivity_norm_ortho_delay.*direction_mode.(cell_type{i}).ipsi_not,1)*projection_multiplier;
        direction_projections.(cell_type{i}).delay.selectivity(k,:) = direction_projections.(cell_type{i}).delay.contra(k,:) - direction_projections.(cell_type{i}).delay.ipsi(k,:)*projection_multiplier;
        
        direction_projections.(cell_type{i}).choice.contra(k,:) = mean(direction_mode.(cell_type{i}).selectivity_norm_ortho_choice.*direction_mode.(cell_type{i}).contra_not,1)*projection_multiplier;
        direction_projections.(cell_type{i}).choice.ipsi(k,:) = mean(direction_mode.(cell_type{i}).selectivity_norm_ortho_choice.*direction_mode.(cell_type{i}).ipsi_not,1)*projection_multiplier;
        direction_projections.(cell_type{i}).choice.selectivity(k,:) = direction_projections.(cell_type{i}).choice.contra(k,:) - direction_projections.(cell_type{i}).choice.ipsi(k,:)*projection_multiplier;
        
        outcome_projections.(cell_type{i}).correct(k,:) = mean(outcome_mode.(cell_type{i}).selectivity_norm_ortho.*outcome_mode.(cell_type{i}).correct,1)*projection_multiplier;
        outcome_projections.(cell_type{i}).error(k,:) = mean(outcome_mode.(cell_type{i}).selectivity_norm_ortho.*outcome_mode.(cell_type{i}).error,1)*projection_multiplier;
        outcome_projections.(cell_type{i}).selectivity(k,:) = outcome_projections.(cell_type{i}).correct(k,:) - outcome_projections.(cell_type{i}).error(k,:)*projection_multiplier;
        
        ramp_projections_lateral.(cell_type{i}).contra(k,:) = mean(ramp_mode_lateral.(cell_type{i}).selectivity_norm.*direction_mode.(cell_type{i}).contra_not,1)*projection_multiplier;
        ramp_projections_lateral.(cell_type{i}).ipsi(k,:) = mean(ramp_mode_lateral.(cell_type{i}).selectivity_norm.*direction_mode.(cell_type{i}).ipsi_not,1)*projection_multiplier;
        ramp_projections_lateral.(cell_type{i}).selectivity(k,:) = ramp_projections_lateral.(cell_type{i}).contra(k,:) - ramp_projections_lateral.(cell_type{i}).ipsi(k,:)*projection_multiplier;
        
        response_projections_lateral.(cell_type{i}).contra(k,:) = mean(response_mode_lateral.(cell_type{i}).selectivity_norm.*direction_mode.(cell_type{i}).contra_not,1)*projection_multiplier;
        response_projections_lateral.(cell_type{i}).ipsi(k,:) = mean(response_mode_lateral.(cell_type{i}).selectivity_norm.*direction_mode.(cell_type{i}).ipsi_not,1)*projection_multiplier;
        response_projections_lateral.(cell_type{i}).selectivity(k,:) = response_projections_lateral.(cell_type{i}).contra(k,:) - response_projections_lateral.(cell_type{i}).ipsi(k,:)*projection_multiplier;
        
        response_projections_central.(cell_type{i}).stim1(k,:) = mean(response_mode_central.(cell_type{i}).selectivity_norm.*stimulus_mode.(cell_type{i}).stim1_not,1)*projection_multiplier;
        response_projections_central.(cell_type{i}).stim2(k,:) = mean(response_mode_central.(cell_type{i}).selectivity_norm.*stimulus_mode.(cell_type{i}).stim2_not,1)*projection_multiplier;
        response_projections_central.(cell_type{i}).selectivity(k,:) = response_projections_central.(cell_type{i}).stim1(k,:) - response_projections_central.(cell_type{i}).stim2(k,:)*projection_multiplier;
    end
    
    all_times.suc_not.d1r(k,:) = mean(d1r_mean_activity.correct_trials_suc_time_not,1);
    all_times.nacl_not.d1r(k,:) = mean(d1r_mean_activity.correct_trials_nacl_time_not,1);
    
    all_times.suc_not.non(k,:) = mean(non_mean_activity.correct_trials_suc_time_not,1);
    all_times.nacl_not.non(k,:) = mean(non_mean_activity.correct_trials_nacl_time_not,1);
    
    all_times.left_not.d1r(k,:) = mean(d1r_mean_activity.left_mean_time_not,1);
    all_times.right_not.d1r(k,:) = mean(d1r_mean_activity.right_mean_time_not,1);
    
    all_times.left_not.non(k,:) = mean(non_mean_activity.left_mean_time_not,1);
    all_times.right_not.non(k,:) = mean(non_mean_activity.right_mean_time_not,1);
%     clearvars d1r_mean_activity non_mean_activity
    disp({'Finished run: ' num2str(k) })
%     toc
end
disp('Finishing calculating modes and projecting activity modes')




%%
CI_multiplier = 1;
linewidth = 1.2;
figure(1000);clf;
subplots_to_use = [4:5;7:8];
for i = 1:2
    
    if i == 1
        
        
        stim_x_time1 = mean(all_times.suc_not.d1r);
        stim_x_time2 = mean(all_times.nacl_not.d1r);
        direction_x_time1 = mean(all_times.right_not.d1r  );
        direction_x_time2 = mean(all_times.left_not.d1r  );
        subplot_locations = 1:8;
    else
        

        stim_x_time1 = mean(all_times.suc_not.non);
        stim_x_time2 = mean(all_times.nacl_not.d1r);
        direction_x_time1 = mean(all_times.right_not.non  );
        direction_x_time2 = mean(all_times.left_not.non  );
        subplot_locations = 9:18;
    end
subplot(3,3,subplots_to_use(1,i))
number_of_resamples = 1;
rectangle('position',[-1.5 -1.5 1.5 5],'edgecolor','none','facecolor',[0.7 0.7 0.7 0.5]);
hold on;
h = boundedline(direction_x_time1,mean(direction_projections.(cell_type{i}).delay.contra),std(direction_projections.(cell_type{1}).delay.contra)/sqrt(number_of_resamples)*CI_multiplier,'alpha','b-','LineWidth',linewidth);
hold on;
h = boundedline(direction_x_time2,mean(direction_projections.(cell_type{i}).delay.ipsi),std(direction_projections.(cell_type{1}).delay.ipsi)/sqrt(number_of_resamples)*CI_multiplier,'alpha','r-','LineWidth',linewidth);
if i == 1
title('Lick direction mode - delay','FontSize',12)

end
line([0 0],[-0.1 0.1],'Color','k','LineStyle','--','Linewidth',1);

    ylim([-0.004 0.0063]);yticks(-0.004:0.002:0.006)

xlim([-7 2])
set(gca,'linewidth',1)

subplot(3,3,subplots_to_use(2,i))
rectangle('position',[0 -3 1.5 5],'edgecolor','none','facecolor',[0.7 0.7 0.7 0.5]);
hold on;
h = boundedline(direction_x_time1,mean(direction_projections.(cell_type{i}).choice.contra),(std(direction_projections.(cell_type{2}).choice.contra)/sqrt(number_of_resamples))*CI_multiplier,'alpha','b-','LineWidth',linewidth);
hold on;
h = boundedline(direction_x_time2,mean(direction_projections.(cell_type{i}).choice.ipsi),(std(direction_projections.(cell_type{2}).choice.ipsi)/sqrt(number_of_resamples))*CI_multiplier,'alpha','r-','LineWidth',linewidth);
line([0 0],[-0.1 0.1],'Color','k','LineStyle','--','Linewidth',1);
if i == 1

    title('Lick direction mode - choice','FontSize',12)

end
set(gca,'linewidth',1)
delay_window = 57:67;
choice_window = 68:80;
BW = 0.0001;
xlim([-7 2])
    ylim([-0.004 0.004]);yticks(-0.004:0.002:0.004)

end
subplot(3,3,6)
rectangle('position',[-1.5 -3 1.5 5],'edgecolor','none','facecolor',[0.7 0.7 0.7 0.5]);

boundedline(direction_x_time1,mean(direction_projections.(cell_type{1}).delay.selectivity),(std(direction_projections.(cell_type{1}).delay.selectivity)/sqrt(number_of_resamples))*CI_multiplier,'alpha','m-','LineWidth',linewidth)
hold on;
boundedline(direction_x_time2,mean(direction_projections.(cell_type{2}).delay.selectivity),(std(direction_projections.(cell_type{2}).delay.selectivity)/sqrt(number_of_resamples))*CI_multiplier,'alpha','g-','LineWidth',linewidth)
hold on;
line([0 0],[-0.1 0.1],'Color','k','LineStyle','--','Linewidth',1);
ylim([-0.0003 0.006]);xlim([-7 2])
yticks([-0.001:0.001:0.006])


subplot(3,3,9)
rectangle('position',[0 -3 1.5 5],'edgecolor','none','facecolor',[0.7 0.7 0.7 0.5]);

boundedline(direction_x_time1,mean(direction_projections.(cell_type{1}).choice.selectivity),(std(direction_projections.(cell_type{1}).choice.selectivity)/sqrt(number_of_resamples))*CI_multiplier,'alpha','m-','LineWidth',linewidth)
hold on;
boundedline(direction_x_time2,mean(direction_projections.(cell_type{2}).choice.selectivity),(std(direction_projections.(cell_type{2}).choice.selectivity)/sqrt(number_of_resamples))*CI_multiplier,'alpha','g-','LineWidth',linewidth)
hold on;
line([0 0],[-0.1 0.1],'Color','k','LineStyle','--','Linewidth',1);
ylim([-0.0006 0.006]);xlim([-7 2])
yticks([-0.001:0.001:0.006])

%
figure(1000);

clearvars y
y = mean(direction_projections.(cell_type{1}).delay.selectivity(:,delay_window),2) - mean(direction_projections.(cell_type{2}).delay.selectivity(:,delay_window),2);
mean(mean(direction_projections.(cell_type{1}).delay.selectivity(:,delay_window),2))
mean(mean(direction_projections.(cell_type{2}).delay.selectivity(:,delay_window),2))
std(mean(direction_projections.(cell_type{1}).delay.selectivity(:,delay_window),2))
std(mean(direction_projections.(cell_type{2}).delay.selectivity(:,delay_window),2))

p_delay = size(find(y<=0),1)/1000
if p_delay < 0.05
    text(0.0015, 0.015, num2str(p_delay));
end
axes('Position',[0.918468468468469 0.407127882599582 0.0674324324324322 0.221383647798742])

clear h
h = histogram(y,'normalization','probability','BinWidth',BW,'FaceColor',[0.2 0.2 0.2],'EdgeColor','none');
h=gca; h.XAxis.TickLength = [0 0]; h.YAxis.TickLength = [0 0];
lb = prctile(y,5);
ub = prctile(y,95);
line([mean(y) mean(y)],[0 100],'Color','k','LineWidth',2);
line([lb lb],[0 100],'Color','k','LineStyle','--');
line([ub ub],[0 100],'Color','k','LineStyle','--');
ylim([0 0.06]);yticks([0:0.02:0.06]);
xlim([-0.001 0.004]);xticks(-0.001:0.001:0.004);
camroll(-90)

clearvars y p_choice

y = mean(direction_projections.(cell_type{1}).choice.selectivity(:,choice_window),2) - mean(direction_projections.(cell_type{2}).choice.selectivity(:,choice_window),2);
p_choice = size(find(y<=0),1)/1000

axes('Position',[0.918468468468469 0.103373616387978 0.0674324324324322 0.221383647798742])

h = histogram(y,'normalization','probability','BinWidth',BW,'FaceColor',[0.2 0.2 0.2],'EdgeColor','none');
h=gca; h.XAxis.TickLength = [0 0]; h.YAxis.TickLength = [0 0];
camroll(-90)
lb = prctile(y,2.5);

ub = prctile(y,97.5);
line([mean(y) mean(y)],[0 100],'Color','k','LineWidth',2);
line([lb lb],[0 100],'Color','k','LineStyle','--');
line([ub ub],[0 100],'Color','k','LineStyle','--');
ylim([0 0.1]);yticks([0:0.02:0.1]);
xlim([-0.003 0.001]);xticks([-0.003:0.001:0.001]);
box off


%%
rng('default')
clc;
cd D:\John_data_figure_materials\Fig4
clearvars -except d1r_neurons non_neurons d1r_resp non_resp
save_figures =0;
delay_window = 57:68;
corr_window = 69:77;% 69:77
choice_window = 70:80;
contra_neurons = d1r_neurons(d1r_resp.all_right_selective);
ipsi_neurons = d1r_neurons(d1r_resp.all_left_selective);
% contra_neurons = d1r_neurons;
% ipsi_neurons = d1r_neurons;

clims = [0 1];
for j =1:size(contra_neurons,2) % for each neuron


    x_td1(j,:) = mean([contra_neurons(j).Right.frames_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.frames_aligned(contra_neurons(j).Left.incorrect==1,:)],1);
    resp_temp1 = [contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.incorrect==1,:)];
    resp_temp2 = [contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.correct==1,:);contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.incorrect==1,:)];

    

    d1r_contra.selectivity(j,:) = abs(mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    d1r_contra.preferred(j,:) = mean(resp_temp1,1);
    d1r_contra.nonpreferred(j,:) = mean(resp_temp2,1);
    all_first_central_lick_time_contra(j) = mean([contra_neurons(j).Left.central_licks(:,1); contra_neurons(j).Right.central_licks(:,1)]);
    clear resp_temp1 resp_temp2
end
%

clear resp_temp1 resp_temp2 d1r_ipsi
for j =1:size(ipsi_neurons,2) % for each neuron


    y(j,:) = mean(ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:),1);
    y_td1(j,:) = mean([ipsi_neurons(j).Left.frames_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.frames_aligned(ipsi_neurons(j).Right.incorrect==1,:)],1);
    resp_temp1 = [ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.incorrect==1,:)];
    resp_temp2 = [ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.correct==1,:);ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.incorrect==1,:)];
    all_first_central_lick_time_ipsi(j) = mean([ipsi_neurons(j).Left.central_licks(:,1); ipsi_neurons(j).Right.central_licks(:,1)]);

    d1r_ipsi.selectivity(j,:) = abs(mean(resp_temp1,1)-mean(resp_temp2,1))/2;
    d1r_ipsi.preferred(j,:) = mean(resp_temp1,1);
    d1r_ipsi.nonpreferred(j,:) = mean(resp_temp2,1);
        clear resp_temp1 resp_temp2

end


clear R2 R2_choice x d1r_correlation
x = d1r_contra.selectivity;
for i = 1:1000
    temp = randsample(size(x,1),size(x,1),'true');
    d1r_correlation_contra1(:,:,i) = corr(x(temp,:));

d1r_correlation(i,:) = corr(mean(x(temp,corr_window),2),x(temp,:));

d1r_correlation_contra(i,:)= d1r_correlation(i,:) - mean(mean(d1r_correlation(i,1:24)));

end
clear R2 R2_choice x d1r_correlation
x = d1r_ipsi.selectivity;
for i = 1:1000
    temp = randsample(size(x,1),size(x,1),'true');
    d1r_correlation_ipsi1(:,:,i) = corr(x(temp,:));

d1r_correlation(i,:) = corr(mean(x(temp,corr_window),2),x(temp,:));

d1r_correlation_ipsi(i,:)= d1r_correlation(i,:) - mean(mean(d1r_correlation(i,1:24)));

end




clear R1 R2

clear contra_neurons ipsi_neurons non_contra non_ipsi x_t y_t
contra_neurons = non_neurons(non_resp.all_right_selective  );
ipsi_neurons = non_neurons(non_resp.all_left_selective  );


for j =1:size(contra_neurons,2) % for each neuron
    x_tnon(j,:) = mean([contra_neurons(j).Right.frames_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.frames_aligned(contra_neurons(j).Left.incorrect==1,:)],1);
    resp_temp1 = [contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.incorrect==1,:)];
    resp_temp2 = [contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.correct==1,:);contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.incorrect==1,:)];

    non_contra.selectivity(j,:) = abs(mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    non_contra.preferred(j,:) = mean(resp_temp1,1);
    non_contra.nonpreferred(j,:) = mean(resp_temp2,1);
    all_first_central_lick_time_contra(j) = mean([contra_neurons(j).Left.central_licks(:,1); contra_neurons(j).Right.central_licks(:,1)]);
    clear resp_temp1 resp_temp2
end
%

clear resp_temp1 resp_temp2 
for j =1:size(ipsi_neurons,2) % for each neuron

    y_tnon(j,:) = mean([ipsi_neurons(j).Left.frames_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.frames_aligned(ipsi_neurons(j).Right.incorrect==1,:)],1);
    resp_temp1 = [ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.incorrect==1,:)];
    resp_temp2 = [ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.correct==1,:);ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.incorrect==1,:)];
    all_first_central_lick_time_ipsi(j) = mean([ipsi_neurons(j).Left.central_licks(:,1); ipsi_neurons(j).Right.central_licks(:,1)]);

    non_ipsi.selectivity(j,:) = abs(mean(resp_temp1,1)-mean(resp_temp2,1))/2;
    non_ipsi.preferred(j,:) = mean(resp_temp1,1);
    non_ipsi.nonpreferred(j,:) = mean(resp_temp2,1);
        clear resp_temp1 resp_temp2

end
%
clear R2 x non_correlation
x = non_contra.selectivity;
for i = 1:1000
temp = randsample(size(x,1),size(x,1),'true');
    non_correlation_contra1(:,:,i) = corr(x(temp,:));

    non_correlation(i,:) = corr(mean(x(temp,corr_window),2),x(temp,:));

non_correlation_contra(i,:)= non_correlation(i,:) - mean(mean(non_correlation(i,1:24)));
% 
end
clear R2 x non_correlation
x = non_ipsi.selectivity;
for i = 1:1000
temp = randsample(size(x,1),size(x,1),'true');
    non_correlation_ipsi1(:,:,i) = corr(x(temp,:));

    non_correlation(i,:) = corr(mean(x(temp,corr_window),2),x(temp,:));

non_correlation_ipsi(i,:)= non_correlation(i,:) - mean(mean(non_correlation(i,1:24)));

end


%
% close all
rng('default')
figure(1);clf;
axes_limits = [-2 1.5];

subplot(2,2,1)

h1 = imagesc(mean(x_td1),mean(x_td1),mean(d1r_correlation_contra1-mean(mean(d1r_correlation_contra1(:,1:24,:),1)),3),clims); hold on;
line([0 0],[-10 4],'color','k','LineWidth',2);hold on
line([-10 4],[0 0],'color','k','LineWidth',2); hold on;
line([-10 4],[mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],'color','k','LineWidth',2); hold on
line([mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],[-10 4],'color','k','LineWidth',2); hold on
% set(gca,'YTick',axes_limits(1):0.5:axes_limits(2));
% set(gca,'XTick',axes_limits(1):0.5:axes_limits(2));
xticks([-2:1.5]);yticks([-2:1.5])
xlim(axes_limits);ylim(axes_limits);ylabel('Time (s)');xlabel('Time (s)')
colorbar
title('D1R+ contra')
% axis square
box off
subplot(2,2,2)
h1 = imagesc(mean(x_tnon),mean(x_tnon),mean(non_correlation_contra1-mean(mean(non_correlation_contra1(:,1:24,:),1)),3),clims); hold on;
line([0 0],[-10 4],'color','k','LineWidth',2);hold on
line([-10 4],[0 0],'color','k','LineWidth',2); hold on;
line([-10 4],[mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],'color','k','LineWidth',2); hold on
line([mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],[-10 4],'color','k','LineWidth',2); hold on
colormap(turbo)
xlim(axes_limits);ylim(axes_limits);ylabel('Time (s)');xlabel('Time (s)')
xticks([-2:1.5]);yticks([-2:1.5])
% colormap(parula)
% axis square
colorbar
title('D1R- contra')
box off

subplot(2,2,3)
h1 = imagesc(mean(y_td1),mean(y_td1),mean(d1r_correlation_ipsi1-mean(mean(d1r_correlation_ipsi1(:,1:24,:),1)),3),clims); hold on;
line([0 0],[-10 4],'color','k','LineWidth',2);hold on
line([-10 4],[0 0],'color','k','LineWidth',2); hold on;
line([-10 4],[mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],'color','k','LineWidth',2); hold on
line([mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],[-10 4],'color','k','LineWidth',2); hold on
colormap(turbo)
xlim(axes_limits);ylim(axes_limits);ylabel('Time (s)');xlabel('Time (s)')
xticks([-2:1.5]);yticks([-2:1.5])
% colormap(parula)
% axis square
colorbar
title('D1R+ ipsi')
box off
subplot(2,2,4)
h1 = imagesc(mean(y_tnon),mean(y_tnon),mean(non_correlation_ipsi1-mean(mean(non_correlation_ipsi1(:,1:24,:),1)),3),clims); hold on;
line([0 0],[-10 4],'color','k','LineWidth',2);hold on
line([-10 4],[0 0],'color','k','LineWidth',2); hold on;
line([-10 4],[mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],'color','k','LineWidth',2); hold on
line([mean(all_first_central_lick_time_contra) mean(all_first_central_lick_time_contra)],[-10 4],'color','k','LineWidth',2); hold on
colormap(turbo)
xlim(axes_limits);ylim(axes_limits);ylabel('Time (s)');xlabel('Time (s)')
xticks([-2:1.5]);yticks([-2:1.5])
% colormap(parula)
% axis square
colorbar
title('D1R- ipsi')
box off


figure(3);clf;
subplot(1,4,1)
rectangle('position',[-1.5 -1 1.5 2],'edgecolor','none','facecolor',[0.7 0.7 0.7 0.5]);
boundedline(mean(x_td1),mean(d1r_correlation_contra,1),std(d1r_correlation_contra),'alpha','m-')
boundedline(mean(x_tnon),mean(non_correlation_contra,1),std(non_correlation_contra),'alpha','g-')
line([0 0],[0 1],'color','k','LineWidth',1,'LineStyle','--');hold on
ylabel('Population selectivity correlation'); xlabel('Time (s)')
xlim([-6 2]);ylim([0 1]);yticks([0:0.2:1]);
subplot(1,4,2)
g1 =  mean(d1r_correlation_contra(:,delay_window,:),2);
g2 = mean(non_correlation_contra(:,delay_window,:),2);
b = bar(1:2,[mean(g1);mean(g2)],'BarWidth', 0.5); hold on
errorbar(1:2,[mean(g1);mean(g2)],[std2(g1) std2(g2)],'linestyle','none','Color','k')
ylim([0 1]);yticks([0:0.2:1]);xticks([1,2]);
b(1,1).FaceColor = 'flat';
b(1).CData(1,:)  = [1 0 1];
b(1).CData(2,:)= [0 1 0];
box off
xtickangle(45);set(gca,'xticklabel',{'D1R+'; 'D1R-'});
xlim([0 3])

subplot(1,4,3)
rectangle('position',[-1.5 -1 1.5 2],'edgecolor','none','facecolor',[0.7 0.7 0.7 0.5]);

boundedline(mean(y_td1),mean(d1r_correlation_ipsi,1),std(d1r_correlation_ipsi),'alpha','m-')
boundedline(mean(y_tnon),mean(non_correlation_ipsi,1),std(non_correlation_ipsi),'alpha','g-')

line([0 0],[0 1],'color','k','LineWidth',1,'LineStyle','--');hold on
ylabel('Population selectivity correlation'); xlabel('Time (s)')
xlim([-6 2]);ylim([0 1]);yticks([0:0.2:1]);
subplot(1,4,4)
g1 =  mean(d1r_correlation_ipsi(:,delay_window,:),2);
g2 = mean(non_correlation_ipsi(:,delay_window,:),2);
b = bar(1:2,[mean(g1);mean(g2)],'BarWidth', 0.5); hold on
errorbar(1:2,[mean(g1);mean(g2)],[std2(g1) std2(g2)],'linestyle','none','Color','k')
b(1,1).FaceColor = 'flat';
b(1).CData(1,:)  = [1 0 1];
b(1).CData(2,:)= [0 1 0];
box off
ylim([0 1])
xtickangle(45);xlim([0 3]);xticks([1,2]);yticks([0:0.2:1]);
set(gca,'xticklabel',{'D1R+'; 'D1R-'});


clearvars p_values x y
g1 =  mean(d1r_correlation_contra(:,delay_window,:),2);
g2 = mean(non_correlation_contra(:,delay_window,:),2);
x = squeeze(g1)';
        y = squeeze(g2)';
        t_stat = mean(x) - mean(y);
        
        pooled = [x,y];
        for b = 1:1000
            temp1 = randsample(1:size([x,y],2),size(x,2),true);
            temp2 = randsample(1:size([x,y],2),size(y,2),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
                p_values.pearson_correlation.d1rcontravsnoncontra = (length(find((x-y)<=0))+1)/1001

    clearvars x y
    
g1 =  mean(d1r_correlation_contra(:,delay_window,:),2);
g2 = mean(d1r_correlation_ipsi(:,delay_window,:),2);
x = squeeze(g1)';
        y = squeeze(g2)';
        t_stat = mean(x) - mean(y);
        
        pooled = [x,y];
        for b = 1:1000
            temp1 = randsample(1:size([x,y],2),size(x,2),true);
            temp2 = randsample(1:size([x,y],2),size(y,2),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
                p_values.pearson_correlation.d1ripsivsd1rcontra = (length(find((x-y)<=0))+1)/1001

clearvars x y
    
g1 =  mean(d1r_correlation_contra(:,delay_window,:),2);
g2 = mean(non_correlation_ipsi(:,delay_window,:),2);
x = squeeze(g1)';
        y = squeeze(g2)';
        t_stat = mean(x) - mean(y);
        
        pooled = [x,y];
        for b = 1:1000
            temp1 = randsample(1:size([x,y],2),size(x,2),true);
            temp2 = randsample(1:size([x,y],2),size(y,2),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
                p_values.pearson_correlation.d1rcontravsnonipsi = (length(find((x-y)<=0))+1)/1001
                
                clearvars x y
    
g1 =  mean(d1r_correlation_ipsi(:,delay_window,:),2);
g2 = mean(non_correlation_ipsi(:,delay_window,:),2);
x = squeeze(g1)';
        y = squeeze(g2)';
        t_stat = mean(x) - mean(y);
        
        pooled = [x,y];
        for b = 1:1000
            temp1 = randsample(1:size([x,y],2),size(x,2),true);
            temp2 = randsample(1:size([x,y],2),size(y,2),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
                p_values.pearson_correlation.d1ripsivsnonipsi = (length(find((x-y)<=0))+1)/1001
                
                  clearvars x y
    
g1 =  mean(non_correlation_contra(:,delay_window,:),2);
g2 = mean(non_correlation_ipsi(:,delay_window,:),2);
x = squeeze(g1)';
        y = squeeze(g2)';
        t_stat = mean(x) - mean(y);
        
        pooled = [x,y];
        for b = 1:1000
            temp1 = randsample(1:size([x,y],2),size(x,2),true);
            temp2 = randsample(1:size([x,y],2),size(y,2),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
                p_values.pearson_correlation.noncontravsnonipsi = (length(find((x-y)<=0))+1)/1001

set(figure(1),'Position',[647 533 591 404]);
set(figure(3),'Position',[653 532 1118 220]);
