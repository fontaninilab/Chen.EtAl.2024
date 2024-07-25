%% Plot example preparatory responses 

clearvars -except d1r_neurons non_neurons d1r_resp non_resp
close all;
clc;

figure_startup
figure(2);clf;
d1r_example_lick_neurons = [36];%% 29 36 d1r_lick_only = 5 6 21 24 27 28 29 35 42 83 32 17 14 20 30 43 54
non_example_lick_neurons = [330];%%  310  169 184 191 192 221 27 87 82 173 248
all_example_neurons = [d1r_example_lick_neurons;non_example_lick_neurons];
J = customcolormap([0 0.5 1], [255,255,255;150,150,150; 0,0,0]./255);
for k = 1:2
    
    cell_type = all_example_neurons(k,:);
    clear population group
    group = 'lick_only_1';
    if k == 1
        population = d1r_neurons(d1r_resp.(group));
    else
        population = non_neurons(non_resp.(group));
    end
    
    for i = 1:size(cell_type,2)
        
        num_t1 = size(population(cell_type(i)).Sucrose.dff,1);  %Get the number of first trial type
        trial_temp = population(cell_type(i)).all_trials.lick_aligned_traces_matrix;
        
        time_temp = population(cell_type(i)).all_trials.lick_aligned_frames_matrix;
        
        clims = [0 6];

        ax(1) = subplot(2,4,k);
        
        imagesc(mean(time_temp),1:size(trial_temp,1),trial_temp,clims); hold on;
        
        line([0 0],[0 size(trial_temp,1)+2],'color','k','LineWidth',1,'LineStyle','--'); hold on;
                colormap(ax(1),J);
        
        xlim([-3 1]);set(gca,'XTick',[]);
        yticks([1 size(trial_temp,1)]); hold on;
        colorbar('Position',[0.95 0.6 0.01 0.38],'YTick', 0:2:6)
        box off
        ylabel('Trial #');
        ax(2) = subplot(2,4,k+4);
        
        boundedline(mean(time_temp),mean(trial_temp),std(trial_temp)./sqrt(size(trial_temp,1)),'alpha','LineWidth',1.5,'color','k'); hold on;
        line([0 0],[min(mean(trial_temp(:,16:43)))-1 max(mean(trial_temp))+2],'color','k','LineWidth',1,'LineStyle','--'); hold on;
        if k == 1
            ylim([-2.5 6]);
             yticks([-2:2:6]); 
        elseif k == 2
            ylim([-0.5 8]);
            yticks([-2:2:8]);
        end
        xlim([-3 1]);
        
        xlabel('Time (s)'); ylabel('\Delta F/F')
        %         axis square
        
    end
    set(gca,'TickDir','out'); set(gca,'linewidth',1)
    
    linkaxes(ax,'x');
end

%
d1r_example_lick_neurons = [68];%d1r_delayonly =58
non_example_lick_neurons = [170];%%220 189 195 223 231 167 169 170 179
all_example_neurons = [d1r_example_lick_neurons;non_example_lick_neurons];
% close all
for k = 1:2
    
    cell_type = all_example_neurons(k,:);
    clear population group
    group = 'delay_resp_only';
    if k == 1
        population = d1r_neurons(d1r_resp.(group));
    else
        population = non_neurons(non_resp.(group));
    end
    
    for i = 1:size(cell_type,2)
        
        num_t1 = size(population(cell_type(i)).Sucrose.dff,1);  %Get the number of first trial type
        trial_temp = [population(cell_type(i)).Right.dff_aligned; population(cell_type(i)).Left.dff_aligned];
        
        time_temp = (population(cell_type(i)).Right.mean_frames_aligned+population(cell_type(i)).Right.mean_frames_aligned)/2;
        clims = [0 8];
        
        %
        ax(1) = subplot(2,4,k+2);
        
        imagesc(time_temp,1:size(trial_temp,1),trial_temp,clims); hold on;
        
        line([0 0],[0 size(trial_temp,1)+2],'color','k','LineWidth',1,'LineStyle','--'); hold on;
                colormap(ax(1),J);
        xlim([-5 2]);set(gca,'XTick',[]);
        yticks([1 size(trial_temp,1)]); hold on;
                colorbar('Position',[0.95 0.15 0.01 0.38],'YTick', 0:2:8)
        box off
        ylabel('Trial #');
        ax(2) = subplot(2,4,k+6);
        
        boundedline(time_temp,mean(trial_temp),std(trial_temp)./sqrt(size(trial_temp,1)),'alpha','LineWidth',1.5,'color','k'); hold on;
        line([0 0],[min(mean(trial_temp(:,16:43)))-1 max(mean(trial_temp))+2],'color','k','LineWidth',1,'LineStyle','--'); hold on;
        if k == 1
            ylim([-2 10])
        elseif k == 2
            ylim([-2 8]);
        end
        xlim([-3 1]);
        xlabel('Time (s)');yticks([-2:2:10]); hold on;
        ylabel('\Delta F/F')        
    end
    set(gca,'TickDir','out'); set(gca,'linewidth',1)
    
    linkaxes(ax,'x');
end
set(gcf,'Position',[417 718 938 260]);


%% Plot mean preparatory responses and quantification of onset and time-to-peak times
clc
% close all
CloseVariableWindows
clearvars -except d1r_neurons non_neurons d1r_resp non_resp
subplot_dimensions = [3,4];
number_of_error_trials = 3;
clear time_of_peak_d1r1 time_of_peak_non1 time_of_peak_non directionAA_non end_time mean_d1r mean_non mean_response
clear peak_d1r peak_non p_values
figure(3);clf;
figure(13);clf;
rng('default')
clims = [0 1];
threshold = 0.25;
FC_preparation = 19:38;%19:38
LL_preparation_early = 36:69;%36:69
lateral_aligned_preparatory_window = 38:71;%38:71
BW = 0.15;
subplot_positions = [9,11,17];
subplot_positions2 = [10,12,18];
J = customcolormap([0 0.5 1], [255,255,217;29,145,192; 8,29,88]./255);

for i = 1:2
    resp_groups = {'lick_only_1','delay_resp_only','lick_delay_resp_only'};
    a = d1r_resp.(resp_groups{i});
    b = non_resp.(resp_groups{i});
    group1 = d1r_neurons;
    group2 = non_neurons;
    
    randomized = 0;
    [d1r_mean_activity, non_mean_activity] = plot_mean_response_for_all_neurons(a,b,group1, group2,number_of_error_trials,randomized);
    
    figure(3);
    
    subplot(subplot_dimensions(1),subplot_dimensions(2),subplot_positions(i))
    
    line([0 0],[min(mean(d1r_mean_activity.all_trials_taste_aligned))-1 max(mean(d1r_mean_activity.all_trials_taste_aligned))+1],'color','k','LineWidth',1.5,'LineStyle','--'); hold on;
    if i == 1
        boundedline(mean(non_mean_activity.all_trials_drylick_aligned_time,'omitnan'),mean(non_mean_activity.all_trials_drylick_aligned_activity ,'omitnan'),std(non_mean_activity.all_trials_drylick_aligned_activity,'omitnan')./sqrt(size(non_mean_activity.all_trials_drylick_aligned_activity,1)),'alpha','g-','LineWidth',1.3)
        boundedline(mean(d1r_mean_activity.all_trials_drylick_aligned_time,'omitnan'),mean(d1r_mean_activity.all_trials_drylick_aligned_activity,'omitnan'),std(d1r_mean_activity.all_trials_drylick_aligned_activity,'omitnan')./sqrt(size(d1r_mean_activity.all_trials_drylick_aligned_activity,1)),'alpha','m-','LineWidth',1.3)
        ylabel('\Delta F/F');xlabel('Time from central lick (s)')
    elseif i == 2
        boundedline(mean(non_mean_activity.all_trials_taste_aligned_time,'omitnan'),mean(non_mean_activity.all_trials_taste_aligned ,'omitnan'),std(non_mean_activity.all_trials_taste_aligned,'omitnan')./sqrt(size(non_mean_activity.all_trials_taste_aligned,1)),'alpha','g-','LineWidth',1.3)
        boundedline(mean(d1r_mean_activity.all_trials_taste_aligned_time,'omitnan'),mean(d1r_mean_activity.all_trials_taste_aligned,'omitnan'),std(d1r_mean_activity.all_trials_taste_aligned,'omitnan')./sqrt(size(d1r_mean_activity.all_trials_taste_aligned,1)),'alpha','m-','LineWidth',1.3)
        ylabel('\Delta F/F');xlabel('Time from central lick (s)')
    end
    xlim([-3 2]);
    if i == 1
        ylim([-0.3 3]);
        yticks([-1:1:3])
    elseif i == 2
        ylim([-1 5.2]);
        yticks([-1:1:5])
    elseif i == 3
        ylim([-0.5 7]);
    end
    figure(3);
    
    subplot(subplot_dimensions(1),subplot_dimensions(2),subplot_positions2(i))
    line([0 0],[min(mean(d1r_mean_activity.all_trials_lateral_aligned))-1 max(mean(d1r_mean_activity.all_trials_lateral_aligned))+5],'color','k','LineWidth',1.5,'LineStyle','--'); hold on;

    
    boundedline(mean(non_mean_activity.all_trials_lateral_aligned_time,'omitnan'),mean(non_mean_activity.all_trials_lateral_aligned,'omitnan'),std(non_mean_activity.all_trials_lateral_aligned,'omitnan')./sqrt(size(non_mean_activity.all_trials_lateral_aligned,1)),'alpha','g-','LineWidth',1.3)
    boundedline(mean(d1r_mean_activity.all_trials_lateral_aligned_time,'omitnan'),mean(d1r_mean_activity.all_trials_lateral_aligned,'omitnan'),std(d1r_mean_activity.all_trials_lateral_aligned,'omitnan')./sqrt(size(d1r_mean_activity.all_trials_lateral_aligned,1)),'alpha','m-','LineWidth',1.3)
    xlim([-6 2]);ylim([-0.5 3.5]);xlabel('Time from lateral lick (s)')
    if i == 1
        ylim([-0.3 2]);
                yticks([-1:1:2])

    elseif i == 2
        ylim([-1 5.2]);
                yticks([-1:1:5])

    end
    mean_response.d1r(i).central_preparatory = mean(d1r_mean_activity.all_trials_drylick_aligned_activity(:,FC_preparation),2);
    mean_response.non(i).central_preparatory = mean(non_mean_activity.all_trials_drylick_aligned_activity(:,FC_preparation),2);
    
    mean_response.d1r(i).lateral_preparatory_late = mean(d1r_mean_activity.all_trials_lateral_aligned(:,lateral_aligned_preparatory_window),2);
    mean_response.non(i).lateral_preparatory_late = mean(non_mean_activity.all_trials_lateral_aligned(:,lateral_aligned_preparatory_window),2);
    
    mean_response.d1r(i).lateral_preparatory_early = mean(d1r_mean_activity.all_trials_taste_aligned(:,LL_preparation_early),2);
    mean_response.non(i).lateral_preparatory_early = mean(non_mean_activity.all_trials_taste_aligned(:,LL_preparation_early),2);
    
    mean_response.d1r(i).central_preparatory_allcells = d1r_mean_activity.all_trials_drylick_aligned_activity;
    mean_response.d1r(i).lateral_preparatory_allcells = d1r_mean_activity.all_trials_taste_aligned;
    mean_response.d1r(i).central_preparatory_allcells_time = d1r_mean_activity.all_trials_drylick_aligned_time;
    mean_response.d1r(i).lateral_preparatory_allcells_time = d1r_mean_activity.all_trials_taste_aligned_time;
    
    mean_response.non(i).central_preparatory_allcells = non_mean_activity.all_trials_drylick_aligned_activity;
    mean_response.non(i).lateral_preparatory_allcells = non_mean_activity.all_trials_taste_aligned;
    mean_response.non(i).central_preparatory_allcells_time = non_mean_activity.all_trials_drylick_aligned_time;
    mean_response.non(i).lateral_preparatory_allcells_time = non_mean_activity.all_trials_taste_aligned_time;
    
    num_of_response.d1r(i) = size(d1r_mean_activity.all_trials_taste_aligned,1);
    num_of_response.non(i) = size(non_mean_activity.all_trials_taste_aligned,1);
    
    dry_lick_aligned_time.d1r = mean(d1r_mean_activity.all_trials_taste_aligned_time,'omitnan');
    lateral_lick_aligned_time.d1r = mean(d1r_mean_activity.all_trials_lateral_aligned_time,'omitnan');
    
    dry_lick_aligned_time.non = mean(non_mean_activity.all_trials_taste_aligned_time,'omitnan');
    lateral_lick_aligned_time.non = mean(non_mean_activity.all_trials_lateral_aligned_time,'omitnan');
    
    
    clear x y
    figure(13);
    if i == 1
        search_window = FC_preparation;
        x = d1r_mean_activity.all_trials_drylick_aligned_activity;
        y = non_mean_activity.all_trials_drylick_aligned_activity;
        x_time = d1r_mean_activity.all_trials_drylick_aligned_time;
        y_time = non_mean_activity.all_trials_drylick_aligned_time;
        
        [d1r_ramp_onset, non_ramp_onset] = ramp_onset(x,y,x_time,y_time,search_window(1),search_window(end));
        clearvars x y bootStrap_sample temp1 temp2 store1 store2
        x = d1r_ramp_onset';
        y = non_ramp_onset';
        t_stat = mean(x) - mean(y);
        
        pooled = [x;y];
        for b = 1:1000
            temp1 = randsample(1:size([x;y],1),size(x,1),true);
            temp2 = randsample(1:size([x;y],1),size(y,1),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
        p_values.bootstrap_centralonset(i) = ((length(find(bootStrap_sample<=t_stat)))+1)/1001;
        
      
        subplot(2,4,i)
        
        h1 = cdfplot(non_ramp_onset);hold on;ylabel('')
        h2 = cdfplot(d1r_ramp_onset);grid off
        set( h2(1), 'LineStyle', '-', 'Color', 'm','LineWidth',1.3);
        grid off
        set( h1(1), 'LineStyle', '-', 'Color', 'g','LineWidth',1.3);
%         text(-1, 0.6,num2str(p_values.bootstrap_centralonset(i)));
        title({'Onset of response' '(central)'});        box off
        
        subplot(2,4,i+1)
        for b = 1:1000
            
            rand_idx1 = randsample(1:size(x,1),size(x,1),true);
            rand_idx2 = randsample(1:size(y,1),size(y,1),true);
            bootstrap_samples.centralonsettime.d1r(b) = mean(x(rand_idx1,1));
            bootstrap_samples.centralonsettime.non(b) = mean(y(rand_idx2,1));
            
        end
        bar(1,mean(bootstrap_samples.centralonsettime.d1r),'FaceColor','none','Edgecolor','m','LineWidth',1.1);
        hold on;
        bar(2,mean(bootstrap_samples.centralonsettime.non),'FaceColor','none','Edgecolor','g','LineWidth',1.1)
        errorbar([1 2], [mean(bootstrap_samples.centralonsettime.d1r) mean(bootstrap_samples.centralonsettime.non)],[std(bootstrap_samples.centralonsettime.d1r) std(bootstrap_samples.centralonsettime.non)],'LineStyle','none','Color','k')
        ylim([-0.8 0])
        box off
%         text(1.25,-0.75,['p = ' num2str(p_values.bootstrap_centralonset(i))])
        set(gca,'YDir','reverse');
        camroll(90)
         
    elseif i == 2
                clearvars x y bootStrap_sample temp1 temp2 store1 store2 t_stat pooled p_values.bootstrap_lateralonset

        search_window = lateral_aligned_preparatory_window;
        x = d1r_mean_activity.all_trials_lateral_aligned;
        y = non_mean_activity.all_trials_lateral_aligned;
        x_time = d1r_mean_activity.all_trials_lateral_aligned_time;
        y_time = non_mean_activity.all_trials_lateral_aligned_time;
        [d1r_ramp_onset, non_ramp_onset] = ramp_onset(x,y,x_time,y_time,search_window(1),search_window(end));
        clearvars x y bootStrap_sample temp1 temp2 store1 store2 t_stat pooled
        x = d1r_ramp_onset';
        y = non_ramp_onset';
        t_stat = mean(x) - mean(y);
        
        pooled = [x;y];
        for b = 1:1000
            temp1 = randsample(1:size([x;y],1),size(x,1),true);
            temp2 = randsample(1:size([x;y],1),size(y,1),true);
            bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        end
        p_values.bootstrap_lateralonset = ((length(find(bootStrap_sample<=t_stat)))+1)/1001;
    
        
        subplot(2,4,i+3)
        
        h1 = cdfplot(non_ramp_onset);hold on;ylabel('')
        h2 = cdfplot(d1r_ramp_onset);grid off
        set( h2(1), 'LineStyle', '-', 'Color', 'm','LineWidth',1.3);
        grid off
        set( h1(1), 'LineStyle', '-', 'Color', 'g','LineWidth',1.3);
%         text(-1, 0.6,num2str(p_values.bootstrap_lateralonset));
        title({'Onset of response' '(lateral)'});
        box off
        subplot(2,4,i+4)
        for b = 1:1000
            
            rand_idx1 = randsample(1:size(x,1),size(x,1),true);
            rand_idx2 = randsample(1:size(y,1),size(y,1),true);
            bootstrap_samples.lateralonsettime.d1r(b) = mean(x(rand_idx1,1));
            bootstrap_samples.lateralonsettime.non(b) = mean(y(rand_idx2,1));
            
        end
        bar(1,mean(bootstrap_samples.lateralonsettime.d1r),'FaceColor','none','Edgecolor','m','LineWidth',1.1);
        hold on;
        bar(2,mean(bootstrap_samples.lateralonsettime.non),'FaceColor','none','Edgecolor','g','LineWidth',1.1)
        errorbar([1 2], [mean(bootstrap_samples.lateralonsettime.d1r) mean(bootstrap_samples.lateralonsettime.non)],[std(bootstrap_samples.lateralonsettime.d1r) std(bootstrap_samples.lateralonsettime.non)],'LineStyle','none','Color','k')
       
        ylim([-2 0])
        box off
%         text(1.5,-1.8,['p = ' num2str(p_values.bootstrap_lateralonset)])
        set(gca,'YDir','reverse');
        camroll(90)
    end
    clearvars d1r_mean_activity non_mean_activity
   
end

figure(3);
colormap parula
clear I_d1r
temp = mean_response.d1r(1).central_preparatory_allcells;
for q = 1:size(temp,1)
    r(q) = length(find(temp(q,[FC_preparation(1)-5:FC_preparation(end)+1])>0.5));
end
    [~,I_d1r.central_sorted] = sort(r);

temp = mean_response.d1r(2).central_preparatory_allcells;
for q = 1:size(temp,1)
    r(q) = length(find(temp(q,[LL_preparation_early(1):LL_preparation_early(end)])>0.3));
end
[~,I_d1r.lateral_sorted] = sort(r);

all_d1r_prep_responses = [mean_response.d1r(1).central_preparatory_allcells(flipud(I_d1r.central_sorted'),:); mean_response.d1r(2).central_preparatory_allcells(flipud(I_d1r.lateral_sorted'),:)];

% subplot(subplot_dimensions(1),subplot_dimensions(2),[1])
% % clearvars I_d1r I_non r
% clearvars y_data
% y_data = mean_response.d1r(1).central_preparatory_allcells(flipud(I_d1r.central_sorted'),:);
% imagesc(dry_lick_aligned_time.d1r,1:size(y_data,1),y_data,clims);
% xlim([-3 2]);yticks([1 size(y_data,1)])
% line([0 0],[0 size(mean_response.d1r(1).central_preparatory_allcells,1)+size(mean_response.d1r(2).central_preparatory_allcells,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% box off
% subplot(subplot_dimensions(1),subplot_dimensions(2),[3])
% clearvars y_data
% y_data = mean_response.d1r(2).central_preparatory_allcells(flipud(I_d1r.lateral_sorted'),:);
% imagesc(dry_lick_aligned_time.d1r,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(mean_response.d1r(1).central_preparatory_allcells,1)+size(mean_response.d1r(2).central_preparatory_allcells,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% box off
% colormap(J)
% hold on;
% line([0 0],[0 size(mean_response.d1r(1).central_preparatory_allcells,1)+size(mean_response.d1r(2).central_preparatory_allcells,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% 
% xlim([-3 2]);yticks([1 size(y_data,1)])

% box off
% clearvars I_non r
% 
% temp = mean_response.non(1).central_preparatory_allcells;
% for q = 1:size(temp,1)
%     r(q) = length(find(temp(q,[FC_preparation(1)-5:FC_preparation(end)+1])>0.5));
% end
% [~,I_non.central_sorted] = sort(r);
% 
% temp = mean_response.non(2).central_preparatory_allcells;
% for q = 1:size(temp,1)
%     r(q) = length(find(temp(q,[LL_preparation_early(1):LL_preparation_early(end)])>0.3));
% end
% [~,I_non.lateral_sorted] = sort(r);
% 
% all_non_prep_responses = [mean_response.non(1).central_preparatory_allcells(flipud(I_non.central_sorted'),:); mean_response.non(2).central_preparatory_allcells(flipud(I_non.lateral_sorted'),:)];
% clearvars y_data 
% subplot(subplot_dimensions(1),subplot_dimensions(2),[5])
% 
% y_data = mean_response.non(1).central_preparatory_allcells(flipud(I_non.central_sorted'),:);
% imagesc(dry_lick_aligned_time.non,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(mean_response.non(1).central_preparatory_allcells,1)+size(mean_response.non(2).central_preparatory_allcells,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% box off

% xlim([-3 2]);yticks([1 size(y_data,1)])
% 
% subplot(subplot_dimensions(1),subplot_dimensions(2),[7])
% clearvars y_data
% y_data = mean_response.non(2).central_preparatory_allcells(flipud(I_non.lateral_sorted'),:);
% imagesc(dry_lick_aligned_time.non,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(y_data,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% 
% xlim([-3 2]);yticks([1 size(y_data,1)])
% box off
% %
% clearvars all_d1r_prep_responses all_non_prep_responses
% 
% all_d1r_prep_responses = [mean_response.d1r(1).lateral_preparatory_allcells(flipud(I_d1r.central_sorted'),:); mean_response.d1r(2).lateral_preparatory_allcells(flipud(I_d1r.lateral_sorted'),:)];
% all_non_prep_responses = [mean_response.non(1).lateral_preparatory_allcells(flipud(I_non.central_sorted'),:); mean_response.non(2).lateral_preparatory_allcells(flipud(I_non.lateral_sorted'),:)];
% subplot(subplot_dimensions(1),subplot_dimensions(2),[2])
% clearvars y_data
% y_data = mean_response.d1r(1).lateral_preparatory_allcells(flipud(I_d1r.central_sorted'),:);
% imagesc(lateral_lick_aligned_time.d1r,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(y_data,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% box off
% xlim([-1.5 2]);set(gca,'ytick',[])


% subplot(subplot_dimensions(1),subplot_dimensions(2),[6])
% clearvars y_data
% y_data =mean_response.non(1).lateral_preparatory_allcells(flipud(I_non.central_sorted'),:);
% imagesc(lateral_lick_aligned_time.d1r,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(y_data,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% box off
% xlim([-1.5 2]);set(gca,'ytick',[])
% 
% subplot(subplot_dimensions(1),subplot_dimensions(2),[4])
% clearvars y_data
% y_data = mean_response.d1r(2).lateral_preparatory_allcells(flipud(I_d1r.lateral_sorted'),:);
% imagesc(lateral_lick_aligned_time.d1r,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(mean_response.d1r(1).central_preparatory_allcells,1)+size(mean_response.d1r(2).central_preparatory_allcells,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% xlim([-1.5 2]);set(gca,'ytick',[])
% box off
% 
% subplot(subplot_dimensions(1),subplot_dimensions(2),[8])
% clearvars y_data
% y_data = mean_response.non(2).lateral_preparatory_allcells(flipud(I_non.lateral_sorted'),:);
% imagesc(lateral_lick_aligned_time.d1r,1:size(y_data,1),y_data,clims);
% line([0 0],[0 size(y_data,1)],'Color','k','LineWidth',1.5,'LineStyle','--')
% % 
% xlim([-1.5 2]);set(gca,'ytick',[])
% 


set(gcf,'Position',[-1690 42 1465 781]);
%
clc;
clear temp temp1 populations population_time_between
BW = 0.25;
min_percentage = 0.2;
max_percentage = 1;
populations = {'d1r','non'};
search_window = FC_preparation;
figure(13)
for j = 1:2
    for i = 1:size(mean_response.(populations{j})(1).central_preparatory_allcells,1)
        
        p = mean_response.(populations{j})(1).central_preparatory_allcells(i,search_window(1):search_window(end));
        p_time = mean_response.(populations{j})(1).central_preparatory_allcells_time(i,search_window(1):search_window(end)+1);
        temp = (p-min(p))./(max(p)-min(p)) ;
        
        n=min_percentage;
        [~,idx]=min(abs(temp(1,:)-n));
        temp1(i,1) = idx;
        n=max_percentage;
        [~,idx1]=min(abs(temp(1,idx:end)-n));
        temp1(i,2) = idx1+idx;
        temp1(i,3) = (temp1(i,2) -temp1(i,1))*0.135;
        temp1(i,4) = ((p_time(1,temp1(i,2)) - p_time(1,temp1(i,1))));
        clear temp p_time p_time clear idx
        
    end
    
    population_time_between.(populations{j}) = temp1;
    clear temp1 temp
end

clearvars x y bootStrap_sample temp1 temp2 store1 store2
x = population_time_between.d1r(:,4);
y = population_time_between.non(:,4);
t_stat = mean(x) - mean(y);
pooled = [x;y];
for i = 1:1000
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(i) = mean(pooled(temp1)) - mean(pooled(temp2));
end
p_values.time_between_central_preparation_bootstrap = ((length(find(bootStrap_sample<=t_stat))+1)/1001);

subplot(2,4,3)

h1 = cdfplot(population_time_between.non(:,4));hold on;ylabel('')
h2 = cdfplot(population_time_between.d1r(:,4));grid off
set( h2(1), 'LineStyle', '-', 'Color', 'm','LineWidth',1.3);
% text(1, 0.6,num2str(p_values.time_between_central_preparation_bootstrap));

grid off
box off 
set( h1(1), 'LineStyle', '-', 'Color', 'g','LineWidth',1.3);hold on;
title('Ramp time (central)'); xlim([0 2.7])
clearvars x y
x = population_time_between.d1r(:,4);
y = population_time_between.non(:,4);
for i = 1:1000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.centralramptime.d1r(i) = mean(x(rand_idx1,1));
    bootstrap_samples.centralramptime.non(i) = mean(y(rand_idx2,1));

end
hold on;


subplot(2,4,4)


bar(1, mean(bootstrap_samples.centralramptime.d1r),'FaceColor','none','Edgecolor','m','LineWidth',1.1); hold on
bar(2, mean(bootstrap_samples.centralramptime.non),'FaceColor','none','Edgecolor','g','LineWidth',1.1); hold on
errorbar([1 2], [mean(bootstrap_samples.centralramptime.d1r) mean(bootstrap_samples.centralramptime.non)],[std(bootstrap_samples.centralramptime.d1r) std(bootstrap_samples.centralramptime.non)],'LineStyle','none','Color','k')
ylabel('Time (s)')
ylim([0 1.5])
% text(1.25,1.1,['p = ' num2str(p_values.time_between_central_preparation_bootstrap)]);

camroll(90)
set(gca,'YDir','reverse');
box off

clc;
% close all
clear temp temp1 populations population_time_between
BW = 0.25;
populations = {'d1r','non'};
search_window = LL_preparation_early;
for j = 1:2
    for i = 1:size(mean_response.(populations{j})(2).lateral_preparatory_allcells,1)
        
        p = mean_response.(populations{j})(2).lateral_preparatory_allcells(i,search_window(1):search_window(end));
        p_time = mean_response.(populations{j})(2).central_preparatory_allcells_time(i,search_window(1):search_window(end)+1);
        temp = (p-min(p))./(max(p)-min(p)) ;
        
        n=min_percentage;
        [~,idx]=min(abs(temp(1,:)-n));
        temp1(i,1) = idx;
        %     clear idx
        n=max_percentage;
        [~,idx1]=min(abs(temp(1,idx:end)-n));
        temp1(i,2) = idx1+idx;
        clear idx
        temp1(i,3) = (temp1(i,2) -temp1(i,1))*0.135;
        temp1(i,4) = (p_time(1,temp1(i,2)) - p_time(1,temp1(i,1)));
        clear temp
    end
    
    population_time_between.(populations{j}) = temp1;
    clear temp1 temp
end
clearvars x y bootStrap_sample temp1 temp2
x = population_time_between.d1r(:,4);
y = population_time_between.non(:,4);
t_stat = mean(x) - mean(y);
pooled = [x;y];
for i = 1:1000
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(i) = mean(pooled(temp1)) - mean(pooled(temp2));
end
p_values.time_between_lateral_preparation_bootstrap = ((length(find(bootStrap_sample<=t_stat))+1)/1001);

box off

subplot(2,4,7)
h1 = cdfplot(population_time_between.non(:,4));hold on;ylabel('')
h2 = cdfplot(population_time_between.d1r(:,4));grid off
set( h2(1), 'LineStyle', '-', 'Color', 'm','LineWidth',1.3);
% text(1, 0.6,num2str(p_values.time_between_lateral_preparation_bootstrap));
grid off
box off 

set( h1(1), 'LineStyle', '-', 'Color', 'g','LineWidth',1.3);
title('Ramp time (lateral)');
subplot(2,4,8)
for i = 1:1000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.lateralramptime.d1r(i) = mean(x(rand_idx1,1));
    bootstrap_samples.lateralramptime.non(i) = mean(y(rand_idx2,1));
end
hold on;
% text(1.25,1.95,['p = ' num2str(p_values.time_between_lateral_preparation_bootstrap)])
bar(1, mean(bootstrap_samples.lateralramptime.d1r),'FaceColor','none','Edgecolor','m','LineWidth',1.1); hold on
bar(2, mean(population_time_between.non(:,4)),'FaceColor','none','Edgecolor','g','LineWidth',1.1); hold on

errorbar([1 2], [mean(bootstrap_samples.lateralramptime.d1r) mean(bootstrap_samples.lateralramptime.non)],[std(bootstrap_samples.lateralramptime.d1r) std(bootstrap_samples.lateralramptime.non)],'LineStyle','none','Color','k')
ylim([0 3]);
ylabel('Time (s)');
camroll(90)
set(gca,'YDir','reverse');
box off

p_values
set(gcf,'Position',[435 85 1038 473]);


