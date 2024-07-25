%% Figure 5 analyses

close all
cd 'D:\John_data_figure_materials\Fig5'
clearvars -except d1r_neurons non_neurons d1r_resp non_resp
save_figures = 0;
clc;
clear resp_proportions
rng('default')

figure_startup
non_resp.all_left_selective = unique([non_resp.sample_left_selective non_resp.delay_left_selective non_resp.choice_left_selective],'sorted');
non_resp.all_right_selective = unique([non_resp.sample_right_selective non_resp.delay_right_selective non_resp.choice_right_selective],'sorted');
non_resp.all_selective = [non_resp.all_left_selective non_resp.all_right_selective];
non_resp.all_nonselective = setdiff(non_resp.all_responsive,non_resp.all_selective);

d1r_resp.all_left_selective = unique([d1r_resp.sample_left_selective d1r_resp.delay_left_selective d1r_resp.choice_left_selective],'sorted');
d1r_resp.all_right_selective = unique([d1r_resp.sample_right_selective d1r_resp.delay_right_selective d1r_resp.choice_right_selective],'sorted');
d1r_resp.all_selective = [d1r_resp.all_left_selective d1r_resp.all_right_selective];
d1r_resp.all_nonselective = setdiff(d1r_resp.all_responsive,d1r_resp.all_selective);

all_d1r_sample = [length(d1r_resp.sample_left_selective) length(d1r_resp.sample_right_selective)];
all_d1r_delay = [length(d1r_resp.delay_left_selective) length(d1r_resp.delay_right_selective)];
all_d1r_choice = [length(d1r_resp.choice_left_selective) length(d1r_resp.choice_right_selective)];

all_non_sample = [length(non_resp.sample_left_selective) length(non_resp.sample_right_selective)];
all_non_delay = [length(non_resp.delay_left_selective) length(non_resp.delay_right_selective)];
all_non_choice = [length(non_resp.choice_left_selective) length(non_resp.choice_right_selective)];

resp_proportions_A(1,1) = sum(all_d1r_sample)/sum(length(d1r_resp.all_selective));
resp_proportions_A(2,1) = sum(all_d1r_delay)/sum(length(d1r_resp.all_selective));
resp_proportions_A(3,1) = sum(all_d1r_choice)/sum(length(d1r_resp.all_selective));

resp_proportions_A(1,2) = sum(all_non_sample)/sum(length(non_resp.all_selective));
resp_proportions_A(2,2) = sum(all_non_delay)/sum(length(non_resp.all_selective));
resp_proportions_A(3,2) = sum(all_non_choice)/sum(length(non_resp.all_selective));

for i = 1:2
    
    resp_proportions(1,i) = all_d1r_sample(i)/sum(all_d1r_sample);
    resp_proportions(2,i) = all_d1r_delay(i)/sum(all_d1r_delay);
    resp_proportions(3,i) = all_d1r_choice(i)/sum(all_d1r_choice);
    
    resp_proportions_non(1,i) = all_non_sample(i)/sum(all_non_sample);
    resp_proportions_non(2,i) = all_non_delay(i)/sum(all_non_delay);
    resp_proportions_non(3,i) = all_non_choice(i)/sum(all_non_choice);
    resp_proportions_1A(1,i) = (all_d1r_sample(i)+all_d1r_delay(i)+all_d1r_choice(i))/(sum(all_d1r_sample)+sum(all_d1r_delay)+sum(all_d1r_choice));
    resp_proportions_1A(2,i) = (all_non_sample(i)+all_non_delay(i)+all_non_choice(i))/(sum(all_non_sample)+sum(all_non_delay)+sum(all_non_choice));
    
end

%
figure(1);clf;
subplot(1,4,1)

p1 = pie(resp_proportions_A(:,1));
patchHand = findobj(p1, 'Type', 'Patch'); 
patchHand(1).FaceColor = [166,206,227]./255;
patchHand(2).FaceColor = [251,128,114]./255;
patchHand(3).FaceColor = [117,112,179]./255;
set(gca,'xticklabel',{'Sample','Delay', 'Choice'},'FontSize',12);
ylabel({'Fraction of'; 'responsive neurons'});
box off


subplot(1,4,2)

p1 = pie(resp_proportions_A(:,2));
patchHand = findobj(p1, 'Type', 'Patch'); 
patchHand(1).FaceColor = [166,206,227]./255;
patchHand(2).FaceColor = [251,128,114]./255;
patchHand(3).FaceColor = [117,112,179]./255;


resp_proportions_1(2,2) = length(non_resp.all_left_selective)/(length(non_resp.all_left_selective)+length(non_resp.all_right_selective));
resp_proportions_1(2,1) = length(non_resp.all_right_selective)/(length(non_resp.all_left_selective)+length(non_resp.all_right_selective));
resp_proportions_1(1,2) = length(d1r_resp.all_left_selective)/(length(d1r_resp.all_left_selective)+length(d1r_resp.all_right_selective));
resp_proportions_1(1,1) = length(d1r_resp.all_right_selective)/(length(d1r_resp.all_left_selective)+length(d1r_resp.all_right_selective));
figure(10);clf;
subplot(1,2,1)

p1 = pie(resp_proportions_1(1,:));
patchHand = findobj(p1, 'Type', 'Patch'); 
patchHand(1).FaceColor = 'b';
patchHand(2).FaceColor = 'r';

legend('Contra', 'Ipsi','Location','northwest')
set(gca,'xticklabel',{'D1R+','D1R-'},'FontSize',12);
ylabel({'Fraction of'; 'selective neurons'});
title('Overall direction preference');
subplot(1,2,2)
p2 = pie(resp_proportions_1(2,:));
patchHand = findobj(p2, 'Type', 'Patch'); 
patchHand(1).FaceColor = 'b';
patchHand(2).FaceColor = 'r';
legend('Contra', 'Ipsi','Location','northwest')
set(gca,'xticklabel',{'D1R+','D1R-'},'FontSize',12);
ylabel({'Fraction of'; 'selective neurons'});
title('Overall direction preference');

box off
figure(1);
subplot(1,4,3)
plot(2:3,resp_proportions(2:3,1),'ro-','LineWidth',2); hold on;
plot(2:3,resp_proportions(2:3,2),'bo-','LineWidth',2); hold on;
xticks([2 3])
set(gca,'xticklabel',{'Delay', 'Choice'},'FontSize',12);

ylabel({'Fraction of'; 'selective neurons'});
ylim([0 1]);xlim([1 4])
title('D1R+ direction tuning');
box off
subplot(1,4,4)

plot(2:3,resp_proportions_non(2:3,1),'ro-','LineWidth',2); hold on;
plot(2:3,resp_proportions_non(2:3,2),'bo-','LineWidth',2); hold on;
xticks([2 3])
set(gca,'xticklabel',{'Delay', 'Choice'},'FontSize',12);
ylabel({'Fraction of'; 'selective neurons'});
ylim([0 1]);xlim([1 4])
box off
set(gcf,'Position',[420 742 1404 236])
title('D1R- direction tuning');

% filename = ['Selectivity_proportions_230425']; 
% if save_figures == 1
%     saveas(gcf,filename,'pdf');
% end
clear n1 n2 N1 N2 pval_proportions
n1 = length(non_resp.all_left_selective); N1 = (length(non_resp.all_left_selective)+length(non_resp.all_right_selective));
n2 = length(non_resp.all_right_selective); N2 = (length(non_resp.all_left_selective)+length(non_resp.all_right_selective));
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(1),pval_overall_selectivity(1)] = crosstab(x1,x2);
subplot(1,4,2)

    text(2.2, 0.75,{'p = ' num2str(pval_overall_selectivity(1))})
    text(2.2, 0.65,{'chistat = ' num2str(chi2stat(1))})
    


clear n1 n2 N1 N2
n1 = length(d1r_resp.all_left_selective); N1 = (length(d1r_resp.all_left_selective)+length(d1r_resp.all_right_selective));
n2 = length(d1r_resp.all_right_selective); N2 = (length(d1r_resp.all_left_selective)+length(d1r_resp.all_right_selective));
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(2),pval_overall_selectivity(2)] = crosstab(x1,x2);


clear n1 n2 N1 N2
n1 = all_d1r_sample(1); N1 = sum(all_d1r_sample);
n2 = all_d1r_sample(2); N2 = sum(all_d1r_sample);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(1,1),pval_proportions(1,1)] = crosstab(x1,x2);

clear n1 n2 N1 N2
n1 = all_d1r_delay(1); N1 = sum(all_d1r_delay);
n2 = all_d1r_delay(2); N2 = sum(all_d1r_delay);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(1,2),pval_proportions(1,2)] = crosstab(x1,x2);

clear n1 n2 N1 N2
n1 = all_d1r_choice(1); N1 = sum(all_d1r_choice);
n2 = all_d1r_choice(2); N2 = sum(all_d1r_choice);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(1,3),pval_proportions(1,3)] = crosstab(x1,x2);


clear n1 n2 N1 N2
n1 = all_non_sample(1); N1 = sum(all_non_sample);
n2 = all_non_sample(2); N2 = sum(all_non_sample);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(2,1),pval_proportions(2,1)] = crosstab(x1,x2);

clear n1 n2 N1 N2
n1 = all_non_delay(1); N1 = sum(all_non_delay);
n2 = all_non_delay(2); N2 = sum(all_non_delay);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(2,2),pval_proportions(2,2)] = crosstab(x1,x2);

clear n1 n2 N1 N2
n1 = all_non_choice(1); N1 = sum(all_non_choice);
n2 = all_non_choice(2); N2 = sum(all_non_choice);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat(2,3),pval_proportions(2,3)] = crosstab(x1,x2);

%

clear x y x_selectivity y_selectivity x1_selectivity y1_selectivity ipsi_neurons contra_neurons x_selectivity_stand resp_temp1 resp_temp2 d1r_contra
contra_neurons = d1r_neurons(d1r_resp.all_right_selective);
ipsi_neurons = d1r_neurons(d1r_resp.all_left_selective  );

for j =1:size(contra_neurons,2) % for each selective neuron
    resp_temp = contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);
    for k = 1:size(resp_temp,1)
        for n = 1:size(resp_temp(k,:),2)
        if(resp_temp(k,n)<0)
            resp_temp(k,n) =rand([1,1]);
        end
        end
    end

    x_t(j,:) = mean([contra_neurons(j).Right.frames_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.frames_aligned(contra_neurons(j).Left.incorrect==1,:)],1);
    resp_temp1 = [contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.incorrect==1,:)];
    resp_temp2 = [contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.correct==1,:);contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.incorrect==1,:)];

    x_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    d1r_contra.preferred(j,:) = mean(resp_temp1,1);
    d1r_contra.nonpreferred(j,:) = mean(resp_temp2,1);

    clear resp_temp1 resp_temp2
end


clear resp_temp1 resp_temp2 d1r_ipsi
for j =1:size(ipsi_neurons,2) % for each selective neuron
    resp_temp = ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);
    for k = 1:size(resp_temp,1)
        for n = 1:size(resp_temp(k,:),2)
        if(resp_temp(k,n)<0)
            resp_temp(k,n) =rand([1,1]);
        end
        end
    end

    y_t(j,:) = mean([ipsi_neurons(j).Left.frames_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.frames_aligned(ipsi_neurons(j).Right.incorrect==1,:)],1);
    resp_temp1 = [ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.incorrect==1,:)];
    resp_temp2 = [ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.correct==1,:);ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.incorrect==1,:)];

    y_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    d1r_ipsi.preferred(j,:) = mean(resp_temp1,1);
    d1r_ipsi.nonpreferred(j,:) = mean(resp_temp2,1);
        clear resp_temp1 resp_temp2

end

%
clear x1 y1 contra_neurons ipsi_neurons resp_temp1 resp_temp2 non_contra non_ipsi
contra_neurons = non_neurons(non_resp.all_right_selective  );
ipsi_neurons = non_neurons(non_resp.all_left_selective  );

for j =1:size(contra_neurons,2) % for each selective neuron
    x1_t(j,:) = mean([contra_neurons(j).Right.frames_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.frames_aligned(contra_neurons(j).Left.incorrect==1,:)],1);
    resp_temp1 = [contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.incorrect==1,:)];
    resp_temp2 = [contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.correct==1,:);contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.incorrect==1,:)];
    x1_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    non_contra.preferred(j,:) =  mean(resp_temp1,1);
    non_contra.nonpreferred(j,:) = mean(resp_temp2,1);
    clear resp_temp1 resp_temp2
end

for j =1:size(ipsi_neurons,2) % for each selective neuron
    y1_t(j,:) = mean([ipsi_neurons(j).Left.frames_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.frames_aligned(ipsi_neurons(j).Right.incorrect==1,:)],1);
    resp_temp1 = [ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.incorrect==1,:)];
    resp_temp2 = [ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.correct==1,:);ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.incorrect==1,:)];
    y1_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    non_ipsi.preferred(j,:) = mean(resp_temp1,1);
    non_ipsi.nonpreferred(j,:) = mean(resp_temp2,1);
    clear resp_temp1 resp_temp2
end


% 
x_selectivity_norm = x_selectivity./max(x_selectivity,[],2);
y_selectivity_norm = y_selectivity./max(y_selectivity,[],2);
x1_selectivity_norm = x1_selectivity./max(x1_selectivity,[],2);
y1_selectivity_norm = y1_selectivity./max(y1_selectivity,[],2);
%
rng('default')

figure(3);clf;
clims = [-1 1];
subplot(3,2,[1 3])
clear r I
for t = 1:size(x_selectivity_norm,1)
    r(t) = length(find(x_selectivity_norm(t,31:end)>0.5));

end
[~,I_contra]=sort(r);

clear r I
for t = 1:size(y_selectivity_norm,1)
    r(t) = length(find(y_selectivity_norm(t,31:end)<0.5));

end
[~,I_ipsi]=sort(r);
imagesc(mean([x_t;y_t],1),1:size([x_selectivity_norm;y_selectivity_norm],1),[x_selectivity_norm(flipud(I_contra'),:);y_selectivity_norm(flipud(I_ipsi'),:)*-1],clims);hold on
line([0 0],[-1.5 size([x_selectivity_norm;y_selectivity_norm],1)],'Color','k','LineWidth',1.5,'LineStyle','--');

xlim([-7 2]); set(gca,'xTick',[]);
yticks([1 size(x_selectivity_norm,1) size([x_selectivity_norm;y_selectivity_norm],1)]); hold on;

title('D1R+ selective neurons')

box off
%
subplot(3,2,[2,4])
clear r I
for t = 1:size(x1_selectivity_norm,1)
    r(t) = length(find(x1_selectivity_norm(t,31:end)>0.5));

end
[~,I_contra]=sort(r);

clear r I
for t = 1:size(y1_selectivity_norm,1)
    r(t) = length(find(y1_selectivity_norm(t,31:end)<0.5));

end
[~,I_ipsi]=sort(r);
imagesc(mean([x1_t;y1_t],1),1:size([x1_selectivity_norm;y1_selectivity_norm],1),[x1_selectivity_norm(flipud(I_contra'),:);y1_selectivity_norm(flipud(I_ipsi'),:)*-1],clims);hold on;
colormap(flipud(redblue));
line([0 0],[-1.5 size([x1_selectivity_norm;y1_selectivity_norm],1)],'Color','k','LineWidth',1.5,'LineStyle','--');

yticks([1 size(x1_selectivity_norm,1) size([x1_selectivity_norm;y1_selectivity_norm],1)]); hold on;
set(gca,'xTick',[]);
box off

xlim([-7 2])
title('D1R- selective neurons')
h2 = colorbar;

newPosition2 = [0.915306042541622 0.407142857142859 0.0226785714285715 0.514285714285716];
newUnits = 'normalized';
set(h2,'Position', newPosition2,'Units', newUnits);


rectangle_start = -1.5;

subplot(3,2,5)
d1r_selectivity = [x_selectivity;y_selectivity];

rectangle('Position',[rectangle_start -1 1.5 10],'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')
rectangle('Position',[0 -1 1.5 10],'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')

boundedline(mean(x1_t,1),mean(x_selectivity),std(x_selectivity)./(sqrt(size(x_selectivity,1))),'alpha','color','b','LineWidth',1.5);hold on
boundedline(mean(y1_t,1),mean(y_selectivity),std(y_selectivity)./(sqrt(size(y_selectivity,1))),'alpha','color','r','LineWidth',1.5);hold on
line([0 0],[-1.5 10],'Color','k','LineWidth',1.5,'LineStyle','--')
ylabel('Selectivity (\DeltaF/F)'); xlabel('Time (s)')
title('D1R+ selective neurons')

xlim([-7 2])
ylim([-0.2 5]);yticks(0:1:5)

clear d1r_selectivity non_selectivity
subplot(3,2,6)


rectangle('Position',[rectangle_start -1 1.5 10],'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')
rectangle('Position',[0 -1 1.5 10],'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')

boundedline(mean(x1_t,1),mean(x1_selectivity),std(x1_selectivity)./(sqrt(size(x1_selectivity,1))),'alpha','color','b','LineWidth',1.5);hold on
boundedline(mean(y1_t,1),mean(y1_selectivity),std(y1_selectivity)./(sqrt(size(y1_selectivity,1))),'alpha','color','r','LineWidth',1.5);hold on
line([0 0],[-1.5 10],'Color','k','LineWidth',1.5,'LineStyle','--');

ylabel('Selectivity (\DeltaF/F)'); xlabel('Time (s)')
title('D1R- selective neurons')
set(gcf,'Position',[680 343 665 635])


xlim([-7 2])
ylim([-0.2 5]); yticks(0:1:5)


%
figure(2);clf;
m_size = 15;
rng('default')

clear mean_selectivity d1r_mean_selectivity p_values1
delay_window = 57:68;
choice_window = 69:80;
mean_selectivity.contraselectivity(:,2) = mean(x1_selectivity(:,delay_window),2);
mean_selectivity.contraselectivity(:,3) = mean(x1_selectivity(:,choice_window),2);
% 
mean_selectivity.ipsiselectivity(:,2) = mean(y1_selectivity(:,delay_window),2);
mean_selectivity.ipsiselectivity(:,3) = mean(y1_selectivity(:,choice_window),2);


d1r_mean_selectivity.contraselectivity(:,2) = mean(x_selectivity(:,delay_window),2);
d1r_mean_selectivity.contraselectivity(:,3) = mean(x_selectivity(:,choice_window),2);
% 
d1r_mean_selectivity.ipsiselectivity(:,2) = mean(y_selectivity(:,delay_window),2);
d1r_mean_selectivity.ipsiselectivity(:,3) = mean(y_selectivity(:,choice_window),2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
subplot(2,2,1)
clearvars x y bootStrap_sample p_values bootstrap_samples pooled t_stat pooled
x = d1r_mean_selectivity.contraselectivity(:,2);
y = d1r_mean_selectivity.ipsiselectivity(:,2);

ranksum(x,y)

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:1000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_delay_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_delay_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        bootstrap_samples.mean_delay_diff(b) = bootstrap_samples.mean_delay_selectivity.d1r(b) - bootstrap_samples.mean_delay_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/1001;
elseif t_stat > 0
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/1001;
end

clearvars x y bootStrap_sample pooled
x = d1r_mean_selectivity.contraselectivity(:,3);
y = d1r_mean_selectivity.ipsiselectivity(:,3);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:1000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_choice_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_choice_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
    
    bootstrap_samples.mean_choice_diff(b) = bootstrap_samples.mean_choice_selectivity.d1r(b) - bootstrap_samples.mean_choice_selectivity.non(b);
end
if t_stat < 0
    
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/1001;
elseif t_stat > 0
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/1001;
end

plot([1, 1.6],[mean(bootstrap_samples.mean_delay_selectivity.d1r),mean(bootstrap_samples.mean_choice_selectivity.d1r)],'b.','MarkerSize',m_size,'LineWidth',1.5);hold on;
plot([1.2, 1.8],[mean(bootstrap_samples.mean_delay_selectivity.non),mean(bootstrap_samples.mean_choice_selectivity.non)],'r.','MarkerSize',m_size,'LineWidth',1.5);hold on;

e = errorbar([1, 1.6],[mean(d1r_mean_selectivity.contraselectivity(:,2),1),mean(d1r_mean_selectivity.contraselectivity(:,3),1)],...
    [std(d1r_mean_selectivity.contraselectivity(:,2))/sqrt(size(d1r_mean_selectivity.contraselectivity(:,2),1)),...
    std(d1r_mean_selectivity.contraselectivity(:,3))/sqrt(size(d1r_mean_selectivity.contraselectivity(:,3),1))],'LineWidth',1); hold on;
e.LineStyle = 'none';e.Color = 'blue';

e = errorbar([1.2, 1.8],[mean(d1r_mean_selectivity.ipsiselectivity(:,2),1),mean(d1r_mean_selectivity.ipsiselectivity(:,3),1)],...
    [std(d1r_mean_selectivity.ipsiselectivity(:,2))/sqrt(size(d1r_mean_selectivity.ipsiselectivity(:,2),1)),...
    std(d1r_mean_selectivity.ipsiselectivity(:,3))/sqrt(size(d1r_mean_selectivity.ipsiselectivity(:,3),1))],'LineWidth',1); hold on;
e.LineStyle = 'none';e.Color = 'red';


text(1.05,3.5, num2str(p_values.bootstrap.delay_selectivity))

text(1.7,3.5, num2str(p_values.bootstrap.choice_selectivity))
if p_values.bootstrap.delay_selectivity < 0.05
    line([1 1.2], [4 4],'LineWidth',1.5,'Color','k');text(1.1,4.2, '*')
end
if p_values.bootstrap.choice_selectivity < 0.05
    line([1.6 1.8], [4.8 4.8],'LineWidth',1.5,'Color','k');text(1.7,4.9, '*')
end
set(gca,'LineWidth',1);
ylim([0.5 5]);xlim([0.8 2]);yticks([0:1:5])
box off
ylabel('Selectivity (\DeltaF/F)');
xticks([1.1 1.7])
xticklabels({'Delay','Choice'})


figure(2)
subplot(2,2,2)
clearvars x y t_stat pooled rand_idx1 rand_idx2 bootstrap_samples bootStrap_sample p_values
x = mean_selectivity.contraselectivity(:,2);
y = mean_selectivity.ipsiselectivity(:,2);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:1000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_delay_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_delay_selectivity.non(b) = mean(y(rand_idx2,1));
    
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
            bootstrap_samples.mean_delay_diff(b) = bootstrap_samples.mean_delay_selectivity.d1r(b) - bootstrap_samples.mean_delay_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/1001;
elseif t_stat > 0
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/1001;
end
clearvars x y t_stat pooled rand_idx1 rand_idx2 bootStrap_sample  
x = mean_selectivity.contraselectivity(:,3);
y = mean_selectivity.ipsiselectivity(:,3);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:1000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_choice_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_choice_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
                bootstrap_samples.mean_choice_diff(b) = bootstrap_samples.mean_choice_selectivity.d1r(b) - bootstrap_samples.mean_choice_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/1001
elseif t_stat > 0
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/(1001)
end

plot([1,1.6],[mean(bootstrap_samples.mean_delay_selectivity.d1r),mean(bootstrap_samples.mean_choice_selectivity.d1r)],'b.','MarkerSize',m_size,'LineWidth',1.5);hold on;
plot([1.2,1.8],[mean(bootstrap_samples.mean_delay_selectivity.non),mean(bootstrap_samples.mean_choice_selectivity.non)],'r.','MarkerSize',m_size,'LineWidth',1.5);hold on;


e = errorbar([1,1.6],[mean(bootstrap_samples.mean_delay_selectivity.d1r),mean(bootstrap_samples.mean_choice_selectivity.d1r)],...
    [std(bootstrap_samples.mean_delay_selectivity.d1r)/1,...
    std(bootstrap_samples.mean_choice_selectivity.d1r)/1],'LineWidth',1); hold on;

e.LineStyle = 'none';e.Color = 'blue';

e = errorbar([1.2,1.8],[mean(bootstrap_samples.mean_delay_selectivity.non),mean(bootstrap_samples.mean_choice_selectivity.non)],...
    [std(bootstrap_samples.mean_delay_selectivity.non)/1,...
    std(bootstrap_samples.mean_choice_selectivity.non)/1],'LineWidth',1,'LineStyle','--'); hold on;
e.LineStyle = 'none';e.Color = 'red';
if p_values.bootstrap.delay_selectivity < 0.05
    line([1 1.2], [4 4],'LineWidth',1.5,'Color','k');text(1.1,4.2, '*')
end
if p_values.bootstrap.choice_selectivity < 0.05
    line([1.6 1.8], [4.8 4.8],'LineWidth',1.5,'Color','k');text(1.7,4.9, '*')
end

text(0.9,3.5, num2str(p_values.bootstrap.delay_selectivity))

text(1.9,3.5, num2str(p_values.bootstrap.choice_selectivity))



set(gca,'LineWidth',1);
ylim([0.5 5]);xlim([0.8 2]);yticks([0:1:5])
box off
ylabel('Selectivity (\DeltaF/F)')
xticks([1.1 1.7])
xticklabels({'Delay','Choice'})
%
rng('default')

numruns = 1000;

figure(2)
subplot(2,2,3)
clearvars x y t_stat pooled rand_idx1 rand_idx2 bootstrap_samples bootStrap_sample p_values
x = d1r_mean_selectivity.contraselectivity(:,2);
y = mean_selectivity.contraselectivity(:,2);

t_stat = mean(x) - mean(y);
pooled = [x;y];

for b = 1:numruns
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_delay_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_delay_selectivity.non(b) = mean(y(rand_idx2,1));
    
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
                bootstrap_samples.mean_delay_diff(b) = bootstrap_samples.mean_delay_selectivity.d1r(b) - bootstrap_samples.mean_delay_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/(numruns+1);
elseif t_stat > 0
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/(numruns+1);
end

%
figure(2)
clearvars x y t_stat pooled rand_idx1 rand_idx2 bootStrap_sample 
x = d1r_mean_selectivity.contraselectivity(:,3);
y = mean_selectivity.contraselectivity(:,3);
t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:numruns
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_choice_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_choice_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
                    bootstrap_samples.mean_choice_diff(b) = bootstrap_samples.mean_choice_selectivity.d1r(b) - bootstrap_samples.mean_choice_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/(numruns+1)
elseif t_stat > 0
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/(numruns+1)
end

plot([1,1.6],[mean(bootstrap_samples.mean_delay_selectivity.d1r),mean(bootstrap_samples.mean_choice_selectivity.d1r)],'m.','MarkerSize',m_size,'LineWidth',1.5);hold on;
plot([1.2,1.8],[mean(bootstrap_samples.mean_delay_selectivity.non),mean(bootstrap_samples.mean_choice_selectivity.non)],'g.','MarkerSize',m_size,'LineWidth',1.5);hold on;
%

e = errorbar([1,1.6],[mean(d1r_mean_selectivity.contraselectivity(:,2)),mean(d1r_mean_selectivity.contraselectivity(:,3))],...
    [std(d1r_mean_selectivity.contraselectivity(:,2))/sqrt(size(d1r_mean_selectivity.contraselectivity(:,2),1)),...
    std(d1r_mean_selectivity.contraselectivity(:,3))/sqrt(size(d1r_mean_selectivity.contraselectivity(:,3),1))],'LineWidth',1); hold on;
% 
e.LineStyle = 'none';e.Color = 'magenta';
% 
e = errorbar([1.2,1.8],[mean(mean_selectivity.contraselectivity(:,2)),mean(mean_selectivity.contraselectivity(:,3))],...
    [std(mean_selectivity.contraselectivity(:,2))/sqrt(size(mean_selectivity.contraselectivity(:,2),1)),...
    std(mean_selectivity.contraselectivity(:,3))/sqrt(size(mean_selectivity.contraselectivity(:,3),1))],'LineWidth',1); hold on;
e.LineStyle = 'none';e.Color = 'green';


if p_values.bootstrap.delay_selectivity < 0.05
    line([1 1.2], [4 4],'LineWidth',1.5,'Color','k');text(1.1,4.2, '*')
end
if p_values.bootstrap.choice_selectivity < 0.05
    line([1.6 1.8], [4.8 4.8],'LineWidth',1.5,'Color','k');text(1.7,4.9, '*')
end


text(0.9,3.5, num2str(p_values.bootstrap.delay_selectivity))

text(1.9,3.5, num2str(p_values.bootstrap.choice_selectivity))

set(gca,'LineWidth',1);
ylim([0.5 5]);xlim([0.8 2]);yticks([0:1:5])
box off
ylabel('Selectivity (\DeltaF/F)');
xticks([1.1 1.7])
xticklabels({'Delay','Choice'})
title('Between populations - Contra preferring','FontSize',7)
%
figure(2)
subplot(2,2,4)
clearvars x y t_stat pooled rand_idx1 rand_idx2 bootstrap_samples bootStrap_sample p_values
x = d1r_mean_selectivity.ipsiselectivity(:,2);
y = mean_selectivity.ipsiselectivity(:,2);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:numruns
    
    rand_idx1 = randsample(size(x,1),size(x,1),true);
    rand_idx2 = randsample(size(y,1),size(y,1),true);
    bootstrap_samples.mean_delay_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_delay_selectivity.non(b) = mean(y(rand_idx2,1));

    
     temp1 = randsample(size([x;y],1),size(x,1),true);
    temp2 = randsample(size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
            bootstrap_samples.mean_delay_diff(b) = bootstrap_samples.mean_delay_selectivity.d1r(b) - bootstrap_samples.mean_delay_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/(numruns+1)
elseif t_stat > 0
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/(numruns+1);
end



figure(2)

clearvars x y t_stat pooled rand_idx1 rand_idx2 bootStrap_sample 
x = d1r_mean_selectivity.ipsiselectivity(:,3);
y = mean_selectivity.ipsiselectivity(:,3);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:numruns
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_choice_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_choice_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
        bootstrap_samples.mean_choice_diff(b) = bootstrap_samples.mean_choice_selectivity.d1r(b) - bootstrap_samples.mean_choice_selectivity.non(b);

end
if t_stat < 0
    
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/(numruns+1)
elseif t_stat > 0
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/(numruns+1)
end
p_values1.ranksum(8) = ranksum(x,y);

plot([1,1.6],[mean(d1r_mean_selectivity.ipsiselectivity(:,2)),mean(d1r_mean_selectivity.ipsiselectivity(:,3))],'m.','MarkerSize',m_size,'LineWidth',1.5);hold on;
plot([1.2,1.8],[mean(mean_selectivity.ipsiselectivity(:,2)),mean(mean_selectivity.ipsiselectivity(:,3))],'g.','MarkerSize',m_size,'LineWidth',1.5);hold on;


e = errorbar([1,1.6],[mean(d1r_mean_selectivity.ipsiselectivity(:,2)),mean(d1r_mean_selectivity.ipsiselectivity(:,3))],...
    [std(d1r_mean_selectivity.ipsiselectivity(:,2))/sqrt(size(d1r_mean_selectivity.ipsiselectivity(:,2),1)),...
    std(d1r_mean_selectivity.ipsiselectivity(:,3))/sqrt(size(d1r_mean_selectivity.ipsiselectivity(:,3),1))],'LineWidth',1); hold on;
% 
e.LineStyle = 'none';e.Color = 'magenta';
% 
e = errorbar([1.2,1.8],[mean(mean_selectivity.ipsiselectivity(:,2)),mean(mean_selectivity.ipsiselectivity(:,3))],...
    [std(mean_selectivity.ipsiselectivity(:,2))/sqrt(size(mean_selectivity.ipsiselectivity(:,2),1)),...
    std(mean_selectivity.ipsiselectivity(:,3))/sqrt(size(mean_selectivity.ipsiselectivity(:,3),1))],'LineWidth',1); hold on;
e.LineStyle = 'none';e.Color = 'green';


text(0.9,3.5, num2str(p_values.bootstrap.delay_selectivity))

text(1.9,3.5, num2str(p_values.bootstrap.choice_selectivity))
if p_values.bootstrap.delay_selectivity < 0.05
    line([1 1.2], [4 4],'LineWidth',1.5,'Color','k');text(1.1,4.2, '*')
end
if p_values.bootstrap.choice_selectivity < 0.05
    line([1.6 1.8], [4.8 4.8],'LineWidth',1.5,'Color','k');text(1.7,4.9, '*')
end
set(gca,'LineWidth',1);
ylim([0.5 5]);xlim([0.8 2]);yticks([0:1:5])
box off
ylabel('Selectivity (\DeltaF/F)');
xticks([1.1 1.7])
box off
ylabel('Selectivity (\DeltaF/F)')
xticklabels({'Delay','Choice'})
title('Between populations - Ipsi preferring','FontSize',7)



set(gcf, 'Renderer', 'painters');
set(gcf,'Position',[460 411 499 506])


%
% close all
number_of_change_points = 1;
start_time = 56;
end_time = 70;
clims = [-2 2];
clc
clear ipt ipt_1_d1r time_of_onset_d1r_contra  p_values time_of_peak_d1r_contra
for i = 1:size(x_selectivity,1)
   
    if isempty(findchangepts(x_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==1
        ipt_1_d1r(i) = NaN;
        time_of_onset_d1r_contra(i)= NaN;

    elseif isempty(findchangepts(x_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==0
        [ipt{i}] = findchangepts(x_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear');
        time_of_onset_d1r_contra(i)= x_t(i,ipt{1,i}(1)+(start_time-1));
    end
        [peak_selec.d1r.contra(i),ind] = max( x_selectivity(i,start_time:end_time));
         time_of_peak_d1r_contra(i) = x_t(i,ind+(start_time-1));
         clearvars ind
end

clear ipt ipt_1_d1r time_of_onset_d1r_ipsi time_of_peak_d1r_ipsi
for i = 1:size(y_selectivity,1)
   
    if isempty(findchangepts(y_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==1
        ipt_1_d1r(i) = NaN;
        time_of_onset_d1r_ipsi(i)= NaN;

    elseif isempty(findchangepts(y_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==0
        [ipt{i}] = findchangepts(y_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear');
        time_of_onset_d1r_ipsi(i)= y_t(i,ipt{1,i}(1)+(start_time-1));
    end
    [peak_selec.d1r.ipsi(i),ind] = max( y_selectivity(i,start_time:end_time));
    time_of_peak_d1r_ipsi(i) = y_t(i,ind+(start_time-1));
    clearvars ind
end

clear ipt ipt_1_d1r time_of_onset_non_contra time_of_peak_non_contra
for i = 1:size(x1_selectivity,1)
   
    if isempty(findchangepts(x1_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==1
        ipt_1_d1r(i) = NaN;
        time_of_onset_non_contra(i)= NaN;

    elseif isempty(findchangepts(x1_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==0
        [ipt{i}] = findchangepts(x1_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear');
        time_of_onset_non_contra(i)= x1_t(i,ipt{1,i}(1)+(start_time-1));
    end
    [peak_selec.non.contra(i),ind] = max( x1_selectivity(i,start_time:end_time));
    time_of_peak_non_contra(i) = x1_t(i,ind+(start_time-1));
    clearvars ind
end

clear ipt ipt_1_d1r time_of_onset_non_ipsi time_of_peak_non_ipsi
for i = 1:size(y1_selectivity,1)
   
    if isempty(findchangepts(y1_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==1
        ipt_1_d1r(i) = NaN;
        time_of_onset_non_ipsi(i)= NaN;

    elseif isempty(findchangepts(y1_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear')) ==0
        [ipt{i}] = findchangepts(y1_selectivity(i,start_time:end_time),'MaxNumChanges',number_of_change_points,'Statistic','linear');
        time_of_onset_non_ipsi(i)= y1_t(i,ipt{1,i}(1)+(start_time-1));
    end
    [peak_selec.non.ipsi(i),ind] = max( y1_selectivity(i,start_time:end_time));
    time_of_peak_non_ipsi(i) = y1_t(i,ind+(start_time-1));
    clearvars ind
end

[fi1,xi1,U1] = ksdensity(time_of_onset_d1r_contra,'Function','pdf');
[fi2,xi2,U2] = ksdensity(time_of_onset_d1r_ipsi,'Function','pdf');
BW = 0.25;

figure(4);
subplot(1,2,1)
plot(xi1,fi1*U1,'m','LineWidth',1.3); hold on;
subplot(1,2,2)

plot(xi2,fi2*U2,'g','LineWidth',1.3); hold on;
subplot(1,2,1)

line([mean(time_of_onset_d1r_contra,'omitnan') mean(time_of_onset_d1r_contra,'omitnan')],[0 0.2],'Color','m','LineStyle','-')
subplot(1,2,2)

line([mean(time_of_onset_d1r_ipsi,'omitnan') mean(time_of_onset_d1r_ipsi,'omitnan')],[0 0.2],'Color','g','LineStyle','-')

%
subplot(1,2,1)

[fi1,xi1,U1] = ksdensity(time_of_onset_non_contra,'Function','pdf');
subplot(1,2,2)

[fi2,xi2,U2] = ksdensity(time_of_onset_non_ipsi,'Function','pdf');
subplot(1,2,1)

plot(xi1,(fi1*U1),'m--','LineWidth',1.3); hold on;
subplot(1,2,2)

plot(xi2,fi2*U2,'g--','LineWidth',1.3); hold on;
subplot(1,2,1)

line([mean(time_of_onset_non_contra,'omitnan') mean(time_of_onset_non_contra,'omitnan')],[0 0.2],'Color','m','LineStyle','--')
box off
ylabel('Onset probability');xlabel('Time (s)')

subplot(1,2,2)

line([mean(time_of_onset_non_ipsi,'omitnan') mean(time_of_onset_non_ipsi,'omitnan')],[0 0.2],'Color','g','LineStyle','--')


BW = 0.25;
box off
ylabel('Onset probability');xlabel('Time (s)')


rng ('default')
clearvars x y bootStrap_sample
x = time_of_onset_d1r_ipsi';
y = time_of_onset_non_ipsi';

        t_stat = mean(x,'omitnan') - mean(y,'omitnan');
        
        
        
        pooled = [x;y];
        for b = 1:1000
            temp1 = randsample(1:size([x;y],1),size(x,1),true);
            temp2 = randsample(1:size([x;y],1),size(y,1),true);
            bootStrap_sample(b) = mean(pooled(temp1),'omitnan') - mean(pooled(temp2),'omitnan');
        end
        if t_stat > 0
        p_values.bootstrap_onset_ipsi = (length(find(bootStrap_sample>=t_stat))+1)/(1000+1)
        elseif t_stat < 0
                    p_values.bootstrap_onset_ipsi = (length(find(bootStrap_sample<=t_stat))+1)/(1000+1)
        end

    clearvars x y bootStrap_sample
x = time_of_onset_d1r_contra';
y = time_of_onset_non_contra';

        t_stat = mean(x,'omitnan') - mean(y,'omitnan');
        
        
        
        pooled = [x;y];
        for b = 1:1000
            temp1 = randsample(1:size([x;y],1),size(x,1),true);
            temp2 = randsample(1:size([x;y],1),size(y,1),true);
            bootStrap_sample(b) = mean(pooled(temp1),'omitnan') - mean(pooled(temp2),'omitnan');
        end
        if t_stat > 0
        p_values.bootstrap_onset_contra = (length(find(bootStrap_sample>=t_stat))+1)/(1000+1)
        elseif t_stat < 0
                    p_values.bootstrap_onset_contra= (length(find(bootStrap_sample<=t_stat))+1)/(1000+1)
        end
        
        
        
        
%%
clearvars -except d1r_neurons non_neurons d1r_resp non_resp

clear x y x_selectivity y_selectivity x1_selectivity y1_selectivity ipsi_neurons contra_ne urons x_selectivity_stand resp_temp1 resp_temp2 d1r_contra
contra_neurons = d1r_neurons(d1r_resp.all_selective  );
ipsi_neurons = d1r_neurons(d1r_resp.all_selective  );

for j =1:size(contra_neurons,2) % for each selective neuron
    resp_temp = contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);
    for k = 1:size(resp_temp,1)
        for n = 1:size(resp_temp(k,:),2)
        if(resp_temp(k,n)<0)
            resp_temp(k,n) =rand([1,1]);
        end
        end
    end

    x_t(j,:) = mean([contra_neurons(j).Right.frames_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.frames_aligned(contra_neurons(j).Left.incorrect==1,:)],1);
    resp_temp1 = [contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.incorrect==1,:)];
    resp_temp2 = [contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.correct==1,:);contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.incorrect==1,:)];

    x_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    d1r_contra.preferred(j,:) = mean(resp_temp1,1);
    d1r_contra.nonpreferred(j,:) = mean(resp_temp2,1);

    clear resp_temp1 resp_temp2
end


clear resp_temp1 resp_temp2 d1r_ipsi
for j =1:size(ipsi_neurons,2) % for each selective neuron
    resp_temp = ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);
    for k = 1:size(resp_temp,1)
        for n = 1:size(resp_temp(k,:),2)
        if(resp_temp(k,n)<0)
            resp_temp(k,n) =rand([1,1]);
        end
        end
    end
%         y(j,:) = mean(resp_temp,1);

    y_t(j,:) = mean([ipsi_neurons(j).Left.frames_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.frames_aligned(ipsi_neurons(j).Right.incorrect==1,:)],1);
    resp_temp1 = [ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.incorrect==1,:)];
    resp_temp2 = [ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.correct==1,:);ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.incorrect==1,:)];

    y_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    d1r_ipsi.preferred(j,:) = mean(resp_temp1,1);
    d1r_ipsi.nonpreferred(j,:) = mean(resp_temp2,1);
        clear resp_temp1 resp_temp2

end

%
clear x1 y1 contra_neurons ipsi_neurons resp_temp1 resp_temp2 non_contra non_ipsi
contra_neurons = non_neurons(non_resp.all_selective  );
ipsi_neurons = non_neurons(non_resp.all_selective  );

for j =1:size(contra_neurons,2) % for each selective neuron
    x1_t(j,:) = mean([contra_neurons(j).Right.frames_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.frames_aligned(contra_neurons(j).Left.incorrect==1,:)],1);
    resp_temp1 = [contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.correct==1,:);contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.incorrect==1,:)];
    resp_temp2 = [contra_neurons(j).Left.dff_aligned(contra_neurons(j).Left.correct==1,:);contra_neurons(j).Right.dff_aligned(contra_neurons(j).Right.incorrect==1,:)];
    x1_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    non_contra.preferred(j,:) =  mean(resp_temp1,1);
    non_contra.nonpreferred(j,:) = mean(resp_temp2,1);
    clear resp_temp1 resp_temp2
end

for j =1:size(ipsi_neurons,2) % for each selective neuron
    y1_t(j,:) = mean([ipsi_neurons(j).Left.frames_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.frames_aligned(ipsi_neurons(j).Right.incorrect==1,:)],1);
    resp_temp1 = [ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.correct==1,:);ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.incorrect==1,:)];
    resp_temp2 = [ipsi_neurons(j).Right.dff_aligned(ipsi_neurons(j).Right.correct==1,:);ipsi_neurons(j).Left.dff_aligned(ipsi_neurons(j).Left.incorrect==1,:)];
    y1_selectivity(j,:) = (mean(resp_temp1,1)-mean(resp_temp2,1))/1;
    non_ipsi.preferred(j,:) = mean(resp_temp1,1);
    non_ipsi.nonpreferred(j,:) = mean(resp_temp2,1);
    clear resp_temp1 resp_temp2
end


%
x_selectivity_norm = x_selectivity./max(abs(x_selectivity),[],2);
y_selectivity_norm = y_selectivity./max(abs(y_selectivity),[],2);
x1_selectivity_norm = x1_selectivity./max(abs(x1_selectivity),[],2);
y1_selectivity_norm = y1_selectivity./max(abs(y1_selectivity),[],2);
%
rng('default')

figure(2);clf;



rectangle_start = -1.5;

subplot(1,2,1)

rectangle('Position',[rectangle_start -1 1.5 10],'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')
rectangle('Position',[0 -1 1.5 10],'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')

boundedline(mean(x_t,1),mean(x_selectivity_norm),std(x_selectivity_norm)./(sqrt(size(x_selectivity_norm,1))),'alpha','color','m','LineWidth',1.5);hold on
boundedline(mean(x1_t,1),mean(x1_selectivity_norm),std(x1_selectivity_norm)./(sqrt(size(x1_selectivity_norm,1))),'alpha','color','g','LineWidth',1.5);hold on

line([0 0],[-1.5 10],'Color','k','LineWidth',1.5,'LineStyle','--')
line([0 0],[-1.5 10],'Color','k','LineWidth',1.5,'LineStyle','--');line([-6 2],[0 0],'Color','k','LineWidth',1.5,'LineStyle','--');

ylabel('Normalized Selectivity (\DeltaF/F)'); xlabel('Time (s)')


xlim([-7 2])
ylim([-0.1 0.3]);yticks([-0.3:0.1:0.3])


%
clear mean_selectivity d1r_mean_selectivity
delay_window = 57:68;
choice_window = 69:80;
mean_selectivity.contraselectivity(:,2) = mean(x1_selectivity_norm(:,delay_window),2);
mean_selectivity.contraselectivity(:,3) = mean(x1_selectivity_norm(:,choice_window),2);
mean_selectivity.ipsiselectivity(:,2) = mean(y1_selectivity(:,delay_window),2);
mean_selectivity.ipsiselectivity(:,3) = mean(y1_selectivity(:,choice_window),2);

d1r_mean_selectivity.contraselectivity(:,2) = mean(x_selectivity_norm(:,delay_window),2);
d1r_mean_selectivity.contraselectivity(:,3) = mean(x_selectivity_norm(:,choice_window),2);
% 
% d1r_mean_selectivity.ipsiselectivity(:,1) = mean(y_selectivity(:,35:43),2);
d1r_mean_selectivity.ipsiselectivity(:,2) = mean(y_selectivity(:,delay_window),2);
d1r_mean_selectivity.ipsiselectivity(:,3) = mean(y_selectivity(:,choice_window),2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
clearvars x y bootStrap_sample p_values bootstrap_samples pooled t_stat pooled
x = d1r_mean_selectivity.contraselectivity(:,2);
y = mean_selectivity.contraselectivity(:,2);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:10000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_delay_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_delay_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
end
if t_stat < 0
    
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/10001;
elseif t_stat > 0
    p_values.bootstrap.delay_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/10001;
end
% p_values.ranksum(i) = ranksum(x,y);

clearvars x y bootStrap_sample pooled
x = d1r_mean_selectivity.contraselectivity(:,3);
y = mean_selectivity.contraselectivity(:,3);

t_stat = mean(x) - mean(y);
pooled = [x;y];
for b = 1:10000
    
    rand_idx1 = randsample(1:size(x,1),size(x,1),true);
    rand_idx2 = randsample(1:size(y,1),size(y,1),true);
    bootstrap_samples.mean_choice_selectivity.d1r(b) = mean(x(rand_idx1,1));
    bootstrap_samples.mean_choice_selectivity.non(b) = mean(y(rand_idx2,1));
    temp1 = randsample(1:size([x;y],1),size(x,1),true);
    temp2 = randsample(1:size([x;y],1),size(y,1),true);
    bootStrap_sample(b) = mean(pooled(temp1)) - mean(pooled(temp2));
end
if t_stat < 0
    
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample<=t_stat))+1)/10001;
elseif t_stat > 0
    p_values.bootstrap.choice_selectivity = (length(find(bootStrap_sample>=t_stat))+1)/10001;
end
line([0,2],[0, 0],'Color','k','LineStyle','--')
hold on;

plot([0.8,1.4],[mean(bootstrap_samples.mean_delay_selectivity.d1r),mean(bootstrap_samples.mean_choice_selectivity.d1r)],'m.','MarkerSize',15,'LineWidth',1.5);hold on;
plot([1,1.6],[mean(bootstrap_samples.mean_delay_selectivity.non),mean(bootstrap_samples.mean_choice_selectivity.non)],'g.','MarkerSize',15,'LineWidth',1.5);hold on;
e = errorbar([0.8,1.4],[mean(bootstrap_samples.mean_delay_selectivity.d1r),mean(bootstrap_samples.mean_choice_selectivity.d1r)],...
    [std(bootstrap_samples.mean_delay_selectivity.d1r)/1,...
    std(bootstrap_samples.mean_choice_selectivity.d1r)/1],'LineWidth',1); hold on;

e.LineStyle = 'none';e.Color = 'magenta';

e = errorbar([1,1.6],[mean(bootstrap_samples.mean_delay_selectivity.non),mean(bootstrap_samples.mean_choice_selectivity.non)],...
    [std(bootstrap_samples.mean_delay_selectivity.non)/1,...
    std(bootstrap_samples.mean_choice_selectivity.non)/1],'LineWidth',1); hold on;

e.LineStyle = 'none';e.Color = 'green';
text(0.7,mean(bootstrap_samples.mean_delay_selectivity.non), num2str(p_values.bootstrap.delay_selectivity))

text(1.7,mean(bootstrap_samples.mean_choice_selectivity.non), num2str(p_values.bootstrap.choice_selectivity))
set(gca,'LineWidth',1);
ylim([-0.05 0.2]);
xlim([0.5 1.9]);
box off
ylabel('Selectivity (\DeltaF/F)');legend('','D1R+','D1R-','Location','northwest')
xticks([0.9 1.5]);yticks([-0.05:0.05:0.2])

xticklabels({'Delay','Choice'})
set(gcf, 'Renderer', 'painters');
set(gcf,'Position',[137 449 1687 529])

%% Plot example neurons
close all
clear temp_mean1 temp_mean2 time_temp1 time_temp2 sem1 sem2 all_trials_sorted trial_temp lick_temp1 lick_temp2 t
clear d1r_example_neurons non_example_neurons group
 
% group = {'sample_left_selective','sample_right_selective','delay_left_selective','delay_right_selective','choice_left_selective','choice_right_selective'};
group = {'sample_taste_selective_sucrose','delay_right_selective','choice_left_selective','choice_right_selective'};

line_width = 1;
figure_startup
for z = 42%:45
for j = 4%:4
    clear time_set
%     if strcmp(group{j},'sample_left_selective')
%         d1r_example_neurons = [1 0];
%         non_example_neurons = [9 10]; % 7 8 
%         all_example_neurons = [d1r_example_neurons;non_example_neurons];
%         time_set = alltimes_d1r.ipsi_pref_activity_time;
%         selec_group = 'all_left_selective';
%         direc = 3;
%         rectangle_start = -1;
%     elseif strcmp(group{j},'sample_right_selective')
%         d1r_example_neurons = [1 0];
%         non_example_neurons = [9 10]; % 7 8 
%         all_example_neurons = [d1r_example_neurons;non_example_neurons];
%         time_set = alltimes_d1r.contra_pref_activity_time;
%         selec_group = 'all_right_selective';
%         direc = 2;
%         rectangle_start = -1;
    if strcmp(group{j},'sample_taste_selective_sucrose')
        d1r_example_neurons = [1 2]; % 11
        non_example_neurons = [1 2]; % 7 8 14
        all_example_neurons = [d1r_example_neurons;non_example_neurons];
%         time_set = alltimes_d1r.ipsi_pref_activity_time;
        selec_group = 'all_left_selective';
        direc = 3;
        rectangle_start = -1;
    elseif strcmp(group{j},'delay_right_selective')
        d1r_example_neurons = [9 4];   %d1r+: 5 4 d1r-:
        non_example_neurons = [96 36];    % %% 11 12
        all_example_neurons = [d1r_example_neurons;non_example_neurons];
%         time_set = alltimes_d1r.contra_pref_activity_time;
        selec_group = 'all_right_selective';
        direc = 1;
        rectangle_start = -1;
    elseif strcmp(group{j},'choice_left_selective')
        d1r_example_neurons = [9 13];   %d1r+: 4 3 9 13
        non_example_neurons = [20 25];    %d1r-: 3 1 7 6
        all_example_neurons = [d1r_example_neurons;non_example_neurons];
%         time_set = alltimes_d1r.ipsi_pref_activity_time;
        selec_group = 'all_left_selective';
        direc = 3;
        rectangle_start = 0;
    elseif strcmp(group{j},'choice_right_selective')
        d1r_example_neurons = [14 8];   %d1r+: 4 7 8 14 16 17 18
        non_example_neurons = [50 51];    %d1r-: 2
        all_example_neurons = [d1r_example_neurons;non_example_neurons];
%         time_set = alltimes_d1r.contra_pref_activity_time;
        selec_group = 'all_right_selective';
        direc = 1;
        rectangle_start = 0;
    end
    for k = 1%1:2 % population type
        clear cell_type
        cell_type = all_example_neurons(k,:);
        clear population
        if k == 1
            population = d1r_neurons(d1r_resp.(group{j}));
        else
            population = non_neurons(non_resp.(group{j}));
        end
        
        
        figure;clf;
        for i = 1:size(cell_type,2)
            clear x_time x_time_temp temp_resp

            
            temp_mean1 = population(cell_type(i)).Left.dff_aligned(population(cell_type(i)).Left.correct==1,:);
            temp_mean2 = population(cell_type(i)).Right.dff_aligned(population(cell_type(i)).Right.correct==1,:);
            
            time_temp1 = population(cell_type(i)).Left.frames_aligned(population(cell_type(i)).Left.correct==1,:);
            time_temp2 = population(cell_type(i)).Right.frames_aligned(population(cell_type(i)).Right.correct==1,:);
            
            lick_temp1 = population(cell_type(i)).Left.central_licks(population(cell_type(i)).Left.correct==1,:);
            lick_temp2 = population(cell_type(i)).Right.central_licks(population(cell_type(i)).Right.correct==1,:);
            
            p1 = population(cell_type(i)).direction_selectivity.delay_p;
            p2 = population(cell_type(i)).direction_selectivity.choice_p;

            sem1 = std(temp_mean1)./sqrt(size(temp_mean1,1));
            sem2 = std(temp_mean2)./sqrt(size(temp_mean2,1));
            
            clear r I t
            for t = 1:size(temp_mean1,1)
                r(t) = length(find(temp_mean1(t,1:end)>0.1));
            end
            [~,I]=sort(r);
            temp_mean1 = temp_mean1(flipud(I'),:);
            lick_temp1 = lick_temp1(flipud(I'),:);
            
            clear r I t
            for t = 1:size(temp_mean2,1)
                r(t) = length(find(temp_mean2(t,1:end)>0.5));
            end
            [~,I]=sort(r);
            temp_mean2 = temp_mean2(flipud(I'),:);
            lick_temp2 = lick_temp2(flipud(I'),:);
            
            all_trials_sorted = [temp_mean1;temp_mean2];
            all_lick_sorted = [lick_temp1(:,1);lick_temp2(:,1)];
            
            trial_temp(1,:) = mean(temp_mean1);
            trial_temp(2,:) = mean(temp_mean2);
            clear ax
            ax(1) = subplot(3,size(cell_type,2),i);
            imagesc(mean([time_temp1;time_temp2],1),1:size(all_trials_sorted,1),all_trials_sorted,[min(min(trial_temp)) max(max(trial_temp))+5]);hold on;
J = customcolormap([0 0.5 1], [255,255,255;150,150,150; 0,0,0]./255);
        yticks([1 size(all_trials_sorted,1)]); hold on;

            colormap(ax(1),J);
            h2 = colorbar;
%                 colorbar('YTick', [round(min(min(trial_temp))-0.3,0) round(max(max(trial_temp))+max(max([sem1;sem2])),0)]);

            if i == 1
                newPosition2 = [0.492557989690722 0.708515283842793 0.0245312500000001 0.215065502183406];
            else 
                newPosition2 = [0.944120489690722 0.708515283842794 0.0245312500000001 0.215065502183406];
            end
            newUnits = 'normalized';
%             
            set(h2,'Position', newPosition2,'Units', newUnits);
%                             colorbar('YTick', [round(min(min(trial_temp))-0.3,0) round(max(max(trial_temp))+max(max([sem1;sem2])),0)]);

            scatter(all_lick_sorted(:,1),1:size(all_lick_sorted,1),15,'g','filled'); hold on;
            line([0 0],[0 size(all_trials_sorted,1)+1],'color',[0 0 0],'LineWidth',2); hold on;
            line([-6.75 -6.75],[0 size(temp_mean1,1)],'color','r','LineWidth',3); hold on;
            line([-6.75 -6.75],[size(temp_mean1,1) size(temp_mean1,1)+size(temp_mean2,1)+1],'color','b','LineWidth',5); hold on;
            xlim([-7 2]);
            title({'P val delay', num2str(p1), 'P val choice', num2str(p2)});
            box off
            
            ax(2) = subplot(3,size(cell_type,2),i+2);
            line([0 0],[min(min(trial_temp))-0.3 max(max(trial_temp))+10],'color',[0.7 0.7 0.7],'LineWidth',line_width,'LineStyle','--'); hold on;
            rectangle('Position',[rectangle_start min(min(trial_temp))-0.3 1 max(max(trial_temp))+10],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
            boundedline(mean(time_temp1),trial_temp(1,:), sem1,'alpha','r','LineWidth',line_width); hold on;
            boundedline(mean(time_temp2),trial_temp(2,:), sem2,'alpha','b','LineWidth',line_width); hold on;
            xlim([-7 2]);
                    yticks([round(min(min(trial_temp))-0.3,0) round(max(max(trial_temp))+max(max([sem1;sem2])),0)]); hold on;
            ylim([round(min(min(trial_temp))-0.8,0) round(max(max(trial_temp))+max(max([sem1;sem2])),0)]);

            ylabel('\Delta F/F')
            xlabel('Time from lateral lick (s)')
            
%             ax(3) = subplot(3,size(cell_type,2),i+4);
%             line([0 0],[min(temp_resp) max(temp_resp)],'color',[0.7 0.7 0.7],'LineWidth',line_width,'LineStyle','--'); hold on;
%             rectangle('Position',[rectangle_start min(temp_resp) 1 max(temp_resp)+0.5],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
%             plot(x_time_temp,temp_resp,'color','k','LineWidth',line_width); hold on;
%             boundedline(mean(time_temp2),trial_temp(2,:), sem2,'alpha','b'); hold on;
%             xlim([-7 3]);
%             box off
%             ylabel('Direction Index')
%             xlabel('Time from lateral lick (s)')
            
            clear temp_mean1 temp_mean2 time_temp1 time_temp2 sem1 sem2 all_trials_sorted trial_temp lick_temp1 lick_temp2 t p1 p2
        end
        clear temp_mean1 temp_mean2 time_temp1 time_temp2 sem1 sem2 all_trials_sorted trial_temp lick_temp1 lick_temp2 t p1 p2
%         linkaxes(ax,'x');
        set(gcf,'Position',[1121 150 694 744]);
        clear filename
            
        filename = [(group{j}) ' ' num2str(k) ' ' num2str(all_example_neurons(k,:))]; 
        if save_figures == 1
            saveas(gcf,filename,'pdf');
        end
        pause(2)
    end
    clear temp_mean1 temp_mean2 time_temp1 time_temp2 sem1 sem2 all_trials_sorted trial_temp lick_temp1 lick_temp2 t p1 p2
end
end



