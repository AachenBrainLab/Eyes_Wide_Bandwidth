%% Here I load all the Directionality results and estimate the Orientation distribution of the total MC generation example
close all; clear all; clc;
%%% Gabor model comparisons from figure 3 
% Step 1: get the processed data from Zenodo [http://doi.org/10.5281/zenodo.14605879]. 
    % The folders low-bandwidth, mid_bandwidth and high_bandiwdth contain frames from an example rendering 
    % of our MC stimuli with the Orientation BW parameters used for this project. 
    % 
    % Optional
    % Orientation_analysis.ijm has been used (can be tested again if desired) to generate the respective .csv containg the orientation
    % distribution analysis for each frame using the Local gradient orientation method in Directionality Fiji plugin.
    % (replace dir1 with the path to the folder you want to analyze and dir2 where you want to save the output)
 
%% Narrow BW example

D = ('C:\Users\balla\Documents\@MC_PAPER\Figure1\panel_b\narrow_bandwidth');
Total_trial_frames = dir(fullfile(D,'*.csv')); 
alpha = numel(Total_trial_frames);
Total_trials = zeros(180,alpha);
for k = 1:alpha
    F = fullfile(D,Total_trial_frames(k).name);
    Total_trial_frames(k).data = readmatrix(F);
    Total_trials(:,k) =  Total_trial_frames(k).data(:,2);
end

n_dist = imgaussfilt(mean(Total_trials,2),3);
n_dist_std = imgaussfilt(std(Total_trials,0,2),3);
narrow_dist = n_dist/max(n_dist);
narrow_dist_sd  = n_dist_std /max(n_dist_std );
narrow_dist_sem = narrow_dist_sd/sqrt(length(narrow_dist_sd));
y = (Total_trial_frames(1).data(:,1) -90)';
clearvars -except narrow_dist narrow_dist_sem alpha y;

%% Mid BW example 

D = ('C:\Users\balla\Documents\@MC_PAPER\Figure1\panel_b\mid_bandwidth');
Total_trial_frames = dir(fullfile(D,'*.csv')); 
Total_trials = zeros(180,alpha);
for k = 1:alpha
    F = fullfile(D,Total_trial_frames(k).name);
    Total_trial_frames(k).data = readmatrix(F);
    Total_trials(:,k) =  Total_trial_frames(k).data(:,2);
end
m_dist = imgaussfilt(mean(Total_trials,2),3);
m_dist_std = imgaussfilt(std(Total_trials,0,2),3);
mid_dist = m_dist/max(m_dist);
mid_dist_sd  = m_dist_std /max(m_dist_std );
mid_dist_sem = mid_dist_sd/sqrt(length(mid_dist_sd));

clearvars -except narrow_dist narrow_dist_sem mid_dist mid_dist_sem alpha y;

%% Broad BW example
D = ('C:\Users\balla\Documents\@MC_PAPER\Figure1\panel_b\broad_bandwidth');
Total_trial_frames = dir(fullfile(D,'*.csv')); 
Total_trials = zeros(180,alpha);
for k =1 : alpha
    F = fullfile(D,Total_trial_frames(k).name);
    Total_trial_frames(k).data = readmatrix(F);
    Total_trials(:,k) =  Total_trial_frames(k).data(:,2);
end
b_dist =  imgaussfilt(mean(Total_trials,2),3);
b_dist_std = imgaussfilt(std(Total_trials,0,2),3);
broad_dist = b_dist/max(b_dist);
broad_dist_sd  = b_dist_std /max(b_dist_std );
broad_dist_sem = broad_dist_sd/sqrt(length(broad_dist_sd));
clearvars -except narrow_dist narrow_dist_sem mid_dist mid_dist_sem broad_dist broad_dist_sem y;


%% Plotting the MC distribution plot with SEM error shading
x2 = [y, fliplr(y)];
narrow_sem = [(narrow_dist + narrow_dist_sem)', fliplr((narrow_dist - narrow_dist_sem)')];
mid_sem = [(mid_dist + mid_dist_sem)', fliplr((mid_dist - mid_dist_sem)')];
broad_sem = [(broad_dist + broad_dist_sem)', fliplr((broad_dist - broad_dist_sem)')];
h1 =fill(x2, narrow_sem, [120/256 ,  197/256   ,  239/256 ]); hold on;   
plot(y,narrow_dist,'LineWidth', 2,'color',[120/256 ,  197/256   ,  239/256 ]); hold on;   
h2 =fill(x2, mid_sem, [81/256 ,  160/256   ,  213/256]); hold on;   
plot(y,mid_dist,'LineWidth', 2, 'color',[81/256 ,  160/256   ,  213/256]); hold on;   
h3 = fill(x2, broad_sem, [44/256 ,  82/256  ,  140/256 ]); hold on;
plot(y,broad_dist,'LineWidth', 2,'color',[44/256 ,  82/256  ,  140/256 ]); hold on;     
box off;
set(h1(1),'Edgecolor',[120/256 ,  197/256   ,  239/256 ],'Facealpha',0.5,'EdgeAlpha',0.5);
set(h2(1),'Edgecolor',[81/256 ,  160/256   ,  213/256],'Facealpha',0.5,'EdgeAlpha',0.5);
set(h3(1),'Edgecolor',[44/256 ,  82/256  ,  140/256 ],'Facealpha',0.5,'EdgeAlpha',0.5);
xlim([-100 100])
ylim([-0.5 1.2])



