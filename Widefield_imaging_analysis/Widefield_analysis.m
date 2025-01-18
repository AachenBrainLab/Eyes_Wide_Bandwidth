% This analysis utilizes the Matlab packages:
% "Image Processing Toolbox" and "Statistics and Machine Learning Toolbox"
% In addition, a number of general function in this repo:
% https://github.com/musallGroup/generalTools
% All other functions, I included in this repo. In case I missed a used
% dependency here, the Gabor-model or the PS-Stats, 
% please contact me at: g.nabbefeld@ucl.ac.uk


%%
clear all; clc; close all;
addpath("functions");


%% Load allen CCF map
load("allenDorsalMapSM.mat", "dorsalMaps"); allenMask = dorsalMaps.allenMask;


%%
% Define session data paths. This code assumes a certain file hierarchy and
% itterates through all available sessions.

% Those 2 dates where the final version of the stimulus protocol:
sessionPaths = [abs_paths_str("Y:\BpodBehavior\2*\MotionMare\Session Data\20230911*\Vc.mat"), ...
                abs_paths_str("Y:\BpodBehavior\2*\MotionMare\Session Data\20230912*\Vc.mat")]; c_lim = 0.025;
folder = fullfile(pwd, "Figures");

% Define stimulus conditions
conditions2plot = [1, 2, 3, 4, 13, 14];
condition_keys = ["Narrow", "Broad SF-BW", "Braod Ori-BW", "Braod SF + Ori-BW", "Narrow 0.16cpd", "Broad Ori-BW 0.16cpd"];
suffix = "";


% This function was used to align the imaging data to the allen CCF and
% safe the transformation in the file "opts2.mat"
% cc_alignBrainToAllen('Y:\BpodBehavior\2673\MotionMare\Session Data\20230911_084849\opts2.mat')


%% Loop over individual sessions
trace_session_cnt = 0;
for cPath = sessionPaths
    % extract parameters from session path
    session_folder = fileparts(cPath);
    [~, session_date] = fileparts(session_folder);
    [~, mouse_id] = fileparts(fileparts(fileparts(fileparts(session_folder))));

    disp(cPath)

    %% Load imaging data and cut into trials
    % define trials window to trialize the continuous recordings.
    baselineDur = 1;  % duration in sec
    responseDur = 3.1;  % duration in sec

    [Vc, U, opts, SessionData] = load_trialized_bpod_wf_data(session_folder, baselineDur, responseDur);
    
    [~, opts.date_time] = fileparts(session_folder);
    opts.date = string(opts.date_time{1}(1:8));

    % Define save folder
    save_folder = fullfile(folder, opts.cAnimal+"_"+opts.date_time+suffix);
    make_dir(save_folder);

    % align recordings
    if isfield(opts, "transParams")
        addpath("C:\Users\nabbefeld\Documents\GitHub\cortexCam\allenCode");  % Force Matlab to use this one!
        U = alignAllenTransIm_for_MotionMule(U, opts.transParams);
    end

    % the way I wrote this script
    responseDur = 3;
    movies = generate_movies_V4(U, Vc, opts, SessionData, responseDur, conditions2plot, condition_keys);

    % store movies over sessions
    if trace_session_cnt == 0
        all_movies = movies;
    else
        all_movies(:, :, :, :, trace_session_cnt+1) = movies;
    end

    area_masks_50 = find_most_responsive_50percent_pixel_masks(movies, opts);

%     % Optinal visualization for individual sessions:
%     %% Visualize average activity
%     baselineDur = 1;
%     responseDur = 3;
%     plot_average_stim_response_V4(U, Vc, opts, SessionData, baselineDur, responseDur, c_lim, conditions2plot, condition_keys);
%     fig_size(gcf, [0 0 10/16, 1]);
%     exportgraphics(gcf, fullfile(save_folder, mouse_id + "_" + session_date + "_average_response_maps_MotionMule_V4.png"), "Resolution", 500);
%     
%     %% Differential Maps
%     responseDur = 2;
%     plot_differential_maps_V4(U, Vc, opts, SessionData, responseDur, conditions2plot, condition_keys, round(c_lim/3,  2));
%     fig_size(gcf, [0 0 1, 1]);
%     exportgraphics(gcf, fullfile(save_folder, mouse_id + "_" + session_date + "_Differential_response_maps_V4.png"), "Resolution", 500);
% 
%     %%
%     temp_c_lim = round(c_lim/6, 2); if temp_c_lim == 0; temp_c_lim = round(2*c_lim/6, 2)/2; end
%     plot_differential_maps_V4(U, Vc, opts, SessionData, responseDur, conditions2plot(5:6), condition_keys(5:6), temp_c_lim);
%     fig_size(gcf, [0 0 1, 1]);
%     exportgraphics(gcf, fullfile(save_folder, mouse_id + "_" + session_date + "_Differential_response_maps_V4_high_SF_against_each_other.png"), "Resolution", 500);
    

    %%
    temp_AUCs = compute_auc_maps(U, Vc, opts, 1:4);

    %%
    if isfield(opts, "transParams")  % case to ensure data was properly aligned first
%         % optional plot to look at traces
%         temp_traces = plot_area_traces_by_stimulus(U, Vc, opts, SessionData, conditions2plot, condition_keys);
%         exportgraphics(gcf, fullfile(save_folder, "Traces.png"), "Resolution", 500);
        
        % This is computes the traces for the most over all response 50% pixels of the area:
        temp_traces_50 = compute_area_traces_by_stimulus_partial_areas(U, Vc, opts, SessionData, conditions2plot, condition_keys, area_masks_50);
        
        % accumulate results across sessions
        trace_session_cnt = trace_session_cnt + 1;
        if trace_session_cnt == 1
%             all_traces = temp_traces;  % [frames, trials, area_id, condition]
            all_traces_50 = temp_traces_50;  % these are the traces for only half the Pixels in a given area (the most responsive ones)
            all_area_masks_50 = area_masks_50;
        else
%             all_traces(:, :, :, :, trace_session_cnt) = temp_traces;
            all_traces_50(:, :, :, :, trace_session_cnt) = temp_traces_50;  % these are the traces for only half the Pixels in a given area (the most responsive ones)
            all_area_masks_50(:, :, :, trace_session_cnt) = area_masks_50;
        end
    end
end

% Save results 
% save(fullfile(fileparts(save_folder), "all_traces.mat"), "all_traces", "all_traces_50", "all_area_masks_50", "-v7.3");
save(fullfile(fileparts(save_folder), "all_traces.mat"), "all_traces_50", "all_area_masks_50", "-v7.3");
save(fullfile(fileparts(save_folder), "all_movies.mat"), "all_movies", "-v7.3");


%% Analysis across sessions and mice
%% Violin plot comparing individual trials for each area
n_conditions = 4;
% plot_violins_over_all_trials_response_integrals(all_traces, n_conditions, condition_keys)
% exportgraphics(gcf, fullfile(fileparts(save_folder), "Integrals_individual_trials.png"), "Resolution", 500);
% savefig(gcf, fullfile(fileparts(save_folder), "Integrals_individual_trials.fig"));

% Since not the entirety of V1 responds to the stimulation we compare the most response 50% of pixels across all stimuli.
plot_violins_over_all_trials_response_integrals(all_traces_50, n_conditions, condition_keys)
exportgraphics(gcf, fullfile(fileparts(save_folder), "Integrals_individual_trials_50perc_px.png"), "Resolution", 500);
savefig(gcf, fullfile(fileparts(save_folder), "Integrals_individual_trials_50perc_px.fig"));



%% Compare orientation modulation of V1 and higher visual areas
targetAreas = ["VISp", "VISpm", "VISam", "VISa", "VISrl", "VISal", "VISl"];
clearvars hndl ax
close all; figure(1); fig_size()

cl = gca().ColorOrder;

bonferroni_corr = 1;
try close(2); catch; end
figure(2); fig_size(gcf, [0.1 0.05 0.3 0.4]);

y_lim = [-5 13];
for condition_id = 1
    for area_id = 1:7
        % obtain the average response towards narrow and Ori stimuli for
        % session and area. 
        temp_data_narrow = squeeze(nanmean(nanmean(all_traces_50(16:45, :, area_id, 1, :), 1), 2));
        temp_data_ori = squeeze(nanmean(nanmean(all_traces_50(16:45, :, area_id, 3, :), 1), 2));

        temp_data = 100*(cat(3, temp_data_narrow', temp_data_ori'));
        
        if area_id == 1
            response_integrals = temp_data;
        else
            response_integrals(area_id, :, :) = temp_data;
        end
    end

    response_integrals_individual = response_integrals;  % keep it for the LME-stats
    response_integrals = response_integrals(:, :, 2) - response_integrals(:, :, 1);  % take the difference for the violins


    ax{condition_id} = subplot(1, 2, condition_id); hold on;

    % Sort areas by effect size:
    [~, order] = sort(median(response_integrals, 2), "descend");
    response_integrals = response_integrals(order, :);
    response_integrals_individual = response_integrals_individual(order, :, :);
    
    plot([0.5 7.5], [0 0], "--k"); hold on;
    v = violinplot(response_integrals');
    for i = 1:7
        v(i).ViolinColor = cl(i, :);
    end
    xlim([0.5 7.5]);


    % LME model
    for i = 1:7
        % compare modulation for each area
        [p, t, ~, ~] = LME_compare(response_integrals_individual(i, :, 1), response_integrals_individual(i, :, 2), [1 2 3 4 1 2 3 4]);

        % bonferroni correction (7 regions)
        p = 7*p;

        % print LME stats
        % This prints the stats for Figure 5c:
        fprintf("%s: p=%0.2e, T=%0.2f\n", targetAreas(order(i)), 7*p, t);

        % indicate significances
        text(i, max(response_integrals(i, :)) + 0.1, get_significance_str(p), "VerticalAlignment", "middle", "HorizontalAlignment", "center");
    end

    xticks(1:7); xticklabels(targetAreas(order)); xtickangle(45);
    ylabel("Ori - narrow (\DeltaF/F in %)")
    title("Ori - narrow session-wise");
    response_difference_by_area = response_integrals;

    ylim(ylim() + [0, 0.15]);
end

base_name = fullfile(fileparts(save_folder), "Ori_vs_narrow_integrals_bonferroni");

exportgraphics(gcf, base_name+".png", "Resolution", 600);
exportgraphics(gcf, base_name+".pdf");
save(base_name+".mat", "response_difference_by_area", "-v7.3");

% save results for visualization
target_areas_in_order = targetAreas(order);
save("Figures\Area_wise_Ori_response_mudulation.mat", "response_integrals", "target_areas_in_order", "-v7.3");




%% Average response plot
try close(1); catch; end
figure(1); fig_size(gcf, [.05 .05 .5 .3]);

c_lim = 3.0;
narrow_resp = squeeze(nanmean(nanmean(all_movies(:, :, opts.preStim+(1:30), 1, :), 5), 3));
for stim_id = 1:4
    subplot(1, 4, stim_id);
    img = squeeze(nanmean(nanmean(all_movies(:, :, opts.preStim+(1:30), stim_id, :), 5), 3));
    mean_img = nanmean(img, 3);
    if isfield(opts, "transParams")
        plot_pretty_allen_brain_GN(100*img, c_lim .* [-1, 1], fireice(6*256), 1.0);
        caxis(c_lim .* [-1, 1]);
    end
    title("[" + condition_keys(stim_id) + "] - [Narrow]");
end
cb = colorbar_without_position_change(); cb.Ticks = c_lim .* [-1 0 1]; ylabel(cb, "Change in \DeltaF/F [%]");
exportgraphics(gcf, fullfile(fileparts(save_folder), "Average_responses.png"), "Resolution", 500);
savefig(gcf, fullfile(fileparts(save_folder), "Average_responses.fig"));


%% Average response difference plot
try close(1); catch; end
figure(1); fig_size(gcf, [.05 .05 .5 .3]);

c_lim = 1.25;
narrow_resp = squeeze(nanmean(nanmean(all_movies(:, :, opts.preStim+(1:30), 1, :), 5), 3));
for stim_id = 2:4
    subplot(1, 3, stim_id-1);
    img = squeeze(nanmean(nanmean(all_movies(:, :, opts.preStim+(1:30), stim_id, :), 5), 3)) - narrow_resp;
    mean_img = nanmean(img, 3);
    if isfield(opts, "transParams")
        plot_pretty_allen_brain_GN(100*img, c_lim .* [-1, 1], fireice(6*256), 1.0);
        caxis(c_lim .* [-1, 1]);
    end
    title(condition_keys(stim_id));
end
cb = colorbar_without_position_change(); cb.Ticks = c_lim .* [-1 0 1]; ylabel(cb, "Change in \DeltaF/F [%]");
exportgraphics(gcf, fullfile(fileparts(save_folder), "Average_response_difference.png"), "Resolution", 500);
savefig(gcf, fullfile(fileparts(save_folder), "Average_response_difference.fig"));


%% Average response for the narrow stimulus alone
try close(1); catch; end
figure(1); fig_size(gcf, [.05 .05 .5 .3]);

c_lim = 2;  % in %
for stim_id = 1
    subplot(1, 3, stim_id);
    img = squeeze(nanmean(nanmean(all_movies(:, :, opts.preStim+(1:30), stim_id, :), 5), 3)) - ...
          squeeze(nanmean(nanmean(nanmean(all_movies(:, :, 1:opts.preStim, :, :), 5), 4), 3));
    mean_img = nanmean(img, 3);
    plot_pretty_allen_brain_GN(100*img, c_lim .* [0, 1], viridis(6*256), 1.0);
    caxis(c_lim .* [0, 1]);
    title(condition_keys(stim_id));
end
cb = colorbar_without_position_change(); cb.Ticks = c_lim .* [-1 0 1]; ylabel(cb, "\DeltaF/F [%]");
exportgraphics(gcf, fullfile(fileparts(save_folder), "Average_Narrow.png"), "Resolution", 500);
savefig(gcf, fullfile(fileparts(save_folder), "Average_Narrow.fig"));


%% Compute the stats Figure 5d:
LME_trials_averaged_over_recordings



