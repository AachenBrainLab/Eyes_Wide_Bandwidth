% Run the portilla simoncelli analysis
n_conditions = 32;
n_frames = 120;

nSamples = n_conditions .* n_frames;

% settings for the texture analysis
Nsc=6;
Nor=4;
Na=7;

% To generate the autocorrelation plots, run this script with these
% settings:
if 0
    Nsc=2;
    Nor=2;
    Na=51;
end


% preallocate outputs for later
all_stats        = nan(6, nSamples);
all_autoCorrReal = nan(Na, Na, Nsc+1, nSamples);  % Not sure if Nsc or Nor
all_autoCorrMag  = nan(Na, Na, Nsc, Nor, nSamples);  % Not sure if Nsc or Nor

figure(1);
Cnt = 0;
for stimulus_folder_id = 1:n_conditions
    disp(stimulus_folder_id + "/" + n_conditions);
    
    filePaths = abs_paths_str("D:\Data\MC_MotionMule_StimFrame_for_analysis\" + sprintf("%03d", stimulus_folder_id) + "\frame*.png");
    for fPath = filePaths(1:n_frames)
        im0 = double(imread(fPath));

        % the algorithm seems to have issues with the original frame size
        % If I undestand this right the frames need to be (2^n) x (2^m)
        % im0 = im0(1:256, 1:512);
        im0 = im0(1:256, 1:256);
        im0 = im0 - nanmean(im0(:));
        im0 = im0 - nanmean(im0(:));  % I do this twice as there as small rounding errors
        im0 = im0 ./ (255./2);
        
        % normalize lum, but doesnt really have an effect
        % im0 = ((im0) ./ (255 ./ 2)) - 1;

        if 0  % Test image
            im0 = zeros(256, 256);
            x = (1:256) - 128.5;
            y = (1:256) - 128.5;
            im0(sqrt(power(x, 2) + power(y', 2)) < 50) = 1;
            im0 = imgaussfilt(im0, 11);
        end
        
        % This runs the PS-texture analysis code on the provided frame:
        % This assumes that you have cloned the repository provided by the
        % Portilla Simoncelli publication and added it to your Matlab path.
        [params] = textureAnalysis(im0, Nsc, Nor, Na);

        Cnt = Cnt + 1;

        all_stats(:, Cnt) = params.pixelStats';
        all_autoCorrReal(:, :, :, Cnt) = params.autoCorrReal;
        all_autoCorrMag(:, :, :, :, Cnt) = params.autoCorrMag;
        
        if Cnt == 1
            all_params = params;
        else
            all_params(1, Cnt) = params;
        end
    end
end


%% Average_pixel parameters
save("IndividualFrames_Resolution_" + Na + "_SpatialScales_" + Nsc + "_Orientations_" + Nor +  ".mat", "n_conditions", "n_frames", "Nsc", "Nor", "Na", "all_stats", "all_autoCorrReal", "all_autoCorrMag", "all_params", "-v7.3");



% %% Keep these for now, but I think this was just me testing things and
%    nothing beyond here used:

% %% 
% n_spatial_scales = size(all_autoCorrMag, 3);
% n_orientations = size(all_autoCorrMag, 4);
% norm_data_keys = ["raw", "normalized"];
% for norm_data = [0 1]
% for condition_id = 1:4
%     figure(condition_id+2); fig_size(gcf, [0.05+0.3*mod(condition_id+1, 2), 0.05+0.45*floor((condition_id-1)/2), 0.3, 0.45])
%     stim_ids = ceil((1:size(all_autoCorrReal, 4))/8) == condition_id;
%     clearvars ax;
%     for i = 1:5
%         if i == 1
%             ax = {subplot(5, 5, i)};
%         else
%             ax{end+1} = subplot(5, 5, i);
%         end
%         imagesc(nanmean(all_autoCorrReal(:, :, i, stim_ids), 4)); colormap(gca, jet); axis image; xticks([]); yticks([]);
%         if i == 3; title(condition_keys(condition_id)); end
%     end
%     if norm_data; match_axis(ax, 1); end
%     
%     clearvars ax;
%     for i = 1:n_spatial_scales
%         for j = 1:n_orientations
%             if i == 1 && j == 1
%                 ax = {subplot(n_spatial_scales+1, n_orientations, n_orientations*i+j)};
%             else
%                 ax{end+1} = subplot(n_spatial_scales+1, n_orientations, n_orientations*i+j);
%             end
%             
%             imagesc(nanmean(all_autoCorrMag(:, :, i, j, stim_ids), 5)); colormap(gca, gray); axis image; xticks([]); yticks([]);
%         end
%     end
%     if norm_data; match_axis(ax, 1); end
%     
%     make_dir(pwd+"\Figures");
%     exportgraphics(gcf, pwd+"\Figures\" + condition_keys(condition_id) + "_" + norm_data_keys(norm_data+1) + ".png", "Resolution", 300);
% end
% drawnow
% end
% 
% 
% 
% %% try to use a PCA to better compare stimuli of diffrent features
% % The main aim here is to show that broad broad is probably further away
% % from the other stimuli then the others are from each other. 
% 
% try close(2); catch; end; figure(2); fig_size(gcf, [.1 .1 .8 .5]);
% n_spatial_scales = size(all_autoCorrMag, 3);
% n_orientations = size(all_autoCorrMag, 4);
% for condition_id = 1:4
%     clearvars ax;
%     for i = 1:n_spatial_scales
%         for j = 1:n_orientations
%             if i == 1 && j == 1
%                 ax = {subplot(n_spatial_scales+1, n_orientations, n_orientations*i+j)};
%             else
%                 ax{end+1} = subplot(n_spatial_scales+1, n_orientations, n_orientations*i+j);
%             end
%             
%             imagesc(nanmean(all_autoCorrMag(:, :, i, j, stim_ids), 5)); colormap(gca, gray); axis image; xticks([]); yticks([]);
%         end
%     end
%     if norm_data; match_axis(ax, 1); end
%     
%     make_dir(pwd+"\Figures");
%     exportgraphics(gcf, pwd+"\Figures\" + condition_keys(condition_id) + "_" + norm_data_keys(norm_data+1) + ".png", "Resolution", 300);
% end
% drawnow
% 
% 
% 
% 
% 
% %% Raw_coefficient_correlation
% close all;
% figure(1); fig_size(gcf, [0.1 0.1 0.4 0.45]);
% for j = 1:2
%     subplot(1, 2, j);
%     for i = 1:4
%         use_ids = (1:10*8) + (i-1) .* 10*8;
%         plot(squeeze(nanmean(all_autoCorrReal(4, 4, :, use_ids), 4)), "LineWidth", 1.5); hold on;
%     end
%     
%     if j == 1
%         legend(["Narrow-SF/Narrow-Ori", "Broad-SF/Narrow-Ori", "Narrow-SF/Broad-Ori", "Broad-SF/Broad-Ori"], "location", "nw")
%     end
%     xlabel("scale of Auto-correlations (i think)");
%     ylabel("Magnituge of Auto-correlations");
% end
% set(gca, "yscale", "log");
% exportgraphics(gcf, pwd+"\Figures\Raw_coefficient_correlation.png", "Resolution", 300);
% 
% 
% %% Coefficient magnitude statistics
% close all;
% figure(1); fig_size(gcf, [0.1 0.1 0.7 0.45]);
% for ori_ids = 1:4
%     subplot(1, 4, ori_ids);
%     for i = 1:4
%         use_ids = (1:10*8) + (i-1) .* 10*8;
%         plot(squeeze(nanmean(all_autoCorrMag(4, 4, ori_ids, :, use_ids), 5)), "LineWidth", 1.5); hold on;
%     end
%     title("Orientation: " + ((ori_ids-1) * 45) + "°")
%     
% %     if j == 1
% %         legend(["Narrow-SF/Narrow-Ori", "Broad-SF/Narrow-Ori", "Narrow-SF/Broad-Ori", "Broad-SF/Broad-Ori"], "location", "nw")
% % %     end
%     xlabel("Spatial scale of Auto-correlations");
%     ylabel("Coefficient magnitudes of Auto-correlations");
% end
% set(gca, "yscale", "log");
% % exportgraphics(gcf, pwd+"\Figures\Coefficient_magnitude_statistics.png", "Resolution", 300);
% 
% 
% 
% 
% 
% %% Look at each frame individually
% close all;
% figure(1);
% n_frames = 120;
% 
% % [parameter_id, n_frames_per_condition, condition]
% params = nan(4, n_frames*8, 4);
% for i = 1:4
%     % use_ids = (1:10*8) + (i-1) .* 10*8;
%     use_ids = (1:n_frames*8) + (i-1) .* n_frames*8;
%     
%     Cnt = 0;
%     for frame_ids = use_ids
%         autoCorrReal = nanmean(cat(4, all_params(frame_ids).autoCorrReal), 4);
%         autoCorrMag  = nanmean(cat(5, all_params(frame_ids).autoCorrMag), 5);
% 
%         autoCorrReal             = sum(vec(abs(autoCorrReal(:))));
% %         autoCorrReal             = nanmean(vec(autoCorrReal(4, 4, :)));
%         cross_position_energy    = nanmean(vec(autoCorrMag(:)));
%         cross_scale_energy       = nanmean(vec(prod(autoCorrMag, 3)));
%         cross_orientation_energy = nanmean(vec(prod(autoCorrMag, 4)));
%         
%         Cnt = Cnt + 1;
%         params(:, Cnt, i) = [autoCorrReal, cross_position_energy, cross_scale_energy, cross_orientation_energy];
%     end
% end
% 
% 
% for params_id = 1:4
% %     params_mean = nanmean(params, 1);
% %     params_error = compute95CI(params, 1);
% 
%     subplot(1, 4, params_id);
%     h = violinplot(vec(params(params_id, :, :)), ceil((1:n_frames*8*4) ./ (n_frames*8)));
%     xticks(1:4);
% 
% %     if i == 1
% %         reference_values = params_mean;
% %     end
% %     
% %     for params_id = 1:4
% %         subplot(1, 4, params_id); bar(params_id, params_mean(params_id) ./ reference_values(params_id)); hold on; xticks(1:4);
% % 
% %         subplot(1, 4, params_id); bar(params_id, params_mean(params_id) ./ reference_values(params_id)); hold on; xticks(1:4);
% %     end
% end
% 
% 
% 
% %% Look at each frame individually
% % close all;
% figure(2);
% n_frames = 120;
% 
% % [parameter_id, n_frames_per_condition, condition]
% params = nan(4, 8, 4);
% for i = 1:4
%     % use_ids = (1:10*8) + (i-1) .* 10*8;
%     use_ids = (1:n_frames:n_frames*8) + (i-1) .* n_frames*8;
%     
%     Cnt = 0;
%     for frame_ids = use_ids
%         autoCorrReal = nanmean(cat(4, all_params(frame_ids + (0:n_frames-1)).autoCorrReal), 4);
%         autoCorrMag  = nanmean(cat(5, all_params(frame_ids + (0:n_frames-1)).autoCorrMag), 5);
% 
% %         autoCorrReal             = nanmean(vec(abs(autoCorrReal(:))));
% %         autoCorrReal             = nanmean(vec(autoCorrReal(:)));
% %         autoCorrReal             = nanmean(vec(autoCorrReal(4, 4, :)));
% 
%         if 1  % Center only
%             autoCorrReal             = nanmean(vec(autoCorrReal(4, 4, 1:end-1)));
%             cross_position_energy    = nanmean(vec(autoCorrMag(4, 4, :, :)));
%             cross_scale_energy       = nanmean(vec(prod(autoCorrMag(4, 4, :, :), 3)));
%             cross_orientation_energy = nanmean(vec(prod(autoCorrMag(4, 4, :, :), 4)));
%         else
%             autoCorrReal             = nanmean(vec(autoCorrReal(:, :, 1:end-1)));
%             cross_position_energy    = nanmean(vec(autoCorrMag));
%             cross_scale_energy       = nanmean(vec(prod(autoCorrMag, 3)));
%             cross_orientation_energy = nanmean(vec(prod(autoCorrMag, 4)));
%         end
%         
% %         cross_scale_energy       = nanmean(vec(prod(autoCorrMag(:, :, 1:2, :), 3)+prod(autoCorrMag(:, :, 2:3, :), 3)+prod(autoCorrMag(:, :, 3:4, :), 3)));
% %         cross_orientation_energy = nanmean(vec(prod(autoCorrMag(:, :, :, 1:2), 4)+prod(autoCorrMag(:, :, :, 2:3), 4)+prod(autoCorrMag(:, :, :, 3:4), 4)));
%         
%         
%         Cnt = Cnt + 1;
%         params(:, Cnt, i) = [autoCorrReal, cross_position_energy, cross_scale_energy, cross_orientation_energy];
%     end
% end
% 
% 
% for params_id = 1:4
%     subplot(1, 4, params_id); hold off;
%     h = violinplot(vec(params(params_id, :, :)), ceil((1:8*4) ./ (8)));
%     xticks(1:4);
% %     if i == 1
% %         reference_values = params_mean;
% %     end
% %     
% %     subplot(1, 4, params_id); bar(params_id, params_mean(params_id) ./ reference_values(params_id)); hold on; xticks(1:4);
% end







