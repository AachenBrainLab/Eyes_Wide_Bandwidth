% Either average these stats for all frames of a given stimulus movie
% realization. In the final submission we showed these stats on the level
% of the individual stimulus frames:
plot_individual_frames = 1;

% I looked at different scales but the final version was version 3, only used fo naming here.
version = 3;
if ~exist("all_params", "var")  % to prevent reloading when running again
    load("IndividualFrames_Resolution_7_SpatialScales_6_Orientations_4.mat", "all_params")
    spatialDimension2use = 6:6;  % use the coarsest scale
end

cl = [214 180 148; 128 110 166; 38 101 168; 146 90 68] ./ 255;


%%
for spatialDimension2use = 6  % use the coarsest scale
condition_keys = ["Narrow", "SF", "Ori", "Mixed"];

close all;
temp = get_average_CrossOrientation(all_params, spatialDimension2use);
figure(1);

temp = reshape(temp, [120*8, 4]);  
data4stats = squeeze(nanmean(reshape(temp, 120, 8, []), 1));  % preallocate data for statistics
temp_mtx = temp;
for i = 1:4
    [~, ids] = rmoutliers(temp(:, i));
    temp_mtx(ids, i) = nan;
end
subplot(1, 2, 1); plot([-1 5], [0 0], "--", "color", 0.5 .* [1 1 1]);
if plot_individual_frames
    v = violinplot(temp_mtx); ylabel("Cross Orientation Energy (A.U.)");
else
    v = violinplot(data4stats); ylabel("Cross Orientation Energy (A.U.)");
end
xticks(1:4); xticklabels(condition_keys);

CrossOrientation_individual_frame_data_table = reshape(temp_mtx, [4*960, 1]); % first all frames of video 1, then 2, 3, 4

%% INSERT FOR LME STATS %
disp("Cross Orientation Energies:")
CrossEnergyStatsMovies
crossOrientationEnergy_p = crossEnergy_p;

%%
vCl = 0.425;
for i = 1:4
    v(i).ViolinColor = cl(i, :);
    if plot_individual_frames
        v(i).ViolinAlpha = 0.15;
        v(i).BoxColor = vCl .* [1 1 1];
        v(i).BoxPlot.Vertices(:, 1) = v(i).BoxPlot.Vertices(:, 1) + [-1, 1, 1, -1]' .* 0.04;
        v(i).EdgeColor = vCl .* [1 1 1];
    
        v(i).ScatterPlot.SizeData = 7;
    end
end

% prepare axis
xlim([0.5 4.5]);

y_lim = ylim();
yRange = 1.6 .* (y_lim(2) - y_lim(1));

% Significance bars etc.
y_offset = yRange/1000;
bar_padding = 0.1;
y_level_size = yRange/13;
y_start = y_lim(2);
draw_significance_with_bar(1, 2, y_start + 0.*y_level_size, bar_padding, y_offset, crossOrientationEnergy_p(1));
draw_significance_with_bar(1, 3, y_start + 1.*y_level_size, bar_padding, y_offset, crossOrientationEnergy_p(2));
draw_significance_with_bar(1, 4, y_start + 3.*y_level_size, bar_padding, y_offset, crossOrientationEnergy_p(3));

draw_significance_with_bar(2, 3, y_start + 0.*y_level_size, bar_padding, y_offset, crossOrientationEnergy_p(4));
draw_significance_with_bar(2, 4, y_start + 2.*y_level_size, bar_padding, y_offset, crossOrientationEnergy_p(5));
draw_significance_with_bar(3, 4, y_start + 0.*y_level_size, bar_padding, y_offset, crossOrientationEnergy_p(6));


% update axis lims
ylim(y_lim + [0, 0.5 .* diff(y_lim)]);
yticks(-10000:1000:10000);

% labels
title("Cross Orientation Energy")
ylabel("Energy (A.U.)")


ax = gca(); ax.YAxis.Exponent = 0;


% %% Cross Scales Energy %% %
temp2 = get_average_CrossScale(all_params);

temp2 = reshape(temp2, [120*8, 4]);

temp_mtx2 = temp2;
data4stats2 = squeeze(nanmean(reshape(temp2, 120, 8, []), 1));
for i = 1:4
    [~, ids] = rmoutliers(temp2(:, i));
    temp_mtx2(ids, i) = nan;
end

CrossScales_individual_frame_data_table = reshape(temp_mtx2, [4*960, 1]); % first all frames of video 1, then 2, 3, 4

% plot
subplot(1, 2, 2); plot([-1 5], [0 0], "--", "color", 0.5 .* [1 1 1]);
if plot_individual_frames
    v = violinplot(temp_mtx2); ylabel("Cross Scale Energy (A.U.)");
else
    v = violinplot(data4stats2); ylabel("Cross Scale Energy (A.U.)");
end

xticks(1:4); xticklabels(condition_keys);

for i = 1:4
    v(i).ViolinColor = cl(i, :);
    if plot_individual_frames
        v(i).ViolinAlpha = 0.15;
        v(i).BoxColor = vCl .* [1 1 1];
        v(i).BoxPlot.Vertices(:, 1) = v(i).BoxPlot.Vertices(:, 1) + [-1, 1, 1, -1]' .* 0.04;
        v(i).EdgeColor = vCl .* [1 1 1];
    
        v(i).ScatterPlot.SizeData = 7;
    end
end
xlim([0.5 4.5]);

y_lim = ylim();
yRange = 1.6 .* (y_lim(2) - y_lim(1));


disp("Cross Scales Energies:")
data4stats = data4stats2;

% compute bonferroni corrected ranksum tests:
CrossEnergyStatsMovies
crossScaleEnergy_p = crossEnergy_p;


% Significance bars etc.
y_offset = yRange/1000;
bar_padding = 0.1;
y_level_size = yRange/13;
y_start = y_lim(2);
draw_significance_with_bar(1, 2, y_start + 0.*y_level_size, bar_padding, y_offset, crossScaleEnergy_p(1));
draw_significance_with_bar(1, 3, y_start + 1.*y_level_size, bar_padding, y_offset, crossScaleEnergy_p(2));
draw_significance_with_bar(1, 4, y_start + 3.*y_level_size, bar_padding, y_offset, crossScaleEnergy_p(3));

draw_significance_with_bar(2, 3, y_start + 0.*y_level_size, bar_padding, y_offset, crossScaleEnergy_p(4));
draw_significance_with_bar(2, 4, y_start + 2.*y_level_size, bar_padding, y_offset, crossScaleEnergy_p(5));
draw_significance_with_bar(3, 4, y_start + 0.*y_level_size, bar_padding, y_offset, crossScaleEnergy_p(6));


%
ylim(y_lim + [0, 0.5 .* diff(y_lim)]);
yticks(-1000:100:1000);

% labels
title("Cross Scale Energy")
ylabel("Energy (A.U.)")


fig_size(gcf, [0.1 0.1 0.275 0.325]);
fontsize_for_every_figure_element(7, "normal");

drawnow; set(gcf, 'Units', 'centimeters', 'Position', [.5, 2, 16 .* (1.1 .* 2 ./ 3), 16 .* (1.1 ./ 3)], 'PaperUnits', 'centimeters', 'PaperSize', [21, 29.7]); drawnow;

folder = "CrossEnergies"; make_dir(folder);
if plot_individual_frames
    exportgraphics(gcf, fullfile(folder, "CrossEnergies_Version"+version+"_Scale"+spatialDimension2use+".png"), "Resolution", 1000);
else
    exportgraphics(gcf, fullfile(folder, "CrossEnergies_Version"+version+"_Scale"+spatialDimension2use+"_movie_averages.png"), "Resolution", 1000);
end

end


%% SAVE THE DATA TABLE
CrossOrientationEnergy = CrossOrientation_individual_frame_data_table;
CrossScalesEnergy = CrossScales_individual_frame_data_table;

StimulusName = string();
ConditionKeys = ["Narrow", "SF", "Ori", "Mixed"];
for i = 1:4
    for movie_id = 1:8
        StimulusName((i-1)*960+(movie_id-1)*120+(1:120), :) = repmat(ConditionKeys(i)+"_StimulusMovie_"+movie_id+"_Frame", [120, 1])+string(reshape(char(sprintf("%04d", (1:120))), [4, 120])');
    end
end
T = table(StimulusName, CrossOrientationEnergy, CrossScalesEnergy);
writetable(T, "Data_table_S7_Portilla_CrossEnergies.xlsx");
clearvars StimulusName CrossOrientationEnergy CrossScalesEnergy CrossOrientation_individual_frame_data_table CrossScales_individual_frame_data_table







%% Functions
function temp = get_average_CrossOrientation(all_params, spatialDimension2use)
    temp = cat(4, all_params(1, :).cousinMagCorr);

    for i = 1:size(temp, 1)
        % exclude the diagonal, as we are not interested in the
        % auto-correlations
        temp(i, i, :, :) = nan;
    end
    
    % in addition we now only want to focus on the 0° vs. +- 45°, so
    % excluding the 0° vs. 90°
    % exclude all comparrision vs. 90°:
    temp(:, 3, :, :) = nan;
    temp(3, :, :, :) = nan;

    % also exclude 45° vs. -45°:
    temp(2, 4, :, :) = nan;
    temp(4, 2, :, :) = nan;

    % note: the last one is always just all zeros
    temp = reshape(nanmean(nanmean(nanmean(temp(:, :, spatialDimension2use, :), 3), 2), 1), [120, 8, 4]);
end


function temp = get_average_CrossScale(all_params)
    temp = cat(4, all_params(1, :).parentMagCorr);

    % Only use the diagonal, as we want to correlation for same 
    % orientations, but between differen scales.
    % Therefore, only use the identity matrix values:
    for i = 1:size(temp, 1)
        for j = 1:size(temp, 2)
            if ~(i == j)
                temp(i, j, :, :) = nan;
            end
        end
    end
    
    % We only want to look at the correlation of the Orientation 1 (0deg) 
    % across scales:
    temp = temp(1, 1, :, :);

    % remove the latest, empty ones:
    temp = temp(:, :, 1:end-1, :);

    % use the last, coarsest scale
    spatialDimension2use = size(temp, 3);
    temp = reshape(nanmean(nanmean(nanmean(temp(:, :, spatialDimension2use, :), 1), 2), 3), [120, 8, 4]);
end



