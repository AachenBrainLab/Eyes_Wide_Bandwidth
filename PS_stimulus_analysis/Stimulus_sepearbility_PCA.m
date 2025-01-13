% define color schema
cl = [214 180 148; 128 110 166; 38 101 168; 146 90 68] ./ 255;


%% load PS-stats
if ~exist("all_params", "var")
    load("IndividualFrames_Resolution_7_SpatialScales_6_Orientations_4.mat", "all_params");
end


%% concat stats to run a single PCA over those
plot_mean = 1;
temp = [reshape(cat(4, all_params.pixelStats), [], 3840); ...
        reshape(cat(4, all_params.pixelLPStats), [], 3840); ...
        reshape(cat(4, all_params.autoCorrReal), [], 3840); ...
        reshape(cat(4, all_params.autoCorrMag), [], 3840); ...
        reshape(cat(4, all_params.magMeans), [], 3840); ...
        reshape(cat(4, all_params.cousinMagCorr), [], 3840); ...
        reshape(cat(4, all_params.parentMagCorr), [], 3840); ... 
        reshape(cat(4, all_params.cousinRealCorr), [], 3840); ...
        reshape(cat(4, all_params.parentRealCorr), [], 3840); ...
        reshape(cat(4, all_params.varianceHPR), [], 3840)];

% Z-Score parameters for PCA
temp = (temp - nanmean(temp, 2)) ./ nanstd(temp, [], 2); 

opts.Display = 1;
coeff = pca(temp, "Options", opts);


%% Plotting
% %% First plot individual frames
% prepare figure
figure(2);
drawnow; set(gcf, 'Units', 'centimeters', 'Position', [.5, 2, 16 .* (1.1 .* 2 ./ 3), 16 .* (1.1 ./ 3)], 'PaperUnits', 'centimeters', 'PaperSize', [21, 29.7]); drawnow;

% loop over conditions to apply color schema
subplot(1, 2, 1); hold off;
for i = 1:4
    frame_ids = (1:960) + (i-1).*960;
    plot(coeff(frame_ids, 1), coeff(frame_ids, 2), "o", "MarkerFaceColor", cl(i, :), "MarkerEdgeColor", cl(i, :), "MarkerSize", 1); hold on;
end
% define axis labels
xlabel("PC1");
ylabel("PC2");
title("Individual frames")

% axis limits
xlim(0.065 .* [-1 1]);
ylim(0.065 .* [-1 1]);



% %% Then do the same for the movie-wise averages
subplot(1, 2, 2); hold off;
data_for_table = zeros(4*8, 2);
StimulusName = string();
ConditionKeys = ["Narrow", "SF", "Ori", "Mixed"];
% loop over conditions to apply color schema
for i = 1:4
    frame_ids = (1:960) + (i-1).*960;
    temp_data = squeeze(nanmean(reshape(coeff(frame_ids, 1:3), 120, 8, []), 1));
    h = plot(temp_data(:, 1), temp_data(:, 2), "o", "MarkerFaceColor", cl(i, :), "MarkerEdgeColor", 1 .* cl(i, :), "MarkerSize", 2); hold on;

    data_for_table((i-1)*8+(1:8), :) = temp_data(:, 1:2);
    StimulusName((i-1)*8+(1:8), :) = repmat(ConditionKeys(i)+"_StimulusMovie_", [8, 1])+int2str((1:8)');
end
% define axis labels
xlabel("PC1");
ylabel("PC2");
title("Grouping of PS-statistics by condition");

% axis limits and ticks
xticks(0.04 .* [-1 -0.5 0 0.5 1]);
yticks(0.04 .* [-1 -0.5 0 0.5 1]);
xlim(0.045 .* [-1 1]);
ylim(0.045 .* [-1 1]);


%% save PC1 and PC2 values for the data-table:
PC1 = data_for_table(:, 1);
PC2 = data_for_table(:, 2);
T = table(StimulusName, PC1, PC2);
writetable(T, "Data_table_S7_Portilla_comp_PCA.xlsx");
clearvars StimulusName PC1 PC2 data_for_table


%% sets the fontsize for all figure elements
fontsize_for_every_figure_element(7, "normal");

%%
exportgraphics(gcf, "S7_Portilla_Simoncelli_stats_PCA.png", "Resolution", 1000);

%% I save a separate version to get the legend
legend(["narrow", "SF", "Ori", "mixed"])
exportgraphics(gcf, "Portilla PCA_legend.png", "Resolution", 1000);

