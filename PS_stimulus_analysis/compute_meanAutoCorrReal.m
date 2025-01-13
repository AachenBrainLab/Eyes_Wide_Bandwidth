%% Either load or recompute the average real-autoCorrelations:
try
    recompute = 0;
    load("meanAutoCorrReal.mat", "meanAutoCorrReal");
catch
    recompute = 1;
    load("IndividualFrames_Resolution_51_SpatialScales_2_Orientations_2.mat", "all_autoCorrReal");

    % fetch the autocorrelations at the coarsest scale and average over frames of a given condition:
    meanAutoCorrReal          = nanmean(all_autoCorrReal(:, :, end, (1-1) .* (120*8) + (1:8*120)), 4);
    meanAutoCorrReal(:, :, 2) = nanmean(all_autoCorrReal(:, :, end, (2-1) .* (120*8) + (1:8*120)), 4);
    meanAutoCorrReal(:, :, 3) = nanmean(all_autoCorrReal(:, :, end, (3-1) .* (120*8) + (1:8*120)), 4);
    meanAutoCorrReal(:, :, 4) = nanmean(all_autoCorrReal(:, :, end, (4-1) .* (120*8) + (1:8*120)), 4);
end


%%
clearvars("ax")
figure();
for i = 1:4
    ax{i} = subplot(2, 2, i);
    imagesc(meanAutoCorrReal(:, :, i)); caxis([-1 1] .* max(abs(caxis()))); axis image; colormap(gca, jet);
end
colorbar_without_position_change();
match_axis(ax, 1);


%%
% Use these data to plot with custom color maps for the main figure
if recompute == 1
    save("meanAutoCorrReal.mat", "meanAutoCorrReal");
end

