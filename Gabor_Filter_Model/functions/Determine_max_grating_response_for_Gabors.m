%% Project: Gabor Model of V1 responses in towards MotionCloud stimuli.
%  This script is used to have a way to normalize gabor-filter responses
%  from gabors with different spatial frequency tuning.
%  The appraoch here is that I "present" each of the
%  Gabors I'm using in the model with the appropriat Sine-wave
%  grating and determine the peak response of each Gabor to this Grating
%  and normalize responses according to this peak-grating-response. 
%  By: Gerion Nabbefeld, 2023.

%%
close all;
use_surroundSupression = 0;


bin_f = 4;
n_ori = 180;
n = round(640/bin_f);
w = n;
h = round(340/bin_f);
n_stims = 5;
results = nan(n_stims, n_ori, n_stims);

SF = sort([0.01 0.02 0.04 0.06 0.08], "descend");
cyc_per_image = SF*130;

n_scales = length(SF);
SpatialAspectRatio = 0.55;
gb = gabor(((640 ./ bin_f) ./ 130) ./ SF, -90+180 .* (0:n_ori-1) ./ n_ori, "SpatialAspectRatio", SpatialAspectRatio, "SpatialFrequencyBandwidth", 2.38);

figure(70);
set(gcf, 'Units', 'centimeters', 'Position', [.5, 2, 16 .* (6/6), 25.2 .* (2/4)], 'PaperUnits', 'centimeters', 'PaperSize', [21, 29.7]); drawnow;
for scale = 1:length(SF)
    subplot(length(SF), 1, scale);
    kernel = imag(gb(scale+n_scales*floor(n_ori/2)).SpatialKernel);
    kernel = kernel(int32(round((size(kernel, 1) ./ 2) + (-200/bin_f:200/bin_f))), ...
                    int32(round((size(kernel, 2) ./ 2) + (-200/bin_f:200/bin_f))));
    imagesc(kernel); colormap("gray"); axis image; colorbar(); caxis([-1 1]);
    xticks([]); yticks([]);
    title(SF(scale));
end
drawnow;


figure(1); fig_size();
for stim_id = 1:n_scales
    y = linspace(0, cyc_per_image(stim_id) .* 2 .* pi, w+1);
    img = repmat(sin(y(1:end-1)), [h, 1, 1]);
    
    % figure(1); subplot(5, n_scales, stim_id);
    figure(1); subplot(5, n_scales, n_scales-stim_id+1);
    imagesc(img); colormap("gray"); axis image; % colorbar();
    xticks([]); yticks([]);
    caxis([-1 1]);
    if stim_id == 1
        cb = colorbar_without_position_change();
    end
    title(sprintf("%0.2fcpd", SF(stim_id)));
    drawnow;

    mag = imgaborfilt(img, gb);
    results(:, :, stim_id) = reshape(nanmean(nanmean(mag, 1), 2), n_scales, n_ori);
    disp("cycles per image: " + cyc_per_image(stim_id));
end

figure(70); 

%%
max_responses = nan(5,1);
for i = 1:5
    max_responses(i, 1) = max(results(i, :, i), [], 2);
end

save_folder = fullfile(pwd, "finding_gabor_max_responses"); make_dir(save_folder);
save("finding_gabor_max_responses\Matched_gabors_max_responses.mat", "max_responses");


%%
figure(1);
plot_peak_resp_over_SFs = 0;
for stim_id = 1:n_stims
    %%
    subplot(4, n_stims, 1*n_stims+n_stims-stim_id+1);
    imagesc(results(:, :, stim_id) ./ max_responses); % axis image;
    xticks(1:45:n_ori+1); xticklabels((180/n_ori) .* ((1:45:n_ori+1)-1)); xlabel("Orientation (deg)")
    if stim_id == n_stims
        yticks(1:length(SF)); yticklabels(SF(1:end)); ylabel("Spatial freq. [cyc/img]");
        ylabel("Gabor-filter"+newline+"spatial frequency");
    else
        yticks([]); yticklabels([]); ylabel("");
    end
    caxis([0 1]);
    if stim_id == 1
        cb = colorbar_without_position_change();
        ylabel(cb, sprintf("%s", "Gabor-response"+newline+"magnitude"));
    end

    %%
    ax = subplot(4, n_stims, 2*n_stims+n_stims-stim_id+1);
    delete(ax);
    if plot_peak_resp_over_SFs
        subplot(4, n_stims, 2*n_stims+n_stims-stim_id+1);
        plot(max(results(:, :, stim_id), [], 2) ./ max_responses);
        xticks(1:length(SF)); xticklabels(SF(1:end)); xlabel("Spatial freq. [cyc/deg]");
        if stim_id == n_stims; ylabel(sprintf("%s", "Gabor-response"+newline+"magnitude")); end
    
        if stim_id == 1
            refResp = max(vec(results(:, :, stim_id)));
        else
        end
    end

    %%
    subplot(4, n_stims, (2+plot_peak_resp_over_SFs)*n_stims+n_stims-stim_id+1);
    plot(results(end:-1:1, :, stim_id)' ./ max_responses(end:-1:1, 1)');
    if stim_id == n_stims; ylabel(sprintf("%s", "Gabor-response"+newline+"magnitude")); end
    xticks(1:45:n_ori+2); xticklabels((180/n_ori) .* ((1:45:n_ori+1)-1)); xlabel("Orientation (deg)")
    xlim([0 181]);
end

fName = fullfile(save_folder, "Matched_gabors_max_resp_"+SpatialAspectRatio+"_aspectRatio");
exportgraphics(gcf, fName+".png", "Resolution", 300);
savefig(gcf, fName+".fig");



%%
drawnow;
fontsize_for_every_figure_element(7, "normal", 0.5);
fontsize_for_every_figure_element(6, "normal", 0.5);
set(findobj(gcf,'Type','line'), "LineWidth", 1.5);


drawnow;
set(gcf, 'Units', 'centimeters', 'Position', [.5, 2, 16 .* (6/6), 25.2 .* (1.6/4)], 'PaperUnits', 'centimeters', 'PaperSize', [21, 29.7]); drawnow;
exportgraphics(gcf, fName+"_TEST.png", "Resolution", 300);


%% save
exportgraphics(gcf, fName + ".pdf", "ContentType", "vector", "BackgroundColor", "none");
exportgraphics(gcf, fName + ".png", "BackgroundColor", "white", "Resolution", 300);
saveas(gcf, fName+".svg");
savefig(gcf, fName+".fig");



%% This is an extra step to save the legend as well
stim_id = 1;
subplot(4, n_stims, 3*n_stims+n_stims-stim_id+1);
plot(results(end:-1:1, :, stim_id)' ./ max_responses(end:-1:1, 1)');
if stim_id == n_stims; ylabel("Norm. Energy"); end
xticks(1:45:n_ori+1); xticklabels((180/n_ori) .* ((1:45:n_ori+1)-1)); xlabel("Orientation (deg)")

if stim_id == 1
    leg = legend(["0.01cpd", "0.02cpd", "0.04cpd", "0.06cpd", "0.08cpd"], "Location", "eastoutside");
    leg.ItemTokenSize = [6, 6];
end
exportgraphics(gcf, fName + "_LEGEND.pdf", "ContentType", "vector", "BackgroundColor", "none");




