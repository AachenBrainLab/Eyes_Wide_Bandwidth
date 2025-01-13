close all;
n_ori = 180;  % compute and plot in 1deg increments
bin_f = 4;  % use the same binning we apply to the stimulus frames
SF = 0.04;
sfb = 2.38;  % Spatial frequency bandwidth for: SpatialAspectRatio == 0.55


%%
all_sar = [0.3, 0.55, 1.];

n = length(all_sar);
results = nan(n, n_ori);
sar_cnt = 0;
px_per_deg = 160./130;
n_px = round(220 ./ bin_f);

x = linspace(0, 0.04*130*2*pi, 641);
img = repmat(sin(x(1:end-1))', 1, 360)';
img = imresize(img, 1/bin_f);

% figure
figure(1); subplot(6, 1, 1); fig_size();
imagesc(img); colormap("gray"); axis image; xticks([]); yticks([]);
cb = colorbar_without_position_change();
x = size(img, 2)./2 + [-n_px, n_px];
y = size(img, 1)./2 + [-n_px, n_px];
caxis([-1 1]);
drawnow;

for sar = all_sar
    sar_cnt = sar_cnt + 1;

    % Construct the full gabor bank once here
    gb = gabor(((640 ./ bin_f) / 130) ./ SF, 180 .* (0:n_ori-1) ./ n_ori, "SpatialAspectRatio", sar, "SpatialFrequencyBandwidth", sfb);

    mag = imgaborfilt(img, gb);

    mag = mag((size(mag, 1)./2+(-18:18)), (size(mag, 2)./2+(-18:18)), :);
    mag = mag ./ (max(vec(mag)));

    % in the model I first scale the magnitudes to match the dF/F
    % I have to do the same here, of the fitted parameters don't properly!
    mag = scale_factor .* mag;

    results(sar_cnt, :) = reshape(mag(19, 19, :), 1, n_ori);
    
    % subplot(3, n, n+sar_cnt);
    subplot(4, n, n+sar_cnt);

    kernel = imag(gb(1).SpatialKernel);

    try
        kernel = kernel(ceil(size(kernel, 1)./2) + (-n_px:n_px), ceil(size(kernel, 2) ./ 2) + (-n_px:n_px));
    catch; end
    imagesc(kernel); caxis([-1 1]); colormap("gray"); axis image; title("Spatial aspect"+newline+"ratio: " + sprintf("%0.2f", sar));
    xticks([]); yticks([]);
    drawnow;
end
results_backup = results;


%%
results = results_backup(:, [ceil(n_ori/2)+1:n_ori, 1:(n_ori/2)])';
results = results ./ nanmax(results, [], 1);

set(gcf, "units", "pixels");
set(gcf, "OuterPosition", [1 1 1920 1080]);
x = linspace(0, 180, n_ori+1); x(end) = [];
for i = 1:n
    subplot(4, n, i+2*n); hold off;

    plot(x, results(:, i), "-k"); hold on;
    temp_ids = results(:, i) > 0.5;
    title("HWHM: " + round(sum(results(:, i) > 0.5, 1) .* (180/n_ori) ./ 2) + "°");

    try
        FWHM_ids = [find(temp_ids, 1, "first"), find(temp_ids, 1, "last")];
        plot(x(FWHM_ids), [0.5, 0.5], ":r");
    catch
    end

    xlim([0 180]); xticks([0:45:180]); 
    xticklabels([-90:45:90]); xlabel("Orientation (°)");
    ylim([0 1]);

    if i == 1
        ylabel("Gabor response"+newline+"magnitude");
    end
end


%% set fontSize, LineWidth, PageSize
drawnow; fontsize_for_every_figure_element(7, "normal", 0.5);
set(findobj(gcf,'Type','line'), "LineWidth", 1.5);
drawnow; set(gcf, 'Units', 'centimeters', 'Position', [.5, 2, 7.7, 14], 'PaperUnits', 'centimeters', 'PaperSize', [21, 29.7]); drawnow;

%% save
folder = fullfile(pwd, "determine_aspect_ratio"); make_dir(folder);
fName = fullfile(folder, "Gabor_AspectRatio_that_gives_FWHM_tuning_width_of_44deg_for_suppl_figure");
exportgraphics(gcf, fName+".png", "Resolution", 600);
savefig(gcf, fName+".fig");
saveas(gcf, fName+".svg");


