use_SurrSupp = 2;
fit_outerSurroundLimit = 1;
model_str = "raw_gabors" + "_sar"+SpatialAspectRatio;

if fit_outerSurroundLimit
    model_str = model_str + "_outerSurroundLim";
end
if use_SurrSupp == 2
    model_str = model_str + "_divisiveSupp";
end


%% Load 2P data
load("MC_cell_data.mat");
n = MC_cell_data.narrow_response{1};
b = MC_cell_data.broad_response{1};
s = MC_cell_data.sf_response{1};
cell_tuning_peak = MC_cell_data.tuning_peak{1};
all_session_ids = MC_cell_data.session_id{1};

sessions2use = 3:13;
response_ratios = nan(3, 5, length(sessions2use));
for session_id = 1:length(sessions2use)
    for i = 1:5
        idx = (cell_tuning_peak == i) & (all_session_ids == sessions2use(session_id));
        response_ratios(:, i, session_id) = [nanmedian(n(idx)), ...
                                             nanmedian(s(idx)), ...
                                             nanmedian(b(idx))];
    end
end
response_ratios_individual = response_ratios;
response_ratios = nanmean(response_ratios, 3);


%%
cl = [214 180 148; 128 110 166; 38 101 168; 146 90 68] ./ 255;
fit_cl = [240 209 1; 103 26 25] ./ 255;

for mode = ["tuning only", "supression"]
load("prediction_files\amplitude_predictions_" + model_str + ".mat", ...
     "amplitude_predictions", ...
     "amplitude_predictions_selected_oris", ...
     "raw_gabor_prediction", ...
     "raw_gabor_prediction_selected_oris", ...
     "scale_factor");

if strcmpi(mode, "tuning only")
    gabor_fit = raw_gabor_prediction;
    gabor_fit_selected_oris = raw_gabor_prediction_selected_oris;

    c_fit_cl = fit_cl(1, :);
elseif strcmpi(mode, "supression")
    gabor_fit = amplitude_predictions;
    gabor_fit_selected_oris = amplitude_predictions_selected_oris;

    c_fit_cl = fit_cl(2, :);
end

if size(gabor_fit, 1) == 2  % I changed the model to only try and fit narrow and Ori, no longer SF!
    % I can only do this once:
    if size(response_ratios, 1) > 2
        response_ratios_individual = response_ratios_individual([1 3], :, :);
        response_ratios = response_ratios([1 3], :);
        cl = cl([1 3], :);
    end
end


if strcmpi(mode, "tuning only")
    close all; fig_size(gcf, [.1 .1 .45 .3]);
end

x = -89:90;
marker_size = 3;


% Compute R2
amp_fit = cat(2, gabor_fit_selected_oris(1, :), gabor_fit_selected_oris(end, :));
data2fit_individual_session = squeeze(cat(2, response_ratios_individual(1, :, :), response_ratios_individual(end, :, :)));
data2fit_mean = nanmean(data2fit_individual_session, 2);
data2fit_median = nanmedian(data2fit_individual_session, 2);
if 0  % mean
    % data2fit = squeeze(cat(2, response_ratios(1, :), response_ratios(end, :)));
    data2fit = data2fit_mean;
else  % median
    data2fit = data2fit_median;
end

data2fit_median = reshape(data2fit_median, [5, 2])';
data2fit_mean = reshape(data2fit_mean, [5, 2])';


if 0
    R2 = 1 - (nansum(nansum(power(data2fit_individual_session - amp_fit', 2))) ./ nansum(nansum(power(data2fit_individual_session - nanmean(nanmean(data2fit_individual_session)), 2))));
else
    temp_data = data2fit; % temp_data = nanmedian(temp_data, 2);
    R2 = 1 - (nansum(nansum(power(temp_data - amp_fit', 2))) ./ nansum(nansum(power(temp_data - nanmean(nanmean(temp_data)), 2))));
end


for stim_id = [1 size(cl, 1)]  % 2:Ori

if stim_id == 1; subplot(1, 2, 1); else; subplot(1, 2, 2); end 
if strcmpi(mode, "tuning only")
    v = violinplot(squeeze(response_ratios_individual(stim_id, :, :))'); hold on;
    for i = 1:5; v(i).ViolinColor = cl(stim_id, :); v(i).ViolinAlpha = 0.45; end
end

h = plot(1:5, gabor_fit_selected_oris(stim_id, :), "o", "color", c_fit_cl, "LineWidth", 2);
h.MarkerFaceColor = c_fit_cl;
h.MarkerSize = marker_size;
plot(1:5, gabor_fit_selected_oris(stim_id, :), "color", c_fit_cl, "LineWidth", 2);

xticks(1:5);
xticklabels(-45:22.5:45);
xlim([0.25 5.75]);
ylim([0 0.045]);
xlabel("Tuning Orientation (°)");
if stim_id == 1; ylabel("Amplitudes (\DeltaF/F)"); end


if stim_id == 1
    text(0.5, 0.040 + (strcmpi(mode, "tuning only") .* 0.0035), "R² = "+sprintf("%0.3f", R2), "color", c_fit_cl);
    title("narrow responses", "FontWeight", "normal");
else
    title("broad Ori responses", "FontWeight", "normal");
end

end

end
save_folder = fullfile(pwd, "Model_fit_plots"); make_dir(save_folder);
exportgraphics(gcf, fullfile(save_folder, "Model_fit_individual_amplitudes_"+model_str+".png"), "Resolution", 600);




