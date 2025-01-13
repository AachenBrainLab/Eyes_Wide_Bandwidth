% In this script we fit the measured amplitudes towards narror and Ori
% stimuli with the responses of Gabor filters towards the stimulus frames
modIdx = @(a, b) squeeze((a-b)./(a+b));

% for reproducibility
rng("default");

SF_ids = 1:5;

% I want to put this in a overarching script running them all in sequence
if ~exist("SpatialAspectRatio", "var")
    SpatialAspectRatio = 0.55;
end


%% load the measured response amplitudes to fit
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
response_ratios = nanmedian(response_ratios_individual, 3);
data2fit = response_ratios;

% For the individual session R2's
temp = table2array(readtable("ORI_amplitude_modulation.csv"));
ORI_amplitude_modulation = temp(2:end, :);


%%
cl = [214 180 148; 128 110 166; 38 101 168; 146 90 68] ./ 255;
if (exist("backup_gabor_for_rerun", "var")) && (currently_loaded_sar == SpatialAspectRatio)  % reload
    % If run multiple times I dont want to load all the data again. Speeds
    % it up a lot
    gabor_responses = backup_gabor_for_rerun;
else
    % if (~exist("gabor_responses", "var")) || (currently_loaded_sar == SpatialAspectRatio)  % reload
    addpath("finding_gabor_max_responses");  % This is where I load the sin-wave responses for normalization from
    load("finding_gabor_max_responses\Matched_gabors_max_responses.mat")

    load("gabor_responses_Matlab_SpatialAspectRatio_"+SpatialAspectRatio+"_all4_stims_SpatialMaps.mat", "gabor_responses")
    gabor_responses = gabor_responses ./ reshape(max_responses, [1 1 length(max_responses) 1 1 1]);
    gabor_responses = gabor_responses(:, :, :, end:-1:1, :, :);
    gabor_responses = gabor_responses(:, :, :, [91:180 1:90], :, :);

    % keep this as backup in case of rerunning:
    backup_gabor_for_rerun = gabor_responses;
    currently_loaded_sar = SpatialAspectRatio;  % Since the 0.37 and the 0.55 have different data, I need to keep track of that!
end

% remove conditions with slightly skewed orientation distributions.
gabor_responses(:, :, :, :, :, vec((1:120)' + [2 4 7] .* 120)) = nan;

% ids of orientations, corresponding to [-45deg, -22.5deg, 0deg, 22.5deg, 45deg]
ori_ids = [45, 67, 90, 112, 135]+1;

if length(SF_ids) < 5
    % subselect spatial frequencies
    gabor_responses = gabor_responses(:, :, SF_ids, :, :, :);
end


% define activation fuction to use
activation_fun = @activation_fun_raw_gabors;
model_str = "raw_gabors";

% construct model_str for saving
model_str = model_str + "_sar"+SpatialAspectRatio;

if fit_outerSurroundLimit
    model_str = model_str + "_outerSurroundLim";
end

if use_SurrSupp == 2
    model_str = model_str + "_divisiveSupp";
end

figName = model_str;



%% Fit the model on [Narrow, Ori]
% This subselects these conditions [1: Narrow, 3: Ori]
SF_amplitude_modulation = nan;
response_ratios = response_ratios([1 3], :);
response_ratios_individual = response_ratios_individual([1 3], :, :);
data2fit = data2fit([1 3], :);
gabor_responses = gabor_responses(:, :, :, :, [1 3], :);
cl = cl([1 3], :);

%%
% For the "Raw Gabor" comparision I need to match the Normalized Gabor
% responses to the dF/F range of our data, otherwise I can't do the 
% surroundSupression fit properly.

% I compute this directly after loading the data.
% Also save the scaled_raw_gabors.
raw_gabor_prediction = squeeze(nanmean(nanmean(gabor_responses(19, 19, :, :, :, :), 6), 3))';

raw_gabor_prediction_selected_oris = raw_gabor_prediction(:, ori_ids);
scale_factor = mean(mean(data2fit)) ./ mean(mean(raw_gabor_prediction_selected_oris));
raw_gabor_prediction_selected_oris = scale_factor .* raw_gabor_prediction_selected_oris;

raw_gabor_prediction = scale_factor .* raw_gabor_prediction;

gabor_responses = scale_factor .* gabor_responses;


% raw center only response
if size(raw_gabor_prediction_selected_oris, 2) == 3
    raw_amplitude_modulation_idx_ori = modIdx(raw_gabor_prediction_selected_oris(:, 3), raw_gabor_prediction_selected_oris(:, 1));
elseif size(raw_gabor_prediction_selected_oris, 2) == 2
    raw_amplitude_modulation_idx_ori = modIdx(raw_gabor_prediction_selected_oris(:, 2), raw_gabor_prediction_selected_oris(:, 1));
end


%%
figure(2);
subplot(1, 2, 1); hold off;
for i = 1:size(cl, 1)
    plot(-90+1:90, squeeze(nanmean(nanmean(gabor_responses(19, 19, :, :, i, :), 3), 6)), "color", cl(i, :), "LineWidth", 2); hold on;
end
xlabel("Orientation");
xticks(-90:45:90); title("Gabor"+newline+"responses");
ylabel("Response"+newline+"magnitude");
xlim([-90 90]);
ylim([0 0.03]);


%%
if rerun  % rerun fit
    max_iter = 100;
    [parameters, fval] = fit_MC_amplitudes(activation_fun, gabor_responses(:, :, :, ori_ids, :, :), data2fit, max_iter, use_SurrSupp, 1, fit_outerSurroundLimit);

    % save fit for later use
    save("amplitude_fit_parameters_"+model_str+".mat", "scale_factor", "parameters", "fval");
else
    % simply load if no refitting is desired
    clearvars parameters
    load("amplitude_fit_parameters_"+model_str+".mat", "scale_factor", "parameters", "fval");
end


%% plot responses with fitted surround suppression
% unpack parameters
surround_threshold = parameters(1);
surround_weight = parameters(2);
surround_offset = parameters(3);
if length(parameters) > 3
    surround_threshold_outer = parameters(4);
else
    surround_threshold_outer = inf;
end
if length(parameters) > 4
    center_weight = parameters(5);
else
    center_weight = 1;
end

if use_SurrSupp == 2
    % this is enforced during the fit, but the overwrite is not saved, so
    % here to match the fit:
    if surround_threshold_outer < surround_threshold; surround_threshold_outer = surround_threshold + 1; end
end


%%
% Save the Amplitude-predictions:
amplitude_predictions = activation_fun(gabor_responses, parameters, use_SurrSupp);
amplitude_predictions_selected_oris = amplitude_predictions(:, ori_ids);
save("prediction_files\amplitude_predictions_" + model_str + ".mat", ...
     "amplitude_predictions", ...
     "amplitude_predictions_selected_oris", ...
     "raw_gabor_prediction", ...
     "raw_gabor_prediction_selected_oris", ...
     "scale_factor", "-v7.3");


%%
center_only_responses = raw_gabor_prediction_selected_oris;
ori_fit_center_only = modIdx(center_only_responses(2, :), center_only_responses(1, :));

% define model associated colors
temp_cl_supp = [103 26 25] ./ 255;
temp_cl_tuning = [240 209 1] ./ 255;

% set up figure
figure(3); fig_size(gcf, [.1 .1 .4 .6]);
set(gcf, "Units", "pixels", "OuterPosition", [50 50 600 550]);

% only retain selected oris
amplitude_predictions = amplitude_predictions(:, ori_ids);

%
subplot(2, 3, 1:2); hold off;
plot([0.5, 5.5], [0 0], "--k", "HandleVisibility", "off"); hold on;

measured_data = modIdx(response_ratios(2, :),  response_ratios(1, :));

% prediction with surround suppression
predicted_data = modIdx(amplitude_predictions(2, :), amplitude_predictions(1, :));

% fetch measure modulation index
session_wise_idx = ORI_amplitude_modulation;

% remove outliers
for i = 1:size(session_wise_idx, 2)
    [~, idx] = rmoutliers(session_wise_idx(:, i), 1);
    session_wise_idx(idx, i) = NaN;
end


% compute R2 with matched mean to focus on the orientation-dependence
predict_R2 = predicted_data - nanmean(predicted_data) + nanmean(nanmean(measured_data));
R2 = 1 - (nansum(power(measured_data - predict_R2, 2)) ./ nansum(power(measured_data - nanmean(measured_data), 2)));

% plot measured data
v = violinplot(session_wise_idx);
for i = 1:length(v); v(i).ViolinColor = cl(2, :); end


predicted_data2 = ori_fit_center_only;
h = plot(predicted_data2, "-", "color", temp_cl_tuning, "LineWidth", 2, "DisplayName", "Tuning only model (raw Gabors)");

% compute R2 with matched mean to focus on the orientation-dependence
predict_R2 = predicted_data2 - nanmean(predicted_data2) + nanmean(nanmean(measured_data));
R2_2 = 1 - (nansum(power(measured_data - predict_R2, 2)) ./ nansum(power(measured_data - nanmean(measured_data), 2)));


% rename for clarity when saving
amplitude_measured_ori = measured_data;
amplitude_prediction_center_only_ori = predicted_data2;
amplitude_prediction_center_surround_ori = predicted_data;


% plot
h(2) = plot(predicted_data, "-", "color", temp_cl_supp, "LineWidth", 2, "DisplayName", "SurroundSupression Model");
ylabel("ORI response modulation");
xlabel("Orientation tuning peak (°)"); xticks(1:5); xticklabels((-2:2) .* 22.5);
xlim([0.5 5.5])
ylim([-0.35 0.39]);
text(1.2, 0.45, "R²=" + sprintf("%0.3f", R2), "color", temp_cl_supp);
text(1.2, 0.55, "R²=" + sprintf("%0.3f", R2_2), "color", temp_cl_tuning);

% set axis lims and ticks
ylim([-0.6 0.6]);
yticks(-0.6:0.3:0.6);

% define legend spoition for better control
legend(h, "Position", [0.707 0.83 0.17 0.09]);

% renamve for clarity when saving
R2_ori = R2;
R2_ori_center_only = R2_2;


%% save results
disp("Saved: "+ "prediction_files\amplitude_modulation_idx_predictions_"+figName+".mat")
save("prediction_files\amplitude_modulation_idx_predictions_"+figName+".mat", ...
     "scale_factor", ...
     "amplitude_measured_ori", ...
     "amplitude_prediction_center_only_ori", ...
     "amplitude_prediction_center_surround_ori", ...
     "R2_ori", "R2_ori_center_only", ...
     "-v7.3");







%% functions
function activation = activation_fun_raw_gabors(data, x, use_SurrSupp)
    surround_threshold = x(1);
    surround_weight    = x(2);
    offset = x(3);
    if length(x) > 3
        surround_threshold_outer = x(4);
    else
        surround_threshold_outer = 15;
    end

    if length(x) > 4
        center_weight = x(5);
    else
        center_weight = 1;
    end
    
    if (size(data, 1) == 37) && (size(data, 2) == 37)
        data = add_surround_supression(data, surround_threshold, surround_weight, offset, surround_threshold_outer, use_SurrSupp, center_weight);
    end

    % redundant but just in case
    % rectify responses - negative responses don't make sense
    data(data < 0) = 0;

    % averaging over dim 1 and 2 is not necessary anymore, I keep it
    % for compatibility. 
    activation = squeeze(nanmean(nanmean(nanmean(nanmean(data, 1), 2), 3), 6))';
end







