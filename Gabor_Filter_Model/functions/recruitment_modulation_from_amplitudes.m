% We now ues the nonlinear function for the amplitude predictions already
% so here we just use a linear-regression
modIdx = @(a, b) squeeze((a-b)./(a+b));

try close(1); catch; end

if raw_gabors
    model_str = "raw_gabors_sar"+SpatialAspectRatio;
end


if fit_outerSurroundLimit
    model_str = model_str + "_outerSurroundLim";
end

if use_SurrSupp == 2
    model_str = model_str + "_divisiveSupp";
end

figName = model_str + "_matchMeans" + match_means_R2_only;


for plot_SurrSupp = 0:1


%% Load the measured recruitment indices
% For the individual session R2's
temp = table2array(readtable("ORI_recruitment_modulation.csv"));
ORI_recruitment_idx = temp(2:end, :);

% remove outliers
for i = 1:size(ORI_recruitment_idx, 2)
    [~, idx] = rmoutliers(ORI_recruitment_idx(:, i), 1);
    ORI_recruitment_idx(idx, i) = NaN;
end


%% %% %% end
load("prediction_files\amplitude_predictions_"+model_str+".mat", "amplitude_predictions", "amplitude_predictions_selected_oris", "raw_gabor_prediction", "raw_gabor_prediction_selected_oris");
if ~plot_SurrSupp
    if exist("raw_gabor_prediction_selected_oris", "var")
        amplitude_predictions_selected_oris = raw_gabor_prediction_selected_oris;
    else
        amplitude_predictions_selected_oris = raw_gabor_prediction;
    end
end
cl = [38 101 168] ./ 255;


% Plot to determine how to fit the data.
figure(1); fig_size(gcf, [.1 .1 .3 .7]); hold off;
x = modIdx(amplitude_predictions_selected_oris(2, :), amplitude_predictions_selected_oris(1, :))';
y = nanmedian(ORI_recruitment_idx, 1)';


% Fit the measured recruitment with the fit amplitude prediction (without intercept)
beta = x\y;
yPredict1 = beta*x;



%% 
try
    % If I have the center already
    if plot_SurrSupp
        load("prediction_files\recruitment_predictions_"+figName+"_center_only.mat", ...
             "ori_recruitment_predicted", "R2_ori");
        ori_recruitment_predicted_center_only = ori_recruitment_predicted;
        R2_ori_center_only = R2_ori;
    end
catch; end

predicted_cell_ratios = reshape(yPredict1, [5, 1])';

% Now also the ratios
figure(3); fig_size(gcf, [0.05 0.05 0.4 0.55])
y = nanmedian(ORI_recruitment_idx, 1);
x = predicted_cell_ratios;

ori_recruitment_measured = y;
ori_recruitment_predicted = x;

% match the means for computing the R2
predict_R2 = ori_recruitment_predicted - nanmean(ori_recruitment_predicted) + nanmean(y);
R2_ori = 1 - (nansum((y - predict_R2).^2) ./ nansum((y - nanmean(y)).^2));


subplot(2, 3, 4:5); hold off;
plot([0 6], [0 0], "--k"); hold on;
if plot_violins
    v = violinplot(ORI_recruitment_idx);
    for i = 1:length(v); v(i).ViolinColor = cl(end, :); end
end
text(1.2, 0.45, "R² = "+sprintf("%0.3f", R2_ori), "color", temp_cl_supp);
% try
if plot_SurrSupp
    h = plot(ori_recruitment_predicted_center_only', "-", "LineWidth", 2, "color", temp_cl_tuning, "DisplayName", "Linear-fit(Tuning only model)");
    text(1.2, 0.55, "R² = "+sprintf("%0.3f", R2_ori_center_only), "color", temp_cl_tuning);
    h(end+1) = plot(ori_recruitment_predicted', "-", "color", temp_cl_supp, "LineWidth", 2, "DisplayName", "Linear-fit(Supression model)");
    xticks(1:5); xticklabels(-45:22.5:45); xlabel("Orientation tuning peak (°)")
    ylabel("Ori recruitment  modulation");
    ylim([-0.6 0.6]); yticks(-0.6:0.3:0.6);
    xlim([0.5 5.5]);
    legend(h, "Position", [0.694 0.354 0.17 0.09]);
end
% catch; end


figure(3)
if ~plot_SurrSupp
    % I now loop over both settings and save the combined one in the other
    % (second) case
    make_dir(fullfile(pwd, "prediction_files"));
    save("prediction_files\recruitment_predictions_"+figName+"_center_only.mat", ...
         "ori_recruitment_measured", ...
         "ori_recruitment_predicted", ...
         "R2_ori", ...
         "beta", ...
         "-v7.3");
else
    make_dir(fullfile(pwd, "prediction_files"));
    save_folder = fullfile(pwd, "Model_fit_plots"); make_dir(save_folder);
    fName = "Recruitment_prediction_"+figName;

    savefig(gcf, fullfile(save_folder, fName+".fig"));
    exportgraphics(gcf, fullfile(save_folder, fName+".png"), "Resolution", 300);
    
    % save the combined data
    save("prediction_files\recruitment_predictions_"+figName+".mat", ...
         "ori_recruitment_measured", ...
         "ori_recruitment_predicted", ...
         "ori_recruitment_predicted_center_only", ...
         "R2_ori", ...
         "R2_ori_center_only", ...
         "beta", ...
         "-v7.3");

    try
        % try clean up and remove the intermediate file
        delete("prediction_files\recruitment_predictions_"+figName+"_center_only.mat");
    catch; end
end

end  % END SurrSupp loop



