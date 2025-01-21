
%%% Gabor model comparisons from figure 3 
% Step 1: get the processed data from Zenodo [http://doi.org/10.5281/zenodo.14605879]. 
load('Gabor_model_comparisons.mat')

%% general variables, this is for pooling accross all cells
narrow = MC_cell_data.narrow_response{1};
broad = MC_cell_data.broad_response{1};
cell_tuning_peak = MC_cell_data.tuning_peak{1};
all_session_ids = MC_cell_data.session_id{1};
sessions = 3:13;
response_ratios = nan(2, 5, length(sessions));
for session_id = 1:length(sessions)
    for i = 1:5
        idx = (cell_tuning_peak == i) & (all_session_ids == sessions(session_id));
        response_ratios(:, i, session_id) = [nanmedian(narrow(idx)), ...
                                             nanmedian(broad(idx))];
    end
end
response_ratios_individual = response_ratios;
narrow_ampl = squeeze(response_ratios_individual(1, :, :))'*100;
broad_ampl = squeeze(response_ratios_individual(2, :, :))'*100;
amp_mod = (broad_ampl - narrow_ampl)./(broad_ampl + narrow_ampl);


%% looking at the pooled accross all cells version of the plots as a first look

%raincloudPlot(data, [color_rgb], ks_density_bandwidth)
raincloudPlot(narrow_ampl,  [ 0.4392    0.6902    0.8392],0.5); hold on;
ylabel('Response amplitude ΔF/F(%)'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''}); ylim([0 4]); xlim([0 6]); 
set(gca, 'XGrid', 'off', 'YGrid', 'off'); 
title('narrow band amplitudes');

% looking at the pooled accross all cells version
raincloudPlot(broad_ampl,  [ 0.0314    0.3176    0.6118],0.5); hold on;
ylabel('Response amplitude ΔF/F(%)'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''}); ylim([0 4.5]); xlim([0 6]); 
set(gca, 'XGrid', 'off', 'YGrid', 'off'); 
title('broad band amplitudes');

% looking at the pooled accross all cells version
raincloudPlot(amp_mod,  [ 0.0314    0.3176    0.6118],0.1); hold on;
ylabel('Response modulation'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''}); ylim([-0.6 0.6]); xlim([0 6]); 
set(gca, 'XGrid', 'off', 'YGrid', 'off'); cardinal_line(true,false, 0,'--', 'k', 1);
title('Amplitude modulation');

%% Plotting the by session data as well as the model predictions as seen in Figure 3 

raincloudPlot(narrowband_amplitudes,  [ 0.4392    0.6902    0.8392],0.7); hold on;
plot(tuning_only_model_amplitude_predictions(1,:), Color ='#F0D101', Marker ='o' , LineWidth=3); ylim([0 4]); xlim([0 6]); 
plot(suppression_model_amplitude_predictions(1,:), Color ='#671A19', Marker ='o', LineWidth=3); ylim([0 4]); xlim([0 6]); 
ylabel('Response amplitude ΔF/F(%)'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''});
set(gca, 'XGrid', 'off', 'YGrid', 'off'); 
title('narrow band amplitudes');
% Statistics: doing the LME model comparison against 0 for each column 
% The first row of p shows the p values and the 2nd row shows the T-statistics
for i = 1:5
 [ p_narrow(1,i),  p_narrow(2,i)] = LME_compare(narrowband_amplitudes(:,i),zeros(length(narrowband_amplitudes),1),session_animal_ids');
end

raincloudPlot(broadband_amplitudes,  [ 0.0314    0.3176    0.6118],0.7); hold on;
plot(tuning_only_model_amplitude_predictions(2,:), Color ='#F0D101', Marker ='o' , LineWidth=3); ylim([0 4]); xlim([0 6]); 
plot(suppression_model_amplitude_predictions(2,:), Color ='#671A19', Marker ='o', LineWidth=3); ylim([0 4]); xlim([0 6]); 
ylabel('Response amplitude ΔF/F(%)'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''});
set(gca, 'XGrid', 'off', 'YGrid', 'off'); 
title('broad band amplitudes');
% Statistics: doing the LME model comparison against 0 for each column 
% The first row of p shows the p values and the 2nd row shows the T-statistics
for i = 1:5
 [ p_broad(1,i),  p_broad(2,i)] = LME_compare(broadband_amplitudes(:,i),zeros(length(broadband_amplitudes),1),session_animal_ids');
end

raincloudPlot(amplitude_modulation,  [ 0.0314    0.3176    0.6118], 0.1); hold on;
plot(modulation_prediction_center_only_model, Color ='#F0D101', Marker ='o' , LineWidth=3); ylim([-0.6 0.6]); xlim([0 6]); 
plot(modulation_prediction_suppression_model, Color ='#671A19', Marker ='o', LineWidth=3); ylim([-0.6 0.6]); xlim([0 6]); 
ylabel('Response modulation'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''});
set(gca, 'XGrid', 'off', 'YGrid', 'off'); 
cardinal_line(true,false, 0,'--', 'k', 1);
title('Response modulation');
% Statistics: doing the LME model comparison against 0 for each column 
% The first row of p shows the p values and the 2nd row shows the T-statistics
for i = 1:5
 [ p_amplitude_modulation(1,i),  p_amplitude_modulation(2,i)] = LME_compare(amplitude_modulation(:,i),zeros(length(amplitude_modulation),1),session_animal_ids');
end

raincloudPlot(recruitment_modulation,  [ 0.0314    0.3176    0.6118],0.1); hold on;
plot(recruitment_tuning_only_model_prediction, Color ='#F0D101', Marker ='o' , LineWidth=3); ylim([-0.6 0.6]); xlim([0 6]); 
plot(recruitment_supression_model_prediction, Color ='#671A19', Marker ='o', LineWidth=3); ylim([-0.6 0.6]); xlim([0 6]); 
ylabel('Recruitment modulation'); xlabel(["Orientation tuning peak [°]" ]); 
xticklabels( {'-45°', '-25°', '0°','25°','45°' ,''});
set(gca, 'XGrid', 'off', 'YGrid', 'off'); 
cardinal_line(true,false, 0,'--', 'k', 1);
title('Recruitment modulation');
% Statistics: doing the LME model comparison against 0 for each column 
% The first row of p shows the p values and the 2nd row shows the T-statistics
for i = 1:5
 [ p_recruitment_modulation(1,i),  p_recruitment_modulation(2,i)] = LME_compare(recruitment_modulation(:,i),zeros(length(recruitment_modulation),1),session_animal_ids');
end


