%% Loading the relevant variables
% Step 1: get the processed data from Zenodo [http://doi.org/10.5281/zenodo.14605879]. 
load('Responses_all_sessions.mat');  % This loads the table of responses for the 13 recording session of the Niell mice (see Figure 4)

%% some further processing steps 
session_nr = 1:size(T_responses,2);  % session number (1:13)
cellNR = cat(1,T_responses(session_nr(:)).cellNR); % number of detected cells
resp_cell_nr = cat(1,T_responses(session_nr(:)).resp_cell_nr); % number of responding cells 
respStrength_all = cat(1,T_responses(session_nr(:)).respStrength_all); % response amplitudes of all cells
respStrength_p = cat(1,T_responses(session_nr(:)).respStrength_p); % Mann Whitney U test against baseline (defines how significantly responsive a cells is
sig_responders = cat(1,T_responses(session_nr(:)).sig_responders); % identities of the significantly responding cells and their response amplitude
resp_AUC = cat(1,T_responses(session_nr(:)).resp_AUC);% AUC of response
absolutes = cat(1,T_responses(session_nr(:)).absolutes); % absolute response (during the stimulus period)
minMaxTuning = cat(1,T_responses(session_nr(:)).minMaxTuning); % amplitude of tuning peak (difference between the peak and the baseline)
tuning_peaks = cat(1,T_responses(session_nr(:)).tuning_peaks); %  used interpolation to get a more fine-grained tuning curve and found the peak area there
Center_Rf = cat(2,T_responses(session_nr(:)).Center_Rf); % logical 0 means receptive field not at the center stimulus and 1 means receptive fields of the neuron is at center stimulus location

%Here we find cells that are responding to the center and any of the center stimuli
fullfieldIdx = 1:4;  % Indices of the FullField stimuli (1 - narrow, 2- SF, 3- ORI, 4 - mixed)
centerIdx = 5:8;     % Indices of the center (15° aperture) stimuli (1 - narrow, 2- SF, 3- ORI, 4 - mixed)
pThresh = 0.025;     % significance threshold 
centerRespones = (respStrength_p(:, centerIdx) < pThresh) & (resp_AUC(:, centerIdx) > 0.5); %logical index to get center responses
center_sig_responders = any(centerRespones,2); %this is for significant center responses in ANY condition
%The following variable holds the identity of each cell pertaining to a certain mouse (1-5)
animal_id_all_cells = vertcat((ones(size(T_responses(1).sig_responders,1) + size(T_responses(2).sig_responders,1),1)),...
    (ones(size(T_responses(3).sig_responders,1)+ size(T_responses(4).sig_responders,1),1)*2), ...
    (ones(size(T_responses(5).sig_responders,1) + size(T_responses(6).sig_responders,1)+ size(T_responses(7).sig_responders,1),1)*3),...
    (ones(size(T_responses(8).sig_responders,1) + size(T_responses(9).sig_responders,1) + size(T_responses(10).sig_responders,1),1)*4),...
    (ones(size(T_responses(11).sig_responders,1) + size(T_responses(12).sig_responders,1)+ size(T_responses(13).sig_responders,1),1)*5));
%picking up the center RF stimulus responsive cells that also respond to (any of) the center stimuli
useIdx = Center_Rf' & center_sig_responders;
animal_id = animal_id_all_cells(useIdx);

% cell responses and surround modulation index
%FF responses 
ff_response = respStrength_all(useIdx,fullfieldIdx);
ff_response_abs = absolutes(useIdx,fullfieldIdx);
ff_response_auc = resp_AUC(useIdx,fullfieldIdx);
%center responses
center_response = respStrength_all(useIdx,centerIdx);
center_response_abs = absolutes(useIdx,centerIdx);
center_response_auc = resp_AUC(useIdx,centerIdx);
%surround modulation index
surroundModulation = (ff_response_abs - center_response_abs) ./ (ff_response_abs + center_response_abs);
%response enrichment index
for i = 1:4
    enrichment_ff(:,i) = (ff_response_abs(:,i) - ff_response_abs(:,1))./ (ff_response_abs(:,i) + ff_response_abs(:,1));
end

%% Control figures (should look similar to Figure 4 d,e,f)
figure(1) % Fullfield responses to all stimuli
violinplot(rmoutliers(ff_response)*100); 
xticklabels({'narrow', 'SF', 'ORI','mixed'}); 
ylabel("Response amplitude ΔF/F(%)");
cardinal_line(true,false, median(rmoutliers(ff_response(:,1)*100)),'--', 'k', 2);
niceFigures; 

figure(2)  % Center responses to all stimuli
violinplot(rmoutliers(center_response)*100); 
xticklabels({'narrow', 'SF', 'ORI','mixed'});
ylabel("Response amplitude ΔF/F(%)");
cardinal_line(true,false, median(rmoutliers(center_response(:,1)*100)),'--', 'k', 2);
niceFigures; 

figure(3)  % Surround modulation of all stimuli
violinplot(rmoutliers(surroundModulation)*100); 
xticklabels({'narrow', 'SF', 'ORI','mixed'});
ylabel("Surround Modulation index");
cardinal_line(true,false, median(rmoutliers(surroundModulation(:,1)*100)),'--', 'k', 2);
niceFigures; 


%% The model stats  
tuneIdx = minMaxTuning(useIdx) > median(minMaxTuning); %get well tuned neurons
norm_tuning = abs(tuning_peaks(useIdx)-200);
useAnimals = animal_id_all_cells(useIdx);
%enrichments
enrichment_ff_sf = enrichment_ff(:,2);
enrichment_ff_ori = enrichment_ff(:,3);
enrichment_ff_mixed = enrichment_ff(:,4);
%surround indices
delta_surround_ori = (surroundModulation(:,1));
delta_surround_sf = (surroundModulation(:,1));
delta_surround_mixed = (surroundModulation(:,1));
% orientation tuning used as regressor in the model
predTuning = norm_tuning(tuneIdx);  %predictor 2 - the orientation tuning of the cell

T_vals = zeros(3,3);
% model for SF
targData_SF = enrichment_ff_sf(tuneIdx); %target data to be predicted is response enrichment
predSurround_SF = surroundModulation(tuneIdx,1); %predictor1 - surround modulation 
predCenter_SF = center_response_auc(tuneIdx,2) - center_response_auc(tuneIdx,1);  %predictor3 - center response enrichment
[pVals, tStats, threeVarModel_SF] = LME_compare_3vars(targData_SF, predSurround_SF, predTuning, predCenter_SF, useAnimals(tuneIdx));
predVar_SF = corr2(targData_SF, fitted(threeVarModel_SF))^2; 
disp('================')
fprintf('3 variable model T-stats for SF:\n t-Stat_surround: %f\n t-Stat_tuning: %f\n t-Stat_center: %f\n', tStats(1),tStats(2),tStats(3))
disp('----------------')
fprintf('3 variable model p-Vals for SF:\n p-Val_surround: %f\n p-Val_tuning: %f\n p-Val_center: %f\n', pVals(1),pVals(2),pVals(3))
fprintf('Explained variance (R^2) full model: %f\n', predVar_SF);
disp('================')
T_vals(:,1) = tStats;

%%model for ORI
targData_ORI = enrichment_ff_ori(tuneIdx);
predSurround_ORI = surroundModulation(tuneIdx,1);
predCenter_ORI = center_response_auc(tuneIdx,3) - center_response_auc(tuneIdx,1);
[pVals, tStats, threeVarModel_ORI] = LME_compare_3vars(targData_ORI, predSurround_ORI, predTuning, predCenter_ORI, useAnimals(tuneIdx));
predVar_ORI = corr2(targData_ORI, fitted(threeVarModel_ORI))^2;
disp('================')
fprintf('3 variable model T-stats for ORI:\n t-Stat_surround: %f\n t-Stat_tuning: %f\n t-Stat_center: %f\n', tStats(1),tStats(2),tStats(3))
disp('----------------')
fprintf('3 variable model p-Vals for ORI:\n p-Val_surround: %f\n p-Val_tuning: %f\n p-Val_center: %f\n', pVals(1),pVals(2),pVals(3))
fprintf('Explained variance (R^2) full model: %f\n', predVar_ORI);
disp('================')
T_vals(:,2) = tStats;

%%model for MIXED
targData_MIX = enrichment_ff_mixed(tuneIdx);
predSurround_MIX = surroundModulation(tuneIdx,1);
predCenter_MIX = center_response_auc(tuneIdx,4) - center_response_auc(tuneIdx,1);
[pVals, tStats, threeVarModel_MIX] = LME_compare_3vars(targData_MIX, predSurround_MIX, predTuning, predCenter_MIX, useAnimals(tuneIdx));
predVar_MIX = corr2(targData_MIX, fitted(threeVarModel_MIX))^2;
disp('================')
fprintf('3 variable model T-stats for MIX:\n t-Stat_surround: %f\n t-Stat_tuning: %f\n t-Stat_center: %f\n', tStats(1),tStats(2),tStats(3))
disp('----------------')
fprintf('3 variable model p-Vals for MIX:\n p-Val_surround: %f\n p-Val_tuning: %f\n p-Val_center: %f\n', pVals(1),pVals(2),pVals(3))
fprintf('Explained variance (R^2) full model: %f\n', predVar_MIX);
disp('================')
T_vals(:,3) = tStats;

%% The polar plot 
figure(5)
T_vals = abs(T_vals)./(max(max(abs(T_vals))));
angles = [pi/2, 7*pi/6, 11*pi/6, pi/2]; % Angles for the innermost triangle
radii1 = vertcat(T_vals(:,1),T_vals(1,1));% Radii for the innermost triangle
radii2 = vertcat(T_vals(:,2),T_vals(1,2)); % Radii for the middle triangle
radii3 = vertcat(T_vals(:,3),T_vals(1,3)); % Radii for the outermost triangle
polarplot(angles, radii1, 'Color',  [0.5020    0.4314    0.6510],'LineWidth',1.5); hold on;
polarplot(angles, radii2, 'Color', [0.1882    0.3373    0.5529],'LineWidth',1.5); hold on;
polarplot(angles, radii3, 'Color', [ 0.5725    0.3529    0.2667],'LineWidth',1.5); hold off;
% the limits and aspect ratio for the polar plot
rmax = 1; % maximum radius
rlim([0 rmax]); rtickangle(-30);
thetaticks(90:120:330); % theta ticks
thetaticklabels({ 'Surround modulation', 'Orientation tuning','center response'});
legend({ 'SF', 'Ori','mixed'});


