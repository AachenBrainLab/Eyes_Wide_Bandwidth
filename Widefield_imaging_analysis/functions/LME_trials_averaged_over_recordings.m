% load data
% load('C:\Users\nabbefeld\Documents\GitHub\Analysis_routines_for_multisensory_task\WF_MotionMule\Figure_V4_clearSkull_mice\all_traces.mat');
% n_conditions = 4;
% condition_keys = ["Narrow", "Broad SF-BW", "Braod Ori-BW", "Braod SF + Ori-BW", "Narrow 0.16cpd", "Broad Ori-BW 0.16cpd"];
% plot_violins_over_all_trials_response_integrals(all_traces_50, n_conditions, condition_keys(1:n_conditions));


load("Figure_V4_clearSkull_mice/all_traces.mat", "all_traces_50");
all_traces = all_traces_50;

%%
response_integrals = squeeze(nanmean(nanmean(all_traces(16:45, :, 1, 1:4, :), 1), 2));

a = response_integrals(1, :)';
b = response_integrals(2, :)';
c = response_integrals(3, :)';
d = response_integrals(4, :)';

temp_responses = vec(response_integrals');


% animal_labels   = repmat(vec(["Mouse1", "Mouse2", "Mouse3", "Mouse4"]), [8, 1]);
animal_labels   = repmat(vec([1, 2, 3, 4]), [8, 1]);
stimulus_labels = vec(repmat(["Narrow", "SF", "Ori", "Mixed"], [8, 1]));
% sf_bw           = vec(repmat(["Narrow", "SF", "Ori", "Mixed"], [4*80, 1]));
% ori_bw          = vec(repmat(["Narrow", "SF", "Ori", "Mixed"], [4*80, 1]));

tbl = cat(2, [a; b; c; d], animal_labels, stimulus_labels);


% 
% try close(2); catch; end
% figure(2);
% subplot(1, 2, 1);
% imagesc(tbl'); caxis([0 4])
% 
% tbl = tbl(~isnan(tbl(:, 1)), :);
% subplot(1, 2, 2);
% imagesc(tbl'); caxis([0 4])
% 
% Build a table
% tbl = table(tbl(:, 1), tbl(:, 2), tbl(:, 3), 'VariableNames', ["Responses", "AnimalID", "StimulusID"]);


% Wit hthe average trial per session this became obsolete
% use_ids = ~isnan(temp_responses);
% temp_responses = temp_responses(use_ids);
% animal_labels = animal_labels(use_ids);
% stimulus_labels = stimulus_labels(use_ids);

% % Match trial-counts
% all_animals = unique(animal_labels);
% n_trials = inf;
% for cAnimal = unique(all_animals)'
%     for cStim = unique(stimulus_labels)'
%         temp = sum(strcmpi(animal_labels, cAnimal) & strcmpi(stimulus_labels, cStim));
%         disp(cAnimal + " - " + cStim + ": " + temp)
%         if temp < n_trials; n_trials = temp; end
%     end
% end
% 
% % now select trials acordingly
% for cAnimal = unique(all_animals)'
%     for cStim = unique(stimulus_labels)'
%         % get all trials from this mouse in this condition
%         temp_ids = find(strcmpi(animal_labels, cAnimal) & strcmpi(stimulus_labels, cStim));
% 
%         if length(temp_ids) > n_trials
%             temp_ids = temp_ids(n_trials+1:end);
%         else
%             % this is important! If I already have the number of trials I
%             % want then skip this.
%             continue
%         end
% 
%         % select only the first n_trials from those
%         plot(isnan(temp_responses)); hold on; drawnow;
%         disp(cAnimal + " - " + cStim)
%         temp_responses(temp_ids) = nan;
%         plot(isnan(temp_responses)); hold off; drawnow;
%     end
% end
% 
% use_ids = ~isnan(temp_responses);
% temp_responses = temp_responses(use_ids);
% animal_labels = animal_labels(use_ids);
% stimulus_labels = stimulus_labels(use_ids);

%
tbl = table(temp_responses, animal_labels, stimulus_labels, 'VariableNames', ["Responses", "AnimalID", "StimulusID"]);
% tbl

% % lme = fitlme(tbl,'Responses~StimulusID+(StimulusID|AnimalID)');
% lme = fitlme(tbl,'Responses ~ 1 + StimulusID + AnimalID + (1 + StimulusID | AnimalID)');
% lme
% 
% 

% sanity check!
% condition1 = "Narrow";

test_comparrisions = ["Narrow", "SF"; ...
                      "Narrow", "Ori"; ...
                      "Narrow", "Mixed"; ...
                      "SF", "Ori"; ...
                      "Ori", "Mixed"];
all_pVals = zeros(size(test_comparrisions, 1), 1);
all_tStat = zeros(size(test_comparrisions, 1), 1);
all_trials = zeros(size(test_comparrisions, 1), 1);
for comparrision_id = 1:size(test_comparrisions, 1)
    condition1 = test_comparrisions(comparrision_id, 1);
    condition2 = test_comparrisions(comparrision_id, 2);
%     if  ~all(strcmpi(tbl(strcmpi(tbl.StimulusID, condition1), :).AnimalID, ...
%                      tbl(strcmpi(tbl.StimulusID, condition2), :).AnimalID))
%         error("Unexpected Missmatch!");
%     end
    
    [pVal_cStim, tStat_cStim, fullmodel, modelCompare] = ...
        LME_compare(tbl(strcmpi(tbl.StimulusID, condition1), :).Responses, ...
                    tbl(strcmpi(tbl.StimulusID, condition2), :).Responses, tbl(strcmpi(tbl.StimulusID, "Narrow"), :).AnimalID);
    all_pVals(comparrision_id) = pVal_cStim;
    all_tStat(comparrision_id) = tStat_cStim;

    all_trials(comparrision_id) = ...
        length(tbl(strcmpi(tbl.StimulusID, condition1), :).Responses) + ...
        length(tbl(strcmpi(tbl.StimulusID, condition2), :).Responses);

    disp(condition1 + " vs. " + condition2 + ": p: " + pVal_cStim);
end

for comparrision_id = 1:size(test_comparrisions, 1)
    condition1 = test_comparrisions(comparrision_id, 1);
    condition2 = test_comparrisions(comparrision_id, 2);
    disp(condition1 + " vs. " + condition2 + ": p: " + all_pVals(comparrision_id));
end


if 1  % I would like to only manually overwrite the results files
    % Here I save the stats for Figure 5d
    fID = fopen("LinearMixedEffectsModel_stats_by_recordings.txt", "w+");
    try
        for comparrision_id = 1:size(test_comparrisions, 1)
            condition1 = test_comparrisions(comparrision_id, 1);
            condition2 = test_comparrisions(comparrision_id, 2);
            fwrite(fID, condition1 + " vs. " + condition2 + ": p: " + all_pVals(comparrision_id) + "; T-Stat: " + all_tStat(comparrision_id) + " (n=" + all_trials(comparrision_id) + ")" + newline);
        end

        fwrite(fID, newline + "Individual Response integrals: mean +- s.e.m." + newline);
        stim_keys = ["Narrow", "SF", "Ori", "Mixed"];
        for i = 1:4
            % fwrite(fID, stim_keys(i) + ": " + (100 .* nanmean(response_integrals(i, :), 2)) + " % +- " + (100 .* sem(response_integrals(i, :), 2)) + " % \DeltaF/F" + newline);
            fprintf(fID, "%s: %0.2f %% +- %0.2f %% %s %s", stim_keys(i), (100 .* nanmean(response_integrals(i, :), 2)), (100 .* sem(response_integrals(i, :), 2)), "\DeltaF/F", newline);
        end

        % Simon decises we need to change everything over to median +- ci
        fwrite(fID, newline + "Individual Response integrals: median +- c.i." + newline);
        stim_keys = ["Narrow", "SF", "Ori", "Mixed"];
        for i = 1:4
            % fwrite(fID, stim_keys(i) + ": " + (100 .* nanmean(response_integrals(i, :), 2)) + " % +- " + (100 .* sem(response_integrals(i, :), 2)) + " % \DeltaF/F" + newline);
            fprintf(fID, "%s: %0.2f %% +- %0.2f %% %s %s", stim_keys(i), (100 .* nanmedian(response_integrals(i, :), 2)), (100 .* compute95CI(response_integrals(i, :), 2)), "\DeltaF/F", newline);
        end
    catch; end
    fclose(fID);
end




