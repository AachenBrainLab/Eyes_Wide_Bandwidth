function AUCs = compute_auc_maps(U, nVc, opts, conditions2plot)

if ~exist("conditions2plot", "var"); conditions2plot = 1:4; end

baselineVc = squeeze(nanmean(nVc(:, 1:opts.preStim, :), 2));
mean_baseline = svdFrameReconstruct(U, nanmean(baselineVc, 2));
clearvars baselineVc;

n_conditions = length(conditions2plot);
n_frames = round(opts.SessionData.stimulusDur * opts.sRate) + opts.preStim;
AUCs = nan(540, 586, n_conditions-1);
for stim_id = 1:n_conditions  % all_stim_ids
    use_ids = ceil((opts.SessionData.stim_id)/8) == conditions2plot(stim_id);
    tempVc = squeeze(nanmean(nVc(:, opts.preStim + (1:n_frames), use_ids), 2));
    individual_trial_movies = svdFrameReconstruct_GN(U, squeeze(tempVc)) - mean_baseline;
    
    if stim_id == 1
        reference_responses = individual_trial_movies;
    else
        labels = cat(2, ones(1, size(reference_responses, 3)), zeros(1, size(individual_trial_movies, 3)));
        AUCs(:, :, stim_id-1) = reshape(colAUC(cat(2, reshape(reference_responses, 540*586, []), reshape(individual_trial_movies, 540*586, []))', labels', "ROC", "plot", 0, "abs", 0), 540, 586);
    end
end


