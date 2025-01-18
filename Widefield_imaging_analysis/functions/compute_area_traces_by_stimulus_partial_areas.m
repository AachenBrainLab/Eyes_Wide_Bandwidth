function all_traces = compute_area_traces_by_stimulus_partial_areas(U, Vc, opts, SessionData, stim_ids2plot, condition_keys, area_masks_50)

%% get area traces
baselineVc = squeeze(nanmean(Vc(:, 1:opts.preStim, :), 2));

baseline = svdFrameReconstruct(U, baselineVc);
mean_baseline = nanmean(baseline, 3);

% stimulusFrames = round(opts.sRate .* opts.SessionData.stimulusDur);

all_traces = nan(opts.preStim + opts.postStim + 1, 8*10, size(area_masks_50, 3), length(stim_ids2plot));
for stim_id = 1:length(stim_ids2plot)
    % This if for V4 now the final version, hopefully!
    cTrialIds  = ceil((SessionData.stim_id)/8) == stim_ids2plot(stim_id);

    img = svdFrameReconstruct_GN(U, Vc(:, :, cTrialIds));
    img = img - mean_baseline;

    for area_id = 1:size(area_masks_50, 3)  % length(targetAreas)
        temp_mask = area_masks_50(:, :, area_id);
        
        % before I just averaged over all pixel right away. Instead I need to compute and compare all the individual responses
        traces = squeeze(nanmean(nanmean(temp_mask .* img, 1), 2));
        all_traces(:, 1:size(traces, 2), area_id, stim_id) = traces;
    end
end




