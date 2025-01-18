function area_masks_50 = find_most_responsive_50percent_pixel_masks(movies, opts, targetAreas)
%% area masks
% This loads a version of the allenCCF map that contains binary maps for
% each area.
load("areaMasks.mat", "labelsSplit", "sidesSplit", "areaMasks");


if ~exist("targetAreas", "var")  % defaults
    targetAreas = ["VISp", "VISpm", "VISam", "VISa", "VISrl", "VISal", "VISl"];
end
area_mask_ids = zeros(1, length(targetAreas), "int32");
area_masks_50 = zeros(540, 586, length(targetAreas), "single");
Cnt = 0;
for areaStr = targetAreas
    Cnt = Cnt + 1;
    area_mask_ids(Cnt) = find(strcmpi(labelsSplit, areaStr) & strcmpi(sidesSplit, "L"));
    area_masks_50(:, :, Cnt) = single(areaMasks(:, :, area_mask_ids(Cnt)));
end
area_masks_50(area_masks_50 < 1) = NaN;


% not average over all stimuli to then find the most responsive pixels (as integrals)
% I compute the average stimulus response
stimulusFrames = round(opts.SessionData.stimulusDur .* opts.sRate);
stim_frames = opts.preStim + (1 : stimulusFrames);
responses = nanmean(nanmean(movies(:, :, stim_frames, :), 3), 4);

for area_id = 1:length(targetAreas)
    % responses
    
    masked_responses = responses .* area_masks_50(:, :, area_id);
%     masked_responses(isnan(masked_responses)) = [];
    
    % to get the 50% pixels I can conveniently use the median to define
    % this threshold
    threshold = nanmedian(vec(masked_responses));
    
    area_masks_50(:, :, area_id) = masked_responses > threshold;
end
area_masks_50(area_masks_50 < 1) = NaN;





