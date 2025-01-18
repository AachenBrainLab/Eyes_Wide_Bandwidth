function movies = generate_movies_V4(U, nVc, opts, SessionData, responseDur, conditions2plot, condition_keys)

if ~exist("conditions2plot", "var"); conditions2plot = 1:4; end
if ~exist("condition_keys", "var"); condition_keys = string(conditions2plot); end

baselineVc = squeeze(nanmean(nVc(:, 1:opts.preStim, :), 2));
% baselineVc = squeeze(nanmean(nVc(:, (-round(baselineDur * opts.sRate):0)+opts.preStim, :), 2));
clearvars ax;
c_lim = 0.05;

baseline = svdFrameReconstruct(U, baselineVc);
mean_baseline = nanmean(baseline, 3);
% std_baseline = nanstd(baseline, [], 3);

% all_stim_ids = sort(unique(SessionData.stim_id));
n_conditions = length(conditions2plot);

n_frames = round(responseDur * opts.sRate) + opts.preStim;
movies = nan(size(U, 1), size(U, 2), n_frames, n_conditions);
for stim_id = 1:n_conditions % all_stim_ids
    use_ids = ceil((SessionData.stim_id)/8) == conditions2plot(stim_id);
    tempVc = nVc(:, 1:n_frames, use_ids);
    img = svdFrameReconstruct_GN(U, squeeze(nanmean(tempVc, 3)));
    movies(:, :, :, stim_id) = img - mean_baseline;
end


