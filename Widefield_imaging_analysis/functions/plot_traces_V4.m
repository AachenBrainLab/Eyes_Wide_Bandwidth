function traces = plot_traces_V4(U, nVc, opts, SessionData, responseDur, conditions2plot, condition_keys, ROI_loc)
if ~exist("conditions2plot", "var"); conditions2plot = 1:4; end
if ~exist("condition_keys", "var"); condition_keys = string(conditions2plot); end

ROI_radius = 10;

% prepare the ROI mask
shape = size(U, 1:2);
mask = nan(shape(1), shape(2), "single");
mask(sqrt(power(((1:shape(1))-ROI_loc(1)), 2) + ...
          power(((1:shape(2))-ROI_loc(2)), 2)') < ROI_radius) = 1;
U = nanmean(nanmean(U .* mask, 1), 2);


baselineVc = squeeze(nanmean(nVc(:, 1:opts.preStim, :), 2));
% baselineVc = squeeze(nanmean(nVc(:, (-round(baselineDur * opts.sRate):0)+opts.preStim, :), 2));
clearvars ax;
% c_lim = 0.05;

baseline = svdFrameReconstruct(U, baselineVc);
mean_baseline = nanmean(baseline, 3);
% std_baseline = nanstd(baseline, [], 3);

% all_stim_ids = sort(unique(SessionData.stim_id));
% n_conditions = length(conditions2plot);
n_conditions = length(unique(SessionData.stim_id));

n_frames = round(responseDur * opts.sRate) + opts.preStim;
% I leave it here, cause I could compute both here at the same time If I 
% have the ROI already!
% movies = nan(size(U, 1), size(U, 2), n_frames, n_conditions);
traces = nan(n_frames, n_conditions);
% traces = nan(n_frames, n_conditions);
for stim_id = 1:n_conditions  % 1:n_conditions % all_stim_ids
%     % use_ids = ceil((SessionData.stim_id)/8) == conditions2plot(stim_id);
%     use_ids = SessionData.stim_id == stim_id;
%     tempVc = nVc(:, 1:n_frames, use_ids);
%     img = svdFrameReconstruct_GN(U, squeeze(nanmean(tempVc, 3)));
%     % movies(:, :, :, stim_id) = img - mean_baseline;
%     temp = img - mean_baseline;
%     traces(:, stim_id) = squeeze(nanmean(nanmean(temp .* mask, 1), 2));
    
    % use_ids = ceil((SessionData.stim_id)/8) == conditions2plot(stim_id);
    use_ids = SessionData.stim_id == stim_id;
    tempVc = nVc(:, 1:n_frames, use_ids);
    img = svdFrameReconstruct_GN(U, squeeze(nanmean(tempVc, 3)));
    % movies(:, :, :, stim_id) = img - mean_baseline;
    % traces2(:, stim_id) = squeeze(nanmean(nanmean(img - mean_baseline, 1), 2));
    traces(:, stim_id) = squeeze(img) - mean_baseline;
end


try close(1); catch; end
figure(1);
% subplot(1, );
% temp_ids = 1:n_conditions;
x = ((1-opts.preStim) + (1:size(traces, 1))) ./ opts.sRate;
cl = gca().ColorOrder;
clearvars h;
for stim_id = 1:4
    use_ids = ceil((1:n_conditions)/8) == stim_id;  % conditions2plot(stim_id);
    
    temp_data = traces(:, use_ids);
    h(stim_id) = plot_with_error_shading_GN(x, 100*nanmean(temp_data, 2), 100*sem(temp_data, 2), 0.2, cl(stim_id, :), cl(stim_id, :), 1); hold on;
%     plot(x, temp_data, "color", cl(stim_id, :)); hold on;
end
legend(h, condition_keys(1:4));
xlabel("Time [s]");
ylabel("\DeltaF/F [%]")

temp_ylim = ylim();
% patch_hndl = patch([0 0 2 2], temp_ylim([1 2 2 1]), "k", "FaceAlpha", 0.2, "EdgeAlpha", 0, "HandleVisibility", "off");
patch_hndl = patch([0 0 2 2], temp_ylim([1 2 2 1]), "k", "FaceAlpha", 0.2, "EdgeAlpha", 0);

h = get(gca, "Children");
set(gca, "Children", [h(2:end); h(1)]);
patch_hndl.HandleVisibility = "off";
ylim(temp_ylim);

%%
%     shape = size(movies, 1:2);
%     mask = zeros(shape(1), shape(2), "single");
%     mask(sqrt(power(((1:shape(1))-ROI_loc(1)), 2) + ...
%               power(((1:shape(2))-ROI_loc(2)), 2)') < 10) = 1;
% 
%     figure(1)
%     subplot(1, 3, 1);
%     imagesc(U(:, :, 1)); colormap(gray); axis image;
%     subplot(1, 3, 2);
%     imagesc(mask); colormap(gray); axis image;
%     subplot(1, 3, 3);
%     imagesc(U(:, :, 1) + 0.05*mask); colormap(gray); axis image;

    