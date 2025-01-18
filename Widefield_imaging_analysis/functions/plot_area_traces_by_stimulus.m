function all_traces = plot_area_traces_by_stimulus(U, Vc, opts, SessionData, stim_ids2plot, condition_keys)
%%
y_lim = [-1.5, 6];


%% area masks
load("areaMasks_with_ALM_MM_ROIs.mat", ...
     "labelsSplit", "sidesSplit", "areaMasks");

targetAreas = ["VISp", "VISpm", "VISam", "VISa", "VISrl", "VISal", "VISl"];
area_mask_ids = zeros(1, length(targetAreas), "int32");
Cnt = 0;
for areaStr = targetAreas
    Cnt = Cnt + 1;
    area_mask_ids(Cnt) = find(strcmpi(labelsSplit, areaStr) & strcmpi(sidesSplit, "L"));
end


%% get area traces
try close(1); catch; end
figure(1); fig_size();
baselineVc = squeeze(nanmean(Vc(:, 1:opts.preStim, :), 2));

baseline = svdFrameReconstruct(U, baselineVc);
mean_baseline = nanmean(baseline, 3);


% stim_ids2plot = sort(unique(SessionData.stim_id));
clearvars ax;
cl = get(gca, 'colororder');
if length(stim_ids2plot) > 8
    cl = cat(1, [0.05 0.05 0.05; 0.75 0.2 0.95], cl);
else
    cl = cat(1, [0 0 0], cl);
end

stimulusDur = opts.SessionData.stimulusDur;

all_traces = nan(opts.preStim + opts.postStim + 1, 8*10, length(targetAreas), length(stim_ids2plot));
for stim_id = 1:length(stim_ids2plot)
    % This was the old version (I think only V1):
    % cTrialIds = SessionData.stim_id == stim_ids2plot(stim_id);

    % This if for V4 now the final version, hopefully!
    cTrialIds  = ceil((SessionData.stim_id)/8) == stim_ids2plot(stim_id);

    img = svdFrameReconstruct_GN(U, Vc(:, :, cTrialIds));
    img = img - mean_baseline;

    for area_id = 1:length(targetAreas)
        ax{stim_id} = subplot(2, 4, area_id); hold on;
        
        if stim_id == 1
            % patch([0 0 1 1], y_lim([1 2 2 1]), "k", "FaceAlpha", 0.10, "EdgeAlpha", 0); hold on;
            patch(stimulusDur .* [0 0 1 1], y_lim([1 2 2 1]), "k", "FaceAlpha", 0.10, "EdgeAlpha", 0); hold on;
            plot([-2 3], [0 0], "k--"); 
        end
        
        temp_mask = areaMasks(:, :, area_mask_ids(area_id));
        temp_mask(temp_mask < 1) = NaN;
        
        traces = squeeze(nanmean(nanmean(temp_mask .* img, 1), 2));
        all_traces(:, 1:size(traces, 2), area_id, stim_id) = traces;
        
        x = (-opts.preStimDur * opts.sRate : opts.postStimDur * opts.sRate) ./ opts.sRate;
        if ~(length(x) == size(traces, 1))
            x = (-opts.preStim : opts.postStim) ./ opts.sRate;
        end
        h = plot_with_error_shading(x, 100*nanmean(traces, 2), 100*sem(traces, 2), 0.15, cl(stim_id, :));
        
        if area_id == 1
            hndl(stim_id) = h;
        end
        
        if stim_id == length(stim_ids2plot)  % last stim
            title(targetAreas(area_id));
            xlabel("Time [s]")
            ylabel("\DeltaF/F [%]")
        end
        
        ylim(y_lim);
        xlim([-opts.preStimDur, opts.postStimDur]);
    end
    drawnow;
end

try
    legend(hndl, condition_keys, "Location", "NW");
catch
    legend(hndl, string(num2str(stim_ids2plot'))', "Location", "NW");
end

% % colormap(inferno(6*256));
% cb = colorbar_without_position_change();
% cb.Ticks = c_lim .* [-1 0 1];
% match_axis(ax, 1);

