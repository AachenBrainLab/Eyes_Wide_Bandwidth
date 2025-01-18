function plot_average_stim_response_V4(U, nVc, opts, SessionData, baselineDur, responseDur, c_lim, conditions2plot, condition_keys)
% Defaults:
% if ~exist("baselineDur", "var"); baselineDur = 1.0; end
if ~exist("responseDur", "var"); responseDur = 2.0; end
if ~exist("c_lim", "var");       c_lim = 0.03; end
% if ~exist("conditions2plot", "var"); conditions2plot = length(SessionData.stimFile) / 8; end
if ~exist("conditions2plot", "var"); conditions2plot = 1:4; end
if ~exist("condition_keys", "var"); condition_keys = string(conditions2plot); end


% 
try close(1); catch; end
figure(1); fig_size();

% n_baseline_frames = round(opts.preStim * opts.sRate);

baselineVc = squeeze(nanmean(nVc(:, 1:opts.preStim, :), 2));
clearvars ax;

baseline = svdFrameReconstruct(U, baselineVc);
mean_baseline = nanmean(baseline, 3);
% std_baseline = nanstd(baseline, [], 3);

% stim_keys = ["Narror BWs", "Broad SF BW", "Broad Ori BW", "Broad SF-BW and Ori-BW"];
% all_stim_ids = sort(unique(SessionData.stim_id));
% all_stim_ids = length(SessionData.stimFile) / 8;  % the unique condition without repetitions
n_response_frames = ceil(responseDur * opts.sRate);
n_conditions = length(conditions2plot);
n_rows = round(sqrt(n_conditions));
n_cols = ceil(n_conditions ./ n_rows);
for stim_id = 1:n_conditions
    ax{stim_id} = subplot(n_rows, n_cols, stim_id);
    
    use_ids = ceil((SessionData.stim_id)/8) == conditions2plot(stim_id);
    tempVc = nVc(:, opts.preStim + (1:n_response_frames), use_ids);
    
    if isfield(opts, "transParams")
        % plot an allen outline
    end
    img = svdFrameReconstruct_GN(U, tempVc);  % This one is for testing
    
    img = svdFrameReconstruct_GN(U, squeeze(nanmean(tempVc, 2)));
    mean_img = nanmean(img, 3);
    
    if 1  % plot mean
% %         % plot_pretty_allen_brain_GN(img, [-0, 0.015], inferno(6*256), 1.0);
% %         plot_pretty_allen_brain_GN(mean_img - mean_baseline, c_lim .* [-1, 1], fireice(6*256), 1.0);
%         imagesc(mean_img - mean_baseline, c_lim .* [-1, 1], fireice(6*256), 1.0);

        if isfield(opts, "transParams")
            plot_pretty_allen_brain_GN(mean_img - mean_baseline, c_lim .* [-1, 1], fireice(6*256), 1.0);
        else
            imagesc(mean_img - mean_baseline);
            axis image;
            caxis(c_lim .* [-1, 1]);
            colormap(fireice(6*256));
        end
        
%         colorbar_without_position_change();
        % title(stim_keys(stim_id));
        title(condition_keys(stim_id));
    else
        % plot d'
        c_lim = 3;
    
        std_img = nanstd(img, [], 3);
        dPrime = abs(mean_img - mean_baseline) ./ ((std_img + std_baseline) ./ 2);

        plot_pretty_allen_brain_GN(dPrime, c_lim .* [0, 1], inferno(6*256), 1.0);
%         title("Stimulus: " + stim_id);
    end
end

cb = colorbar_without_position_change();
cb.Ticks = c_lim .* [-1 0 1];
match_axis(ax, 1);

xticks([]); yticks([]);

% %% Manual overwrite:
% caxis([-1 1] .* 0.05)
% match_axis(ax, 1);


