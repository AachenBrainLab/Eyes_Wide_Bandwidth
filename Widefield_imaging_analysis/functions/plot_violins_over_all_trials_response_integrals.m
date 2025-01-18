function plot_violins_over_all_trials_response_integrals(all_traces, n_conditions, condition_keys)

targetAreas = ["VISp", "VISpm", "VISam", "VISa", "VISrl", "VISal", "VISl"];

if ~exist("n_conditions", "var")
    n_conditions = 4;
end

% cl = gca().ColorOrder;
% cl = cat(1, [0.05 0.05 0.05], cl);

cl = [214 180 148; 128 110 166; 38 101 168; 146 90 68] ./ 255;

clearvars hndl ax
close all; figure(1); fig_size()

try close(2); catch; end; figure(2);
if n_conditions > 3
    fig_size(gcf, [0.1 0.05 0.6 0.9]);
else
    fig_size(gcf, [0.1 0.05 0.45 0.9]);
end

% y_lim = [-1.2 4.5];
y_lim = [-5 13];
for area_id = 1:7
%     response_integrals = nan(4, 8);
    for condition_id = 1:n_conditions
        temp_data = all_traces(:, :, area_id, condition_id, :);
        temp_data = reshape(temp_data, size(temp_data, 1), []);
        
%         stim_ids = (1:8) + 8.*(condition_id-1);
        if condition_id == 1
            response_integrals = 100*nanmean(temp_data(16:45, :), 1);  % 100* for percent
        else
            response_integrals(condition_id, :) = 100*nanmean(temp_data(16:45, :), 1);  % 100* for percent
        end
    end

    ax{condition_id} = subplot(2, 4, area_id); hold on;

    v = violinplot(response_integrals', vec(repmat(condition_keys, [8, 1])));
    for i = 1:n_conditions
        v(i).ViolinColor = cl(i, :);
        v(i).ScatterPlot.SizeData = 8;
    end

%     ylim(y_lim);
    temp_max = max(response_integrals(:));

    % Insert for new stats:
    LME_trials_averaged_over_recordings

    levels = 0.065;
%     f = 1+1*levels; p = ranksum(response_integrals(1, :), response_integrals(2, :)); disp(p); plot([1.1, 1.9], temp_max .* f .* [1 1], "-k"); text(1.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
%     f = 1+1*levels; p = ranksum(response_integrals(2, :), response_integrals(3, :)); disp(p); plot([2.1, 2.9], temp_max .* f .* [1 1], "-k"); text(2.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
%     f = 1+3*levels; p = ranksum(response_integrals(1, :), response_integrals(3, :)); disp(p); plot([1.1, 2.9], temp_max .* f .* [1 1], "-k"); text(2.0, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
%     if n_conditions > 3
%         f = 1+5*levels; p = ranksum(response_integrals(1, :), response_integrals(4, :)); disp(p); plot([1.1, 3.9], temp_max .* f .* [1 1], "-k"); text(2.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
%         f = 1+1*levels; p = ranksum(response_integrals(3, :), response_integrals(4, :)); disp(p); plot([3.1, 3.9], temp_max .* f .* [1 1], "-k"); text(3.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");

    f = 1+1*levels; p = all_pVals(1); disp(p); plot([1.1, 1.9], temp_max .* f .* [1 1], "-k"); text(1.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
    f = 1+1*levels; p = all_pVals(4); disp(p); plot([2.1, 2.9], temp_max .* f .* [1 1], "-k"); text(2.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
    f = 1+3*levels; p = all_pVals(2); disp(p); plot([1.1, 2.9], temp_max .* f .* [1 1], "-k"); text(2.0, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
    if n_conditions > 3
        f = 1+5*levels; p = all_pVals(3); disp(p); plot([1.1, 3.9], temp_max .* f .* [1 1], "-k"); text(2.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
        f = 1+1*levels; p = all_pVals(5); disp(p); plot([3.1, 3.9], temp_max .* f .* [1 1], "-k"); text(3.5, f*temp_max, get_significance_str(p), "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
        xticks(1:4); xticklabels(condition_keys(1:4)); xtickangle(45);
        xlim([0.5 4.5]);
    else
        xticks(1:3); xticklabels(condition_keys(1:3)); xtickangle(45);
        xlim([0.5 3.5]);
    end
    
    if (area_id == 1) || (area_id == 5)
        ylabel("Response Integrals" + newline + "\DeltaF/F [%]");
    end
    
    title(targetAreas(area_id));
    ylim(ylim() + [0, 0.6]);

    % INSERT FOR TESTING
%     test_script_for_LME
    % 

    disp(targetAreas(area_id) + ":");
    for i = 1:n_conditions
        fprintf("%0.2f +- %0.2f\n", nanmean(response_integrals(i, :), 2), sem(response_integrals(i, :), 2));
    end
    disp(newline)
end
