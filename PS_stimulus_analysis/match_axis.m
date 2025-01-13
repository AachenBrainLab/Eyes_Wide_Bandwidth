function ax_lim = match_axis(ax, img)

ax = ax(:);
ax(cellfun(@(x) isempty(x), ax)) = [];
ax_lim = [inf, -inf];
for i = 1:size(ax, 1)
    if img
        temp_lim = caxis(ax{i});
    else
        temp_lim = get(ax{i}, 'YLim');
    end
    
    if ax_lim(1) > temp_lim(1); ax_lim(1) = temp_lim(1); end
    if ax_lim(2) < temp_lim(2); ax_lim(2) = temp_lim(2); end
end

for i = 1:size(ax, 1)
    if img
        caxis(ax{i}, ax_lim);
    else
        set(ax{i}, 'YLim', ax_lim);
    end
end

