function cb = colorbar_without_position_change(ax)
if ~exist("ax", "var")
    ax = gca;
end

drawnow();  % Not sure if this is required, just to make sure
originalPosition = get(ax, 'Position');
cb = colorbar;

% Reset axes to original size.
set(ax, 'Position', originalPosition);
