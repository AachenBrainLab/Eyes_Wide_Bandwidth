function niceFigures(ax)
% This function adjusts various axis and child properties for better
% visualization, excluding intermediate grid lines and ticks.

% if no axis is specified, use the current axis
if nargin == 0 || isempty(ax)
    ax = gca;
end

% General axis properties
try
    ax.Parent.Renderer = 'painters'; % High-quality rendering for vector graphics
end
ax.LineWidth = 2;          % set axis line width
ax.FontSize = 14;          % set axis font size
ax.TickDir = 'out';        % set ticks to point outward
ax.Box = 'off';            % remove top and right box lines
ax.TickLength = [0.01 0.01]; % set tick length for better scaling
ax.XMinorTick = 'off';     % disable minor ticks for X-axis
ax.YMinorTick = 'off';     % disable minor ticks for Y-axis

% remove grid lines
ax.XGrid = 'off';
ax.YGrid = 'off';

% adjust axis labels and titles
if isprop(ax, 'Title')
    ax.Title.FontSize = 16;
    ax.Title.FontWeight = 'bold';
end
if isprop(ax, 'XLabel')
    ax.XLabel.FontSize = 14;
end
if isprop(ax, 'YLabel')
    ax.YLabel.FontSize = 14;
end

% children properties
children = ax.Children;
for idx = 1:numel(children)
    child = children(idx);
    switch class(child)
        case 'matlab.graphics.chart.primitive.Line'
            child.LineWidth = 2;
        case 'matlab.graphics.chart.primitive.ErrorBar'
            child.LineWidth = 2;
            child.MarkerSize = 6;
            child.MarkerFaceColor = 'w';
        case 'matlab.graphics.chart.primitive.Scatter'
            child.SizeData = max(child.SizeData, 36); % minimum marker size
        case 'matlab.graphics.chart.primitive.Bar'
            child.LineWidth = 1.5;
        case 'matlab.graphics.chart.primitive.Histogram'
            child.EdgeColor = 'none'; % remove edges for cleaner look
        case 'matlab.graphics.chart.primitive.BoxChart'
            child.LineWidth = 2;
        otherwise
            try
                % applying general settings to other objects
                if isprop(child, 'LineWidth')
                    child.LineWidth = 2;
                end
                if isprop(child, 'MarkerSize')
                    child.MarkerSize = max(child.MarkerSize, 6);
                end
            catch
                
            end
    end
end

% ensure hidden lines and other elements are updated
hiddenLines = findall(ax, 'Type', 'line');
for idx = 1:numel(hiddenLines)
    hiddenLines(idx).LineWidth = 2;
end

% update legends if present
legendObj = findobj(ax.Parent, 'Type', 'legend');
for idx = 1:numel(legendObj)
    legendObj(idx).FontSize = 12;
    legendObj(idx).Box = 'off'; % remove legend box for cleaner appearance
end

end
