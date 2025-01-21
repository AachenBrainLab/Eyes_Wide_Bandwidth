function cardinal_line(horizontal, vertical, cardinal, lineStyle, lineColor, lineWidth)
    % @E.B.2023
    % Draws a horizontal or a vertical line in the current plot.
    %
    % Inputs:
    %   horizontal  - this is to identify the direcetion of the cardinal line 
    %   horizontal  - this is to identify the direcetion of the cardinal line
    %   cardinal    - The x or y-coordinate of the cardinal line.
    %   lineStyle   - (Optional) Line style, e.g., '--', ':', '-'. Default is '-'.
    %   lineColor   - (Optional) Line color, e.g., 'r', 'b', [0, 0.5, 0]. Default is 'k' (black).
    %   lineWidth   - (Optional) Line width. Default is 1.
    %
    % Examples: 
    %   cardinal_line(0,1, 0.5); % Draws a black verical line at y = 0.5
    %   cardinal_line(1,0, 0.5,'--', 'r', 2); % Draws a red dashed horizontal line with width 2 at y = 0.5

    if nargin < 4 || isempty(lineStyle)
        lineStyle = '-';
    end
    if nargin < 5 || isempty(lineColor)
        lineColor = 'k';
    end
    if nargin < 6 || isempty(lineWidth)
        lineWidth = 1;
    end

    % Get the current axis limits
    xLimits = xlim;
    yLimits = ylim;

    % Draw the horizontal line
    if horizontal
        hold on; % Ensure the current plot is not cleared
        plot(xLimits, [cardinal, cardinal], 'LineStyle', lineStyle, 'Color', lineColor, 'LineWidth', lineWidth);
        hold off;
            % Draw the horizontal line
    elseif vertical
        hold on; % Ensure the current plot is not cleared
        plot([cardinal, cardinal], yLimits, 'LineStyle', lineStyle, 'Color', lineColor, 'LineWidth', lineWidth);
        hold off;
end
