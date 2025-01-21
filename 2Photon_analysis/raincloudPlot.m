function raincloudPlot(data,color,bw)
    % raincloudPlot(data, [color_rgb], ks_density_bandwidth)
    % @E.B.2023
    % plots half-violin plots with customizable colors and optional raincloud plots
    %
    % Inputs:
    %   data        - a m x n matrix where each column is a violin
    %   color       - rgb of preferred color, where each row is the color that belongs to given column/violin
    %                 can also be a single color for all violins e.g.  [ 0.4392    0.6902    0.8392]
    %   bw          - bandwidth of the ks density, default is 0.5
    
    % Examples:
    %   raincloudPlot(data, [0.4392    0.6902    0.8392], 0.5)
    
if nargin < 3 || isempty(bw)
    bw = 0.5;
end

%   plot (raincloud plot) for each column
figure; niceFigures;
hold on;

% parameters
numCols = size(data, 2); % number of columns
jitterRange =-0.25; % range for scatter jitter
colors = zeros(numCols,3); %in case all columns need to have the same color
if size(color,1) == 1
    colors = repmat(color,numCols,1);
end
% looping through each column to create raincloud plots
for i = 1:numCols
    % extracting the column data
    colData = data(:, i);

    %kernel density estimate
    [f(i,:), xi(i,:)] = ksdensity(colData, linspace(prctile(colData,1),prctile(colData,99),100),"Bandwidth",bw );

end

%rescaling f so that all violins max their belly at the same point
new_range = 0.5;
f_min = min(min(f));
f_max = max(max(f));
scaled_f = (f / (f_max - f_min)) * new_range;

for i = 1:numCols
    % extracting the column data
    colData = data(:, i);

    % density kernel plot
    fill( [scaled_f(i,:), zeros(size(scaled_f(i,:)))] + i, [xi(i,:), fliplr(xi(i,:))], colors(i, :), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9);

    % add scatter points (with jitter)
    jitter = (rand(size(colData))) * jitterRange;
    scatter(jitter + i, colData, 30, colors(i, :),'LineWidth',1);

    % add boxplot
    quartiles = prctile(colData, [25 50 75]);
    boxWidth = 0.10; % Width of the box
    fill([i - boxWidth, i + boxWidth, i + boxWidth, i - boxWidth], ...
        [quartiles(1), quartiles(1), quartiles(3), quartiles(3)], ...
        colors(i, :), 'FaceAlpha', 0.8, 'EdgeColor', 'k');
    line([i - boxWidth, i + boxWidth], [quartiles(2), quartiles(2)], ...
        'Color', 'k', 'LineWidth', 1);

    % whiskers
    whiskerLow = min(colData);
    whiskerHigh = max(colData);
    line([i, i], [whiskerLow, quartiles(1)], 'Color', 'k', 'LineWidth', 0.5); % Lower whisker
    line([i, i], [quartiles(3), whiskerHigh], 'Color', 'k', 'LineWidth', 0.5); % Upper whisker
end

% plot adaptation
xlim([0, numCols + 1]);
ylim([min(data(:)) - 1, max(data(:)) + 1]);
xticks(1:numCols);

% xlabel('Group');
% ylabel('Data Value');
% xticklabels({'Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group5'});
% title('Raincloud Plot with Vertical Violins and Boxplots');
% set(gca, 'XGrid', 'off', 'YGrid', 'off');
% hold off;