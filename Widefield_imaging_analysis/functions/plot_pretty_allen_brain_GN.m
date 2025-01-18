function plot_pretty_allen_brain_GN(img, cRange, cMap, contourWidth, areaLineWidth, extended_allenMask, limit_plotted_region)
% limit_plotted_region: I added this because the aligned images have a
% large empty margine around them. I can pack plots much closer if I remove
% this. 

if ~exist("contourWidth", "var")
    % contourWidth = 2.0;
    contourWidth = 0.5;
end

if ~exist("cMap", "var")
    cMap = colormap_blueblackred(6*256);
end

if ~exist("areaLineWidth", "var")
%     areaLineWidth = 0.1;
    areaLineWidth = 0.3;
end

if ~exist("extended_allenMask", "var"); extended_allenMask = 0; end

if ~exist("limit_plotted_region", "var"); limit_plotted_region = 1; end

addpath("C:\Users\nabbefeld\Documents\GitHub\rateDisc");
addpath("C:\Users\nabbefeld\Documents\GitHub\generalTools");

load('allenDorsalMapSM.mat', 'dorsalMaps'); allenMask = dorsalMaps.allenMask;

if extended_allenMask  % this is an optional overwrite to enlarge the allenMask
    allenMask = imerode(allenMask, strel("disk", 40), "same");
end

% figure(1);
if all(size(img) == size(allenMask))
    try
        img = arrayCrop(img, allenMask);
    catch; end
end

if exist("cRange", "var")
    if ~isempty(cRange)
        if length(cRange) == 2
            imageScale(img, cRange(2), cRange(1) ~= 0);
            % imageScale has weird behavior and doesnt let me set both lims
            caxis(cRange);
        else
            imageScale(img, cRange);
        end
    else
        imageScale(img);
    end
else
    imageScale(img);
end
colormap(gca, cMap);
addpath("C:\Users\nabbefeld\Documents\GitHub\rateDisc");  % for: rateDisc_plotAllenOutline.m
rateDisc_plotAllenOutline(gca, "Both", 0, areaLineWidth);  % area outlines
hold on;
% contour(allenMask, 'k', 'LineWidth', contourWidth);  % black outline


%% New part where I try to plot a prettier contour
[M, c] = contour(allenMask, "k", "Levels", 1, "LineWidth", contourWidth);
try
    contourTable = getContourLineCoordinates(M);
    delete(c); clearvars M;
    
    % It took me forever to figure out why the outline always came out shitty.
    % turns out Matlab draws line properly (smooth corners) only aslong as you
    % use < 2001 values to plot. The Allen Outline usually hast 2059 which
    % would result in artefacts at every line segment. So ...
    % I only draw every other vertice and that fixes it. Stupid, but works.
    x = contourTable.X;
    y = contourTable.Y;
    plot(x([1:2:end, 3]), y([1:2:end, 3]), "k", "LineWidth", contourWidth); axis image;
catch ME
    disp(ME)
end


% If desired remove the empty/NaN margines aournd the brain
if limit_plotted_region
    lim_padding = 1;
    x1 = find(any(allenMask == 1), 1, 'first') - lim_padding;  x2 = find(any(allenMask == 1), 1, 'last') + lim_padding;
    y1 = find(any(allenMask' == 1), 1, 'first') - lim_padding; y2 = find(any(allenMask' == 1), 1, 'last') + lim_padding;
    xlim(gca, [x1, x2]); ylim(gca, [y1, y2]);
end




