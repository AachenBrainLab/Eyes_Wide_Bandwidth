%% stim explanation
% Folder 01 - 08 narrow orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 09 - 16 narrow orienation band and high sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 17 - 24 broad orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 25 - 32 broad orienation band and broad sf band ( central sf 0.04 cpd and 0° central ori)
%
% Center_mask versions (15°)
% Folder 33 - 40 narrow orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 41 - 48 narrow orienation band and high sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 49 - 56 broad orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 57 - 64 broad orienation band and broad sf band ( central sf 0.04 cpd and 0° central ori)
%
% Center aperture versions (15°)
% Folder 65 - 72 narrow orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 73 - 79 narrow orienation band and high sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 80 - 88 broad orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 89 - 96 broad orienation band and broad sf band ( central sf 0.04 cpd and 0° central ori)
%
% Folder 97 - 104 high central SF (0.16 cpd) narrow orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)
% Folder 105 - 112 high central SF (0.16 cpd) broad orienation band and narrow sf band ( central sf 0.04 cpd and 0° central ori)

%% basic variables
opts.Location = {'V1'};
opts.baseDur = 0.5;
opts.postStim = 2.5;
opts.brainRange = 3000; %max. depth in micrometers
opts.verbose = false; %flag to supress some of the text outputs
opts.brainThresh = 1;
opts.groupColors = {[0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]};
opts.optoGenetics = true; %flag to isolate sensory responses during optogenetics
opts.trialTriggerChan = 2;
opts.stimGroups = {1:8, 9:16, 17:24, 25:32; 33:40, 41:48, 49:56, 57:64; 65:72, 73:80,  81:88, 89:96}';
opts.stimGroups = opts.stimGroups(1 : 4, 1); %only use first 4 stimgroups for now
opts.stimNames = {'narrow ORI, narrow SF', 'narrow ORI, broad SF', 'broad ORI, narrow SF', 'broad ORI, broad SF'};
opts.selectivityRange = 1; %time after stimulus onset to compute selectivity index
opts.basePath = 'F:\\MotionClouds\';
opts.savePath = 'D:\\MotionClouds\';
opts.reload = true;

%% run over groups
recs = dir(opts.basePath);
recs = recs(~contains({recs.name}, '.'));

nrGroups = length(recs);
allLFP = cell(1, nrGroups);
allStimStd = cell(1, nrGroups);
animalIDs = nan(1, nrGroups);
scIdx = false(1, nrGroups);

for iRecs = 1 : size(recs, 1)
    try
        
        % get LFP data
        cPath = fullfile(recs(iRecs).folder, recs(iRecs).name);
        temp = textscan(recs(iRecs).name, '%f%s', 'delimiter', '_');
        animalIDs(iRecs) = temp{1};
        [allLFP{iRecs}, allStimStd{iRecs}, lfMeta, scIdx(iRecs)] = pC_lfpResponse_MotionMare(cPath, opts);
        
    catch ME
        disp(['Aborted: ' ME.message]);
    end
end

allMice = unique(animalIDs);
animalNrs = nan(1, length(animalIDs));
for iMice = 1 : length(allMice)
    animalNrs(ismember(animalIDs, allMice(iMice))) = iMice;
end

%% make brain surface figure - meanStd is the index to identify the brain surface
figure
% show alignment to check for inconsistencies
hold on;
cStd = cat(2,allStimStd{:});
plot(cStd)
title('Brain surface for each recording');
axis square
nhline(opts.brainThresh);

% check for response alignment
mergeStim = cat(4,allLFP{:});
[~,idx] = min(squeeze(mergeStim(:, 1450, 3, :)));

%% combined figures - LFP for each group
lfpRange = 25; %for visual response plot
lfpTraceRange = [-50 50]; %plot range for traces

stimLabel = opts.stimNames;
respLabels = {'LFP absolute'};
nrStims = length(opts.stimGroups);
layerDepths = [100 300 500 700 900];
cortexRange = [0, 1000];
nrLayers = length(layerDepths);
layerLabels = {'0-100', '100-300', '300-500', '500-700', '700-900'};
timeRange = [-.2 1]; %used time range in seconds

clear lines
h1 = figure('name', (sprintf('%s - %s', respLabels{1}, 'original response')));
h2 = figure('name', (sprintf('%s - %s', respLabels{1}, 'differential response')));
h3 = figure('name', respLabels{1});
h4 = figure('name', respLabels{1});
refTrace = cell(1, nrLayers);
selIdx = cell(nrStims, nrLayers);
meanDiff = cell(nrStims, nrLayers);
meanTrace = cell(nrStims, nrLayers);
for iStims = 1 : nrStims
    
    mergeStim = cat(4,allLFP{:});
    cRange = lfpRange;
    depthRange = (1 : 20 : size(mergeStim,1)*20)-1;
    traceRange = lfpTraceRange;
    
    
    % get current stimulus and difference to first stimulus
    diffStim = squeeze(mergeStim(:,:,iStims,:)) - squeeze(mergeStim(:,:,1,:));
    mergeStim = squeeze(mergeStim(:,:,iStims,:));
    
    % get traces from top, middle and bottom layers
    useLayers = nan(1, nrLayers);
    mergeTrace = cell(1, nrLayers);
    diffTrace = cell(1, nrLayers);
    for iLayers = 1 : nrLayers
        [~, useLayers(iLayers)] = min(abs(depthRange - layerDepths(iLayers)));
        if iLayers == 1
            cIdx = 1 : useLayers(iLayers);
        else
            cIdx = useLayers(iLayers-1) : useLayers(iLayers);
        end
        
        mergeTrace{iLayers} = -squeeze(nanmean(mergeStim(cIdx, :, :), 1)) + (max(traceRange)*2);
        if iStims == 1
            refTrace{iLayers} = mergeTrace{iLayers};
        end
        
        diffTrace{iLayers} = mergeTrace{iLayers} - refTrace{iLayers};
        sumTrace = mergeTrace{iLayers} + refTrace{iLayers} - (max(traceRange)*4); % to compute selectivity index
        
        % compute selectivity index for current versus first stimulus
        selRange = round((opts.baseDur + [0.1 opts.selectivityRange]) .*lfMeta.sRateHz);
        meanTrace{iStims, iLayers} = mean((mergeTrace{iLayers}(selRange(1):selRange(2),:)));
        meanDiff{iStims, iLayers} = mean((diffTrace{iLayers}(selRange(1):selRange(2),:)));
        selIdx{iStims, iLayers} = mean((diffTrace{iLayers}(selRange(1):selRange(2),:))) ./ mean(abs(sumTrace(selRange(1):selRange(2),:)));
    end
    
    % show average images
    for xx = 1 : 2
        depthIdx = depthRange > cortexRange(1) & depthRange < cortexRange(2);
        if xx == 1
            % show regular responses
            figure(h1);
            cStim = nanmean(mergeStim(depthIdx, :, :), 3);
            currRange = cRange;
        elseif xx == 2
            % show differential response
            figure(h2);
            cStim = nanmean(diffStim(depthIdx, :, :), 3);
            currRange = cRange/2;
        end
        
        subplot(1, nrStims, iStims)
        cStim = smoothCol(cStim, 2, 100);
        cImg = imagesc(cStim);
        ax = cImg.Parent;
        ax.XTick = 1 : lfMeta.sRateHz*0.5 : size(cStim,2);
        useTime = (((0 : lfMeta.sRateHz*0.5 : size(cStim,2)-1)./lfMeta.sRateHz) - opts.baseDur)*1000;
        ax.XTickLabel = useTime;
        
        xRange = (opts.baseDur + timeRange) .*lfMeta.sRateHz;
        xlim(xRange);
        colormap(ax, flipud(viridis(256)));
        caxis([-currRange currRange]);
        
        title(stimLabel(iStims))
        xlabel('time(ms)')
        nvline(opts.baseDur*lfMeta.sRateHz, 'w', 'linewidth', 2);
        
        depthSteps = round(200 / mean(diff(depthRange)));
        ax.YTick = 1 : depthSteps : size(cStim,1);
        ax.YTickLabels = depthRange(1 : depthSteps : length(depthRange));
        ylabel('Depth (um)');
        axis square;
        nhline(useLayers, 'w', 'linewidth', 2);
        drawnow;
    end
    
    % plot traces
    figure(h3);
    Cnt = 0;
    for xx = 1 : 2
        for iLayers = 1 : nrLayers
            
            if xx == 1
                % show regular responses
                cStim = mergeTrace{iLayers} - (max(traceRange)*2);
            elseif xx == 2
                % show differential response
                cStim = diffTrace{iLayers};
            end
            
            Cnt = Cnt + 1;
            subplot(2, nrLayers, Cnt); hold on;
            t = (0 : size(mergeTrace{1},1)-1)./ lfMeta.sRateHz - opts.baseDur;
            xlim(timeRange); axis square;
            %                 ylim(traceRange)
            nvline(0, '--k');
            nhline(0, '--k');
            %                 cStim = cStim - cStim(abs(t) == min(abs(t)),:);
            lines(iStims) = stdshade(cStim', 0.1, opts.groupColors{iStims}, t);
            if iLayers == 1
                cRange = [0 layerDepths(1)];
            else
                cRange = [layerDepths(iLayers-1) layerDepths(iLayers)];
            end
            title(sprintf('Depth range: %d-%d um', cRange(1), cRange(2)));
            ylabel('LFP deflection (uV)');
            xlabel('time (s)');
            drawnow;
        end
    end
end
legend(lines, opts.stimNames);

% plot differences as errorbar plot
figure(h4);
%     subplot(1,2,1);
errorbar(cellfun(@mean,meanDiff)', cellfun(@sem, meanDiff)', '-o', 'linewidth', 2)
xlim([0 6]); axis square;
title('Mean response difference');
view(90,90);
ax = gca;
ax.XTick = 1 : nrLayers;
ax.XTickLabel = layerLabels;
legend(opts.stimNames, 'location', 'southeast');
ax.YTick = 0:3:20;
grid on;
niceFigure

%% do statistical test for each condition
narrowData = nanmean(cat(1, meanTrace{1, :}),1);
pVal = nan(1, size(meanDiff,1)-1);
tStat = nan(1, size(meanDiff,1)-1);
disp('====================')
for iCond = 2 : size(meanDiff,1)
    
    fprintf(opts.stimNames{iCond})
    fprintf('\n');
    cData = nanmean(cat(1, meanTrace{iCond, 2:end}));
    diffData = nanmean(cat(1, meanDiff{iCond, 2:end}));
    %     diffData = cat(2, meanDiff{iCond, :});
    
    [pVal(iCond-1), tStat(iCond-1)] = LME_compare(narrowData, cData, animalNrs);
    
    fprintf('Mean response difference: %.2f +- %.2f\n',  mean(diffData), sem(diffData))
    fprintf('pVal for difference to narrow: %.10f\n',  pVal(iCond-1))
    fprintf('tStat for difference to narrow: %.10f\n',  tStat(iCond-1))
    fprintf('n = %.0f recordings\n',length(cData));
    disp('====================')
    
end

%% do statistical test for each condition and layer
animalNumbers = repmat(animalNrs, size(meanDiff,2), 1);
animalNumbers = animalNumbers(:);
narrowData = nanmean(cat(1, meanTrace{1, :}),1);
pVal = nan(1, size(meanDiff,1)-1);
tStat = nan(1, size(meanDiff,1)-1);
disp('====================')
for iCond = 2 : size(meanDiff,1)
    
    disp(opts.stimNames{iCond})
    for iLayer = 1 : size(meanDiff,2)
        cData = meanTrace{iCond, iLayer};
        diffData = meanDiff{iCond, iLayer};
        
        [pVal(iCond-1), tStat(iCond-1)] = LME_compare(narrowData, cData, animalNrs);
        
        disp('------------')
        fprintf('Layer %i\n', iLayer)
        fprintf('Mean response difference: %.2f +- %.2f\n',  mean(diffData), sem(diffData))
        fprintf('pVal for difference to narrow: %.10f\n',  pVal(iCond-1))
        fprintf('tStat for difference to narrow: %.10f\n',  tStat(iCond-1))
        
    end
    fprintf('n = %.0f recordings\n',length(cData));
    disp('====================')
end

%% same analysis for SC
scSurface = [80, 103, 80];
useRecs = find(scIdx);
lfpRange = 25; %for vision
stimLabel = opts.stimNames;
respLabels = {'LFP absolute'};
baseDur = round(opts.baseDur .*lfMeta.sRateHz);
nrStims = length(opts.stimGroups);
layerDepths = 150 : 150 : 900;
scRange = 900;
nrLayers = length(layerDepths);
timeRange = [-.2 1]; %used time range in seconds
selRange = [1500, 2000];
clear lines

h1 = figure('name', (sprintf('%s - %s', respLabels{1}, 'original response')));
h2 = figure('name', (sprintf('%s - %s', respLabels{1}, 'differential response')));
refTrace = cell(1, nrLayers);
selIdx = cell(nrStims, nrLayers);
meanDiff = cell(nrStims, nrLayers);
meanTrace = cell(nrStims, nrLayers);
for iStims = 1 : nrStims
    
    useRecs = useRecs(~isnan(scSurface));
    brainRange = (1 : 20 : size(allLFP{1},1)*20)-1;
    for iRecs = 1 : length(useRecs)
        
        scStart = brainRange(scSurface(iRecs)); %sc start for current recording
        depthIdx = brainRange >= scStart & brainRange < scStart + scRange;
        cStim = allLFP{useRecs(iRecs)}(depthIdx, :, :);
        
        if iRecs == 1
            mergeStim = cStim;
        else
            mergeStim = cat(4, mergeStim, cStim);
        end
    end
    cRange = lfpRange;
    depthRange = (1 : 20 : size(mergeStim,1)*20)-1;
    
    % get current stimulus and difference to first stimulus
    diffStim = squeeze(mergeStim(:,:,iStims,:)) - squeeze(mergeStim(:,:,1,:));
    mergeStim = squeeze(mergeStim(:,:,iStims,:));
    
    % get traces from top, middle and bottom layers
    useLayers = nan(1, nrLayers);
    mergeTrace = cell(1, nrLayers);
    diffTrace = cell(1, nrLayers);
    for iLayers = 1 : nrLayers
        
        
        [~, useLayers(iLayers)] = min(abs(depthRange - layerDepths(iLayers)));
        if iLayers == 1
            cIdx = 1 : useLayers(iLayers);
        else
            cIdx = useLayers(iLayers-1) : useLayers(iLayers);
        end
        
        mergeTrace{iLayers} = -squeeze(nanmean(mergeStim(cIdx, :, :), 1)) + (max(traceRange)*2);
        if iStims == 1
            refTrace{iLayers} = mergeTrace{iLayers};
        end
        
        diffTrace{iLayers} = mergeTrace{iLayers} - refTrace{iLayers};
        sumTrace = mergeTrace{iLayers} + refTrace{iLayers} - (max(traceRange)*4); % to compute selectivity index
        
        % compute selectivity index for current versus first stimulus
        selRange = round((opts.baseDur + [0.1 0.4]) .*lfMeta.sRateHz);
        meanTrace{iStims, iLayers} = mean((mergeTrace{iLayers}(selRange(1):selRange(2),:)));
        meanDiff{iStims, iLayers} = mean((diffTrace{iLayers}(selRange(1):selRange(2),:)));
        meanDiff{iStims, iLayers} = [mean(meanDiff{iStims, iLayers}), meanDiff{iStims, iLayers}];
        selIdx{iStims, iLayers} = mean((diffTrace{iLayers}(selRange(1):selRange(2),:))) ./ mean(abs(sumTrace(selRange(1):selRange(2),:)));
    end
    
    % show average images
    for xx = 1 : 2
        if xx == 1
            % show regular responses
            figure(h1);
            cStim = nanmean(mergeStim, 3);
            currRange = cRange;
        elseif xx == 2
            % show differential response
            figure(h2);
            cStim = nanmean(diffStim, 3);
            currRange = cRange/2;
        end
        
        subplot(1, nrStims, iStims)
        cStim = smoothCol(cStim, 2, 200);
        cImg = imagesc(cStim);
        ax = cImg.Parent;
        ax.XTick = 1 : lfMeta.sRateHz*0.5 : size(cStim,2);
        useTime = (((0 : lfMeta.sRateHz*0.5 : size(cStim,2)-1)./lfMeta.sRateHz) - opts.baseDur)*1000;
        ax.XTickLabel = useTime;
        
        xRange = (opts.baseDur + timeRange) .*lfMeta.sRateHz;
        xlim(xRange);
        colormap(ax, flipud(viridis(256)));
        caxis([-currRange currRange]);
        
        title(stimLabel(iStims))
        xlabel('time(ms)')
        nvline(opts.baseDur*lfMeta.sRateHz, 'w', 'linewidth', 2);
        
        depthSteps = round(200 / mean(diff(depthRange)));
        ax.YTick = 1 : depthSteps : size(cStim,1);
        ax.YTickLabels = depthRange(1 : depthSteps : length(depthRange));
        ylabel('Depth (um)');
        axis square;
        nhline(useLayers, 'w', 'linewidth', 2);
        drawnow;
    end
end

%% plot differences as errorbar plot
layerLabels = cell(1, length(layerDepths));
layerLabels{1} = ['0 - ' num2str(layerDepths(1))];

% Generate labels for each layer
for i = 1:length(layerDepths) - 1
    layerLabels{i+1} = sprintf('%d-%d', layerDepths(i), layerDepths(i+1));
end

% Display the result
figure;
errorbar(cellfun(@mean,meanDiff)', cellfun(@sem, meanDiff)', '-o', 'linewidth', 2)
xlim([0 6]); axis square;
title('Mean response difference');
view(90,90);
ax = gca;
ax.XTick = 1 : nrLayers;
ax.XTickLabel = layerLabels;
legend(opts.stimNames, 'location', 'southeast');
ax.YTick = 0:5:30;
xlabel('Depth from SC surface');
ylabel('Response difference (uV)');
grid on;
niceFigure
