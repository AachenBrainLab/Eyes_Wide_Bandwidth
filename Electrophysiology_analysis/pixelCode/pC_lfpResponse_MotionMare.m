function [lfpResp, meanStd, lfMeta, hasSC] = pC_lfpResponse_MotionMare(basePath, opts)

%% set paths for where to find the data
% get paths
fprintf('Current path: %s. ', basePath);
lfMetaName = dir([basePath, '*' filesep '*' filesep '*_g0_imec0' filesep '*.lf.meta']);
lfpDir = lfMetaName.folder;
    
[~, recName] = fileparts(basePath);
saveFile = fullfile(basePath, 'lfpData.mat');

opts.reload = false; %set this flag to false since binary data could not be included in the online repo
if ~opts.reload && exist(saveFile, 'file')
    try
        % try to get data from local save file
        disp('Loading local data ... ');
        load(saveFile, 'lfpResp', 'meanStd', 'lfMeta', 'hasSC');
        if  ~exist('hasSC', 'var')
            hasSC = false;
        end
            
    catch ME
        opts.reload = true;
        disp('Error when loading saved file occurred. Loading raw data instead.')
        disp(ME.message);
    end
end

if opts.reload || ~exist(saveFile, 'file')    
    %% get meta files and digital triggers
    fprintf('Loading raw data ... ');
    lfMeta = gt_readSpikeGLXmeta(fullfile(lfMetaName.folder,lfMetaName.name));
    
    %% load settings from stimulation software
    % get bpod behavior
    [~, trialOnTimes, bhv] = getBpodTriggers(basePath, opts.trialTriggerChan);
    
    % get stim times
    for x = 1 : length(bhv.stim_id)
        barTime(x) = bhv.RawEvents.Trial{x}.States.trialCode1(1);
        stimTime(x) = bhv.RawEvents.Trial{x}.States.PlayStimulus(1);
    end
    realStimTime = trialOnTimes + stimTime - barTime; %real time of the stimulus osnet in the recording
    
    %% get analog data
    lfFile = dir(fullfile(lfpDir,'*.lf.bin'));
    lfFile = fullfile(lfpDir, lfFile.name);
    baseSize = round(opts.baseDur * lfMeta.sRateHz);
    winSize = round(opts.postStim * lfMeta.sRateHz);
    
    %% get mean standard deviations of different parts of the recording to determine the brain surface.
    % check for bad channels results - tends to be more accurate if saline was applied
    badChansFile = fullfile(lfpDir, 'bad_channel_ids_list.npy');
    
    if exist(badChansFile, 'file')
        badChans = readNPY(badChansFile);
        badChans(badChans == 191) = []; %ignore the reference channel
        
        % find brain surface based on last consecutive block of bad channels
        badChans = flipud(badChans+1);
        cIdx = min([find(diff(badChans) < -1, 1), length(badChans)]);
        cBrainIdx = true(1, lfMeta.nSites);
        if ~isempty(badChans) && badChans(1) >= 384
            cBrainIdx(badChans(cIdx):end) = false;
        end
        cBrainIdx = fliplr(cBrainIdx);
        %     cBrainIdx = cBrainIdx(chanSteps);
        meanStd = cBrainIdx';
        
    else
        % if badchannels are not available use LFP method
        meanStd = medfilt1(pC_findBrainSurface_lfp(lfFile, lfMeta));
        meanStd = flipud(meanStd);
        meanStd = meanStd ./ std(meanStd);
        meanStd = meanStd - meanStd(1);
        cBrainIdx = meanStd' > opts.brainThresh; %contacts in the brain
    end
    
    % get depth of individual contacts
    depthRange = (1 : length(cBrainIdx)/2) * 20;
    depthRange = sort(repmat(depthRange,1,2));
    brainStart = find(cBrainIdx, 1); %cortical surface
    depthRange = depthRange - depthRange(brainStart);
    
    depthIdx = depthRange >= 0 & depthRange < opts.brainRange;
    if sum(cBrainIdx) < sum(depthIdx)
        error('Not enough contacts in brain to cover requested depth');
    end
    cBrainIdx = cBrainIdx & depthIdx; %contacts at desired depth range
    
    %% load analog data for diffrent stim groups
    lfpResp = zeros(round(sum(cBrainIdx)/2), winSize+baseSize, numel(opts.stimGroups), 'single');
    winIdx = -baseSize+1:winSize;
    for iStim = 1 : 4
        
        Cnt = 0;
        cStim = realStimTime(ismember(bhv.stim_id, opts.stimGroups{iStim}));
        cData = zeros(sum(cBrainIdx), length(winIdx));
        nrChans = length(find(cBrainIdx));
        disp(['Stim group = ' num2str(iStim)]);
        for iChans = length(cBrainIdx) - find(cBrainIdx)
            Cnt = Cnt + 1;
            if ismember(iChans, badChans)
                cData(Cnt, :) = NaN;
            else
                temp = pC_extractAnalogChannel(lfFile, lfMeta.nChans, iChans, cStim, winIdx, opts); %get current channel
                cData(Cnt, :) = nanmean(temp,2); %get trial average
            end
            if rem(iChans, round(nrChans/3)) == 0
                fprintf('Current channel %i / %i\n', Cnt, nrChans);
            end
        end
        
        % average columns together. If first depth is only one contact, add a row of NaNs here
        if rem(Cnt, 2) ~= 0
            cData = [nan(1, size(cData,2)); cData];
        end
        cData = reshape(cData', size(cData,2), 2, []); %combine columns
        
        % get LFP and CSD
        lfpResp(:, :, iStim) = smoothCol(fillgaps(squeeze(nanmean(cData,2))'), 1, 5); %do some smoothing
        lfpResp(:, :, iStim) = lfpResp(:, :, iStim) - mean(lfpResp(:, 1:baseSize, iStim),2);
    end
    
    %save result file for faster loading
    save(saveFile, 'lfpResp', 'meanStd', 'lfMeta');
    fprintf('done\n');
end