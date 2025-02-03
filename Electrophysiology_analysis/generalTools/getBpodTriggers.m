function [trialNumbers, trialOnTimes, bhv] = getBpodTriggers(fPath, trialTriggerChan)
% function to get trial triggers from a digital line that was recorderded
% with the NI-DAQ with spikeGLX.

% TTL interpretation (to read barcodes)
shortInt = 2.5;
longInt = 5.5;

%% get behavioral data and compare trial counts between SessionData and trial triggers in SpikeGLX file
bhvFiles = findBhvFile(fPath);
load(fullfile(bhvFiles(1).folder, bhvFiles(1).name));

%% get digital triggers from niDAQ data
nidaqDir = dir(fullfile(fPath,'*_g0'));
nidaqDir = fullfile(fPath,nidaqDir.name);
metaName = dir(fullfile(nidaqDir,'*.nidq.meta'));
niMeta = gt_readSpikeGLXmeta(fullfile(metaName.folder, metaName.name));

trigDat = pC_extractDigitalChannel(nidaqDir, niMeta.nChans, niMeta.nChans, 'nidq'); %get digital channel from ni-daq. usually the lat channel
niSyncDat = trigDat{1}; %sync trrigers should be in the first digital line

% make sure first event is a trigger onset. otherwise delete first sample from trialEvents.
trialEvents = floor(trigDat{1}{trialTriggerChan}{1} *  1000); %get trigger onset times
if trigDat{1}{trialTriggerChan}{2}(1) > trigDat{1}{trialTriggerChan}{3}(1) 
    trialEvents(1) = [];
end

trialStartTrace = false(1, round(niMeta.fileTimeSecs * 1000));
for iEvents = 1 : 2 : length(trialEvents)
    if iEvents + 1 > length(trialEvents)
        trialStartTrace(trialEvents(iEvents) : end) = true;
    else
        trialStartTrace(trialEvents(iEvents) : trialEvents(iEvents+1)-1) = true;
    end
end

trialOnTimes = trigDat{1}{trialTriggerChan}{2}; %get trigger onset times
if length(trialOnTimes) > SessionData.nTrials * 1.5 %check for barcode triggers. This would lead to much more trial triggers compared to SessionData structure.
    [trialNumbers, trialOnTimes] = segmentVoltageAndReadBarcodes(trialStartTrace, shortInt, longInt, 0.1, 0.5);
    trialOnTimes = trialOnTimes(trialNumbers > -1) / 1000;
    trialNumbers = trialNumbers(trialNumbers > -1);
else
    trialNumbers = 1 : length(trialOnTimes); %raw onset triggers. Assume that each trigger is a new trial.
end

% compare number of trial triggers
bhv = SessionData;
if length(trialNumbers) < SessionData.nTrials
    fprintf('!!! Bpod data has %d more trials as trial triggers in the SpikeGLX data. Make sure this is OK !!!\n', SessionData.nTrials - length(trialNumbers));
    bhv = selectBehaviorTrials(SessionData, trialNumbers);
elseif any(trialNumbers > SessionData.nTrials)
    cIdx = trialNumbers > SessionData.nTrials;
    fprintf('!!! SpikeGLX has %d more trial triggers as trials in bpod data. Make sure this is OK !!!\n', sum(cIdx));
    trialNumbers(cIdx) = [];
    trialOnTimes(cIdx) = [];
end