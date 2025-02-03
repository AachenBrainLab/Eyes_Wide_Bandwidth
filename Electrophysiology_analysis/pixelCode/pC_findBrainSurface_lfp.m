function meanStd = pC_findBrainSurface_lfp(filename, metaInfo, eventTimes)

testRange = 0.01; %amount of data for testing
winSize = 5; %testsize in seconds
nrWins = floor(metaInfo.fileTimeSecs * testRange / winSize);
timeSteps = metaInfo.fileTimeSecs / 2 / nrWins;
if ~exist('eventTimes', 'var') || isempty(eventTimes)
    eventTimes = (metaInfo.fileTimeSecs / 2  - winSize*2) + (1 : timeSteps : metaInfo.fileTimeSecs/2); %get evenly sampled events from second half of recording
end
winIdx = 1:round(winSize * metaInfo.sRateHz); 

%% get datafile and nr of samples
fprintf(1,' loading analog data to identify brain surface %s\n', filename);
d = dir(filename);
fid = fopen(filename, 'r');

Cnt = 0;
meanStd = zeros(metaInfo.nChans-1, 1);
for x = 1 : length(eventTimes)
    
    cOn = eventTimes(x) * metaInfo.sRateHz + winIdx(1) ; %onset for current window
    
    fseek(fid, (metaInfo.nChans*2)* round(cOn - Cnt), 'cof'); %move pointer to current window onset
    tempDat = fread(fid, [1, length(winIdx)*metaInfo.nChans], 'int16=>int16'); %read data according to window size
    tempDat = reshape(tempDat, metaInfo.nChans, []);
    
    [b, a] = butter(2,1000/(metaInfo.sRateHz), 'high');
    tempDat = filtfilt(b,a,double(tempDat(1:end-1,:)'))';
    meanStd = meanStd + std(double(tempDat), [], 2);

%     meanStd = meanStd + std(double(tempDat(1:end-1,:)), [], 2);
    Cnt = cOn + length(winIdx);

end
fclose(fid);

