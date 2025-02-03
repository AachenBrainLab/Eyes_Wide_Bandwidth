function analogDat = pC_extractAnalogChannel(filename, numChans, chanIdx, eventTimes, winIdx, opts)
% usage: analogDat = pC_extractAnalogChannel(filename, numChans, chanIdx, eventTimes, winIdx)
% function to extract analog channel from binary data file.
% fileName is the name of the binary file, numChans says how many channels
% there are in total, chanIdx identifies the channel of interest, eventTimes
% are the timpoints around which data should be returned and winIdx is a
% timeIndex that defines how much data is returned around each event in
% eventTimes.

maxReadSize = 1e9; %maximum number of bytes to be read in one step

if ~exist('eventTimes', 'var')
    eventTimes = []; %no events, get the whole trace
end

if ~exist('winIdx', 'var')
    winIdx = []; %no events, get the whole trace
end

if ~exist('opts', 'var')
    opts.verbose = true; %show feedback
end

% get datafile and nr of samples
if opts.verbose
    fprintf(1,' loading analog channel %d from %s\n', chanIdx, filename);
end
d = dir(filename);
nSamp = d.bytes/2/numChans;
if isempty(eventTimes)
    analogDat = zeros(1, nSamp, 'int16');
else
    analogDat = zeros(length(winIdx), length(eventTimes), 'int16');
end

% skip over the first samples of the other channels
fid = fopen(filename, 'r');
fseek(fid, (chanIdx-1)*2, 'bof');

if isempty(eventTimes)
    % whole trace
    nBatch = floor(nSamp/maxReadSize);
    for b = 1:nBatch
        analogDat(1, (b-1)*maxReadSize+1:b*maxReadSize) = fread(fid, [1, maxReadSize], 'int16=>int16', (numChans-1)*2); % load data, skipping other channels
    end
    analogDat(1, nBatch*maxReadSize+1:end) = fread(fid, [1, Inf], 'int16=>int16', (numChans-1)*2);
    
else
    metaName = strrep(filename, '.bin', '.meta');
    metaVars = pC_readSpikeGLXmeta(metaName);
    
    % get specific times
    Cnt = 0;
    for x = 1 : length(eventTimes)
        
        cOn = eventTimes(x) * metaVars.sRateHz + winIdx(1) ; %onset for current window        
        fseek(fid, (numChans*2)* round(cOn - Cnt), 'cof'); %move pointer to current window onset
        analogDat(:,x) = fread(fid, [1, length(winIdx)], 'int16=>int16', (numChans-1)*2); %read data according to window size
        
        Cnt = cOn + length(winIdx);
        
    end
end
fclose(fid);


