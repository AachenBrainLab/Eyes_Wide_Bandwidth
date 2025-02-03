function syncDat = pC_digitalFromRaw(filename, numChans, syncChanIndex, loadRaw)
% usage: syncDat = pC_digitalFromRaw(filename, numChans, syncChanIndex, loadRaw)
% extract digital channels from a binary file, like the nidq.bin file from spikeGLX. 
% extraChanIndices are 1-indexed

if ~exist('loadRaw', 'var')
    loadRaw = false; %use this if triggers should be loaded from raw file
end

maxReadSize = 1e9; %maximum number of bytes to be read in one step
trigThresh = 32; %threshold to detect logical events

d = dir(filename); 
[folder,fn] = fileparts(filename);
syncFname =  fullfile(folder, [fn '_di' num2str(syncChanIndex) '.mat']);

% check if sync file exist and get data from there if not empty
checker = false;
if exist(syncFname, 'file')
    syncSize = dir(syncFname);
    
    if syncSize.bytes > 0
        checker = true;
    end
end

if checker && ~loadRaw
    % get syncData from dedicated file
    load(syncFname, 'syncDat');
    
else
    
    fprintf(1,' loading raw file %s\n', filename);
    
    %check size and make syncData array
    nSamp = d.bytes/2/numChans;
    syncDat = zeros(1, nSamp, 'int16');
    
    % skip over the first samples of the other channels
    fid = fopen(filename, 'r');
    fseek(fid, (syncChanIndex-1)*2, 'bof');

    % loop trough data in batches
    nBatch = floor(nSamp/maxReadSize);
    for b = 1:nBatch
        syncDat((b-1)*maxReadSize+1:b*maxReadSize) = fread(fid, [1, maxReadSize], 'int16=>int16', (numChans-1)*2); % skipping other channels
    end

    % all the other samples
    syncDat(nBatch*maxReadSize+1:end) = fread(fid, [1, Inf], 'int16=>int16', (numChans-1)*2); % skipping other channels

    fclose(fid);
    
    % translate digital channel bits and then to events
    metaName = strrep(filename, '.bin', '.meta');
    metaVars = pC_readSpikeGLXmeta(metaName);
    syncDat = pC_digitalParse(syncDat, metaVars.sRateHz);

    save(syncFname, 'syncDat'); %save for later
end