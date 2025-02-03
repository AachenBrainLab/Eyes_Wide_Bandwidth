function syncDat = pC_extractDigitalChannel(folder, numChans, syncChanIndex, datFile)
% usage: syncDat = pC_extractDigitalChannel(folder, numChans, syncChanIndex, datFile)
% extract digital channels from a int16 binary file. Every digital channel
% contains 16 bits and is returned as a cell array. Each cell contains
% event times from a given channel (cell 1 is channel 0).
% This can be used to get digital channels from the 'nidq' binary file from spikeGLX.

if ~exist('datFile', 'var') || isempty(datFile)
    datFile = 'nidq'; %load from NI-DAQ by default
end

dataFiles = dir(fullfile(folder,['*.' datFile '.bin']));
syncDat = cell(1, length(dataFiles));
for d = 1:length(dataFiles)
    
    filename = fullfile(folder, dataFiles(d).name);
    syncDat{d} = pC_digitalFromRaw(filename, numChans, syncChanIndex);
    
end