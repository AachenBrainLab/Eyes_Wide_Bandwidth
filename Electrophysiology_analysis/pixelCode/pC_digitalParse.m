function eventTimes = pC_digitalParse(digitalChannel, Fs)
% function eventTimes = pC_digitalParse(digitalChannel, Fs)
%
% returns the event times of all 16 digital inputs recorded by SpikeGLX
% They are returned in a length=16 cell array. Each cell has three cells,
% which contain the times of all events (sorted), the times of onsets, and
% the times of offsets, respectively. 
%
% Fs is optional input argument - if supplied, output will be in units of
% seconds, otherwise samples

digitalChannel = digitalChannel(:); % column

if isa(digitalChannel, 'int16')
    % change to uint16 but don't change the underlying bits
    digitalChannel = typecast(digitalChannel, 'uint16');
end

bitRepresentation = pC_dec2bin(single(digitalChannel), 16);

Cnt = 0;
eventTimes = cell(1, 16);
for b = 16:-1:1
    Cnt = Cnt + 1;
    eventTimes{Cnt} = pC_digitalToTimestamp(bitRepresentation(:,b), Fs);
end


