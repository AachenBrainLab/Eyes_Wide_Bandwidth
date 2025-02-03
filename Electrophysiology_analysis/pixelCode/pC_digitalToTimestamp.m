function cEvents = pC_digitalToTimestamp(binData, Fs)
% convert binary data to timestamps

db = diff([0; single(binData(:))]);
bOn = find(db>0);
bOff = find(db<0);

if nargin>1
    cEvents{1} = sort([bOn; bOff])/Fs;
    cEvents{2} = bOn/Fs;
    cEvents{3} = bOff/Fs;
else
    cEvents{1} = sort([bOn; bOff]);
    cEvents{2} = bOn;
    cEvents{3} = bOff;
end