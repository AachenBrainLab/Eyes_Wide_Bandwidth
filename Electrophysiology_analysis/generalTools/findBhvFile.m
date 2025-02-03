function bhvFiles = findBhvFile(bhvPath)
% function to find behavioral files from a Bpod system. This checks the
% content of mat files in a folder and checks for the variable SessionData
% which usually contains behavioral data.

bhvFiles = dir(fullfile(bhvPath,'*.mat'));
fileSelect = false(1, length(bhvFiles));
for iFiles = 1 : length(bhvFiles)
    try
        cFile = fullfile(bhvPath, bhvFiles(iFiles).name);
        a = whos('-file', cFile);
        fileSelect(iFiles) = any(strcmpi({a(:).name}, 'SessionData'));
    end
end

if sum(fileSelect) == 0
    disp('No files containing a SessionData variable found. Check behavioral path.');
elseif sum(fileSelect) > 1
    disp('Found more than 1 file, containing a SessionData variable. Check behavioral path.');
end
bhvFiles = bhvFiles(fileSelect);
