function [nVc, U, opts, SessionData] = load_trialized_bpod_wf_data(session_folder, preStimDur, postStimDur)
%%
% load(fullfile(session_folder, "Vc.mat"), "Vc", "trigOn", "U", "ledInfo", "frameInfo", "trialOn");
load(fullfile(session_folder, "Vc.mat"), "Vc", "trigOn", "U", "ledInfo", "frameInfo");
load(fullfile(session_folder, "opts2.mat"), "opts");

if ~isfield(opts, "cAnimal")
    [~, opts.cAnimal] = fileparts(fileparts(fileparts(fileparts(session_folder))));
end

[~, dateStr] = fileparts(session_folder);
bhvPath = abs_paths_str(fullfile(session_folder, opts.cAnimal + "*" + dateStr + ".mat"));
load(bhvPath, "SessionData");
opts.SessionData = SessionData;  % pass this on for other analysis that might need it

% settings
opts.sRate = 1 / mean(diff(frameInfo(:,2)));            % sampling rate in Hz
% opts.preStimDur = 2.0;                                  % baseline duration in seconds
% opts.postStimDur = 3.0;                                 % postStimulus duration in seconds
if exist("preStimDur", "var");  opts.preStimDur  = preStimDur;  else; opts.preStimDur  = 2; end
if exist("postStimDur", "var"); opts.postStimDur = postStimDur; else; opts.postStimDur = 3; end
opts.preStim = round(opts.sRate .* opts.preStimDur);    % nr. baseline frames
opts.postStim  = round(opts.sRate .* opts.postStimDur); % nr. postStimulus frames

%% reshape Vc into trials
% find frames that are closes to trialstarts
nrTrials = size(trigOn, 1);
startFrames = zeros(1, nrTrials);
for iTrials = 1 : size(trigOn, 1)
    % find next blue frame after trigger
    startFrames(iTrials) = find((ledInfo(:,3) - trigOn(iTrials, 3)) > 0, 1);
end
% I observed an offset in the responses
offset = 4;

trialIdx = (-opts.preStim:opts.postStim)' + startFrames + offset;
invalidFrameIds = trialIdx > size(Vc, 2);

% replace those frame ids with placeholders
trialIdx(invalidFrameIds) = 1;

nVc = Vc(:, trialIdx(:));  % select trials frames from continuous Vc
nVc(:, vec(invalidFrameIds)) = NaN;  % set invalid frames to NaN

% reshape into trials
nVc = reshape(nVc, [size(Vc,1), size(trialIdx)]);


