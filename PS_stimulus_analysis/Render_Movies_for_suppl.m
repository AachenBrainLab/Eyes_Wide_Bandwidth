%% NOTE: If we include this script, we need to add the write_avi function to this repo!!!


%% This is a quick script to render out the movies for suppl. data
%  By: Gerion Nabbefeld, 2023.
close all;

% provide the base-folder containing the stimulus frames:
StimulusFolder = "D:\Data\MC_MotionMule_StimFrame_for_analysis";


%
n_frames = 120;  % 120;  % 120 frames per stimulus realization
n_folders = 8;   % individual realizations of MotionClouds per condition

%%
all_frames = zeros(360, 640, 120, 4);
for stim_id = 1:4
    frame_Cnt = 0;
    for folder_id = (1:n_folders) + (stim_id-1)*8
        for frame_id = 1:n_frames
            img_path = StimulusFolder + "\" + sprintf("%03d%sframe_%03d.png", folder_id, filesep, frame_id);
            if ~exist(img_path, "file")
                img_path = StimulusFolder + "\" + sprintf("%03d%sframe%06d.png", folder_id, filesep, frame_id);
            end

            % convert to range: [-1 1]
            img = 2*(single(imread(img_path)) ./ 255)-1;

            frame_Cnt = frame_Cnt + 1;
            all_frames(:, :, frame_id, stim_id) = img;
            disp([stim_id, frame_id]);
        end
    end
end
all_frames = permute(all_frames, [2, 1, 3, 4]);


%%
folder = fullfile(pwd, "ExampleMovie"); make_dir(folder);
MovieTitles = ["Narrow", "SF", "Ori", "Mixed"];
for stim_id = 1:4
    write_avi(all_frames(:, :, :, stim_id), [], folder, MovieTitles(stim_id), 60, [-1 1], gray(256));
end


