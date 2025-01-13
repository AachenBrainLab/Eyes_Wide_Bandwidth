%% Main-script to compute Gabor-responses to the MotionClouds
%  I compute the response magnituges of the individual gabors in response
%  to the Narror- and Broad Orientation Bandwidth MotionClouds. This is
%  performed for 5 center spatial frequencies: [.01 .02 .04 .06 .08] cpd
%  180 orientation in 1deg steps. And for all possible phases as I'm
%  computing the energy of the complex-Gabor kernel response. 
%  Based on these Response curves I subsequently compute the amplitude
%  prediction and fit the cell-recruitment prediction. 
%  By: Gerion Nabbefeld, 2023.
% addpath("...\GitHub\MotionCloud_Stimuli_PC_analysis\finding_gabor_max_responses");

close all;

% stimulus settings
n_frames = 120;  % 120 frames (2s) per stimulus realization
n_ori = 180;     % 180 orientation -> 1deg orientation steps (don't need the other 180 as they are identical to the opposite ortientation with a shifter phase)
n_folders = 8;   % individual realizations of MotionClouds per condition (8=all)
bin_f = 4;       % spatial binning of stimulus frames

% stimulus/RF settings
spatial_magnitudes = 1;  % this one frlags if I want to keep the spatial map of the gabor responses
img_size = [90 160];  % after binning, assuming bin_f=4

fullfield_coverage = 130;  % the stimulus covers 130deg on the monitor in the horizontal
rf_radius = img_size(2) * 15 ./ fullfield_coverage;  % I need to cut down on what I save to save memory
center_region_x = round((img_size(1) ./ 2)-rf_radius): round((img_size(1) ./ 2) +rf_radius);
center_region_y = round((img_size(2) ./ 2)-rf_radius): round((img_size(2) ./ 2) +rf_radius);


% Spatial-Aspect-Ratio of Gabors. Derived from simulating sinewave tuning curves
SpatialAspectRatio = 0.55;
SF = [0.01 0.02 0.04 0.06 0.08];  % center spatial frequency of the gabors
n_scales = length(SF);

% Gabor settings
Gabor_str = "_Matlab_SpatialAspectRatio_"+SpatialAspectRatio;
gabor_radius = 200/bin_f;  % For visualization only!


% Visualization of spatial filter kernels
figure(70); fig_size();
for scale = 1:length(SF)
    subplot(length(SF), 1, scale);

    wavelength = (640 ./ (bin_f * 130)) ./ SF(scale);  % determine [px/cyc]
    theta = (0/180)*pi; % ori

    gb = gabor(wavelength, 0, "SpatialAspectRatio", SpatialAspectRatio, "SpatialFrequencyBandwidth", 2.38);
    gb = gb.SpatialKernel;
    n_padding = 200;
    temp = zeros(size(gb) + 2*n_padding);
    temp(n_padding+(1:size(gb,1)), n_padding+(1:size(gb,1))) = gb;
    gb = temp;
    gb = gb(int32((-gabor_radius:gabor_radius)+size(gb,1)./2), int32((-gabor_radius:gabor_radius)+size(gb,1)./2));

    imagesc(imag(gb)); colormap("gray"); axis image; colorbar(); caxis([-1 1]);

    title(SF(scale) + " cpd");
end
exportgraphics(gcf, "Gabors"+Gabor_str+".png", "Resolution", 300);
savefig(gcf, "Gabors"+Gabor_str+".fig");


%%
% Construct the full gabor bank once here (sfb hardcoded for 0.55 sar Gabors!)
gb = gabor(((640 ./ bin_f) / 130) ./ SF, 180 .* (0:n_ori-1) ./ n_ori, "SpatialAspectRatio", SpatialAspectRatio, "SpatialFrequencyBandwidth", 2.38);

% preallocate Gabor filter responses
results = nan(length(center_region_x), length(center_region_y), n_scales, n_ori, 4, n_frames*n_folders, "single");

for stim_id = 1:4
frame_Cnt = 0;
for folder_id = (1:n_folders) + (stim_id-1)*8
    for frame_id = 1:n_frames
        stimulus_folder = "D:\Data\MC_MotionMule_StimFrame_for_analysis";
        img_path = stimulus_folder + "\" + sprintf("%03d%sframe_%03d.png", folder_id, filesep, frame_id);
        % A subset of frames had more leading zeros in the name. This accounts for these:
        if ~exist(img_path, "file")
            img_path = "D:\Data\MC_MotionMule_StimFrame_for_analysis\" + sprintf("%03d%sframe%06d.png", folder_id, filesep, frame_id);
        end

        % convert to range: [-1 1]
        img = 2*(single(imread(img_path)) ./ 255)-1;

        % apply spatial downsampling
        img = imresize(img, 1/bin_f);
        
        % visualization
        figure(1); subplot(4, 4, stim_id);
        imagesc(img); colormap("gray"); axis image; colorbar();
        caxis([-1 1]); title(folder_id + " - " + frame_id);
        drawnow;
    
        % store results
        frame_Cnt = frame_Cnt + 1;
        mag = imgaborfilt(img, gb);

        mag = reshape(mag(center_region_x, center_region_y, :), length(center_region_x), length(center_region_y), 5, n_ori);  % unpack individual gabor responses
        results(1:size(mag, 1), 1:size(mag, 2), :, :, stim_id, frame_Cnt) = single(mag);
        disp([stim_id, frame_Cnt]);
    end
end
end

if n_ori == 180
    % all 4 stims
    gabor_responses = results;
    save("gabor_responses"+Gabor_str+"_all4_stims_SpatialMaps.mat", "gabor_responses", "-v7.3");
else
    gabor_responses = results;
    save("gabor_responses"+Gabor_str+"_all4_stims_SpatialMaps_"+n_ori+"oris.mat", "gabor_responses", "-v7.3");
end

% I shift the orientations by 90deg, it's nicer for visualization
results = cat(2, results(:, 91:180, :), results(:, 1:90, :));


%% Visualization of individual Gabor responses
figure(1); fig_size(gcf, [0 0 0.3 1]);
bin_f = 4;
scaleIds = 1:length(SF);
for stim_id = 1:2  % 1:4
    subplot(4, size(results, 3), size(results, 3)+stim_id);
    imagesc(results(:, :, stim_id)); % axis image;
    colormap("gray")
    xticks(0:n_ori/8:n_ori); xticklabels((180/n_ori) .* ((0:n_ori/8:n_ori))); xlabel("Orientation [deg]")
    
    yticks(1:length(SF)); yticklabels(SF); ylabel("Spatial freq. [cyc/img]")
    cb = colorbar(); ylabel(cb, "Energy");
    
    
    %%
    subplot(4, size(results, 3), 2*size(results, 3)+stim_id);
    plot(sum(results(:, :, stim_id), 2));
    xticks(1:length(SF)); xticklabels(SF); xlabel("Spatial freq. [cyc/img]")
    
    if stim_id == 1
        refResp = sum(sum(results(:, :, stim_id)));
    else
    end
    cResp = sum(sum(results(:, :, stim_id)));
    title("sum: " + sprintf("%0.1f %%", 100 * (cResp / refResp)));
    
    
    %%
    subplot(4, size(results, 3), 3*size(results, 3)+stim_id);
    plot(results(scaleIds, :, stim_id)');
    xticks(0:n_ori/8:n_ori); xticklabels((180/n_ori) .* ((0:n_ori/8:n_ori))); xlabel("Orientation [deg]")
    
    if stim_id == 1
        legend(num2str(SF') + " cpd");
    end
end

exportgraphics(gcf, "GaborResponses"+Gabor_str+".png", "Resolution", 300);
savefig(gcf, "GaborResponses"+Gabor_str+".fig");






