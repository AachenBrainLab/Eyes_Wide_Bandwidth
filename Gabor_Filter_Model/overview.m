%% overview:
% general setting concerning the tuning-width of the Gabors
% 0.55: ~44° FWHM
SpatialAspectRatio = 0.55;

% what model to use:
rerun = 0;  % This runs the fit again
raw_gabors = 1;
fit_outerSurroundLimit = 1;  % Sets if we want to have an outer limit for the surround area (this was good to test but didnt improve it)
use_SurrSupp = 2;  % 1: subtractive surrSupp and 2: divisive surround suppression!
use8oris = 0;  % This was just a test to use a larger surround, but it didnt help!
use_orientation_unspecific_surround = 0;  % This was used during testing

compute_individual_session_R2 = 0;

%% 1. Determine SpatialAspectRatio:
if 0  % set 0 to skip this if you dont want to rerun this script
    % Computes the tuning curves for different Gabors in response to
    % sine-wave grating. 
    determine_gabor_aspect_ratio
end


%% 2. Determine MaxResponses to sinewaves for normalization:
% DEPRICATED: Determine_max_grating_response_for_Gabors: 
%     this determins the max response of the imgaborfilt(), gabor kernels.
% INSTEAD:
if 0
    Determine_max_grating_response_for_Gabors
end


%% 3. Compute the Gabor responses for every frame in the stimulus set.
% NOTE: this is the most time intensive step! Depending on hardware, can
% take 2 days to run!
if 0
    % Here I'm computing the gabor energy-responses for each of our
    % Stimulus frames. These response curves are what the subsequent
    % analysis are based on.
    Compute_Gabor_responses_for_stimulus_frames.m
end


%% 4. to render out the individual Gabor images for visualization
if 0
    Save_gabors_individually
end


%% 5. Fit the amplitudes based on the Gabor-responses to the MotionClouds
match_means = 0;  % rather to adjust the fit itself to only focus on the orientation-dependence. Not done in the final version anymore.
match_means_R2_only = 1;  % used for computing the R2

% This script fits the raw Gabor responses to the measured amplitudes with
% narror and Ori stimuli. It also plots the derived response modulation.
amplitudes_raw_gabors_with_surroundSupression


%% 6. Based on these fitted amplitudes I now predict the recruitment. 
% Assuming a linear relationship. The idea here is that I'm describing the 
% average response amplitudes of the population above and a higher response
% amplitude should related to a higher probability of the cells to register
% as responsive.
use_intercept = 0;
recruitment_modulation_from_amplitudes


%% 7. also plot the individual amplitude fits
% The models are fit individually to the amplitude responses
% towards narror and Ori stimuli. This fit is visualized here for
% Figure 3
plot_Amplitude_fit_for_Figure_3





