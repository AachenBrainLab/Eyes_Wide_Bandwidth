# Eyes_Wide_Bandwidth
This repository contains imaging analysis scripts from Balla et al. Nature Communications, 2025.


# Collection of analysis routines part of the MotionClouds project
PS_stimulus_analysis: Contains scrips for the Portilla Simoncelli Analysis in S7
PS_stimulus_analysis/PC_stats_individual_frames.m: computes the PS-stats for each stimulus frame individually. Used by scripts below.
PS_stimulus_analysis/Stimulus_sepearbility_PCA.m: generates PCA plot for S7
PS_stimulus_analysis/PS_CrossEnergy_statistics_for_Figure_S7.m: Plots and compares the CrossScales and CrossOrientation Energies.
PS_stimulus_analysis/compute_meanAutoCorrReal.m: Fetches the individual frame auto-correlations from the PS-analysis and averages by conditions for visualization


# Gabor filter model
Gabor_Filter_Model/overview.m: contains an overview over the different script used to characterize the Gabor filters and fit the model to the measured responses. 