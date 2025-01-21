# Eyes_Wide_Bandwidth 
This repository contains imaging analysis scripts from Balla et al. Nature Communications, 2025. 
Step 1: get the processed data from Zenodo [http://doi.org/doi 10.5281/zenodo.14605879]. 


# 2Photon_analysis 
Gabor_model_comparisons.m: in reference to Figure 3, plots the raincloud figures as control,
runs the 	linear mixed effects model function for statistics and plots the model predictions.
- requires ('Gabor_model_comparisons_data.mat') from Zenodo
  
Multiple_regressors_polar_plot.m: in reference to Figure 4, plots the raincloud figures as control
and runs the linear mixed effects model with multiple regressors, plots the polar plot for the effect 
of orientation tuning, FF input and surround suppresion on orientation bandwidth reponses enrichment.	
- requires ('Responses_all_sessions.mat') from Zenodo

motion_clouds_disributions.m picks up the data from the excel sheets in the respective stimulus folders 
and plots the orientation distributions as seen in figure 1 panel b.
- requires (narrow_bandwidth, mid_bandwidth and broad bandwidth files ) from Zenodo
  
orientation_analysis.ijm: a batch file for running directionality on the example stimuli

***helper functions***

cardinal_line.m: plots horizontal or vertical lines at desired postitions/ color specs etc.

LME_compare.m: does a linear mixed effects model comparison of variables, while controlling for the impact of the random variable 'randomVar'.

LME_compare_3vars.m: compares the role of three variables to explain target data , while controlling for the impact of the random variable 'randomVar'. 

niceFigures.m:f function that changes axes and other figure parameters for better aesthetics 

raincloudPlot.m: plots vertical raincloud plots (i.e. half violins + boxplots + scatterplots) 


# Scripts for Portilla and Simoncelli analysis, presented in S7
PS_stimulus_analysis: Contains scrips for the Portilla Simoncelli Analysis in S7.  

PS_stimulus_analysis/PC_stats_individual_frames.m: computes the PS-stats for each stimulus frame individually. Used by scripts below.  

PS_stimulus_analysis/Stimulus_sepearbility_PCA.m: generates PCA plot for S7.  

PS_stimulus_analysis/PS_CrossEnergy_statistics_for_Figure_S7.m: Plots and compares the CrossScales and CrossOrientation Energies.  

PS_stimulus_analysis/compute_meanAutoCorrReal.m: Fetches the individual frame auto-correlations from the PS-analysis and averages by conditions for visualization.  


# Gabor filter model
Gabor_Filter_Model/overview.m: contains an overview over the different script used to characterize the Gabor filters and fit the model to the measured responses. 
