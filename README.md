# Maier_Hare_Emotion_and_dietary_selfregulation
Code for Maier, S.U. & Hare, T.A. (2020): Emotion and dietary self-regulation



Whenever using this code please cite the accompanying manuscript:

Maier, S.U. &amp; Hare, T.A. (2020), Greater BOLD signal during successful emotional stimulus reappraisal 
is associated with better dietary self-control, bioRxiv, preprint current as of May 27, 2020 (v2)
DOI: https://doi.org/10.1101/542712



1) Experiment code

The code that presented the experiment in the fMRI session is contained in two functions:

esc1_er_v4.m presented the emotion regulation paradigm
esc1_fc_v4.m presented the dietary self-control paradigm

Our combination of the IAPS pictures in arousal-equated sets is documented in the function

esc1_er_allot.m

The ratings for the food images were collected with the rating scale that was executed by the function

ratings_escf1.m 
(which calls the draw function drawscale_escf1.m that draws the actual elements of the scale on the screen)

The 180 food stimuli are contained in the folder “food_images”, and the phase-scrambled versions 
of these stimuli live in “food_images_scr”. The health symbol that was displayed during the inter trial interval 
lives in the folder “fixation_item”.

For a test run of the fMRI code, we provide the design files “imagelist_er_test.mat” for the emotion paradigm 
and “imagelist_test.mat” for the dietary choice paradigm.

The function “randInt.m” is a legacy version of a Matlab function that our code relied upon 
(https://www.mathworks.com/matlabcentral/answers/376258-what-does-x-randint-1-1-1-n) 
and thus also included here for reference.


2) Behavioral analyses presented in the paper

The main analysis script that performs all behavioral analyses presented in the paper is
Maier_Hare_Emotion_and_dietary_self-regulation_models.R

To perform the Bayesian correlation and T-test-like analyses that are based on John Kruschke’s code 
described in his book

	Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: A Tutorial with R, JAGS, 	and Stan. Academic Press / Elsevier.

	and its accompanying webpage (http://doingbayesiandataanalysis.blogspot.com/2017/06/bayesian-estimation-of-correlations-and.html), 

the scripts “bayesCorrelation.R” and “DBDA2E-utilities.R” need to live in the same folder as the main analysis script. 
Note that we also added one modified script: “bayesRankCorrelation.R” performs a rank correlation.

The dataframe for the emotion regulation behavioral analyses is contained in:
Maier_Hare_dataframe_behavior_emotion_ratings.txt

The dataframe for the dietary choice analyses is called:
Maier_Hare_dataframe_behavior_food_choice.txt

Variable keys for both dataframes are provided in
Maier_Hare_dataframes_variable_keys.rtf


3) Sanity check of the emotion ratings for the IAPS stimuli in the Swiss/German sample

In the Behavioral analyses folder, we provide an additional analysis script for validating 
the valence ratings our Swiss/German sample made in comparison to valence ratings provided 
by young adults from Germany: This is done in Maier_Hare_CH_vs_DE_young_adults_IAPS_ratings.R 

The ratings from the German sample are based on the publication of 

	Gruehn, D. and Scheibe, S. (2008). Age-related differences in valence and arousal ratings of pictures from the International 		Affective Picture System (IAPS): Do ratings become more extreme with age? Behavior Research Methods, 40 (2), 512-521.

and the accompanying dataset. The link to the original publication of the dataset together 
with the article can be retrieved here:

	https://link.springer.com/article/10.3758%2FBRM.40.2.512
	

4) fMRI analyses

The code for the preprocessing of the fMRI data lives in the folder fMRI_preprocessing:

Data for the emotion regulation paradigm and dietary self-control paradigm were preprocessed 
with SPM using the batch-generating scripts that are documented in 

preprocessing_emotion_regulation.m 
and 
preprocessing_food_choice.m

Afterwards, additional parameters for physiological noise correction (heartbeat, breathing) 
were generated by applying the RETROICOR algorithms as implemented in the physIO toolbox 
from the TAPAS collection (https://www.tnu.ethz.ch/en/software/tapas/documentations/physio-toolbox.html)

The batch-generating scripts for these are: 
RETROICOR_with_physIO_tapas_spm12_emotion_regulation_task.m 
and 
RETROICOR_with_physIO_tapas_spm12_food_choice_task.m


The fMRI modeling scripts live in fMRI_models. For a detailed description please see the preprint referenced above.

The emotion regulation data were modeled with GLM_ER.m. We also ran a region of interest (ROI) analysis 
with GLM_ER_ROI_analyses_master.m that calls the script roi_contrast_meanLocal_RegSuccessViewEmo.m 
to pull out the betas for the specified contrast. The first-level data for this analysis are provided in 
the folder GLM_ER_first_level and the anatomical and functional masks that defined the regions of interest 
can be found in the folder anatomical_masks.

GLM_FC.m modeled the subjective food value representation and relies on the combination of 
taste and health weighted food ratings that were calculated with the regression of food choices 
on taste and health aspects that is specified ini weighted_food_value.m 
(raw data for this regression live in the folder “behavior”).

GLM_TH.m modeled health- and taste-related signals.

GLM_SCS.m modeled the self-control success and failure.

GLM_ST.m modeled the self-control stakes.


For drift diffusion model (DDM) analyses in the revised manuscript v2 (May 27, 2020), the script ESC_DDM_Bins_alltrials.R fits the DDM using the C++ code in 2ddm_r_cpp_2.cpp and requires the dataframe dataframe_escf1_basic_ddm.txt as input. The results for the current study can be found in ESC_fits_alltrials.csv, and the comparison fits from a previous study by Hare, Malmaud & Rangel 2011 live in the files IAC_fits_Cue1_all_3000.RData and IAC_fits_Cue3_all_3000.RData.

When using these pieces of code, please cite the accompanying paper:

Silvia U. Maier*, Anjali Raja Beharelle*, Rafael Polania, Christian C. Ruff**, Todd A. Hare** (2020) Dissociable mechanisms govern when and how strongly reward attributes affect decisions, Nature Human Behaviour, accessible online: https://rdcu.be/b4ysL
