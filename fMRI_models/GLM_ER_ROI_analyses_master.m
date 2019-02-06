% Scatter plots to illustrate relationship with self-control
% between RegSuccess - View Emo ME contrast and dietary self-control level

% Load behavioral dataframe to retrieve information 
% on self-control success rate
df = importdata('../Behavior_models_R_code/Maier_Hare_dataframe_behavior_food_choice.txt',' ',1);

% Initialize trial counter (always read out last position of dataframe for
% each participant, participant-level variables are always repeated on 
% all 100 trials)
c = 0;

% Column 1 of the dataframe df.data contains all participant identifiers
for i = 1:length(unique(df.data(:,1))) 
    
    c = c+100;
    
    % Dataframe column 13 already contains the self-control success rate for
    % each participant (in percent)
    perc_sc_all(i) = df.data(c,13);
    
    % Dataframe column 19 tells whether this participant was excluded from 
    % food choice MRI analyses (1 = yes, 0 = no)
    excludeFCmri(i) = df.data(c,19);
    
    % Dataframe column 20 tells whether this participant was excluded from
    % emotion regulation MRI analyses (1 = yes, 0 = no)
    excludeERmri(i) = df.data(c,20);
end

% Read out dietary self-control scores from 31 participants who are included
% in both the emotion reappraisal and dietary self-control fMRI analyses
perc_sc = perc_sc_all(excludeFCmri == 0 & excludeERmri == 0);



% Amygdala
mask = 'bilateral_Amygdala_HOA_50perc';

[means, Sbetas] = roi_contrast_meanLocal_RegSuccessViewEmo(mask);

means_amy = means;
Sbetas_amy = Sbetas;

clear means Sbetas



% dlPFC (functionally active regions in Hare et al 2009 or Maier et al 2015)
mask = 'Hare2009_Maier2015_SelfControl_functional_regions';
[means, Sbetas] = roi_contrast_meanLocal_RegSuccessViewEmo(mask);

means_dlPFC = means;
Sbetas_dlPFC = Sbetas;

clear means Sbetas


% Hippocampus
mask = 'bilateral_HC_HOA_50perc';

[means, Sbetas] = roi_contrast_meanLocal_RegSuccessViewEmo(mask);

means_hc = means;
Sbetas_hc = Sbetas;

clear means Sbetas


% vmPFC
mask = 'vmPFC_anat_HOA';

[means, Sbetas] = roi_contrast_meanLocal_RegSuccessViewEmo(mask);

means_vmPFC = means;
Sbetas_vmPFC = Sbetas;

clear means Sbetas


% Striatum
mask = 'bilateral_striatum_HOA';

[means, Sbetas] = roi_contrast_meanLocal_RegSuccessViewEmo(mask);

means_str = means;
Sbetas_str = Sbetas;

clear means Sbetas






axis_fontsize = 32;

fh1 = figure(1);
set(gca, 'fontsize', axis_fontsize)
box('off')
scatter(Sbetas_hc, perc_sc, 200, 'k', 'filled');
hold on

brob = robustfit(Sbetas_hc, perc_sc);
plot(Sbetas_hc, brob(1)+brob(2)*Sbetas_hc, 'k', 'LineWidth', 4)
xlabel('Betas Emotion Regulation Success > View', 'FontSize', axis_fontsize)
ylabel('Dietary Self-Control Success', 'FontSize', axis_fontsize)
ylim([0 100])
title('Hippocampus', 'FontSize', axis_fontsize)

set(fh1, 'position', [0,0,750,650])
ah1 = get(fh1, 'CurrentAxes');
set(ah1,'LineWidth',2);
set(ah1,'TickLength',[0,0]);


fh2 = figure(2);
set(gca, 'fontsize', axis_fontsize)
box('off')
scatter(Sbetas_vmPFC, perc_sc, 200, 'k', 'filled');
hold on

brob = robustfit(Sbetas_vmPFC, perc_sc);
plot(Sbetas_vmPFC, brob(1)+brob(2)*Sbetas_vmPFC, 'k', 'LineWidth', 4)
xlabel('Betas Emotion Regulation Success > View', 'FontSize', axis_fontsize)
ylabel('Dietary Self-Control Success', 'FontSize', axis_fontsize)
ylim([0 100])
title('vmPFC', 'FontSize', axis_fontsize)

set(fh2, 'position', [0,0,750,650])
ah2 = get(fh2, 'CurrentAxes');
set(ah2,'LineWidth',2);
set(ah2,'TickLength',[0,0]);


fh3 = figure(3);
set(gca, 'fontsize', axis_fontsize)
box('off')
scatter(Sbetas_str, perc_sc, 200, 'k', 'filled');
hold on

brob = robustfit(Sbetas_str, perc_sc);
plot(Sbetas_str, brob(1)+brob(2)*Sbetas_str, 'k', 'LineWidth', 4)
xlabel('Betas Emotion Regulation Success > View', 'FontSize', axis_fontsize)
ylabel('Dietary Self-Control Success', 'FontSize', axis_fontsize)
ylim([0 100])
title('Striatum', 'FontSize', axis_fontsize)

set(fh3, 'position', [0,0,750,650])
ah3 = get(fh3, 'CurrentAxes');
set(ah3,'LineWidth',2);
set(ah3,'TickLength',[0,0]);



fh4 = figure(4);
set(gca, 'fontsize', axis_fontsize)
box('off')
scatter(Sbetas_amy, perc_sc, 200, 'k', 'filled');
hold on

brob = robustfit(Sbetas_amy, perc_sc);
plot(Sbetas_amy, brob(1)+brob(2)*Sbetas_amy, 'k', 'LineWidth', 4)
xlabel('Betas Emotion Regulation Success > View', 'FontSize', axis_fontsize)
ylabel('Dietary Self-Control Success', 'FontSize', axis_fontsize)
ylim([0 100])
title('Amygdala', 'FontSize', axis_fontsize)

set(fh4, 'position', [0,0,750,650])
ah4 = get(fh4, 'CurrentAxes');
set(ah4,'LineWidth',2);
set(ah4,'TickLength',[0,0]);



fh5 = figure(5);
set(gca, 'fontsize', axis_fontsize)
box('off')
scatter(Sbetas_dlPFC, perc_sc, 200, 'k', 'filled');
hold on

brob = robustfit(Sbetas_dlPFC, perc_sc);
plot(Sbetas_dlPFC, brob(1)+brob(2)*Sbetas_dlPFC, 'k', 'LineWidth', 4)
xlabel('Betas Emotion Regulation Success > View', 'FontSize', axis_fontsize)
ylabel('Dietary Self-Control Success', 'FontSize', axis_fontsize)
ylim([0 100])
title('dlPFC', 'FontSize', axis_fontsize)

set(fh5, 'position', [0,0,750,650])
ah5 = get(fh5, 'CurrentAxes');
set(ah5,'LineWidth',2);
set(ah5,'TickLength',[0,0]);


[r,p,lo,hi] = corrcoef(Sbetas_amy, perc_sc)
[r,p,lo,hi] = corrcoef(Sbetas_dlPFC, perc_sc)
[r,p,lo,hi] = corrcoef(Sbetas_vmPFC, perc_sc)
[r,p,lo,hi] = corrcoef(Sbetas_hc, perc_sc)
[r,p,lo,hi] = corrcoef(Sbetas_str, perc_sc)
