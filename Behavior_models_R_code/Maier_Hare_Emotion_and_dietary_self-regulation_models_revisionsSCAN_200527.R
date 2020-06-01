# Behavioral models for Maier & Hare, "Emotion and dietary self-regulation"

# Version 1.1, 27 May 2020: added revisions and analyses in response to reviewer questions

# Author: Silvia Maier, Zurich Center for Neuroeconomics 
# Contact: silvia.maier@econ.uzh.ch

# Abbreviations:
# ER = emotion regulation
# SC = dietary self-control

# neg = negative
# pos = positive
# neut = neutral
# Reg = regulate condition
# View = view condition


### Load libraries

#Pirateplots
library("yarrr") 
#Bayesian regression models
library("brms") 
#Bayesian T-like tests
library("BEST") 


library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(matrixStats)


### PART 1: Emotion regulation experiment

# Load and prepare emotion ratings dataframe
escf1.emo.df=read.table("Maier_Hare_dataframe_behavior_emotion_ratings.txt", 
                        header=T)
escf1.emo.df = as.data.frame(escf1.emo.df)

# Restrict to participants included in emotion reappraisal task

# Note we have 36 participants here because 1062 was only excluded 
# from fMRI due to excess motion; the exclude_er variable only flags 
# participants who have to be excluded from the behavioral analyses as well
escf1.emo.df = escf1.emo.df[escf1.emo.df$exclude_er==0,]

# Factorize trial type: 
# 1 = negView, 2 = negReg, 3 = neutView, 4 = posReg, 5 = posView
escf1.emo.df$type = as.factor(escf1.emo.df$type) 

# Save separate dataframe for plotting 
escf1.emo.df.plot = escf1.emo.df

# Relevel type information for interpreting the behavioral model 
# with regard to neutral view as a baseline
escf1.emo.df$type = relevel(escf1.emo.df$type, 3)

# For the behavioral model:
# Exclude rows with NAs (when participant missed rating in scanner)
escf1.emo.df = escf1.emo.df[!is.na(escf1.emo.df$scan_rating),] 

# trial-wise success
escf1.emo.df$regSuccessNeg[escf1.emo.df$type==2] = escf1.emo.df$scan_rating[escf1.emo.df$type==2] - escf1.emo.df$post_rating[escf1.emo.df$type==2]
escf1.emo.df$regSuccessPos[escf1.emo.df$type==4] = escf1.emo.df$post_rating[escf1.emo.df$type==4] - escf1.emo.df$scan_rating[escf1.emo.df$type==4]


### DESCRIPTIVES

### CREATE Figure 2

# Re-name type information for plotting
levels(escf1.emo.df.plot$type) <- c( "Negative View", 
                               "Negative Regulate",
                               "Neutral View",
                               "Positive Regulate", 
                               "Positive View")

# Re-name variables for plotting
names(escf1.emo.df.plot) <- c("subID", 
                         "iapsID", 
                         "setID", 
                         "Mean Rating in Block", 
                         "Mean View Rating", 
                         "Block Type", 
                         "exclude_er")

# To report mean ratings per block type (during block):
# Aggregate ratings
meanViewRatings = aggregate(`Mean Rating in Block` ~ subID 
                            + `Block Type`, 
                            data = escf1.emo.df.plot, 
                            FUN = "mean")

# Plot mean ratings by type
pirateplot(`Mean Rating in Block` ~ `Block Type`, 
           data = meanViewRatings, 
           cex.names = 1.5, 
           cex.axis = 1.5, 
           cex.lab = 1.5, 
           inf.method = "se", 
           inf.f.col = "gray", 
           inf.b.col = "black",
           bean.b.col = "white", 
           bean.f.col = "white") 

#for pupil paper: with pupil colors
pirateplot(`Mean Rating in Block` ~ `Block Type`, 
           data = meanViewRatings, 
           cex.names = 1.5, 
           cex.axis = 1.5, 
           cex.lab = 1.5, 
           inf.method = "se", 
           inf.f.col = c("purple","light pink","light blue","light green", "dark green"), 
           inf.b.col = "black",
           bean.b.col = "white", 
           bean.f.col = "white") 

# Read out descriptives: mean and SD of ratings in different conditions / block types
mean(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Negative View"])
sd(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Negative View"])

mean(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Negative Regulate"])
sd(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Negative Regulate"])

mean(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Neutral View"])
sd(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Neutral View"])

mean(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Positive Regulate"])
sd(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Positive Regulate"])

mean(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Positive View"])
sd(meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Positive View"])




### FIT the behavioral model 
  
# Test differences in ratings between conditions
fit.eq1 = brm(scan_rating ~ type 
                            + (1 + type|subID), 
                            data = escf1.emo.df)
add_waic(fit.eq1)
fit.eq1$R2 = bayes_R2(fit.eq1)

# Test whether the IAPS pictures worked at all: 

# Note that the neutral view rating after releveling takes the value of 0 
# (to serve as baseline for comparisons between factors in the model)

# Test: negative view rating is smaller than neutral rating 
# (Estimate should be negative)
hypothesis(fit.eq1, "type1 < 0")
# Test: neutral view rating smaller than positive view rating 
# (Estimate should be negative)
hypothesis(fit.eq1, "0 < type5")

# See whether rappraisal worked:
# Test: negative regulate rating bigger than negative view rating
# (Estimate should be positive)
hypothesis(fit.eq1, "type2 > type1")
# Test: positive regulate rating smaller than positive view rating
# (Estimate should be negative)
hypothesis(fit.eq1, "type4 < type5")

# Test posterior probability (PP) directly: extract the samples for each 
# regressor, count how often comparison of regressors is greater or less 
# than zero (depending on hypotheses above) and put in variable called psamps 
# -> evaluates the two against each other

s.pSamps=posterior_samples(fit.eq1, "^b")
ppRatingXCond=c("negView<neutView"=mean(s.pSamps$`b_type1` < 0), 
                "neutView<posView"=mean(0 < s.pSamps$`b_type5`), 
                "negReg>negView"=mean(s.pSamps$`b_type2`>s.pSamps$`b_type1`), 
                "posReg<posView"=mean(s.pSamps$`b_type4`<s.pSamps$`b_type5`))



##### 27.05.20 - In response to reviewer questions:

# In response to reviewer question: R3 Q4
# check for alternative definition of regulation success per block that controls for potential spillover effect
meanRegSuccessNeg = (meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Negative Regulate"]) - 
  (meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Negative View"])

meanRegSuccessPos = (meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Positive View"]) -
  (meanViewRatings$`Mean Rating in Block`[meanViewRatings$`Block Type`=="Positive Regulate"])

meanRegSuccessBlock = (meanRegSuccessNeg + meanRegSuccessPos) / 2

meanTWregSuccessNeg = aggregate (regSuccessNeg ~ subID, data = escf1.emo.df, FUN = "mean")
meanTWregSuccessPos = aggregate (regSuccessPos ~ subID, data = escf1.emo.df, FUN = "mean")
meanTWregSuccess = (meanTWregSuccessNeg$regSuccessNeg + meanTWregSuccessPos$regSuccessPos) / 2

# correlation of regSuccess measures trial-wise and per block across all emotion participants:
cor.test(meanRegSuccessBlock, meanTWregSuccess, alternative = "two.sided", method = "pearson")

# In response to reviewer question: R1 Q3: 
# Are participants better at regulating positive compared to negative stimulus content?
# collapse positive and negative regulation success in one variable
escf1.emo.df$regSuccess = escf1.emo.df$regSuccessPos
escf1.emo.df$regSuccess[escf1.emo.df$type==2] = escf1.emo.df$regSuccessNeg[escf1.emo.df$type==2]

# make binary type variable for positive and negative regulation trials
escf1.emo.df$typePosNeg[escf1.emo.df$type==2] = 1 #negative trial
escf1.emo.df$typePosNeg[escf1.emo.df$type==4] = 2 #positive trial
escf1.emo.df$typePosNeg = as.factor(escf1.emo.df$typePosNeg)

fit.regSuccessByTrialType = brm(regSuccess ~ typePosNeg + (1 + typePosNeg |subID), 
                                 data = escf1.emo.df)
add_waic(fit.regSuccessByTrialType)
fit.regSuccessByTrialType$R2 = bayes_R2(fit.regSuccessByTrialType)


# Plot Supplementary Figure 1:
meanRegSuccessPosNeg = aggregate(regSuccess ~ subID + typePosNeg, data = escf1.emo.df, FUN = "mean")
meanRegSuccessPosNeg = as.data.frame(meanRegSuccessPosNeg)
names(meanRegSuccessPosNeg) = c("subID", "Reappraisal Type", "Reappraisal Success")
levels(meanRegSuccessPosNeg$`Reappraisal Type`) = c("Negative","Positive")

pirateplot(`Reappraisal Success` ~ `Reappraisal Type`, 
           data = meanRegSuccessPosNeg, 
           cex.names = 1.5, 
           cex.axis = 1.5, 
           cex.lab = 1.5, 
           inf.method = "se", 
           inf.f.col = "gray", 
           inf.b.col = "black",
           bean.b.col = "white", 
           bean.f.col = "white") 

##### 




### PART 2: Food choice experiment

# Read in and prepare dataframe 

# This dataframe only contains participants included in the dietary 
# self-control analyses. Note that here we leave in participant 1036 
# who was excluded from fMRI analyses due to excess movement, 
# but we can still analyze behavior.

escf1.df=read.table("Maier_Hare_dataframe_behavior_food_choice.txt", 
                    header=T)
#1026 never accepted any food during challenges at all -> no interpretable behavior
escf1.df=escf1.df[escf1.df$subID != 1026,] 

# Exclude missed trials 
escf1.df=escf1.df[!is.na(escf1.df$yesToOnscreen),] 

# Factorize variables
escf1.df$male = as.factor(escf1.df$male)
escf1.df$scanOrder = as.factor(escf1.df$scanOrder)

# Standardize and mean-center taste and health ratings between participants
escf1.df$tasteOnscreenZ = scale(escf1.df$tasteOnscreen) 
escf1.df$healthOnscreenZ = scale(escf1.df$healthOnscreen) 

#  convert hunger ratings to percent of maximal hunger 
# (hunger rating scale ranges from 0 to maximally 426 points/pixels)
escf1.df$vasH = (escf1.df$vasH / 4.26)

# read out mean + SD for hunger rating across food choice participants
HungerPerSub = aggregate(escf1.df$vasH, 
                       by = list(escf1.df$subID), 
                       FUN="mean")
mean(HungerPerSub$x) # 67.37089
sd(HungerPerSub$x) # 17.52199


# Standardize and mean-center BMI, hunger ratings and restrained eating 
# scores between participants

escf1.df$vasHZ = scale(escf1.df$vasH) 
escf1.df$bmiZ = scale(escf1.df$bmi) 
escf1.df$fevrZ = scale(escf1.df$fevr) 

# Make "yes/no to onscreen" response a (0,1) binary variable
escf1.df$yesToOnscreen[escf1.df$yesToOnscreen==2]=0

# Take log of rtFood for reaction time model
escf1.df$logRTfood = log(escf1.df$rtFood)

# Calculate overall emotion regulation success score 
# (=mean of positive and negative emotion regulation success)
escf1.df$overallERsuccess = (escf1.df$posERsuccess + escf1.df$negERsuccess)/2

# Label self-control success and failures: 
escf1.df$scSuccess = as.numeric(escf1.df$cMTLH==0 | escf1.df$cLTMH == 1) 
escf1.df$scFail = as.numeric(escf1.df$cMTLH==1 | escf1.df$cLTMH == 0)

# Upper and lower end of the neutral zone on the rating scale:
neutraltop = 213 + 21
neutralbottom = 213 - 21

# Create index for trial type
# LTLH = untasty-unhealthy
type1_idx = (escf1.df$healthOnscreen < neutralbottom & escf1.df$tasteOnscreen < neutralbottom) 
# LTMH = healty-untasty 
type2_idx = (escf1.df$healthOnscreen > neutraltop & escf1.df$tasteOnscreen < neutralbottom) 
# MTLH = tasty-unhealthy
type3_idx = (escf1.df$healthOnscreen < neutralbottom & escf1.df$tasteOnscreen > neutraltop) 
# MTMH = tasty-healthy
type4_idx = (escf1.df$healthOnscreen > neutraltop & escf1.df$tasteOnscreen > neutraltop)

# Bin items (as presented, so excluding neutral zone on the rating scale)
# Vector for high-taste high-health items 
escf1.df$HTHH = as.numeric(type4_idx)
# Vector for high-taste low-health items
escf1.df$HTLH = as.numeric(type3_idx)
# Vector for low-taste high-health items
escf1.df$LTHH = as.numeric(type2_idx)
# Vector for low-taste low-health items
escf1.df$LTLH = as.numeric(type1_idx)
# Vector for items neutral on either taste or health
escf1.df$neutH = as.numeric(escf1.df$healthOnscreen>=neutralbottom & escf1.df$healthOnscreen<=neutraltop)
escf1.df$neutT = as.numeric(escf1.df$tasteOnscreen>=neutralbottom & escf1.df$tasteOnscreen<=neutraltop)

# Summarize the no challenge bin: Members of HTHH, LTLH, 
# neutH (neutral on health), or neutT (neutral on taste)
escf1.df$nochall = as.numeric(escf1.df$neutT==1 | escf1.df$neutH==1 | escf1.df$HTHH | escf1.df$LTLH)

# Label challenge  trials
escf1.df$ischall = as.numeric(!escf1.df$nochall)


# Separate palatable-unhealthy, healthy-unpalatable and no-challenge trials
# for model 4 (reaction times)
# 1 is MTLH, tasty-unhealthy
escf1.df$type[type3_idx] = 1 
# 2 is LTMH, healthy-untasty
escf1.df$type[type2_idx] = 2 
# 0 for no challenge
escf1.df$type[is.na(escf1.df$type)] = 0
escf1.df$type = as.factor(escf1.df$type)

# Separate healthy-unpalatable from all other trials
# for model 3 (self-control success model)
# LTMH binary variable
escf1.df$type2[type2_idx] = 1 
# 0 for all other trials
escf1.df$type2[is.na(escf1.df$type2)] = 0 
escf1.df$type2 = as.factor(escf1.df$type2) 

# Separate all four trial types for plotting responses and reaction times
# 1 is LTLH
escf1.df$type4[type1_idx] = 1 
# 2 is LTMH
escf1.df$type4[type2_idx] = 2 
# 3 is MTLH
escf1.df$type4[type3_idx] = 3 
# 4 is MTMH
escf1.df$type4[type4_idx] = 4 
escf1.df$type4 = as.factor(escf1.df$type4)

# (|tr| + |hr|) - stakes
escf1.df$stakes = abs(escf1.df$tasteOnscreen) + abs(escf1.df$healthOnscreen) 

# (| weighted tr + weighted hr| * -1)  - subjective difficulty
# weight for subjective difficulty comes from regression per participant
x= unlist(by(escf1.df, escf1.df$subID, function(x) coef(lm(yesToOnscreen ~ tasteOnscreen + healthOnscreen, data = x))))
taste.coef=x[seq(2, 114, by =3)] #taste is 2nd coef., nsubs * 3 coefs, by = 3 coefs
health.coef=x[seq(3, 114, by =3)] #health is 3rd coef.
subs = unique(escf1.df$subID)
for (i in 1:length(subs)) {
  escf1.df$sdiffic[escf1.df$subID==subs[i]] = abs(escf1.df$tasteOnscreen[escf1.df$subID==subs[i]] * taste.coef[i] 
                                                  + escf1.df$healthOnscreen[escf1.df$subID==subs[i]] * health.coef[i])
}



# Subset to dataframe with just challenge trials 
# for self-control success analyses
escf1.chall.df = escf1.df[!is.na(escf1.df$scSuccess) | !is.na(escf1.df$scFail),] 

# In the subset challenge dataframe, replace NA in scSuccess / scFail 
# variable with zeros to make variables binary
escf1.chall.df$scSuccess[is.na(escf1.chall.df$scSuccess)] = 0
escf1.chall.df$scFail[is.na(escf1.chall.df$scFail)] = 0


### DESCRIPTIVES 

# COUNT the number of presented cases for the four trial types

# Sum over participants: 
# Mean number of high-taste, low-health items presented
HTLHperSub = aggregate(escf1.df$HTLH, 
                       by = list(escf1.df$subID), 
                       FUN="sum")
mean(HTLHperSub$x)
sd(HTLHperSub$x)

min(HTLHperSub$x)
median(HTLHperSub$x)
max(HTLHperSub$x)


# Mean number of low-taste, high-health items presented
LTHHperSub = aggregate(escf1.df$LTHH, 
                       by = list(escf1.df$subID), 
                       FUN="sum")
mean(LTHHperSub$x)
sd(LTHHperSub$x)

min(LTHHperSub$x)
median(LTHHperSub$x)
max(LTHHperSub$x)

# Mean number of non-challenging trials presented
NoChallperSub = aggregate(escf1.df$nochall, 
                          by = list(escf1.df$subID), 
                          FUN="sum")
mean(NoChallperSub$x)
sd(NoChallperSub$x)


# SELF-CONTROL SUCCESS RATE and failures by challenge type

# Aggregate mean self-control scores per participant
meanSCperSub = aggregate(escf1.chall.df$scSuccess, 
                        by = list(escf1.chall.df$subID), 
                        FUN="mean")
mean(meanSCperSub$x)
sd(meanSCperSub$x)

# Split successes per challenge type (say no to high-taste, low-health food  
# OR say yes to low-taste, high-health food)
escf1.chall.df$noHTLH = as.numeric(escf1.chall.df$HTLH==1 & escf1.chall.df$yesToOnscreen==0)
escf1.chall.df$yesLTHH = as.numeric(escf1.chall.df$LTHH==1 & escf1.chall.df$yesToOnscreen==1)

# Sum up the number of self-control successes by challenge type
sumSCnoHTLH =  aggregate(escf1.chall.df$noHTLH, 
                         by = list(escf1.chall.df$subID), 
                         FUN="sum")
sumSCyesLTHH =  aggregate(escf1.chall.df$yesLTHH, 
                          by = list(escf1.chall.df$subID), 
                          FUN="sum")

# Calculate proportion of trials out of all self-control challenges 
# of this type in which the participant was successful
percSCnoHTLH = sumSCnoHTLH$x / HTLHperSub$x
percSCyesLTHH = sumSCyesLTHH$x / LTHHperSub$x

# Proportion of self-control success by refusing high-taste/low-health items
mean(percSCnoHTLH[!is.na(percSCnoHTLH)])
sd(percSCnoHTLH[!is.na(percSCnoHTLH)])

median(percSCnoHTLH[!is.na(percSCnoHTLH)])
mad(percSCnoHTLH[!is.na(percSCnoHTLH)])

# Proportion of self-control success by accepting low-taste/high-health items
mean(percSCyesLTHH[!is.na(percSCyesLTHH)])
sd(percSCyesLTHH[!is.na(percSCyesLTHH)])

median(percSCyesLTHH[!is.na(percSCyesLTHH)])
mad(percSCyesLTHH[!is.na(percSCyesLTHH)])


# PARTICIPANT DEMOGRAPHICS (by self-reported gender)
allBMI = aggregate(bmi ~ subID 
                         + male, 
                         data = escf1.df, 
                         FUN="mean")
maleBMImean = mean(allBMI$bmi[allBMI$male==1])
maleBMIsd = sd(allBMI$bmi[allBMI$male==1])
femaleBMImean = mean(allBMI$bmi[allBMI$male==0]) 
femaleBMIsd = sd(allBMI$bmi[allBMI$male==0])

allAGE = aggregate(age ~ subID 
                         + male, 
                         data = escf1.df, 
                         FUN="mean")
maleAGEmean = mean(allAGE$age[allAGE$male==1])
maleAGEsd = sd(allAGE$age[allAGE$male==1])
femaleAGEmean = mean(allAGE$age[allAGE$male==0])
femaleAGEsd = sd(allAGE$age[allAGE$male==0])


### CREATE Figure 3

# Part 3a
# Plot out proportion of yes in 4 trial categories 
meanYESs = aggregate(yesToOnscreen ~ subID 
                                     + type4, 
                                     data = escf1.df, 
                                     FUN = "mean")
levels(meanYESs$type4) <- c("Unhealthy-Unpalatable", 
                            "Healthy-Unpalatable", 
                            "Palatable-Unhealthy", 
                            "Palatable-Healthy")


# Plot shows mean number of "yes" answers by participant and type
pirateplot(yesToOnscreen ~  type4, 
                            data = meanYESs, 
                            inf.method = "se", 
                            cex.names = 1.5, 
                            cex.axis = 1.5, 
                            cex.lab = 1.5,
                            inf.f.col = "gray", 
                            inf.b.col = "black",
                            bean.b.col = "white", 
                            bean.f.col = "white",
                            xlab = "Choice Type",
                            ylab = "Proportion of Yes Responses")


# Part 3b 

# Aggregate reaction times (RTs)
meanRTs = aggregate(rtFood ~ subID 
                             + yesToOnscreen 
                             + type4, 
                             data = escf1.df, 
                             FUN = "mean")
as.data.frame(meanRTs)
names(meanRTs) <- c("subID", 
                    "Choice", 
                    "Type", 
                    "rtFood")
meanRTs$Choice = as.factor(meanRTs$Choice)
levels(meanRTs$Choice) <- c("No", "Yes")
levels(meanRTs$Type) <- c("Unhealthy-Unpalatable", 
                          "Healthy-Unpalatable", 
                          "Palatable-Unhealthy", 
                          "Palatable-Healthy")

# Plot mean RTs by answer and type
pirateplot(rtFood ~ `Choice` 
                    + `Type`, 
                    data = meanRTs, 
                    inf.method = "se",
                    cex.names = 1.5, 
                    cex.axis = 1.5, 
                    cex.lab = 1.5,
                    inf.f.col = "gray", 
                    inf.b.col = "black",
                    bean.b.col = "white", 
                    bean.f.col = "white",
                    xlab = "Choice Type",
                    ylab = "Reaction Time in Seconds")  


### MEAN REACTION TIMES

# Read out reaction times

# Mean RT for non-challenging trials
mean(escf1.df$rtFood[escf1.df$nochall==1]) 
sd(escf1.df$rtFood[escf1.df$nochall==1])/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying no to MTMH
mean(escf1.df$rtFood[escf1.df$cMTMH==0], na.rm=TRUE)
sd(escf1.df$rtFood[escf1.df$cMTMH==0],na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying no to LTLH
mean(escf1.df$rtFood[escf1.df$cLTLH==0], na.rm=TRUE) 
sd(escf1.df$rtFood[escf1.df$cLTLH==0],na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying yes to MTMH
mean(escf1.df$rtFood[escf1.df$cMTMH==1], na.rm=TRUE) 
sd(escf1.df$rtFood[escf1.df$cMTMH==1],na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying yes to LTLH
mean(escf1.df$rtFood[escf1.df$cLTLH==1], na.rm=TRUE) 
sd(escf1.df$rtFood[escf1.df$cLTLH==1],na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying no to MTLH
mean(escf1.df$rtFood[escf1.df$cMTLH==0], na.rm=TRUE) 
sd(escf1.df$rtFood[escf1.df$cMTLH==0],na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying yes to MTLH
mean(escf1.df$rtFood[escf1.df$cMTLH==1], na.rm=TRUE) 
sd(escf1.df$rtFood[escf1.df$cMTLH==1],na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying no to LTMH
mean(escf1.df$rtFood[escf1.df$cLTMH==0], na.rm=TRUE) 
sd(escf1.df$rtFood[escf1.df$cLTMH==0], na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))

# Mean RT for saying yes to LTMH
mean(escf1.df$rtFood[escf1.df$cLTMH==1], na.rm=TRUE)
sd(escf1.df$rtFood[escf1.df$cLTMH==1], na.rm=TRUE)/sqrt(length(unique(escf1.df$subID)))



### FIT Behavioral models for the food choice task

# basic food choice model
fit.eq2 = brm(yesToOnscreen ~ (tasteOnscreenZ + healthOnscreenZ) * scanOrder
                              + male 
                              + vasHZ 
                              + bmiZ 
                              + fevrZ 
                              + (1 + (tasteOnscreenZ + healthOnscreenZ) * scanOrder | subID),
                              data = escf1.df, 
                              family = "bernoulli",
                              warmup = 3000, 
                              iter = 6000, 
                              chains = 4,
                              control = list(adapt_delta = 0.95))
# add information criteria (WAIC and Bayesian R-squared)
add_waic(fit.eq2) 
fit.eq2$R2 = bayes_R2(fit.eq2)


# self-control success model (accounting for challenge type)
fit.eq3 = brm(scSuccess ~ type2 * ( tasteOnscreenZ + healthOnscreenZ) + 
                          (1 + type2 * ( tasteOnscreenZ + healthOnscreenZ)  | subID), 
                          data = escf1.chall.df, 
                          family = "bernoulli")
# add information criteria (WAIC and Bayesian R-squared)
add_waic(fit.eq3)
fit.eq3$R2 = bayes_R2(fit.eq3)


# # (previous) reaction time model
# fit.eq4 = brm(logRTfood ~ yesToOnscreen * type 
#                           + (1 + yesToOnscreen * type |subID), 
#                           data = escf1.df)
# # add information criteria (WAIC and Bayesian R-squared)
# add_waic(fit.eq4)
# fit.eq4$R2 = bayes_R2(fit.eq4)


##### 27.05.20 - In response to reviewer questions:

# R1 last question: Updated our RT model and presented the model below in the paper.

# The RT model should control for strength of preference, which in our case is equal to difficulty. 
# We can construct a model with both difficulty and stakes as regressors: 
# log(rt) = β0+ β1 Yes + β2 Type + β3 (Yes x Type) + (|tr| + |hr|) + (|weighted tr + weighted hr| * -1)  + ε
# So we are adding in stakes + diffic on each trial to RT model but keep Type as defined in orig RT model (fit.eq4)

fit.eq4.sdiffic.stakes.slopes = brm(logRTfood ~ yesToOnscreen * type 
                            + scale(stakes) + scale(sdiffic) 
                            + (1 + (yesToOnscreen * type) + scale(stakes) + scale(sdiffic) |subID), 
                            data = escf1.df)

fit.eq4.sdiffic.stakes.slopes$R2 = bayes_R2(fit.eq4.sdiffic.stakes.slopes)



# R2 Q3: Test fo correlations of hunger level with...
# ...dietary self-control level 
cor.test(HungerPerSub$x, meanSCperSub$x)
# ... and emotion control level
# (to compare between ER and FC task, remove participants excluded from food choice analyses from ER dataset)
meanTWregSuccessNegFC = meanTWregSuccessNeg[meanTWregSuccessNeg$subID != 1013 & 
                                              meanTWregSuccessNeg$subID != 1017 &
                                              meanTWregSuccessNeg$subID != 1026 & 
                                              meanTWregSuccessNeg$subID != 1033,]
meanTWregSuccessPosFC = meanTWregSuccessPos[meanTWregSuccessPos$subID != 1013 & 
                                              meanTWregSuccessPos$subID != 1017 &
                                              meanTWregSuccessPos$subID != 1026 & 
                                              meanTWregSuccessPos$subID != 1033,]
meanTWregSuccessFC = (meanTWregSuccessNegFC$regSuccessNeg + meanTWregSuccessPosFC$regSuccessPos) / 2
# (vice versa, remove participants excluded from emotion analyses from FC dataset)
HungerPerSubER = HungerPerSub[HungerPerSub$Group.1 != 1012 &
                                HungerPerSub$Group.1 != 1015 &
                                HungerPerSub$Group.1 != 1018 &
                                HungerPerSub$Group.1 != 1033 &
                                HungerPerSub$Group.1 != 1036 &
                                HungerPerSub$Group.1 != 1038 &
                                HungerPerSub$Group.1 != 1060,]
cor.test(HungerPerSubER$x, meanTWregSuccessFC)



# R3 Q12: If one used a DDM approach, can se wee that participants in the current study
# have a starting-point bias towards saying "No"?

# We take DDM fits from Maier, Raja Beharelle et al. (2020) Nature Human Behaviour 
# and test against the "natural" and "health" cue fits from Hare, Malmaud & Rangel (2011)) 
# (in the following named study "IAC" - instructed attention cue). 
# We compare the results to this study because participants chose also between 
# eating the proffered foods or eating nothing, as in the ESC study.

# DDM fits from the current ESC study (fits are across all trials)
ESCfits = read.csv("ESC_fits_alltrials.csv", header = TRUE)
ESCfits = ESCfits[,2:7]

#read out means and SDs for Table S1
round(colMeans(ESCfits),2)
round(colSds(as.matrix(ESCfits)),2)

# tDDM fits from the IAC study (condition 1 = health, condition 3 = natural cue)
load("IAC_fits_Cue1_all_3000.RData")
IACfits_health = fitsF
rm(fitsF)
IACfits_health = IACfits_health[,1:6]

load("IAC_fits_Cue3_all_3000.RData")
IACfits_natural = fitsF
rm(fitsF)
IACfits_natural = IACfits_natural[,1:6]


twosamplett_bias_ESC_vs_IAC_health = BESTmcmc(ESCfits$bias, IACfits_health$bias)
plot(twosamplett_bias_ESC_vs_IAC_health)

# mean difference bias term according to BEST: ESC bias term is more pronounced
# towards no (more negative) than IAC bias term in the health condition
mean(twosamplett_bias_ESC_vs_IAC_health$mu1 - twosamplett_bias_ESC_vs_IAC_health$mu2)
#posterior probability that ESC bias term is smaller than IAC bias term
post_prob = mean(twosamplett_bias_ESC_vs_IAC_health$mu1 - twosamplett_bias_ESC_vs_IAC_health$mu2 < 0)
#hdi
hdi(twosamplett_bias_ESC_vs_IAC_health$mu1 - twosamplett_bias_ESC_vs_IAC_health$mu2)

# compare to the natural condition in Hare et al. (2011)
twosamplett_bias_ESC_vs_IAC_natural = BESTmcmc(ESCfits$bias, IACfits_natural$bias)
plot(twosamplett_bias_ESC_vs_IAC_natural)

# mean difference bias term according to BEST: ESC bias term is more pronounced
# towards no (more negative) than IAC bias term in the natural condition
mean(twosamplett_bias_ESC_vs_IAC_natural$mu1 - twosamplett_bias_ESC_vs_IAC_natural$mu2)
#posterior probability that ESC bias term is smaller than IAC bias term
post_prob = mean(twosamplett_bias_ESC_vs_IAC_natural$mu1 - twosamplett_bias_ESC_vs_IAC_natural$mu2 < 0)
#hdi
hdi(twosamplett_bias_ESC_vs_IAC_natural$mu1 - twosamplett_bias_ESC_vs_IAC_natural$mu2)

#####



### BEHAVIORAL LINK between dietary and emotion self-regulation success

# subset food choice dataframe to participants with proper emotion regulation data
# note that we have n = 32 datasets here because 1062 was just excluded from the
# emotion fMRI analyses due to excess motion, but behavior is interpretable
# -> hence the exclusion vector excludeERbev (behavior only)
escf1.er.df = escf1.df[escf1.df$excludeERbev==0,]

# Bayesian correlations between emotion-regulation successand dietary self-control success

source("bayesRankCorrelation.R")

# create aggregate dataframe
overallERsuccess = aggregate(overallERsuccess ~ subID, 
                                                data = escf1.er.df, 
                                                FUN = "mean")
posERsuccess = aggregate(posERsuccess ~ subID, 
                                        data = escf1.er.df, 
                                        FUN = "mean")
negERsuccess = aggregate(negERsuccess ~ subID, 
                                        data = escf1.er.df, 
                                        FUN = "mean")
SCsuccess = aggregate(percSCsuccess ~ subID, 
                                      data = escf1.er.df, 
                                      FUN = "mean")

# overall emotion regulation (ER) success
y=cbind(overallERsuccess$overallERsuccess, SCsuccess$percSCsuccess)
bayesCORRoverall=bayesRankCor(y)

rhoOverall = mean(bayesCORRoverall)
hdiOverall = hdi(bayesCORRoverall)
# Test posterior probability greater than 0 
# because we expect it to be positive (one-sided test)
ppOverall = mean(bayesCORRoverall > 0) 

# positive ER success
y=cbind(posERsuccess$posERsuccess, SCsuccess$percSCsuccess)
bayesCORRpos=bayesRankCor(y)

rhoPos = mean(bayesCORRpos)
hdiPos = hdi(bayesCORRpos)
# Test posterior probability greater than 0 
# because we expect it to be positive (one-sided test)
ppPos = mean(bayesCORRpos > 0) 

# negative ER success 
y=cbind(negERsuccess$negERsuccess, SCsuccess$percSCsuccess)
bayesCORRneg=bayesRankCor(y)

rhoNeg = mean(bayesCORRneg)
hdiNeg = hdi(bayesCORRneg)
# Test posterior probability greater than 0 
# because we expect it to be positive (one-sided test)
ppNeg = mean(bayesCORRneg > 0) 


##### 27.05.20 - In response to reviewer questions:

# R3 Q6: "The dietary choice success measure is calculated simply as the percentage 
# of “self-control success” trials, without taking into account their difficulty. 
# It seems to me that a more comparable measure in the food choice task might be to do 
# something like measuring how difficult each trial that resulted in self-control success is, 
# or how much tastiness is “sacrificed” by restraint, with trials without self-control success 
# counting as 0 or “no restraint”. It would be nice to know whether, if more comparable 
# self-control measures were used, one might observe more of a correlation across the behaviors."

# Option 1 - to parallel the emotion success measure more closely:
# 0) sign taste ratings (subtract the middle of the scale)
# 1) take the taste ratings on successful self-control trials as "price of the sacrifice": 
# (multiply by -1: sacrifice made)
# 2) taste rating on the failed self-control trials (multiply by 1: "taste obtained / sacrifice not made") 
# 3) take the mean across 1) and 2) to calculate the overall score

escf1.chall.eronly.df$tasteOnscreenSigned = escf1.chall.eronly.df$tasteOnscreen - 213

escf1.chall.eronly.df$sacrifice[escf1.chall.eronly.df$yesToOnscreen==1] = -1
escf1.chall.eronly.df$sacrifice[escf1.chall.eronly.df$yesToOnscreen==0] = 1

# flip the sign on the failure trials according to 1) and 2)
escf1.chall.eronly.df$tasteSacrificed = escf1.chall.eronly.df$tasteOnscreenSigned * escf1.chall.eronly.df$sacrifice

# take the mean subjective difficulty across all challenge trials to calculate 3):
meanSdiffScore = aggregate (tasteSacrificed ~ subID, data = escf1.chall.eronly.df, FUN = mean)

# correlate meanSdiffScore with emotion regulation success score
cor.test(meanSdiffScore$tasteSacrificed,overallERsuccess$overallERsuccess)
# separate by valence
cor.test(meanSdiffScore$tasteSacrificed,posERsuccess$posERsuccess)
cor.test(meanSdiffScore$tasteSacrificed,negERsuccess$negERsuccess)


# Option 2 - make emotion regulation measure binary (success/fail) to make it more comparable to the
# overall dietary self-control success rate

# binary success ratings for negative regulate trials:
# success if inscan rating during negative regulate is higher than post rating
escf1.emo.df$binsuccess[escf1.emo.df$type==2 & (escf1.emo.df$scan_rating > escf1.emo.df$post_rating)] = 1 
# no success if inscan rating during negative regulate is equal to or smaller than post rating
escf1.emo.df$binsuccess[escf1.emo.df$type==2 & (escf1.emo.df$scan_rating <= escf1.emo.df$post_rating)] = 0

# binary success ratings for positive regulate trials:
# success if inscan rating during positive regulate is smaller than post rating
escf1.emo.df$binsuccess[escf1.emo.df$type==4 & (escf1.emo.df$scan_rating < escf1.emo.df$post_rating)] = 1
# no success if inscan rating during positive regulate is equal to or higher than post rating
escf1.emo.df$binsuccess[escf1.emo.df$type==4 & (escf1.emo.df$scan_rating >= escf1.emo.df$post_rating)] = 0

# aggregate binary success measures and take the mean across pos / neg conditions
# subsetting only pulls out the regulation trials (type 4 = pos reg, type 2 = neg reg)
# taking the mean across the binsuccess logical gives us the proportion of successful regulation trials
binSuccessERpos = aggregate(binsuccess[escf1.emo.df$type==4] ~ subID[escf1.emo.df$type==4],
                            data = escf1.emo.df, FUN = mean) 

binSuccessERneg = aggregate(binsuccess[escf1.emo.df$type==2] ~ subID[escf1.emo.df$type==2],
                            data = escf1.emo.df, FUN = mean)

# reduce to participants included in the food analyses
binSuccessERpos = binSuccessERpos[binSuccessERpos$`subID[escf1.emo.df$type == 4]` != 1013 & 
                                    binSuccessERpos$`subID[escf1.emo.df$type == 4]` != 1017 &
                                    binSuccessERpos$`subID[escf1.emo.df$type == 4]` != 1026 &
                                    binSuccessERpos$`subID[escf1.emo.df$type == 4]` != 1033 ,]
binSuccessERneg = binSuccessERneg[binSuccessERneg$`subID[escf1.emo.df$type == 2]` != 1013 & 
                                    binSuccessERneg$`subID[escf1.emo.df$type == 2]` != 1017 &
                                    binSuccessERneg$`subID[escf1.emo.df$type == 2]` != 1026 &
                                    binSuccessERneg$`subID[escf1.emo.df$type == 2]` != 1033 ,]

binSuccessERoverall = (binSuccessERpos$`binsuccess[escf1.emo.df$type == 4]` + 
                         binSuccessERneg$`binsuccess[escf1.emo.df$type == 2]`)/2  

# calculate percent self-control success per participant 
# for participants who are included in both ER and FC analyses
escf1.chall.eronly.df = escf1.chall.df[escf1.chall.df$excludeERbev == 0,]
overallSCsuccessRate = aggregate(scSuccess ~subID, data = escf1.chall.eronly.df, FUN = "mean")

# run correlation between overall emotion regulation and self-control success rates

source("bayesRankCorrelation.R")

# run correlations separatly by positive and negative domains
y=cbind(overallSCsuccessRate$scSuccess,binSuccessERoverall)

bayesCORRbinER_FC=bayesRankCor(y)
rhoCORRbinER_FC = mean(bayesCORRbinER_FC)
hdiCORRbinER_FC = hdi(bayesCORRbinER_FC)
# Test posterior probability greater than 0 
# because we expect it to be positive (one-sided test)
ppCORRbinER_FC = mean(bayesCORRbinER_FC > 0) 




# R3 Q7: "One might speculate that the ability to down-regulate positive affect 
# might be more relevant only for the tasty but unhealthy foods, 
# while the ability to down-regulate negative affect might be more relevant 
# for the unappealing but healthy foods. If self-control success for these two distinct 
# dietary choice types is used, and correlated with the distinct types of emotion regulation, 
# does a stronger correlation emerge?"

# succeed by refusing tasty-unhealthy
type1_success_rate = aggregate(scSuccess[escf1.chall.eronly.df$type==1] 
                               ~ subID[escf1.chall.eronly.df$type==1], 
                               data = escf1.chall.eronly.df, FUN = "mean")
type1_success_per_sub = type1_success_rate$`scSuccess[escf1.chall.eronly.df$type == 1]`
# succeed by accepting healthy-untasty (note that here we miss 1007 and 1016)
type2_success_rate = aggregate(scSuccess[escf1.chall.eronly.df$type==2] 
                               ~ subID[escf1.chall.eronly.df$type==2], 
                               data = escf1.chall.eronly.df, FUN = "mean")
type2_success_per_sub = type2_success_rate$`scSuccess[escf1.chall.eronly.df$type == 2]`

# emotion regulation success means
er_pos_reg_success_fc_corr = aggregate(posERsuccess ~ subID,
                                       data = escf1.chall.eronly.df, FUN = "mean")
er_pos_reg_success_fc_per_sub = er_pos_reg_success_fc_corr$posERsuccess
er_neg_reg_success_fc_corr = aggregate(negERsuccess[escf1.chall.eronly.df$subID != 1007 & escf1.chall.eronly.df$subID != 1016] 
                                       ~ subID[escf1.chall.eronly.df$subID != 1007 & escf1.chall.eronly.df$subID != 1016],
                                       data = escf1.chall.eronly.df, FUN = "mean")
er_neg_reg_success_fc_per_sub = er_neg_reg_success_fc_corr$`negERsuccess[escf1.chall.eronly.df$subID != 1007 & escf1.chall.eronly.df$subID != 1016]`

source("bayesRankCorrelation.R")

# run correlations separatly by positive and negative domains
y=cbind(type1_success_per_sub,er_pos_reg_success_fc_per_sub)

bayesCORRposERfcType1=bayesRankCor(y)
rhoCORRposERfcType1 = mean(bayesCORRposERfcType1)
hdiCORRposERfcType1 = hdi(bayesCORRposERfcType1)
# Test posterior probability greater than 0 
# because we expect it to be positive (one-sided test)
ppCORRposERfcType1 = mean(bayesCORRposERfcType1 > 0) 


y=cbind(type2_success_per_sub,er_neg_reg_success_fc_per_sub)

bayesCORRnegERfcType2=bayesRankCor(y)
rhoCORRnegERfcType2 = mean(bayesCORRnegERfcType2)
hdiCORRnegERfcType2 = hdi(bayesCORRnegERfcType2)
# Test posterior probability greater than 0 
# because we expect it to be positive (one-sided test)
ppCORRnegERfcType2 = mean(bayesCORRnegERfcType2 > 0) 

#####