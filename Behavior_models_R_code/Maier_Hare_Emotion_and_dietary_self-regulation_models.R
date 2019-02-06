# Behavioral models for Maier & Hare, "Emotion and dietary self-regulation"

# Version 1.0, Preprint, 18 December 2018

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



### PART 2: Food choice experiment

# Read in and prepare dataframe 

# This dataframe only contains participants included in the dietary 
# self-control analyses. Note that here we leave in participant 1036 
# who was excluded from fMRI analyses due to excess movement, 
# but we can still analyze behavior.

escf1.df=read.table("Maier_Hare_dataframe_behavior_food_choice.txt", 
                    header=T)

# Exclude participant 1026 who never chose to eat during the self-control
# challenges -> no interpretable data 
escf1.df=escf1.df[escf1.df$subID != 1026,] 

# Exclude missed trials 
escf1.df=escf1.df[!is.na(escf1.df$yesToOnscreen),] 

# Factorize variables
escf1.df$male = as.factor(escf1.df$male)
escf1.df$scanOrder = as.factor(escf1.df$scanOrder)

# Standardize and mean-center taste and health ratings between participants
escf1.df$tasteOnscreenZ = scale(escf1.df$tasteOnscreen) 
escf1.df$healthOnscreenZ = scale(escf1.df$healthOnscreen) 

# Standardize and mean-center BMI, hunger ratings and restrained eating 
# scores between participants
# and convert hunger ratings to percent of maximal hunger 
# (hunger rating scale ranges from 0 to maximally 426 points/pixels)
escf1.df$vasHZ = scale(escf1.df$vasH / 4.26) 
escf1.df$bmiZ = scale(escf1.df$bmi) 
escf1.df$fevrZ = scale(escf1.df$fevr) 

# Make "yes/no to onscreen" respones a binary variable
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


# reaction time model
fit.eq4 = brm(logRTfood ~ yesToOnscreen * type 
                          + (1 + yesToOnscreen * type |subID), 
                          data = escf1.df)
# add information criteria (WAIC and Bayesian R-squared)
add_waic(fit.eq4)
fit.eq4$R2 = bayes_R2(fit.eq4)


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