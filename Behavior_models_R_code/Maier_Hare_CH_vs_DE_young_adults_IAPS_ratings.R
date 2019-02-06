# Compare Swiss (CH) and German (DE) valence ratings (from Grühn & Scheibe sample) 
# for stimuli of the International Affective Picture System (IAPS)

# Version 1.0, Preprint, 18 December 2018
# Author: Silvia Maier, silvia.maier@econ.uzh.ch


# load CH ratings dataframe
escf1.df=read.table("Maier_Hare_dataframe_behavior_emotion_ratings.txt", 
                    header=T)
escf1.df = as.data.frame(escf1.df)

# restrict to participants included in emotion reappraisal task
# note we have 36 participants here because 1032 was only 
# excluded from fMRI due to excess motion
escf1.df = escf1.df[escf1.df$exclude_er==0,]

# variable keys:
# iapsID      IAPS identification number
# scan_rating valence rating of young adults 
#             (during an emotion reappraisal task)
# type        trial type in task: 
#             1 = negative View, 
#             2 = negative Regulate, 
#             3 = neutral View, 
#             4 = positive Regulate, 
#             5 = positive View 
# post_rating valence rating of young adults in CH 
#             (view rating for the reappraised pictures)

# for consistency with DE ratings, re-name image 6570.1 to 6570 
# (it is the same image)
escf1.df$iapsID[escf1.df$iapsID==6570.1]=6570

# load DE ratings dataframe from their original publication
gs.df = read.table("Gruhn-Scheibe_2008_PictureData.txt", header=T)
gs.df = as.data.frame(gs.df)

# variable keys in Gruehn & Scheibe 2008 dataset:
# IAPS	  IAPS identification number
# val_y	  Valence rating of Young adults
# aro_y	  Arousal rating of Young adults



# calculate mean + sd (standard deviation), 
# median + mad (median absolute deviation)
# for CH ratings for each iaps image (valence ratings)
ch.mean.val.neg = aggregate(scan_rating ~ iapsID, 
                            data = escf1.df[escf1.df$type==1,], 
                            FUN = "mean")
ch.mean.val.pos = aggregate(scan_rating ~ iapsID, 
                            data = escf1.df[escf1.df$type==5,], 
                            FUN = "mean")
ch.mean.val.neu = aggregate(scan_rating ~ iapsID, 
                            data = escf1.df[escf1.df$type==3,], 
                            FUN = "mean")

ch.sd.val.neg = aggregate(scan_rating ~ iapsID, 
                          data = escf1.df[escf1.df$type==1,], 
                          FUN = "sd")
ch.sd.val.pos = aggregate(scan_rating ~ iapsID, 
                          data = escf1.df[escf1.df$type==5,], 
                          FUN = "sd")
ch.sd.val.neu = aggregate(scan_rating ~ iapsID, 
                          data = escf1.df[escf1.df$type==3,], 
                          FUN = "sd")

ch.median.val.neg = aggregate(scan_rating ~ iapsID, 
                              data = escf1.df[escf1.df$type==1,], 
                              FUN = "median")
ch.median.val.pos = aggregate(scan_rating ~ iapsID, 
                              data = escf1.df[escf1.df$type==5,], 
                              FUN = "median")
ch.median.val.neu = aggregate(scan_rating ~ iapsID, 
                              data = escf1.df[escf1.df$type==3,], 
                              FUN = "median")

ch.mad.val.neg = aggregate(scan_rating ~ iapsID, 
                           data = escf1.df[escf1.df$type==1,], 
                           FUN = "mad")
ch.mad.val.pos = aggregate(scan_rating ~ iapsID, 
                           data = escf1.df[escf1.df$type==5,], 
                           FUN = "mad")
ch.mad.val.neu = aggregate(scan_rating ~ iapsID, 
                           data = escf1.df[escf1.df$type==3,], 
                           FUN = "mad")


# subset DE dataframe to the stimuli used 
# in our positive and negative picture sets

# DE ratings for our negative sets: note that our original ID is 6570.1, 
# theirs 6570 -  both refer to the same image! 
# corrected above in our dataframe
gs.df.val.neg = gs.df[gs.df$IAPS== 1300 | gs.df$IAPS== 1525 | gs.df$IAPS== 2055.1 |
                        gs.df$IAPS== 2095 | gs.df$IAPS== 2352.2 | gs.df$IAPS== 2800 |
                        gs.df$IAPS== 3015 | gs.df$IAPS== 3181 | gs.df$IAPS== 3301 |
                        gs.df$IAPS== 3550 | gs.df$IAPS== 6020 | gs.df$IAPS== 6370 |
                        gs.df$IAPS== 6415 | gs.df$IAPS== 6540 | gs.df$IAPS== 6838 |
                        gs.df$IAPS== 9040 | gs.df$IAPS== 9180 | gs.df$IAPS== 9265 |
                        gs.df$IAPS== 9435 | gs.df$IAPS== 9520 | gs.df$IAPS== 2683 |
                        gs.df$IAPS== 2981 | gs.df$IAPS== 3051 | gs.df$IAPS== 3250 |
                        gs.df$IAPS== 3530 | gs.df$IAPS== 6212 | gs.df$IAPS== 6312 |
                        gs.df$IAPS== 6570 | gs.df$IAPS== 9140 | gs.df$IAPS== 9181 |
                        gs.df$IAPS== 9250 | gs.df$IAPS== 9252 | gs.df$IAPS== 9253 |
                        gs.df$IAPS== 9300 | gs.df$IAPS== 9430 | gs.df$IAPS== 9561 |
                        gs.df$IAPS== 9570 | gs.df$IAPS== 9571 | gs.df$IAPS== 9635.1 |
                        gs.df$IAPS== 9800,]

# DE ratings for our positive sets
gs.df.val.pos = gs.df[gs.df$IAPS== 1460 | gs.df$IAPS== 1710 | gs.df$IAPS== 1721 |
                        gs.df$IAPS== 1750 | gs.df$IAPS== 1810 | gs.df$IAPS== 1920 |
                        gs.df$IAPS== 2050 | gs.df$IAPS== 2080 | gs.df$IAPS== 2091 |
                        gs.df$IAPS== 2224 | gs.df$IAPS== 2260 | gs.df$IAPS== 2311 |
                        gs.df$IAPS== 2351 | gs.df$IAPS== 2375.2 | gs.df$IAPS== 2550 |
                        gs.df$IAPS== 5600 | gs.df$IAPS== 5626 | gs.df$IAPS== 5831 |
                        gs.df$IAPS== 5890 | gs.df$IAPS== 8190 | gs.df$IAPS== 1440 |
                        gs.df$IAPS== 1463 | gs.df$IAPS== 1720 | gs.df$IAPS== 1731 |
                        gs.df$IAPS== 2058 | gs.df$IAPS== 2071 | gs.df$IAPS== 2150 |
                        gs.df$IAPS== 2170 | gs.df$IAPS== 2303 | gs.df$IAPS== 2345 |
                        gs.df$IAPS== 2620 | gs.df$IAPS== 2655 | gs.df$IAPS== 2660 |
                        gs.df$IAPS== 5390 | gs.df$IAPS== 5594 | gs.df$IAPS== 5628 |
                        gs.df$IAPS== 5830 | gs.df$IAPS== 7580 | gs.df$IAPS== 8461 |
                        gs.df$IAPS== 8497,]

# DE ratings for our neutral set
gs.df.val.neu = gs.df[gs.df$IAPS== 2020 | gs.df$IAPS== 2200 | gs.df$IAPS== 2214 |
                        gs.df$IAPS== 2357 | gs.df$IAPS== 2480 | gs.df$IAPS== 2493 |
                        gs.df$IAPS== 2570 | gs.df$IAPS== 2880 | gs.df$IAPS== 2890 |
                        gs.df$IAPS== 7030 | gs.df$IAPS== 7036 | gs.df$IAPS== 7150 |
                        gs.df$IAPS== 7161 | gs.df$IAPS== 7170 | gs.df$IAPS== 7186 |
                        gs.df$IAPS== 7224 | gs.df$IAPS== 7590 | gs.df$IAPS== 7705 |
                        gs.df$IAPS== 7830 | gs.df$IAPS== 8030,]


# calculate correlation between stimulus sets - valence ratings CH and DE

allCHratings = c(ch.mean.val.neg$scan_rating, 
                 ch.mean.val.pos$scan_rating, 
                 ch.mean.val.neu$scan_rating)

allDEratings = c(gs.df.val.neg$val_y, 
                 gs.df.val.pos$val_y, 
                 gs.df.val.neu$val_y)

# calculate overall correlation
cor.test(allCHratings,allDEratings)
plot(allCHratings,allDEratings)

# calculate deviance
dev = abs(allCHratings-allDEratings)
mean(dev) 
# mean deviance of ratings is 0.47 points (on a 9-point Likert-type scale)
dev_signed = allCHratings-allDEratings 
# signed deviance is positive: CH ratings tend to be slightly higher 
# than DE ratings