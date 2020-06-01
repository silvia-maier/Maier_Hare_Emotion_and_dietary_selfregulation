#setwd("~/Dropbox/TSC_shared/delayDDM_DE/ESC1fMRI")

library(Rcpp)
#library(lme4)
library(plyr)
#install.packages('DEoptim')
library(DEoptim)
library(pracma)

# load the c++ file containing the functions to simulate the delayed DDM
sourceCpp("2ddm_r_cpp_2.cpp")

# use Ncores - 1
# on Silvia's iMac, just runs with 1 Thread, otherwise crashes R
RcppParallel::setThreadOptions(numThreads = 1)

tab = read.table("dataframe_escf1_basic_ddm.txt", header=T)
dataBeh = as.data.frame(tab)

## remove bad trials
idx = is.na(dataBeh$yesToOnscreen)
dataBeh = dataBeh[idx==0,] #idx is false when it is not NA 

## bin the taste and health ratings 
# in this dataset, we used a rating scale with a "neutral zone" in the region around zero (+- 5% of the whole scale length)
# rating scale was from 0 to 426 points, bins -5 to -1 go from 0 to (213-21.3), bins 1 to 5 go from (213+21.3) to 426 points)
# items in the neutral zone are coded as bin "0"

# bin width = 191.7 / 5 = 38.34
edges <- c(0,38.34,76.68,115.02,153.36,191.7, 234.3,272.64,310.98,349.32,387.66,426)


tmp = histc(as.numeric(dataBeh$health),edges)
tmp$bin = tmp$bin-6
tmp$bin[tmp$bin==6] = 5 #assign the maximum to the top bin
dataBeh$hbin = tmp$bin #assign to dataframe variable

tmp2 = histc(as.numeric(dataBeh$taste),edges)
tmp2$bin = tmp2$bin-6
tmp2$bin[tmp2$bin==6] = 5 #assign the maximum to the top bin
dataBeh$tbin = tmp2$bin #assign to dataframe variable

# save the taste and health differences in a column - BINNED VERSION
dataBeh$td = dataBeh$tbin
dataBeh$hd = dataBeh$hbin


###---------------------------------------------
# fit the model 

dataBeh2 = dataBeh

#m1 = glmer(cCURR ~ td + hd + (1+td+hd|subID), data=dataBeh2, family=binomial(link = "logit"))


# number of subjects
ns = length(unique(dataBeh2$subID))

# reindex the subjects in continuous order
auxVec = dataBeh2$subID
dataBeh2$idxP = mapvalues(auxVec, unique(dataBeh2$subID), to=1:ns)


# assign negative RTs to left item
# SUM, 29.11.2016: NOTE that this is why the coefficients in the SSM output turned out negative:
# I assigned negative RTs to left item, but then also calculated taste and health differences as L(eft) minus R(ight) item.
# For the purpose of easier reading, in the manuscript the signs for health and taste drift were flipped again!

# SUM, 04.01.2017: for better correspondence to other studies, negative RTs are now assigned to choosing the right item!
dataBeh2$RTddm = dataBeh2$rt
dataBeh2$RTddm[dataBeh2$yesToOnscreen==0] = dataBeh2$RTddm[dataBeh2$yesToOnscreen==0] * -1
ntrials = length(dataBeh2$rt)


# bin the prob. density space of RTs
xpos = seq(-5,5,length.out=1024) #one trial was 3 seconds + on average 2 seconds ITI (1-3 seconds)
dt = xpos[2] - xpos[1]
dataBeh2$RTddm_pos = 0
for (i in 1:ntrials) {
	dataBeh2$RTddm_pos[i] = which.min(abs(xpos - dataBeh2$RTddm[i]))
}

# save the results of the model fit
fits=matrix(0,ns,9)


# subject's loop
for (i in 1:ns) {

# take the trials for the current subject_i
idx = which(dataBeh2$idxP==i)
dataBeh3 = dataBeh2[idx,]

# save each combination of hd,td
data1 = ddply(dataBeh3, .(td, hd), summarize, acc= mean(yesToOnscreen))
vd = data1$td
hd = data1$hd


ll_ddm2 <- function(x) {

	d_v = x[1] # drift rate for taste
	d_h = x[2] # drift rate for health
	thres = x[3] # DDM threshold
	nDT = x[4] # non decision time
	tHin = x[5] # time of health-taste in
	bias = x[6] # bias towards one of the alternatives (e.g. left or right, CURR or ALT)

	probs = NULL
	for (tr in 1:length(vd)) {
	
		# last parameter of "ddm2_parallel" is the number of simulations 
		rts <- ddm2_parallel(d_v,d_h,thres,nDT,tHin,bias,vd[tr],hd[tr],1500) # simulate the DDM for a given x
		rts = rts[rts!=0]
    	
    	# generate the densitty function of the simulated DDM
    	xdens = density(rts, from=-5, to=5, n=1024, bw=0.11)
    	
    	# get the index of the current trial in the PDF space
    	idx = which(dataBeh3$td==vd[tr] & dataBeh3$hd==hd[tr])
    	
    	# save the probability value of the simulated DDM for a given trial(i.e. position in the x-axis of the RT sapace)
    	probs = c(probs, dt*xdens$y[dataBeh3$RTddm_pos[idx]])
	}

	probs[probs==0] = 1e-100
	return (-sum(log(probs)))

} #end of creating function ll_ddm2

# set the bounds of the parameters
lower <- c(-2,-2,0.6,0.01,-1,-1)
upper <- c(2,2,3,1,1,1)

# run the optimizer
fit_s = DEoptim(ll_ddm2, lower, upper, DEoptim.control(itermax = 150))


fits[i,1:6] = fit_s$optim$bestmem # readout the fitted params
fits[i,7] = fit_s$optim$bestval # -LL
fits[i,8] = 2*fit_s$optim$bestval + length(lower)*log(length(vd)) # AIC
fits[i,9] = 2*fit_s$optim$bestval + 2*length(lower) #BIC

} # Close subjects loop


save(fits, file=paste("fits_ESCfMRI_binned_allTrials.RData", sep=""))

#to write out fits to a csv: write.csv(fits, 'ESC_fits_alltrials.csv')


###ll_ddm2(fit_s$optim$bestmem)
###pars = c(0.72155599, 1.4614106, 1.1893310, 0.7428796, 0.0013423594, -0.61846340)
###ll_ddm2(pars)

