#   Validation of SAS SCORECI macro against R functions
#   PropCIs::diffscoreci  (for MN interval)
#   and ratesci::scoreci  (for SCAS interval)
#

install.packages('PropCIs')
library(PropCIs)
install.packages('ratesci')
library(ratesci)

SASCIs<-read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2.csv")
nsamp<-dim(SASCIs)[[1]]
head(SASCIs)

RCIs<-array(NA,dim=c(nsamp,2))
dimnames(RCIs)[[2]]<-c("lclR","uclR")
for (i in 1:nsamp) {
	try(RCIs[i,] <- unlist(diffscoreci(x1=SASCIs[i,"e1"],n1=SASCIs[i,"n1"],
	                                   x2=SASCIs[i,"e0"],n2=SASCIs[i,"n0"],
	                          conf.level=SASCIs[i,"CONFLEV"])))
	#('try' function allows the loop to continue in the event of an error)
}

allCIs <- cbind(SASCIs,RCIs)
head(allCIs)
attach(allCIs)

#identify any instances of missing intervals from diffscoreci
allCIs[is.na(lclR),]
allCIs[is.na(uclR),]
lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

summary(lcld)
summary(ucld)

allCIs[is.na(lcld),]

# hist(lcld)
# hist(ucld)

max(abs(lcld[!is.na(lclR)])) #maximum discrepancy between SAS macro & diffscoreci version
max(abs(ucld[!is.na(uclR)])) #maximum discrepancy between SAS macro & diffscoreci version
#NB precision in the diffscoreci code is +/-(1e-07) so discrepancies of that order are to be expected

detach(allCIs)


#Now do the same with skewness corrected intervals
SASCIs2 <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2skew.csv")
nsamp <- dim(SASCIs2)[[1]]
head(SASCIs2)

RCIs2 <- ratesci::scoreci(x1=SASCIs2[,"e1"], n1=SASCIs2[,"n1"],
                          x2=SASCIs2[,"e0"], n2=SASCIs2[,"n0"],
                          level=SASCIs2[,"CONFLEV"],
                          skew=TRUE, precis=10)$estimates[,c(1,3)]

allCIs2 <- cbind(SASCIs2, RCIs2)
head(allCIs2)

attach(allCIs2)

#identify any instances of missing intervals from diffscoreci
allCIs2[is.na(Lower),]
allCIs2[is.na(Upper),]
lcld <- L_BOUND-Lower
ucld <- U_BOUND-Upper

summary(lcld)
summary(ucld)

allCIs2[is.na(lcld),]
allCIs2[ucld>0.00001,]

# hist(lcld)
# hist(ucld)

max(abs(lcld[!is.na(Lower)])) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld[!is.na(Upper)])) #maximum discrepancy between SAS macro & scoreci version

detach(allCIs2)


