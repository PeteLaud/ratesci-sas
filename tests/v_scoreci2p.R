#   Validation of SAS SCORECI macro against R functions
#   ratesci::scoreci  (for MN & SCAS interval for Poisson RD)
#

install.packages('ratesci')
library(ratesci)

SASCIs<-read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2p.csv")
nsamp<-dim(SASCIs)[[1]]
head(SASCIs)

RCIs <- ratesci::scoreci(x1=SASCIs[,"e1"], n1=SASCIs[,"n1"],
                          x2=SASCIs[,"e0"], n2=SASCIs[,"n0"],
                          level=SASCIs[,"CONFLEV"],
                          skew=FALSE, precis=10, distrib='poi')$estimates[,c(1,3)]

allCIs <- cbind(SASCIs,RCIs)
head(allCIs)
attach(allCIs)

#identify any instances of missing intervals from ratesci
allCIs[is.na(Lower),]
allCIs[is.na(Upper),]
lcld <- L_BOUND - Lower
ucld <- U_BOUND - Upper

summary(lcld)
summary(ucld)


allCIs[is.na(lcld),]

# hist(lcld)
# hist(ucld)

max(abs(lcld[!is.na(Lower)])) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld[!is.na(Upper)])) #maximum discrepancy between SAS macro & scoreci version

detach(allCIs)


#Now do the same with skewness corrected intervals
SASCIs2 <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2pskew.csv")
nsamp <- dim(SASCIs2)[[1]]
head(SASCIs2)

RCIs2 <- ratesci::scoreci(x1=SASCIs2[,"e1"], n1=SASCIs2[,"n1"],
                          x2=SASCIs2[,"e0"], n2=SASCIs2[,"n0"],
                          level=SASCIs2[,"CONFLEV"],
                          skew=TRUE, precis=10, distrib='poi')$estimates[,c(1,3)]

allCIs2 <- cbind(SASCIs2, RCIs2)
head(allCIs2)

attach(allCIs2)

#identify any instances of missing intervals from ratesci
allCIs2[is.na(Lower),]
allCIs2[is.na(Upper),]
lcld <- L_BOUND - Lower
ucld <- U_BOUND - Upper

summary(lcld)
summary(ucld)

allCIs2[is.na(lcld),]
allCIs2[ucld>0.00001,]

# hist(lcld)
# hist(ucld)

max(abs(lcld[!is.na(Lower)])) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld[!is.na(Upper)])) #maximum discrepancy between SAS macro & scoreci version

detach(allCIs2)


