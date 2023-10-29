#   Validation of SAS SCORECI macro against R function
#    ratesci::scoreci  (for stratified MN and SCAS intervals)
#

install.packages('ratesci')
library(ratesci)

SASCIs <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval3.csv")
nsamp <- max(SASCIs$sample)
head(SASCIs)

RCIs <- array(NA, dim = c(nsamp, 2))
dimnames(RCIs)[[2]] <- c("lclR", "uclR")
for (i in 1:nsamp) {
  onesample <- SASCIs[SASCIs$sample == i, ]
  try(RCIs[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
                            x2 = onesample[,"e0"], n2 = onesample[,"n0"],
                            level = onesample[,"CONFLEV"][1], stratified = TRUE,
                            skew = FALSE, precis=10)$estimates[,c(1,3)])
  #('try' function allows the loop to continue in the event of an error)
}

allCIs <- cbind(SASCIs[SASCIs$stratum == 1,],RCIs)
head(allCIs)
attach(allCIs)

# Summarise differences between SAS and R
lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

summary(lcld)
summary(ucld)

hist(lcld)
hist(ucld)

max(abs(lcld)) #maximum discrepancy between SAS macro & diffscoreci version
max(abs(ucld)) #maximum discrepancy between SAS macro & diffscoreci version

detach(allCIs)


# Now do the same with skewness corrected intervals
SASCIs2<-read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval3skew.csv")
nsamp <- max(SASCIs2$sample)
head(SASCIs2)

RCIs2 <- array(NA, dim = c(nsamp, 2))
dimnames(RCIs2)[[2]] <- c("lclR", "uclR")
for (i in 1:nsamp) {
  onesample <- SASCIs2[SASCIs2$sample == i, ]
  try(RCIs2[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
                                    x2 = onesample[,"e0"], n2 = onesample[,"n0"],
                                    level = onesample[,"CONFLEV"][1], stratified = TRUE,
                                    skew = TRUE, precis=10)$estimates[,c(1,3)])
  #('try' function allows the loop to continue in the event of an error)
}

allCIs2 <- cbind(SASCIs2[SASCIs2$stratum == 1,], RCIs2)
head(allCIs2)

attach(allCIs2)

lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

summary(lcld)
summary(ucld)

hist(lcld)
hist(ucld)

max(abs(lcld)) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld)) #maximum discrepancy between SAS macro & scoreci version

detach(allCIs2)


