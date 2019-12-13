# Quick analysis of digestibility data from Nick Santoro's group looking at
# glucose and pentose release under standardized conditions. 
# A single tube of ground material was sent for each individual tissue source
# (i.e. plant), each tube received 3 technical replicates and appears to have 
# been measured 4 times. Using a mixed effects model to try to account for
# that. 
#
# David Duncan
# 2015-09-01

library(lme4)
library(lsmeans)
library(reshape)

dig.dat <- read.table("../data/zip_digestibility.txt", header=TRUE,sep="\t")
dig.dat$clone <- factor(dig.dat$clone, levels=c("Control","5","6","7","8",""))


# Ugh, there are two levels of nesting going on here:
# - 3 technical replicaes
# - 4 measurements, maybe of the same technical rep?
# Second level is, unfortnately, unstacked, and so needs to be reshaped

glu.mlt <- melt(dig.dat[,8:12], id="glu.out")
glu.full <- merge(dig.dat[,c(3:5,8)],glu.mlt[,c(1,3)])
names(glu.full)[5] <- "glu"

glu.lmr <- lmer(glu ~ clone*tissue + (1|plant/glu.out), 
                data=subset(glu.full, clone!=""))

lsmeans(glu.lmr, trt.vs.ctrl ~ clone|tissue)
lsmeans(glu.lmr, pairwise ~ tissue|clone)
