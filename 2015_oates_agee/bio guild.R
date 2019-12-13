#####
# Preliminary analysis of split-depth study at AARS
# Soils sampled in 2013 as a follow-up of Chao et al., 2013 (Agricultural Sciences)
# The objective is to explore changes to microbial community, as defined
# by PLFA profiles, that have occurred at AARS BSCE since establishment
# in 2008. 
# Project under fearless leadership of LG Oates.
# David Duncan
# 14.02.19

library(nlme)
library(lsmeans)
library(AICcmodavg)
library(vegan)
setwd("C:/Users/AndyJ/Dropbox/Articles/AARS split depth study 2013/Analysis")
mass.data <- read.table("../Data/split bio guild.txt", sep = "\t", header = TRUE)
mass.data$plot <- factor(mass.data$plot)
mass.data$rep <- factor(mass.data$rep)
mass.data <- mass.data[mass.data$SampleNo != 98,]

## Biomass ####
# Units of umol lipid/g dry soil

bm.lme0 <- lme(biomass ~ depth, random = ~1|rep/plot/subplot, data = mass.data)
bm.lme0.1 <- update (bm.lme0, weights = varIdent(form = ~1|depth))
anova(bm.lme0, bm.lme0.1)
# Variances are significantly different between depths, and depths have
# significantly different biomasses. Recommend splitting depths for further
# analysis, since direct comparison across depths is unlikely to be
# fruitful

# Just the top 10 cm then
AIC(bm.10.lme0 <- lme(biomass ~ system*subplot, random = ~1|rep/plot/subplot, data = mass.data[mass.data$depth == "0-10",]))
AIC(bm.10.lme0.1 <- update(bm.10.lme0, weights = varIdent(form = ~ 1|system), control = lmeControl(msMaxIter=200)))
AIC(bm.10.lme0.2 <- update(bm.10.lme0, weights = varIdent(form = ~ 1|system*subplot),control = lmeControl(msMaxIter=200)))
anova(bm.10.lme0, bm.10.lme0.1, bm.10.lme0.2)
# Big difference in variances, seems like we need to separate systems and
# subplots, but should also see how sensitive results are to this
anova(bm.10.lme0.2)
lsmeans(bm.10.lme0.2, pairwise ~ system*subplot)

# 10-25 cm now
AIC(bm.25.lme0 <- lme(biomass ~ system*subplot, random = ~1|rep/plot/subplot, data = mass.data[mass.data$depth == "10-25",]))
AIC(bm.25.lme0.1 <- update(bm.25.lme0, weights = varIdent(form = ~ 1|system), control = lmeControl(msMaxIter=200)))
AIC(bm.25.lme0.2 <- update(bm.25.lme0, weights = varIdent(form = ~ 1|system*subplot),control = lmeControl(msMaxIter=200)))
anova(bm.25.lme0, bm.25.lme0.1, bm.25.lme0.2)
# Once agan, separate systems appear to be justified
anova(bm.25.lme0.2)
anova(bm.25.lme0.2)
lsmeans(bm.25.lme0.2, pairwise ~ system*subplot)

### AMF ####
AICc(amf.10.lme0 <- lme(AMF ~ system*subplot, random = ~1|rep/plot/subplot, data = mass.data[mass.data$depth == "0-10",]))
AICc(amf.10.lme0.1 <- update(amf.10.lme0, weights = varIdent(form = ~ 1|system), control = lmeControl(msMaxIter=200)))
AICc(amf.10.lme0.2 <- update(amf.10.lme0, weights = varIdent(form = ~ 1|system*subplot),control = lmeControl(msMaxIter=200)))
anova(amf.10.lme0, amf.10.lme0.1, amf.10.lme0.2)
# Once agan, separate systems appear to be justified
anova(amf.10.lme0.2)
lsmeans(amf.10.lme0.1, pairwise ~ subplot|system, data = mass.data)

AIC(amf.25.lme0 <- lme(AMF ~ system*subplot, random = ~1|rep/plot/subplot, data = mass.data[mass.data$depth == "10-25",]))
AIC(amf.25.lme0.1 <- update(amf.25.lme0, weights = varIdent(form = ~ 1|system), control = lmeControl(msMaxIter=200)))
AIC(amf.25.lme0.2 <- update(amf.25.lme0, weights = varIdent(form = ~ 1|system*subplot),control = lmeControl(msMaxIter=200)))
anova(amf.25.lme0, amf.25.lme0.1, amf.25.lme0.2)
# Once agan, separate systems appear to be justified
anova(amf.25.lme0.2)
lsmeans(amf.25.lme0.2, pairwise ~ system*subplot)

#### GM Neg. ####
AICc(gmn.10.lme0 <- lme(Gm.neg~system*subplot, random = ~1|rep/plot/subplot, data=mass.data[mass.data$depth=="0-10",]))
AICc(gmn.10.lme0.1 <- update(gmn.10.lme0, weights = varIdent(form = ~1|system), , control = lmeControl(msMaxIter=200)))
AICc(gmn.10.lme0.2 <- update(gmn.10.lme0, weights = varIdent(form = ~1|system*subplot), , control = lmeControl(msMaxIter=200)))
# separate variance by system but not by subplot
lsmeans(gmn.10.lme0.2, pairwise ~ subplot|system, data = mass.data)
