# Revisiting AARS PLFA study
#
# Chao Liang took PLFA samples from soils in auxilliary plots in late summer 2008
# to get a baseline measure of soil microbial community composition and test whether
# the different land use history between blocks 1-3 and 4-5 were detectable in the
# community. We are revisiting that study with year 5 data to see how the treatments
# have impacted the microbial community, with a particular interest in observing the 
# effects of N fertilization. The entire analysis for the paper is included in this
# script.
#
# Version 15.08.26
# David Duncan

### Libraries ###
library(vegan)
library(ggplot2)
library(lme4)
library(lsmeans)
library(grid)

cbbPalette <- c("#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#E69F00", "#999999")
#c("gray92", "gray55", "gray30")
shapes <- c(21,22,23,24,25)


### Data import ###
plfa.data <- read.table("../Data/AARS split depth data for R.txt", header = TRUE, sep = "\t")

# Remove samples that don't have measurements for both depths
temp.drop <- c("Liang BF-1", "Liang BF-2", "Liang BF-3", "Liang BF-5")
plfa.data <- subset(plfa.data, !(Sample.ID %in% temp.drop))

# Isolate plfa concentrations, replace blanks with 0
plfa.markers <- plfa.data[,15:83]
plfa.markers[is.na(plfa.markers)] <- 0

# arcsine-square root transform PLFA markers to rein in outliers somewhat
plfa.markers <- asin(sqrt(plfa.markers))/(pi/2)
# remove lipids 20 carbons and longer (likely to be plant contaminants)
plfa.markers <- plfa.markers[,c(1:51,60:69)]

plfa.anc <- plfa.data[,1:14]
plfa.anc$OldRep <- plfa.anc$Rep
levels(plfa.anc$OldRep) <- c("W", "W", "W", "E", "E")
plfa.anc$Treatment <- factor(plfa.anc$Treatment, 
                             levels = levels(plfa.anc$Treatment)[c(1:7,10:11,8:9)])
plfa.anc$base <- factor(plfa.anc$Treatment=="Baseline")

### 0-10 cm ####
## Splitting data by depth because we really aren't interested in comparisons
## between depths when doing this. 
surf.markers <- plfa.markers[plfa.anc$Depth == "0 to 10",]
surf.anc <- plfa.anc[plfa.anc$Depth == "0 to 10",]

### Surface ordination ####
surf.mds1 <- metaMDS(surf.markers, autotransform = FALSE, distance = "euclidian", trymax = 200)
surf.mds2 <- metaMDS(surf.markers, autotransform = FALSE, distance = "euclidian", trymax = 200, previous.best = surf.mds1)
surf.mds3 <- metaMDS(surf.markers, autotransform = FALSE, distance = "euclidian", trymax = 200, previous.best = surf.mds2)
surf.mds4 <- metaMDS(surf.markers, autotransform = FALSE, distance = "euclidian", trymax = 200, previous.best = surf.mds3)
surf.mds5 <- metaMDS(surf.markers, autotransform = FALSE, distance = "euclidian", trymax = 200, previous.best = surf.mds4)
surf.mds6 <- metaMDS(surf.markers, autotransform = FALSE, distance = "euclidian", trymax = 200, previous.best = surf.mds5)
# final stress: 0.05921554 
stressplot(surf.mds6)
# Non-metric fit R2 = 0.996

surf.mds.scores <- scores(surf.mds6)
surf.mds.scores <- cbind(surf.mds.scores, surf.anc, 
                         trt=with(surf.anc, paste(System, Nfert, Stover)))

surf.mds.points <- merge(aggregate(cbind(avg.x=surf.mds.scores$NMDS1, avg.y=surf.mds.scores$NMDS2)~System*Nfert*Stover, surf.mds.scores, mean),
                         aggregate(cbind(sd.x=surf.mds.scores$NMDS1, sd.y=surf.mds.scores$NMDS2)~System*Nfert, surf.mds.scores, sd))
surf.mds.species <- data.frame(wascores(surf.mds6$points, surf.mds.scores[,11:16]))


surf.plot <- ggplot(surf.mds.points, aes(x=avg.x, y = avg.y)) +
  scale_size(range=c(0.5,6)) + scale_color_manual(values=cbbPalette)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill="white"),
        legend.position="none") +
  geom_hline(aes(size=0.5, yintercept=0)) +
  geom_vline(aes(size=0.5, xintercept=0)) +
  geom_errorbar(aes(x=avg.x, ymax=avg.y+sd.y, ymin=avg.y-sd.y, size=0.5, width = 0.005)) +
  geom_errorbarh(aes(y=avg.y, xmax=avg.x+sd.x, xmin=avg.x-sd.x, size=0.5, height = 0.003)) +
  geom_point(aes(size = 6, shape = Nfert, color = System)) +
  geom_segment(data=surf.mds.species[4:6,], aes(x=0, xend=MDS1, y=0, yend=MDS2, size=0.75), arrow = arrow(length=unit(0.3, "cm")))

ggsave("../Figures_results andTable S1/surf_mds.pdf", width=17, height=17, units="cm")

write.table(surf.mds.scores, "../Results & figures/surface NMDS.txt", sep = "\t", quote=FALSE, row.names=FALSE)

### Soil surface adonis ###
surf.dist <- vegdist(surf.markers, distance = "euclidian")
(surf.adon.plot <- adonis(surf.dist ~ OldRep, data = surf.anc, permutation = 9999))
(surf.adon.trt <- adonis(surf.dist ~ base + System * (Nfert + Stover), data = surf.anc, permutation = 9999))


### Surface functional groups ####
## Using functional groups to ask several questions:
## 1. Did biomasses change from 2008 to 2013
## 2. What are the differences among 2013 samples (all vs all pairwise)
## 3. Was there a consistent effect of nitrogen fertilization (4 systems only)
## 4. Are there differences between field halves in the baseline and new data?

## Biomass ##
# treatments vs 2008
surf.bio.lmer <- lmer(Biomass ~ Treatment -1 + (1|Rep/Plot), data = surf.anc)
test.lm <- lm(Biomass ~ Treatment -1, data=surf.anc)
surf.bio.lme <- lme(Biomass ~ Treatment, random=~1 | Rep / Plot, data=surf.anc,
                    weights=varIdent(form=~1 | Treatment))
cld(lsmeans(surf.bio.lme, trt.vs.ctrl1 ~ Treatment), sort=FALSE)

# all vs all treatments
surf.bio.lmer2 <- lmer(Biomass ~ Treatment  -1 +(1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
cld(lsmeans(surf.bio.lmer2, pairwise ~ Treatment, adjust=p.adjust.methods[7]))

# N fertilizer effect
surf.bio.lmer3 <- lmer(Biomass ~ Nfert + (1|Rep/Plot), data = surf.anc[!(surf.anc$System %in% c("Baseline", "Corn")),])
lsmeans(surf.bio.lmer3, pairwise ~ Nfert)

# 2013 land use history effect
surf.bio.lmer4 <- lmer(Biomass ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.bio.lmer4, pairwise ~ OldRep)

# 2008 land use history effect
surf.bio.lmer5 <- lmer(Biomass ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment == "Baseline",])
lsmeans(surf.bio.lmer5, pairwise ~ OldRep)

## All fungi ##
# treatments vs 2008
surf.fun.lmer <- lmer(Fungi ~ Treatment + (1|Rep/Plot), data = surf.anc)
surf.fun.lme <- lme(Fungi ~ Treatment, random=~1 | Rep / Plot, data=surf.anc,
                    weights=varIdent(form=~1 | Treatment), method="REML")

lsmeans(surf.fun.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])
lsmeans(surf.fun.lme, trt.vs.ctrl1 ~ Treatment)

# all vs all treatments
surf.fun.lmer2 <- lmer(Fungi ~ Treatment + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.fun.lmer2, pairwise ~ Treatment)

# N fertilizer effect
surf.fun.lmer3 <- lmer(Fungi ~ Nfert + (1|Rep/Plot), data = surf.anc[!(surf.anc$System %in% c("Baseline", "Corn")),])
lsmeans(surf.fun.lmer3, pairwise ~ Nfert)

# 2013 land use history effect
surf.fun.lmer4 <- lmer(Fungi ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.fun.lmer4, pairwise ~ OldRep)

# 2008 land use history effect
surf.fun.lmer5 <- lmer(Fungi ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment == "Baseline",])
lsmeans(surf.fun.lmer5, pairwise ~ OldRep)

## All bacteria ##
surf.bac.lmer <- lmer(Bacteria ~ Treatment + (1|Rep/Plot), data = surf.anc)
lsmeans(surf.bac.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

surf.bac.lmer2 <- lmer(Bacteria ~ Treatment + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.bac.lmer2, pairwise ~ Treatment)

surf.bac.lmer3 <- lmer(Bacteria ~ Nfert + (1|Rep/Plot), data = surf.anc[!(surf.anc$System %in% c("Baseline", "Corn")),])
lsmeans(surf.bac.lmer3, pairwise ~ Nfert)

surf.bac.lmer4 <- lmer(Bacteria ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.bac.lmer4, pairwise ~ OldRep)

surf.bac.lmer5 <- lmer(Bacteria ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment == "Baseline",])
lsmeans(surf.bac.lmer5, pairwise ~ OldRep)

## AMF ##
surf.amf.lmer <- lmer(AMF ~ Treatment + (1|Rep/Plot), data = surf.anc)
lsmeans(surf.amf.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

surf.amf.lmer2 <- lmer(AMF ~ Treatment + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
cld(lsmeans(surf.amf.lmer2, pairwise ~ Treatment, , adjust=p.adjust.methods[7]))

surf.amf.lmer3 <- lmer(AMF ~ Nfert + (1|Rep/Plot), data = surf.anc[!(surf.anc$System %in% c("Baseline", "Corn")),])
lsmeans(surf.amf.lmer3, pairwise ~ Nfert)

surf.amf.lmer4 <- lmer(AMF ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.amf.lmer4, pairwise ~ OldRep)

surf.amf.lmer5 <- lmer(AMF ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment == "Baseline",])
lsmeans(surf.amf.lmer5, pairwise ~ OldRep)

## Gram negative bacteria ##
surf.gmn.lmer <- lmer(GMn ~ Treatment + (1|Rep/Plot), data = surf.anc)
summary(surf.gmn.lmer)
lsmeans(surf.gmn.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

surf.gmn.lmer2 <- lmer(GMn ~ Treatment + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
cld(lsmeans(surf.gmn.lmer2, pairwise ~ Treatment, adjust=p.adjust.methods[7]))

surf.gmn.lmer3 <- lmer(GMn ~ Nfert + (1|Rep/Plot), data = surf.anc[!(surf.anc$System %in% c("Baseline", "Corn")),])
lsmeans(surf.gmn.lmer3, pairwise ~ Nfert)

surf.gmn.lmer4 <- lmer(GMn ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.gmn.lmer4, pairwise ~ OldRep)

surf.gmn.lmer5 <- lmer(GMn ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment == "Baseline",])
lsmeans(surf.gmn.lmer5, pairwise ~ OldRep)

## Gram-positive bacteria ##
surf.gmp.lmer <- lmer(GMp ~ Treatment + (1|Rep/Plot), data = surf.anc)
summary(surf.gmp.lmer)
lsmeans(surf.gmp.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])


surf.gmp.lmer2 <- lmer(GMp ~ Treatment + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.gmp.lmer2, pairwise ~ Treatment)

surf.gmp.lmer3 <- lmer(GMp ~ Nfert + (1|Rep/Plot), data = surf.anc[!(surf.anc$System %in% c("Baseline", "Corn")),])
lsmeans(surf.gmp.lmer3, pairwise ~ Nfert)

surf.gmp.lmer4 <- lmer(GMp ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment != "Baseline",])
lsmeans(surf.gmp.lmer4, pairwise ~ OldRep)

surf.gmp.lmer5 <- lmer(GMp ~ OldRep + (1|Rep/Plot), data = surf.anc[surf.anc$Treatment == "Baseline",])
lsmeans(surf.gmp.lmer5, pairwise ~ OldRep)

### 10-25 cm ####
## Splitting data by depth because we really aren't interested in comparisons
## between depths when doing this. 
sub.markers <- plfa.markers[plfa.anc$Depth == "10 to 25",]
sub.anc <- plfa.anc[plfa.anc$Depth == "10 to 25",]

### Subsurface ordination ####
sub.mds1 <- metaMDS(sub.markers, distance = "euclidian", trymax = 200)
sub.mds2 <- metaMDS(sub.markers, distance = "euclidian", trymax = 200, previous.best = sub.mds1)
sub.mds3 <- metaMDS(sub.markers, distance = "euclidian", trymax = 200, previous.best = sub.mds2)
sub.mds4 <- metaMDS(sub.markers, distance = "euclidian", trymax = 200, previous.best = sub.mds3)
sub.mds5 <- metaMDS(sub.markers, distance = "euclidian", trymax = 200, previous.best = sub.mds4)
sub.mds6 <- metaMDS(sub.markers, distance = "euclidian", trymax = 200, previous.best = sub.mds5)
# stress was 0.0615413 
stressplot(sub.mds6)

sub.mds.scores <- scores(sub.mds6)
sub.mds.scores <- cbind(sub.mds.scores, sub.anc, 
                        trt=with(sub.anc, paste(System, Nfert, Stover)))
sub.mds.points <- merge(aggregate(cbind(avg.x=sub.mds.scores$NMDS1, avg.y=sub.mds.scores$NMDS2)~System*Nfert*Stover, sub.mds.scores, mean),
                         aggregate(cbind(sd.x=sub.mds.scores$NMDS1, sd.y=sub.mds.scores$NMDS2)~System*Nfert, sub.mds.scores, sd))
sub.mds.species <- data.frame(wascores(sub.mds6$points, sub.mds.scores[,11:16]))


sub.plot <- ggplot(sub.mds.points, aes(x=avg.x, y = -avg.y)) +
  scale_size(range=c(0.5,6)) + scale_color_manual(values=cbbPalette)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill="white"),
        legend.position="none") +
  geom_hline(aes(size=0.5, yintercept=0)) +
  geom_vline(aes(size=0.5, xintercept=0)) +
  geom_errorbar(aes(x=avg.x, ymax=-avg.y+sd.y, ymin=-avg.y-sd.y, size=0.5, width = 0.005)) +
  geom_errorbarh(aes(y=-avg.y, xmax=avg.x+sd.x, xmin=avg.x-sd.x, size=0.5, height = 0.003)) +
  geom_point(aes(size = 6, shape = Nfert, color = System)) +
  geom_segment(data=sub.mds.species[4:6,], aes(x=0, xend=MDS1, y=0, yend=-MDS2, size=0.75), arrow = arrow(length=unit(0.3, "cm")))

ggsave("../Figures_results andTable S1/sub_mds.pdf", width=17, height=17, units="cm")


write.table(sub.mds.scores, "../Results & figures/subsurface NMDS.txt", sep = "\t", quote=FALSE, row.names=FALSE)

### Soil subsurface adonis ###
sub.dist <- vegdist(sub.markers, distance = "euclidian")
(sub.adon.plot <- adonis(sub.dist ~ OldRep, data = sub.anc, permutation = 9999))
(sub.adon.trt <- adonis(sub.dist ~ base + System * (Nfert + Stover), data = sub.anc, permutation = 9999))

### Subsurface functional groups ####
## Using functional groups to ask several questions:
## 1. Did biomasses change from 2008 to 2013
## 2. What are the differences among 2013 samples (all vs all pairwise)
## 3. Was there a consistent effect of nitrogen fertilization (4 systems only)
## 4. Are there differences between field halves in the baseline and new data?
## Biomass
sub.bio.lmer <- lmer(Biomass ~ Treatment + (1|Rep/Plot), data = sub.anc)
summary(sub.bio.lmer)
lsmeans(sub.bio.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

sub.bio.lmer2 <- lmer(Biomass ~ Treatment + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
cld(lsmeans(sub.bio.lmer2, pairwise ~ Treatment, adjust=p.adjust.methods[7]))

sub.bio.lmer3 <- lmer(Biomass ~ Nfert + (1|Rep/Plot), data = sub.anc[!(sub.anc$System %in% c("Baseline", "Corn")),])
lsmeans(sub.bio.lmer3, pairwise ~ Nfert)

sub.bio.lmer4 <- lmer(Biomass ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.bio.lmer4, pairwise ~ OldRep)

sub.bio.lmer5 <- lmer(Biomass ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment == "Baseline",])
lsmeans(sub.bio.lmer5, pairwise ~ OldRep)

## Fungi
sub.fun.lmer <- lmer(Fungi ~ Treatment + (1|Rep/Plot), data = sub.anc)
lsmeans(sub.fun.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

sub.fun.lmer2 <- lmer(Fungi ~ Treatment + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.fun.lmer2, pairwise ~ Treatment)

sub.fun.lmer3 <- lmer(Fungi ~ Nfert + (1|Rep/Plot), data = sub.anc[!(sub.anc$System %in% c("Baseline", "Corn")),])
lsmeans(sub.fun.lmer3, pairwise ~ Nfert)

sub.fun.lmer4 <- lmer(Fungi ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.fun.lmer4, pairwise ~ OldRep)

sub.fun.lmer5 <- lmer(Fungi ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment == "Baseline",])
lsmeans(sub.fun.lmer5, pairwise ~ OldRep)

## Bacteria
sub.bac.lmer <- lmer(Bacteria ~ Treatment + (1|Rep/Plot), data = sub.anc)
lsmeans(sub.bac.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

sub.bac.lmer2 <- lmer(Bacteria ~ Treatment + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.bac.lmer2, pairwise ~ Treatment)

sub.bac.lmer3 <- lmer(Bacteria ~ Nfert + (1|Rep/Plot), data = sub.anc[!(sub.anc$System %in% c("Baseline", "Corn")),])
lsmeans(sub.bac.lmer3, pairwise ~ Nfert)

sub.bac.lmer4 <- lmer(Bacteria ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.bac.lmer4, pairwise ~ OldRep)

sub.bac.lmer5 <- lmer(Bacteria ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment == "Baseline",])
lsmeans(sub.bac.lmer5, pairwise ~ OldRep)

## AMF
sub.amf.lmer <- lmer(AMF ~ Treatment + (1|Rep/Plot), data = sub.anc)
lsmeans(sub.amf.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

sub.amf.lmer2 <- lmer(AMF ~ Treatment + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
cld(lsmeans(sub.amf.lmer2, pairwise ~ Treatment, , adjust=p.adjust.methods[7]))

sub.amf.lmer3 <- lmer(AMF ~ Nfert + (1|Rep/Plot), data = sub.anc[!(sub.anc$System %in% c("Baseline", "Corn")),])
lsmeans(sub.amf.lmer3, pairwise ~ Nfert)

sub.amf.lmer4 <- lmer(AMF ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.amf.lmer4, pairwise ~ OldRep)

sub.amf.lmer5 <- lmer(AMF ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment == "Baseline",])
lsmeans(sub.amf.lmer5, pairwise ~ OldRep)

## GMn
sub.gmn.lmer <- lmer(GMn ~ Treatment + (1|Rep/Plot), data = sub.anc)
lsmeans(sub.gmn.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])

sub.gmn.lmer2 <- lmer(GMn ~ Treatment + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
cld(lsmeans(sub.gmn.lmer2, pairwise ~ Treatment, adjust=p.adjust.methods[7]))

sub.gmn.lmer3 <- lmer(GMn ~ Nfert + (1|Rep/Plot), data = sub.anc[!(sub.anc$System %in% c("Baseline", "Corn")),])
lsmeans(sub.gmn.lmer3, pairwise ~ Nfert)

sub.gmn.lmer4 <- lmer(GMn ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.gmn.lmer4, pairwise ~ OldRep)

sub.gmn.lmer5 <- lmer(GMn ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment == "Baseline",])
lsmeans(sub.gmn.lmer5, pairwise ~ OldRep)

## GMp
sub.gmp.lmer <- lmer(GMp ~ Treatment + (1|Rep/Plot), data = sub.anc)
lsmeans(sub.gmp.lmer, trt.vs.ctrl1 ~ Treatment, adjust=p.adjust.methods[7])


sub.gmp.lmer2 <- lmer(GMp ~ Treatment + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.gmp.lmer2, pairwise ~ Treatment)

sub.gmp.lmer3 <- lmer(GMp ~ Nfert + (1|Rep/Plot), data = sub.anc[!(sub.anc$System %in% c("Baseline", "Corn")),])
lsmeans(sub.gmp.lmer3, pairwise ~ Nfert)

sub.gmp.lmer4 <- lmer(GMp ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment != "Baseline",])
lsmeans(sub.gmp.lmer4, pairwise ~ OldRep)

sub.gmp.lmer5 <- lmer(GMp ~ OldRep + (1|Rep/Plot), data = sub.anc[sub.anc$Treatment == "Baseline",])
lsmeans(sub.gmp.lmer5, pairwise ~ OldRep)

#### Surface vs subsurface ####

protest(sub.mds6, surf.mds6, permutations=9999)

#### Stress indicators ####

# Some samples are missing one or more of the stress indicator lipids.
# Rather than drop the data, we'll assume the values are just very small with
# "small" defined as "half the smallest value > 0"

# Convenience function to generate our definition of "small"
min0.5 <- function(x){
  x <- as.vector(t(x))
  x <- x[x>0]
  return(min(x)/2)
}

# Organizing the indicator lipids
gmp_num <- c("X15.0.ISO","X17.0.ISO")
gmp_den <- c("X15.0.ANTEISO","X17.0.ANTEISO")
gmn_num <- c("X17.0.CYCLO")
gmn_den <- "X16.1.w7c"

# Subsurface
sub_str_mark <- sub.markers[,c(gmp_num,gmp_den,gmn_num,gmn_den)]
# Make a dataframe of the "small values 
sub_min <- apply(sub_str_mark,2, min0.5)
sub_min_mat <- data.frame(matrix(sub_min,ncol=ncol(sub_str_mark),
                                 nrow=nrow(sub_str_mark),byrow=TRUE))


sub_str <- sub_str_mark
sub_str[sub_str<=0] <- sub_min_mat[sub_str<=0]

sub.anc$gmn_str <- log(sub_str[,gmn_num]/sub_str[,gmn_den])
sub.anc$gmp_str <- log(rowSums(sub_str[,gmp_num])/
                         rowSums(sub_str[,gmp_den]))


sub.gmn_str <- lme(gmn_str~Treatment, random=~1|Rep/Plot,
                   weights=varIdent(form=~1|Treatment),data=sub.anc)
lsmeans(sub.gmn_str, trt.vs.ctrl1 ~ Treatment)

sub.gmp_str <- lme(gmp_str~Treatment, random=~1|Rep/Plot,
                   weights=varIdent(form=~1|Treatment),data=sub.anc)
lsmeans(sub.gmp_str, trt.vs.ctrl1 ~ Treatment)

# Surface
surf_str_mark <- surf.markers[,c(gmp_num,gmp_den,gmn_num,gmn_den)]
# Make a dataframe of the "small values 
surf_min <- apply(surf_str_mark,2, min0.5)
surf_min_mat <- data.frame(matrix(surf_min,ncol=ncol(surf_str_mark),
                                 nrow=nrow(surf_str_mark),byrow=TRUE))


surf_str <- surf_str_mark
surf_str[surf_str<=0] <- surf_min_mat[surf_str<=0]

surf.anc$gmn_str <- log(surf_str[,gmn_num]/surf_str[,gmn_den])
surf.anc$gmp_str <- log(rowSums(surf_str[,gmp_num])/
                         rowSums(surf_str[,gmp_den]))


surf.gmn_str <- lme(gmn_str~Treatment, random=~1|Rep/Plot,
                   weights=varIdent(form=~1|Treatment),data=surf.anc)
lsmeans(surf.gmn_str, trt.vs.ctrl1 ~ Treatment)

surf.gmp_str <- lme(gmp_str~Treatment, random=~1|Rep/Plot,
                   weights=varIdent(form=~1|Treatment),data=surf.anc)
lsmeans(surf.gmp_str, trt.vs.ctrl1 ~ Treatment)
