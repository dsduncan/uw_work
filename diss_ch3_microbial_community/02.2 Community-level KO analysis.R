# Community level microbial analysis ----
#
# Let's try this again, this time using KO annotation instead of COGs.
# "But David, you love COGs" I can hear you say and it's true, I'm very fond of
# them. After working with them for a while though, I am increasingly worried 
# about the level of support and curation given to their annotation. There is
# also the difficulty of placing them into a context, all of which the KEGG 
# framework does much better. Plus it seems to be curated more frequently and
# carefully, AND it is much easier to look into pathways and reactions.
# Ok, now the bad news: I'm going to be missing 33 samples because of an 
# annotation issue at IMG. 

# Libraries ----
library(dplyr)
library(ggplot2)
library(lme4)
library(lsmeans)
library(vegan)

# Data import ----
usic <- read.delim("../data/starting/elite KO USiCGs.txt", 
                   stringsAsFactors = FALSE, header=FALSE) %>% unlist

# Import data
# Remove Old Field treatment (too noisy and hard to characterize)
# Relevel crop, turn year into a factor, create trt, plot, and site-year
# Strip out USiCs (they should all be ~1 anyway)
# 
ko <- read.delim("../data/R outputs/KO copy numbers.txt") %>%
  select(-one_of(usic)) %>%
  filter(crop != "Old.Field") %>% droplevels %>%
  mutate(crop=factor(crop, levels=c("Corn", "Switchgrass", "Miscanthus",
                                    "Native.Grass", "Poplar", "Prairie")),
         year=factor(year),
         trt=interaction(crop, is.micro, drop=TRUE),
         plot=interaction(rep, crop, drop=TRUE),
         site.yr=interaction(site, year, drop=TRUE)) %>%
  select(site:m.seqs, trt:site.yr, K00001:K16881)

# 3.1 overview and effort ----
ko %>% select(year, site, trt) %>% table

# Sequencing effort spread
max(ko$gbp); min(ko$gbp); median(ko$gbp)

# Sequencing effort differences
summary(aov(gbp ~ crop * year * site, data = ko))
cld(lsmeans(glm(gbp ~ crop, data=ko), pairwise ~ crop))
cld(lsmeans(glm(gbp ~ year, data=ko), pairwise ~ year))


# 3.2 Site-Year effects ----
# Note that this may be *heavily* changed if we're able to rescue the 
# KBS samples from KO hell
sityr.dist <- ko %>% select(-site:-site.yr) %>% vegdist
adonis(sityr.dist ~ gbp + site + year, data=ko, permutations=9999)

sityr.mds <- metaMDS(sityr.dist)
for(t.i in 1:9){
  sityr.mds <- metaMDS(sityr.dist, previous.best=sityr.mds)
}
sityr.mds$stress #0.1267

# Plotting
sityr.cent <- aggregate(scores(sityr.mds) ~ site * year, data = ko, mean)
sityr.plot <- merge(cbind(select(ko, site:site.yr), scores(sityr.mds)),
                    sityr.cent, 
                    by =c("site", "year"), suffixes=c("", ".cent"))
ggplot(sityr.plot) +
  geom_point(aes(x=NMDS1, y=NMDS2, color=year, shape=site), size=4)

# 3.3 Cropping system effects (within site-year)
# Holders for adonis and mds outputs, so we can loop them all at once first
crop.adon <- list()
crop.adon2 <- list()
crop.plot <- c()
crop.plot2 <- c()

for(t.i in seq_along(levels(ko$site.yr))){
  # Subset
  t.dat <- ko %>% filter(site.yr==levels(site.yr)[t.i])
  # Distance
  t.dist <- t.dat %>% select(-site:-site.yr) %>% vegdist
  # ADONIS (with size effect)
  crop.adon[[t.i]] <- adonis(t.dist ~ gbp + crop, data=t.dat, permutations=9999)
  # MDS
  t.mds <- metaMDS(t.dist, trymax=1000)
  for(t.i in 1:9){
    t.mds <- metaMDS(t.dist, trymax=1000, previous.best=t.mds)
  }
  t.cent <- aggregate(scores(t.mds) ~ trt, data = t.dat, mean)
  t.plot <- merge(cbind(select(t.dat, site:rep, gbp, trt), scores(t.mds),
                        site.yr=t.dat$site.yr[1], size=TRUE),
                      t.cent, 
                      by =c("trt"), suffixes=c("", ".cent"))
  crop.plot <- rbind(crop.plot, t.plot)
  
  # Remove size effect via capscale
  t.dist2 <- capscale(t.dist ~ t.dat$gbp) %>% residuals
  # ADONIS (without size effect)
  crop.adon2[[t.i]] <- adonis(t.dist2 ~ crop, data = t.dat, permutatins=9999)
  # MDS without size effects
  t.mds2 <- metaMDS(t.dist2, trymax=1000)
  for(t.i in 1:9){
    t.mds2 <- metaMDS(t.dist2, trymax=1000, previous.best=t.mds2)
  }
  t.cent2 <- aggregate(scores(t.mds2) ~ trt, data = t.dat, mean)
  t.plot2 <- merge(cbind(select(t.dat, site:rep, gbp, trt), scores(t.mds2),
                        site.yr=t.dat$site.yr[1], size=FALSE),
                  t.cent2, 
                  by =c("trt"), suffixes=c("", ".cent"))
  crop.plot2 <- rbind(crop.plot2, t.plot2)
}

crop.plot2[crop.plot2$site.yr=="ARL.2011", "NMDS1"] <- 
  crop.plot2[crop.plot2$site.yr=="ARL.2011", "NMDS1"] * -1
crop.plot2[crop.plot2$site.yr=="ARL.2011", "NMDS1.cent"] <- 
  crop.plot2[crop.plot2$site.yr=="ARL.2011", "NMDS1.cent"] * -1

ggplot(crop.plot2, aes(color=trt)) +
  geom_point(aes(x=NMDS1, y=NMDS2, shape = site), size=4) +
  geom_point(data=unique(select(crop.plot2, trt, site, site.yr, NMDS1.cent:NMDS2.cent)),
             aes(x=NMDS1.cent, y=NMDS2.cent, shape=site), size = 7) +
  facet_wrap(~site.yr, nrow=2)


