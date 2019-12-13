# Analyzing COG profiles----
#
# Let's get real. we need to analyze these functional genes somehow. All well 
# and good, but how to do this? Initial tests suggest things are pretty messy
# in COGland. Fortunately, COGs are placed into functional groups, things like
# "carbohydrate metabolism". This might provide an interesting alternative to 
# throwing everything against the wall (a.k.a. "Plan A"). At this time, my two
# primary approaches are ADONIS (overview) and NMDS (value-added). Other things
# I want to try should probably go in their own script.
#
# 2016-05-10
# David Duncan

# Libraries ----
library(dplyr)
library(ggplot2)
library(lme4)
library(lsmeans)
library(phytools)
library(reshape2)
library(vegan)

# Plotting ----
ord.theme <- theme(panel.grid=element_blank(),
                   panel.background=element_rect(color="black", 
                                                 fill="transparent"),
                   plot.background=element_rect(fill="transparent"),
                   axis.text=element_text(color="black", size=12),
                   axis.title=element_text(size=18),
                   legend.key=element_blank(),
                   legend.text=element_text(size=14),
                   aspect.ratio=1,
                   strip.background=element_blank(),
                   strip.text=element_text(size=18),
                   aspect.ratio=1)

bar.theme <- theme(panel.grid=element_blank(),
                   panel.background=element_rect(color="black", 
                                                 fill="transparent",
                                                 size=1),
                   plot.background=element_rect(fill="transparent"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black"),
                   axis.text.x=element_text(size=11),
                   axis.text.y=element_text(size=10),
                   axis.title=element_text(size=14),
                   strip.background=element_blank(),
                   strip.text=element_text(size=14),
                   legend.key=element_blank())

# Data import ----
# I'm using COG counts relativized by length and by average read depth for a
# set of 37 single-copy housekeeping genes. The file I'm using was generated
# by script 01 in this project on 2016-04-02.

# Information on COG groupings
# Exclude single-copy COGs (which should be invariant post-relativization) and
# eukaryotic genes removed during 2014 update
cog.meta <- read.delim("../data/starting/COG meta.txt") %>% 
  filter(!outdated, !sngl_cpy)

# Relativized COG data
# Exclude the Old Field (G9) treatment
# Relevel factors to be in order of ecological intensification (ish)
# Turn Year into a factor
# Create terms for crop (accounting for fertilization), plot and site-year
# (trust me, we'll want these later)
# Drop single copy COGs and place in the same order as cogs.meta
cog.rel <- read.csv("../data/R outputs/cog_relativized.csv") %>%
  filter(crop != "Old.Field") %>% droplevels() %>%
  mutate(year=factor(year),
         crop=interaction(crop, is.micro, drop = TRUE),
         crop=factor(crop, levels=levels(crop)[c(1,6,2:4,7,8,5)],
                     labels=c("Corn", "Switchgrass", "Miscanthus", 
                              "Native grass","Poplar", "Prairie", 
                              "Unfertilized swichgrass", "Unfertilized prairie"),
                ordered=TRUE),
         plot = interaction(rep, crop, drop = TRUE),
         site.yr=factor(interaction(site, year, drop = TRUE),
                        labels=c("ARL 2010", "ARL 2011", "ARL 2012", 
                                 "KBS 2012"))) %>% 
  select(c(site:is.micro, gbp, sequencing_batch, plot:site.yr, 
           match(cog.meta$cog, names(.))))

# 3.1 Sequencing depth ----
# How much does depth of sequencing vary among samples?

# Effects of metagenome size
# Metagenome size (gbp) accounts for ~4% of intersample distances, small but
# significant. Thus, we need to account for it in our analyses.
summary(cog.rel$gbp)
summary(aov(gbp ~ sequencing_batch + crop * site * year, data = cog.rel))
cld(lsmeans(glm(gbp ~ sequencing_batch, data = cog.rel), 
            pairwise ~ sequencing_batch))
cld(lsmeans(glm(gbp ~ year, data = cog.rel), pairwise ~ year))

### Figure S1 ###
# Samples sequenced in different batches are systematically differnt, although
# samples are still closely linked:

# Identify samples that are duplicates
fs1.dat1 <- cog.rel %>% filter(is.redo) %>% 
  select(site, crop, year, rep, is.micro) 
# Identify their counterparts, then merge everything
fs1.dat2 <- cog.rel %>% filter(!is.redo, site!="KBS") %>%
  merge(fs1.dat1) %>% 
  select(site, crop, year, rep, is.micro) %>% 
  intersect(fs1.dat1) %>%
  merge(cog.rel)

# Run NMDS, extract scores, add sample data back in, sort so sample pairs
# are adjacent
fs1.mds <- select(fs1.dat2, -site:-site.yr) %>%
  metaMDS(trymax=1000) %>% scores %>%
  cbind(select(fs1.dat2, site:sequencing_batch)) %>%
  arrange(crop, year, rep, is.micro)

# Make start/end points for line segments that connect sample pairs
fs1.seg <- cbind(filter(fs1.mds, sequencing_batch=="B2") %>%
                   rename(x1=NMDS1, y1=NMDS2),
                 filter(fs1.mds, sequencing_batch=="B3") %>%
                   select(NMDS1, NMDS2) %>% rename(x2=NMDS1, y2=NMDS2))

# Acutal plot
ggplot(fs1.mds, aes(x=NMDS1, y=NMDS2, color=crop, shape=sequencing_batch)) +
  geom_segment(data=fs1.seg, aes(x=x1, xend=x2, y=y1, yend=y2)) +
  geom_point(size=6) +
  labs(color="Crop", shape="Batch") +
  scale_color_brewer(palette="Dark2") +
  coord_cartesian(xlim = c(-0.035, 0.065), ylim=c(-0.05, 0.05)) +
  ord.theme +
  theme(aspect.ratio=1)
ggsave("../figures/Figure S1.tiff", width=7.5, units="in",
       dpi=600, compression="lzw")

# 3.2 Site-year effects ----
# Are site-years (remember, 3 years at ARL, 1 at KBS) very different?

# Overall adonis (in-text)
adonis(select(cog.rel, -site:-site.yr) ~ gbp + sequencing_batch + (site + year) * crop,
       data=cog.rel, permutations=9999)

### Figure 1/Figure S2 ###
# Size correction
s2.dist <- cog.rel %>% select(-site:-site.yr) %>% vegdist
f1.dist <- capscale(s2.dist ~ cog.rel$gbp + cog.rel$sequencing_batch) %>% residuals

# Figure 1 (corrected)
f1.nmds <- metaMDS(f1.dist, trymax=1000)
for(t.i in 1:9){
  f1.nmds <- metaMDS(f1.dist, trymax=1000, previous.best = f1.nmds)
}
rm(t.i)
f1.nmds$stress #0.134

# Figure S1 (uncorrected)
s2.nmds <- metaMDS(s2.dist, trace=0, trymax=1000)
for(t.i in 1:9){
  s2.nmds <- metaMDS(s2.dist, trace=0, trymax=1000, previous.best = s2.nmds)
}
rm(t.i)
protest(f1.nmds, s2.nmds)

### Figure 1 ###
# Centroids
f1.cent <- aggregate(scores(f1.nmds) ~ site*year, data=cog.rel, mean)
f1.plot <- cbind(cog.rel[ , 1:9], scores(f1.nmds)) %>%
  merge(f1.cent, by = c("site", "year"), suffixes=c("", ".cent"))
# Plot
ggplot(f1.plot, aes(x=NMDS1, y=NMDS2, color=year, shape=site)) +
  geom_segment(aes(xend=NMDS1.cent, yend=NMDS2.cent), alpha=0.5) +
  geom_point(aes(x=NMDS1, y=NMDS2), size=4, alpha=0.6) +
  geom_point(aes(x=NMDS1, y=NMDS2), data=f1.cent, size=7) +
  geom_text(aes(x=-0.06, y=0.045, 
                label=paste0("Stress: ", round(f1.nmds$stress, 3))),
            vjust=0, hjust=0, size=6, color="black") +
  ord.theme +
  labs(color="Year", shape="Site") +
  scale_color_manual(values = c("#a6611a", "#008837", "#7b3294")) +
  coord_cartesian(xlim = c(-0.06, 0.04), ylim=c(-0.055, 0.045))
ggsave("../figures/Figure 1.tiff", width=7.5, units="in",
       dpi=600, compression="lzw")

### Figure S2 ###
# Centroids
s2.cent <- aggregate(scores(s2.nmds) ~ site*year, data=cog.rel, mean)
s2.plot <- cbind(cog.rel[ , 1:9], scores(s2.nmds)) %>%
  merge(s2.cent, by = c("site", "year"), suffixes=c("", ".cent"))
# Plot
ggplot(s2.plot, aes(x=NMDS1, y=NMDS2, color=year, shape=site)) +
  geom_segment(aes(xend=NMDS1.cent, yend=NMDS2.cent), alpha=0.5) +
  geom_point(aes(x=NMDS1, y=NMDS2), size=4, alpha=0.5) +
  geom_point(aes(x=NMDS1, y=NMDS2), data=s2.cent, size=7) +
  geom_text(aes(x=-0.075, y=0.04, 
                label=paste0("Stress: ", round(s2.nmds$stress, 3))),
            vjust=0, hjust=0, size=6, color="black") +
  ord.theme +
  labs(color="Year", shape="Site") +
  scale_color_manual(values = c("#a6611a", "#008837", "#7b3294")) +
  coord_cartesian(xlim = c(-0.075, 0.04), ylim=c(-0.075, 0.04))
ggsave("../figures/Figure S2.tiff", width=7.5, units="in",
       dpi=600, compression="lzw")

# Are differences among (sub)plots measured on different years correlated?
# Subset out the (sub)plots that were sampled in both years
# Aggregate so that replicates (if any) are averaged
# Distance matrices, then mantel them

# 2010 vs 2011
mantel(subset(cog.rel, year == "2010" & 
                plot %in% subset(cog.rel, year == "2011")$plot) %>%
         select(-(1:10), plot) %>% aggregate(. ~ plot, ., mean) %>% 
         arrange(plot) %>% select(-plot) %>% vegdist(),
       subset(cog.rel, year == "2011" & 
                plot %in% subset(cog.rel, year == "2010")$plot) %>%
         select(-(1:10), plot) %>% aggregate(. ~ plot, ., mean) %>% 
         arrange(plot) %>% select(-plot) %>% vegdist(),
       permutations = 9999)

# 2010 vs 2012
mantel(subset(cog.rel, year == "2010" & 
                plot %in% subset(cog.rel, year == "2012")$plot) %>%
         select(-(1:10), plot) %>% aggregate(. ~ plot, ., mean) %>% arrange(plot) %>%
         select(-plot) %>% vegdist(),
       subset(cog.rel, year == "2012" & 
                plot %in% subset(cog.rel, year == "2010")$plot) %>%
         select(-(1:10), plot) %>% aggregate(. ~ plot, ., mean) %>% arrange(plot) %>%
         select(-plot) %>% vegdist(),
       permutations = 9999)

# 2011 vs 2012
mantel(subset(cog.rel, year == "2011" &
                plot %in% subset(cog.rel, year == "2012")$plot) %>%
         select(-(1:10), plot) %>% aggregate(. ~ plot, ., mean) %>% arrange(plot) %>%
         select(-plot) %>% vegdist(),
       subset(cog.rel, year == "2012" & 
                plot %in% subset(cog.rel, year == "2011")$plot) %>%
         select(-(1:10), plot) %>% aggregate(. ~ plot, ., mean) %>% arrange(plot) %>%
         select(-plot) %>% vegdist(),
       permutations = 9999)

# 3.3 Cropping system effects ----
# Effects by cropping system within site-year

# Make subsets by site-year
# Get adonis tables out of it for our trouble
sy.rel <- list()
for(t.i in 1:length(levels(cog.rel$site.yr))){
  sy.rel[[t.i]] <- cog.rel %>% filter(site.yr == levels(site.yr)[t.i])
  print(adonis(sy.rel[[t.i]][, -1:-10] ~ sequencing_batch + gbp + crop, 
               data = sy.rel[[t.i]], 
         permutations=9999))
}
rm(t.i)

# Loop site-years and calculate:
# - Original vs size-corrected data (Procrustes)
# - Spit out size corrected NMDS scores
# - Spit out size corrected, no corn NMDS scores
sy.scores <- c()
sy.stress <- data.frame()
for(t.i in seq_along(sy.rel)){
  # Distances
  t.dist <- vegdist(sy.rel[[t.i]][ , -1:-10])
  t.ns.dist <- capscale(t.dist ~ sy.rel[[t.i]]$gbp + 
                          sy.rel[[t.i]]$sequencing_batch) %>% residuals

  # Full data NMDS
  t.mds <- metaMDS(t.dist, trymax=1000, trace=0)
  for(t.j in 1:9){
    t.mds <- metaMDS(t.dist, trymax=1000, trace=0, previous.best=t.mds)
  }
  
  # Size-corrected NMDS
  t.ns.mds <- metaMDS(t.ns.dist, trymax=1000, trace=0)
  for(t.ns.j in 1:9){
    t.ns.mds <- metaMDS(t.ns.dist, trymax=1000, trace=0, previous.best=t.ns.mds)
  }

  # Output procrustes test
  print(sy.rel[[t.i]]$site.yr[1])
  print(protest(t.mds, t.ns.mds))
  
  # Track stress, for entering in the figures
  sy.stress <- rbind(sy.stress,
                     cbind(select(sy.rel[[t.i]], site.yr) %>% unique, 
                           stress=t.ns.mds$stress))

  # NMDS scores (all crops)
  sy.scores <- rbind(sy.scores,
                     cbind(sy.rel[[t.i]][ , 1:10], 
                           scores(t.ns.mds)))
}

# Calculate centroids all at once, 'cause that's how I roll now
sy.scores <- aggregate(cbind(NMDS1, NMDS2) ~ crop*site.yr, 
                       data=sy.scores, mean) %>% 
  merge(sy.scores, by=c("crop", "site.yr"), suffixes=c(".cent", ""))

sy.stress <- mutate(sy.stress,
                    x=c(-0.036, -0.025, -0.023, -0.026),
                    y=c(0.018, 0.025, 0.02, 0.03))

### Figure 2 ###
sy.scores <- mutate(sy.scores,
                    fert=factor(!grepl("fertilized", crop), 
                                levels=c("TRUE", "FALSE"), 
                                labels=c("Fertilized", "Unfertilized")),
                    avp=factor(crop, 
                               levels=levels(crop)[c(1:3, 5, 4, 6:8)]))
levels(sy.scores$avp) <- levels(sy.scores$avp)[c(1:6, 2, 6)]
sy.scores <- mutate(sy.scores,
                    invis=factor((avp=="Corn" | (avp=="Prairie" & fert=="Unfertilized"))))
r2.lab <- as_labeller(c("ARL 2010"="Crop~R^2*':'~0.229~'***'",
                        "ARL 2011"="Crop~R^2*':'~0.113~italic(n.s.)",
                        "ARL 2012"="Crop~R^2*':'~0.121~'*'",
                        "KBS 2012"="Crop~R^2*':'~0.284~'***'"),
                      default=label_parsed)

sy.scores %>%
  filter(site.yr=="ARL 2010") %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=avp, alpha=invis, shape=fert)) +
  geom_segment(aes(xend=NMDS1.cent, yend=NMDS2.cent), size=1) +
  geom_point(color="white", alpha=1, size=4) +
  geom_point(size=4) +
#  geom_text(data=filter(sy.stress), color="black",
#                        aes(label=paste0("Stress ", round(stress, 3)),
#                            x=x, y=y), hjust=0) +
  scale_color_manual(
    values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d")) +
  scale_alpha_discrete(range=c(0,1), guide=FALSE) +
  scale_shape_manual(values=c(19, 17)) +
  labs(color="Cropping system", shape="Fertilization") +
  facet_wrap(~site.yr, ncol=2, scales="free") +
  ord.theme +
  theme(legend.position="right",
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="transparent"),
        axis.text=element_blank(),
        axis.ticks=element_blank())
ggsave("../figures/Figure 2.tiff", height=10, width=10, units="in",
       dpi=600, compression="lzw")
ggsave("../../images/site-year ndms a10.png", height=7, width=9.5, units="in",
       bg="transparent")

# 3.4 COG function sets ----

# Define COG categories
# Function groups are given single letter codes in the COG database.
# A given COG can have more than one code (i.e. fall into more than one 
# functional group). I place the COG into *all* groups to which it belongs.
cog.grp <- list(all = "[A-Z]",         # Everything with a function    
                def = "[V]",           # Defense mechanisms
                sig = "[T]",           # Signal transduction
                nrg = "[C]",           # Energy production & conversion
                crb = "[G]",           # Carbohydrate transport & metabolism
                npk = "[P]",           # Inorganic ion transport & metabolism
                scd = "[Q]")           # 2o metabolite biosynth, transp, catab

# Generate a list of COG subsets and a table comparing the subsets to the full
# dataset
# The list of subsets can be used in subsequent analyses so I don't have to 
# keep subsetting all the time

grp.sub <- vector("list", length(cog.grp))   # Subsets
grp.dist <- vector("list", length(cog.grp))  # Inter-sample distances
grp.table <- data.frame(grp = names(cog.grp),
                        ncogs = 0,          # Number of COGs in group
                        copies.05 = 0,      # 5th, 50th, 95th percentile of
                        copies.50 = 0,      # copies per cell
                        copies.95 = 0,     
                        mantel.r= 0,        # Mantel r between group and all
                        mantel.p = 1,       # Mantel p-value
                        adonis.batch = 0,   # R2 of sequencing batch
                        adonis.site = 0,    # R2 of site
                        adonis.year = 0,    # R2 of year
                        adonis.crop = 0,    # R2 of cropping system
                        adonis.batch.p = 0, # Batch factor p-val
                        adonis.site.p = 0,  # Site factor p-val
                        adonis.year.p = 0,  # Year factor p-val
                        adonis.crop.p = 0)  # Crop factor p-val

### Table 2 ###
for(t.i in seq_along(cog.grp)){
  grp.sub[[t.i]] <- cog.rel %>% select(c(1:10, match(
    (cog.meta %>% filter(grepl(cog.grp[t.i], funct)))$cog, names(.))))
  
  # Number of COGs in the group + average copy number per cell
  grp.table[t.i, 2] <- length(grp.sub[[t.i]]) - 10
  grp.table[t.i, 3:5] <- grp.sub[[t.i]] %>% select(-(1:10)) %>% 
    rowMeans %>% quantile(c(0.05, 0.5, 0.95))
  
  grp.dist[[t.i]] <- grp.sub[[t.i]] %>% select(-(1:10)) %>% vegdist()
  
  # Mantel statistic for group vs full set of COGs
  grp.table[t.i, 6:7] <- mantel(grp.dist[[t.i]], grp.dist[[1]])[
    c("statistic", "signif")]
  
  # ADONIS R2 values for year and site
  grp.table[t.i, 8:15] <- 
    adonis(grp.dist[[t.i]] ~ sequencing_batch + site + year + crop, 
           data = grp.sub[[t.i]], permutations = 9999)$aov.tab[1:4, 5:6] %>% 
    unlist()
}

# Is the number of COGs correlated to mantel correlation?
grp.table %>% filter(grp != "all") %>%
  lm(mantel.r ~ ncogs, data = .) %>% summary()
grp.table %>% filter(grp != "all") %>% select(ncogs, mantel.r) %>% cor()

rm(list = ls(pattern = "t\\."))
save.image("04 Step 02 complete.RData")

### Table 3 ###
# Make lists and dataframes to hold everything, just like before
crop.table <- cbind(merge(levels(cog.rel$site.yr), names(cog.grp)),
                    mantel.r = 0, 
                    mantel.p = 0, 
                    adonis.batch = 0,
                    adonis.crop = 0, 
                    adonis.batch.p = 0,
                    adonis.crop.p = 0) %>% rename(site.yr = x, grp = y)
crop.sub <- list()
crop.dist <- list()

# Loop to make the big table
for(t.i in cog.grp %>% names %>% seq_along()){
  for(t.j in cog.rel$site.yr %>% levels %>% seq_along()){
    # Useful index
    t.k <- (t.i - 1) * (cog.rel$site.yr %>% levels() %>% length()) + t.j
    # Make subset
    crop.sub[[t.k]] <- grp.sub[[t.i]] %>% 
      filter(site.yr == levels(site.yr)[t.j])
    # Distance matrix
    crop.dist[[t.k]] <- crop.sub[[t.k]] %>% select(-(1:10)) %>% vegdist()
    # Mantel of function vs full dataset
    crop.table[t.k, 3:4] <- mantel(crop.dist[[t.k]], crop.dist[[t.j]])[
      c("statistic", "signif")]
    # Adonis values
    crop.table[t.k, 5:8] <- 
      adonis(crop.dist[[t.k]] ~ sequencing_batch + crop, data = crop.sub[[t.k]], 
             permutations = 9999)$aov.tab[1:2, 5:6] %>% unlist()
  }
}
rm(list = ls(pattern = "t\\."))

### Fig S3 ###
crop.scores <- c()
crop.stress <- data.frame()

for(t.i in seq_along(names(cog.grp))){
  for(t.j in seq_along(levels(cog.rel$site.yr))){
    
    # Subset
    t.sub <- grp.sub[[t.i]] %>% filter(site.yr == levels(site.yr)[t.j])
    # correct distance by metagenomesize
    t.dist <- t.sub %>% select(-1:-10) %>% vegdist
    t.dist <- capscale(t.dist ~ t.sub$gbp + t.sub$sequencing_batch) %>% residuals
    
    # NMDS
    t.mds <- metaMDS(t.dist, trace=0, trymax = 1000)
    for(t.l in 1:9){
      t.mds <- metaMDS(t.dist, trace=0, trymax = 1000, previous.best = t.mds)
    }
    
    # Outputs
    crop.stress <- rbind(crop.stress,
                         cbind(select(t.sub, site.yr) %>% unique, 
                               grp=names(cog.grp)[t.i],
                               stress=t.mds$stress))
    crop.scores <- rbind(crop.scores,
                         cbind(select(t.sub, crop, site, site.yr),
                               grp=names(cog.grp)[t.i],
                               scores(t.mds)))
    }
}

crop.scores <- aggregate(cbind(NMDS1, NMDS2) ~ crop*site.yr*grp,
                         data=crop.scores, mean) %>%
  merge(crop.scores, by=c("crop", "site.yr", "grp"), suffixes=c(".cent", ""))

# Plotting
crop.scores %>% filter(grp=="sig") %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=crop)) +
  geom_point(size=4) +
  geom_segment(aes(xend=NMDS1.cent, yend=NMDS2.cent), size=1) +
  scale_color_brewer(palette="Dark2") +
  labs(color="") +
  facet_wrap(~site.yr, ncol=2, scales="free") +
  ord.theme +
  theme(legend.position="bottom")
ggsave("../figures/Figure S3.tiff", height=10, width=10, units="in",
       dpi=600, compression="lzw")


### Figure 3 ###
# Ok, we've seen differences in overall profiles, now are certain function sets
# more abundant in particular systems? 
# Does this change among site-years?

# Unstack data so each COGs abundance is a single row
cog.unstk <- cog.rel %>% melt(id.vars = c(1:10))
cog.vals <- c()

for(t.i in seq_along(cog.grp[-1])){
  t.lmer <- cog.unstk %>% 
    subset(variable %in% (cog.meta %>% filter(grepl(cog.grp[t.i+1], funct)))$cog) %>% 
    lmer(asinh(value) ~ crop * site.yr + gbp + 
           (1|rep/plot) + (1 + gbp|variable), data = .) %>% 
    lsmeans(pairwise ~ crop|site.yr) %>% 
    cld(Letters=letters) %>%
    select(crop, site.yr, lsmean, SE, .group) %>%
    mutate(.group=gsub(" ", "", .group),
           t.mean=sinh(lsmean),
           t.min=sinh(lsmean-SE),
           t.max=sinh(lsmean+SE))
  cog.vals <- rbind(cog.vals,
                    cbind(t.lmer, grp=names(cog.grp)[t.i+1]))
}

cog.vals <- filter(cog.vals, !is.na(lsmean))
### Fig 3 ###
cog.labs <- as_labeller(c(def="Defense (V)", sig="Signalling (T)", 
                          nrg="Energy (C)", crb="Carbohydrates (G)",
                          npk="Inorganic ions (P)", 
                          scd="Secondary metabolites (Q)"))
cog.vals %>% 
  ggplot(aes(x=site.yr, fill=crop, y=t.mean)) + 
  geom_crossbar(aes(ymin=t.min, ymax=t.max),
                position="dodge") +
  geom_text(aes(label=.group, y=t.max), vjust=0.2, hjust=-0.3,
            fontface="italic", position=position_dodge(width=0.9), angle=90) +
  scale_fill_brewer(palette="Dark2") +
  facet_wrap(~grp, ncol=2, scales="free_y", labeller=cog.labs) +
  labs(y=expression(Estimated~copies~cell^{-1}),
       x="", fill="") +
  scale_y_continuous(expand=c(0.1, 0)) +
  bar.theme +
  theme(legend.position = "bottom")
ggsave("../figures/Figure 3.tiff", width=7.5, height=9, units="in",
       dpi=600, compression="lzw")

# 3.5 Denitrification pathway COGs ----

den.cog <- read.delim("../data/starting/denit COGs.txt")

### Figure 4 ###
den.vals <- c()

cog.unstk <- cog.unstk %>% mutate(crop.g=crop)
levels(cog.unstk$crop.g) <- c("Corn", rep("Fertilized perennials", 5), 
                              rep("Unfertilized perennials", 2))

cog.denit <- merge(cog.unstk, den.cog, by.x="variable", by.y="cog") %>%
  mutate(fold=log2(value)) %>% droplevels

for(t.i in seq_along(levels(cog.denit$step))){
  t.sub <- cog.denit %>% filter(step==levels(step)[t.i])
  t.means <- t.sub %>% group_by(variable) %>%
    summarise(m.fold=mean(fold)) %>% ungroup %>% data.frame
  t.sub <- merge(t.sub, t.means) %>%
    mutate(fold.change=fold-m.fold)
  t.lmer <- lmer(fold.change ~ crop.g*site.yr*gbp + (1|plot/site.yr) + (1|variable/site.yr), 
                 data=t.sub) %>%
    lsmeans(pairwise ~ crop.g|site.yr) %>% cld(Letters=LETTERS) %>%
    select(crop.g, site.yr, lsmean, SE, .group)
  
  den.vals <- rbind(den.vals,
                     cbind(t.lmer, step=levels(cog.denit$step)[t.i])) %>%
    mutate(.group=gsub(" ", "", .group))
}

den.vals %>% filter(!is.na(lsmean)) %>%
  mutate(step=factor(step, levels=levels(step)[4:1])) %>%
  ggplot(aes(x=step, y=lsmean, color=crop.g)) +
  geom_hline(yintercept=0) +
  geom_point(size=3, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.5,
                position=position_dodge(width=0.5)) +
  geom_text(aes(label=.group, y=lsmean+SE+0.05, group=crop.g), color="black",
            position=position_dodge(width=0.5)) +
  facet_wrap(~site.yr) +
  labs(x="", color="",
       y="Fold-change from global mean copy number") +
  scale_color_brewer(palette="Dark2") +
  bar.theme +
  theme(axis.text.x = element_text(color="white"),
        legend.position="bottom")
ggsave("../../images/denit.png", width=9.5, height=5.6, units="in")
ggsave("../figures/Figure 4.tiff", width=10.5, units="in",
       dpi=600, compression="lzw")
