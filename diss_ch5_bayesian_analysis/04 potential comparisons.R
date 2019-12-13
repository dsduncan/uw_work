# The rest of the story ----
# 
# We have three ways of trying to get at the potential N2O productino of a
# cropping system:
# - Cumulative fluxes are time-weighted over an entire season
#   Benefit from being inherently intersting, sensitive to weather
# - Peak fluxes are the best indicator we have of potential
#   Still sensitive to weather, less inherently interesting, should be tighter
# - Potential de/nitrification parameters from pseudo NOE model
#   Does this even work?
#
# Step 1 is to compare everything to everything. Are these ranking plots
# similarly? Esp in re: within variability. 
# Step 2 is elastic net modeling of the quantities with microbial functional
# gene profiles. Deal with that in due time though...
#
# 


# Libraries ----
library(dplyr)
library(ggplot2)
library(glmnet)
library(reshape2)

# Plotting ----
main.theme <- theme(panel.grid=element_blank(),
                    panel.background=element_rect(color="black", fill="white"),
                    axis.text=element_text(color="black", size=8),
                    axis.title=element_text(size=12),
                    legend.key=element_blank(),
                    legend.position="bottom",
                    strip.background=element_blank())

# 15-level allegedly dichromat-friendly palette
cf.pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d")

# Functions ----

# Data entry ----
agg.dat <- read.delim("../data/R outputs/n2o agg.txt") %>%
  rename(agg.flux=n2o.agg) %>%
  select(site:avp, agg.flux, in.microb)

peak.dat <- read.delim("../data/R outputs/n2o peaks.txt") %>%
  filter(status=="peak") %>% 
  rename(peak.flux=n2o.ha.day) %>%
  mutate(plot.year=factor(interaction(plot, year))) %>%
  select(site:block, plot:in.microb, peak.flux)

pot.dat <- read.delim("../data/R outputs/N2O potentials from RStan.txt")

microb.dat <- read.delim("../data/R outputs/microbe measurements.txt") %>%
  mutate(year=factor(year))

# Parameter placeholder
# Need PDR, PNR values
# Merge needs to allow NAs (some plots won't have enough data for NOE)

all.dat <- merge(agg.dat, peak.dat) %>% merge(pot.dat, all=TRUE) %>%
  mutate(avp=factor(avp, levels=levels(agg.dat$avp)[c(1, 7, 2, 3, 5, 4, 6)]))

# 3.5 Correlationsville ----
### Figure 9 ###
# Correlating other terms to aggregate emissions
all.dat %>%
  melt(measure.vars=c("peak.flux", "PDR", "PNR")) %>%
  filter(!is.na(value)) %>%
  mutate(variable=factor(variable, labels=c("Peak flux", "PDR", "PNR"))) %>%
  ggplot(aes(x=agg.flux, y=value)) +
  geom_point(aes(color=avp), size=2) + 
  stat_smooth(method="lm", se=FALSE, color="black") +
  labs(x=expression(Aggregate~N[2]*O~flux~(g~N~ha^{-1}~year^{-1})),
       y=expression(N[2]*O~flux~(g~N~ha^{-1}~day^{-1})),
       color="Cropping system")+
  scale_x_log10() + scale_y_log10() +
  scale_color_brewer(palette="Dark2") +
  facet_grid(variable~site, scale="free") +
  main.theme +
  theme(aspect.ratio=1,
        legend.text=element_text(size=10),
        axis.text=element_text(size=10),
        strip.text=element_text(size=12))
  
ggsave("../figures/Figure 9.tiff", height=9.5, units="in",
       dpi=600, compression="lzw")

# Write out ----
all.dat %>% filter(in.microb) %>%
  merge(microb.dat) %>% droplevels %>%
  write.table("../data/R outputs/Pot and microb.txt", sep="\t",
              row.names=FALSE, quote=FALSE)
