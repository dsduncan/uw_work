# Assembling the data
#
# Here is where we'll pull all of the data together for the first time.
# We have flux/ancillary data and then microbial community data for a subset
# of observations. The basic unit of measurement here is the plot-year.
# In this section, I want to:
# - Note which field observations have associated microbial data
# - Figure out the number of observations (flux and full data) per plot-year
# - Compare distributions of field variables with and without microbial data
#
# 2016-06-05
# David Duncan

# Libraries ----
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(lubridate)
library(reshape2)
library(scales)

# Plotting ----
main.theme <- theme(panel.grid=element_blank(),
                    panel.background=element_rect(color="black", fill="white"),
                    axis.text=element_text(color="black", size=8),
                    axis.title=element_text(size=12),
                    legend.key=element_blank(),
                    legend.position="bottom",
                    strip.background=element_blank())

# Functions for asinh scaling
asinh_trans <- function(){trans_new("asinh", asinh, sinh, asinh_breaks())}
asinh_breaks <- function(){
  function(x){
    rng <- c(-1*10^(100:0), 0, 10^(0:100))
    min <- max(which(rng < min(x)))
    max <- min(which(rng > max(x)))
    return(rng[min:max])
  }
}

# Data import ----
# Microbial data (identifying info + normalized COG abundances)
# Need to create a plot.year variable to match up the field data
microb <- read.delim("../data/R inputs/crop to trt.txt") %>%
  merge(read.csv("../data/R inputs/cog_relativized.csv")) %>%
  rename(block=rep) %>%
  mutate(plot.year=factor(interaction(block, trt, year))) %>%
  select(plot.year, crop:is.redo, COG0001:COG5665)

# Field data (fluxes, temp, wfps, nitrogen)
# Create plot.year variable
field <- read.delim("../data/R inputs/crop to trt.txt") %>%
  merge(read.delim("../data/R inputs/full flux dat.txt")) %>%
  droplevels %>%
  mutate(date=as.Date(date),
         year=year(date),
         plot=factor(interaction(block, trt)),
         plot.year=factor(interaction(plot, year)),
         system=factor(sub("M", "", trt)), 
         avp=system, fert=trt,
         in.microb=plot.year %in% microb$plot.year) %>%
  merge(read.delim("../data/R inputs/fert dates.txt")) 

# Were treatments fertilized (note this doesn't mean fertilized in any given
# year, only fertilized at all)
levels(field$fert) <- list(fertilized=levels(field$trt)[c(1:9,11,13,15,16,19)],
                           unfertilized=levels(field$trt)[c(10,12,14,17,18)])
# Lumping annual vs perennial because I can't even with the rotations
levels(field$avp) <- c(rep("annual", 4), levels(field$system)[5:10])
field$avp <- factor(field$avp, 
                      labels=c("Annual", "Switchgrass", "Miscanthus", 
                               "Native grass", "Poplar", "Old field",
                               "Prairie"))
# Useful to compare fertilized/unfertilized by site
field <- field %>% 
  mutate(site.fert=factor(interaction(fert, site, sep=" "), 
                          labels=c("ARL fertilized", "ARL unfertilized", 
                                   "KBS fertilized", "KBS unfertilized")),
         avp.fert=factor(interaction(avp, fert, sep=" "), 
                         levels=c("Annual fertilized", "Switchgrass fertilized", 
                                  "Miscanthus fertilized", "Native grass fertilized",
                                  "Poplar fertilized", "Old field fertilized", 
                                  "Prairie fertilized", "Switchgrass unfertilized",
                                  "Miscanthus unfertilized", "Native grass unfertilized",
                                  "Old field unfertilized", "Prairie unfertilized")))

# Correct fertilization dates for unfertilized treatments
field[field$fert=="unfrt", "fert.time"] <- "unfrt"

# How many observations? ----
# Non-negative flux observations for which we have microbial data and all 
# environmental data
field %>% 
  filter(in.microb, !is.na(soilT), !is.na(wfps), !is.na(no3), n2o.ha.day>0) %>%
  select(plot.year) %>% droplevels %>% table
# Median: 13.5, min: 6, max: 19

### Figure 1 ###
# How are things distributed?
# Labeler to give proper names to things
f1.lab <- 
  as_labeller(c(n2o.ha.day="N[2]*O~flux~'('*g~N~ha^{-1}~day^{-1}*')'",
                no3="NO[3]^'-'~'('*mu*g~N~g^{-1}*')'",
                nh4="NH[4]^'+'~'('*mu*g~N~g^{-1}*')'",
                soilT="Soil~temp~(degree~C)",
                wfps="WFPS~('%')"),
              default=label_parsed)

# Top panels (need asinh transformation)
f1 <- field %>% 
  select(n2o.ha.day, no3, nh4, soilT, wfps, site, in.microb) %>%
  melt %>% filter(!is.na(value)) %>%
  ggplot(aes(y=value, x=site, fill=in.microb)) +
  geom_hline(yintercept=0) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free", switch="y", labeller=f1.lab) +
  scale_fill_brewer(palette="PRGn") +
  scale_y_continuous(trans="asinh") +
  labs(fill="Environmental and microbial data", xlab="Site") +
  main.theme + 
  theme(legend.position=c(0.85, 0.4),
        axis.title.y=element_blank())

# Lower panels (don't need transformation)
f1.a <- f1 +
  scale_y_continuous()
  
f1 <- ggplotGrob(f1)
f1.a <- ggplotGrob(f1.a)

f1$grobs[c(5,6,15,16)] <- f1.a$grobs[c(5,6,15,16)]
grid::grid.draw(f1)

ggsave("../figures/Figure 1.tiff", f1,
       dpi=600, compression="lzw")

# Data export ----
write.table(field, "../data/R outputs/field measurements.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)

write.table(microb, "../data/R outputs/microbe measurements.txt", sep="\t",
            quote=FALSE, row.names=FALSE)