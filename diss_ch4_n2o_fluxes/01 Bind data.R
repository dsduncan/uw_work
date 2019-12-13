# Combining data
#
# Here goes. I need to combine N2O flux data with whatever anicillary data I
# can pull together. This may take a few iterations, so it makes more sense to
# modularize things at this juncture rather than try to do it all manually.
# Let's keep things simple people!
#
# 2016-05-31
# David Duncan

# Libraries ----
library(dplyr)
library(ggplot2)
library(lubridate)

# N2O flux data ----
# N2O fluxes (g N2O-N per ha per day) on a daily time step, 2009-2014
# Need to filter out bad negative fluxes
# Rules: 1) Drop all fluxes smaller than -0.85 g/ha/day
#           This retains 76% of negative fluxes
#        2) Drop all negative fluxes fit with HMR, coded as no flux, or fit
#           with only 3 points, given the weakness of evidence there
# These cut off 500 values (out of ~11,000), so the damange isn't too bad...
n2o <- read.delim("../data/R inputs/N2O_fluxes.txt") %>%
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(n2o.ha.day > -0.85, 
         !(n2o.ha.day <= 0 & (Method =="HMR" | 
                                Warning %in% levels(Warning)[3:5]))) %>%
  select(site:date, n2o.ha.day) %>% 
  aggregate(n2o.ha.day ~ site + trt + block + date, ., mean)

# Inorganic N ----
# (ug N per g dry soil)
#====> NOTE <====#
# I'm presently removing all of 2010-10-20, since nearly all of those numbers
# are negative, but JOS is looking into it to see if he can salvage anything.
# Sub-zero values replaced with machine's detection limit (0.05 ug g)
minN <- read.delim("../data/R inputs/soil_inorganic_N.txt") %>%
  select(-date_real) %>% mutate(date=as.Date(date, "%m/%d/%Y")) %>%
  filter(date > "2009-12-31", date != "2010-10-20")
minN[minN$no3 <= 0, "no3"] <- 0.05 # Detection limit

# Bulk density ----
# Download directly from data catalog (Table 221)
# We'll use deep core sampling campaigns from 2008 & 2013, discarding the 
# surface sampling campaigns. We'll only keep the top 10 cm for G1-G10.
bd <- read.csv("../data/R inputs/soil_bulk_density.csv") %>%
  rename(trt=treatment, block=replicate, date=sample_date, 
         bulk_density=bulk_density_gravel_free) %>%
  filter(campaign == "deep core sampling", !(trt %in% c("G11", "G12")),
         horizon_top_depth == 0, bulk_density > 0) %>% droplevels() %>%
  select(year, site:station, bulk_density) %>%
  mutate(block=interaction(site, block), 
         trt=factor(trt, levels=levels(trt)[c(1,3:10,2)], ordered=TRUE,
                    labels=c(paste0("G0", 1:9), "G10")),
         site=factor(site, labels=c("ARL", "KBS")))

# Group treatments as annuals, perennials, poplar
# Group blocks 1&3, 2, 4&5 at ARL, group everything from KBS since it seems
# pretty homogeneous (see histograms below)
bd$trt.grp <- bd$trt
levels(bd$trt.grp) <- list(annual=levels(bd$trt)[1:4], poplar="G08",
                           perren=levels(bd$trt)[c(5:7,9:10)])
levels(bd$block) <- sub("([AK]).*R", "\\1", levels(bd$block))
bd$plot <- interaction(bd$trt, bd$block)
bd$block.grp <- bd$block
levels(bd$block.grp) <- list(A13=c("A1", "A3"), K=paste0("K",1:5), A2="A2", 
                             A45=c("A4", "A5"))

# Plotting bulk density, don't need to run every time
#ggplot(filter(bd, site=="ARL"), aes(x=bulk_density, fill=factor(year))) + 
#  geom_histogram() +  facet_grid(block.grp~trt.grp)
#ggplot(filter(bd, site=="KBS"), aes(x=bulk_density, fill=factor(year))) +
#  geom_histogram() +  facet_grid(.~trt.grp)

# Get means for treatments, blocks, get year slope
# Drop top 9 outliers
# Predict bulk density 
bd.lm <- lm(bulk_density~year*block.grp*trt.grp, 
            data=bd[-c(261, 418, 420, 463, 476, 499, 515, 561, 568),])
bd.avg <- expand.grid(year=2008:2015, trt.grp=levels(bd$trt.grp),
                      block.grp=levels(bd$block.grp)) %>%
  mutate(bd_avg=as.numeric(predict(bd.lm, .)))

bd.out <- bd %>% select(trt, trt.grp, block, block.grp, plot) %>% unique() %>%
  merge(2009:2015) %>% rename(year=y) %>% merge(bd.avg) %>% 
  select(year, plot, bd_avg)

# Temp and moisture ----
# Ultimately, want to merge these with KBS data on data catalog, and should
# probably import everything at once.
arl.tm <- read.delim("../data/R inputs/soil_Tgwc_arl.txt") %>%
  mutate(date = as.Date(date, "%m/%d/%Y"), 
         date_real = as.Date(date_real, "%m/%d/%Y"))

# KBS data downloaded from catalog 2016-05-12

# Temp data just need some relabeling of factors
kbs.t <- read.csv("../data/R inputs/soil_T_kbs.csv") %>%
  mutate(date_real=as.Date(sample_date, "%m/%d/%Y")) %>% 
  select(-sample_date) %>%
  filter(date_real > "2008-12-31", date_real < "2015-01-31") %>% 
  rename(soilT=soil_temperature, trt=treatment, block=replicate) %>%
  mutate(site=factor("KBS"),
         date=date_real,
         trt=factor(trt, labels=paste0("G", c("0", "", rep("0", 8)), 
                                       c(1, 10, 2:9))),
         block=factor(block, labels=gsub("R", "K", levels(block))))

# Moisture data
# Some dates don't match N2O sampling, so we need to fix that
# Treatment and block need to be turned into factors
# Non-mainsite data should be dropped

# First, we need the list of dates when we sampled N2O at KBS
samp.dates <- n2o %>% filter(site == "KBS") %>% select(date) %>% 
  unique %>% as.vector

kbs.m <- read.csv("../data/R inputs/soil_gwc_kbs.csv") %>%
  mutate(date_real=as.Date(sample_date, "%m/%d/%Y"),
         date=samp.dates[sapply(date_real, function(x){
           which.min(abs(x-samp.dates[,1]))}), 1],
         trt=factor(treatment, labels=paste0("G", c(rep("0",9), ""), 1:10)),
         block=factor(replicate, labels=paste0("K", 1:4)),
         vwc=moisture*100,
         site=factor("KBS")) %>%
  filter(grepl("Main", Site), 
         date_real < "2015-01-31", date_real > "2008-12-31") %>%
  select(date_real:site)

# Merge the two KBS datasets (remember to keep missing data!) 
kbs.tm <- merge(kbs.t, kbs.m, all=TRUE)
vwc_t <- rbind(arl.tm, kbs.tm) %>%
  mutate(plot=interaction(substr(trt, 1, 3), block),
         year=year(date))

# Now, merge vwc_t and bd to link plot-level bulk densities
# Note that these are main plot measurements only, but that's the
# best we can do for the microplots
# We need to convert KBS GWC values (just... don't even start) to VWC
# values. At that point, we might as well also calculate WFPS, although at
# this point, I have serious doubts about this whole "soil moisture" metric.
# This is also our last chance to deal with duplicate values before the big
# merge.
wfps.dat <- merge(vwc_t, bd.out) 
wfps.dat <- rbind(filter(wfps.dat, site == "ARL"),
                  filter(wfps.dat, site == "KBS") %>% 
                    mutate(vwc = vwc * bd_avg)) %>%
  mutate(wfps = vwc / (1 - bd_avg / 2.65)) %>%
  select(-plot, -bd_avg, -date_real, -year)
wfps.dat <- merge(
  aggregate(soilT ~ site * trt * block * date, data = wfps.dat, mean),
  aggregate(wfps ~ site * trt * block * date, data = wfps.dat, mean), all=TRUE)

# Merge everything, then write it out
full.dat <- merge(n2o, wfps.dat, all.x=TRUE, all.y=FALSE) %>%
  merge(minN, all.x=TRUE)

write.table(full.dat, "../data/R outputs/full flux dat.txt", sep = "\t", 
            row.names=FALSE, quote=FALSE)

# Error-check AARS ancillary data ----
# Saw some strange numbers, so I had to investigate. It went pretty well.
# This doesn't need to be rerun, I don't think.
vwc.lm <- lm(vwc ~ factor(date), data = subset(vwc_t, site == "ARL"))
st.lm <- lm(soilT ~ factor(date), data = subset(vwc_t, site == "ARL"))
vwc.check <- cbind(vwc_t %>% filter(!is.na(vwc), site == "ARL"),
              vwc_resid = abs(residuals(vwc.lm)))
temp.check <- cbind(vwc_t %>% filter(!is.na(soilT), site == "ARL"),
                  temp_resid = abs(residuals(st.lm)))
write.table(vwc.check, "../data/vwc check.txt", row.names = FALSE, quote = FALSE)
write.table(temp.check, "../data/temp check.txt", row.names = FALSE, quote = FALSE)
