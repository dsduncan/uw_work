# Data entry ----

# Fix date format
# Filter date range
# Remove 2009 ARL Miscanthus (remember the replant!)
# Create some additional factors for useful comparisons
n2o.dat <- read.delim("../data/R outputs/full flux dat.txt") %>%
  mutate(date=as.Date(date)) %>% 
  filter(date>"2008-12-31", !(trt=="G06" & site=="ARL" & date<"2010-01-01")) %>%
  mutate(year=year(date), plot=factor(interaction(trt, block)),
         system=factor(sub("M", "", trt)),
         fert=trt, avp=system) %>%
  merge(read.delim("../data/R inputs/fert dates.txt")) 

# Were treatments fertilized (note this doesn't mean fertilized in any given
# year, only fertilized at all)
levels(n2o.dat$fert) <- list(fert=levels(n2o.dat$trt)[c(1:9,11,13,15,16,19)],
                             unfrt=levels(n2o.dat$trt)[c(10,12,14,17,18)])
# Lumping annual vs perennial because I can't even with the rotations
levels(n2o.dat$avp) <- c(rep("annual", 4), levels(n2o.dat$system)[5:10])
n2o.dat$avp <- factor(n2o.dat$avp, 
                      labels=c("Annual", "Switchgrass", "Miscanthus", 
                               "Native grass mix", "Poplar", "Old field",
                                  "Prairie"))
# Useful to compare fertilized/unfertilized by site
n2o.dat <- n2o.dat %>% 
  mutate(site.fert=factor(interaction(site, fert), 
                          labels=c("ARL fertilized", "KBS fertilized", 
                                   "ARL unfertilized", "KBS unfertilized")),
         avp.fert=factor(interaction(avp, fert)))

# Correct fertilization dates for unfertilized treatments
n2o.dat[n2o.dat$fert=="unfrt", "fert.time"] <- "unfrt"

# Buddhist dataset
full.dat <- n2o.dat %>% 
  filter(!is.na(soilT), !is.na(wfps), !is.na(no3), !is.na(nh4))

# 3.2 cumulative annual N2O fluxes ----

# Assume that fluxes are 0 while soil is frozen (incorrect, but what to do?)
# Read in a list of first freeze dates before/after sampling season 
frz_days <- read.delim("../data/R inputs/Freeze dates.txt") %>%
  mutate(date=as.Date(date), 
         n2o.ha.day=0)

# Aggregate fluxes at a yearly basis
# Key assumptions:
# - Linear interpolation between observations (no gapfill)
# - Fluxes are 0 when soils are frozen (deny reality)
# - Calendar year (Jan 1 - Jan 1)
n2o.agg.dat <- c()

for(t.i in seq_along(levels(n2o.dat$plot))){
  # Subset out a single plot
  t.dat <- filter(n2o.dat, plot==levels(n2o.dat$plot)[t.i])
  # We only really care about fluxes and dates
  t.min_dat <- select(t.dat, date, n2o.ha.day)
  # Add in the freeze dates
  t.min_dat <- rbind(t.min_dat,
                     frz_days %>% 
                       filter(site==t.dat$site[1], 
                              year(date) %in% unique(t.dat$year)) %>%
                       select(-site)) %>%
    arrange(date)
  # Make the ends of years
  t.strt <- data.frame(date= as.Date(paste0(min(t.dat$year), "-01-01")) +
                         years(0:(max(t.dat$year)-min(t.dat$year))),
                       n2o.ha.day=NA)
  # Merge in ends of years
  t.min_dat <- rbind(t.min_dat, t.strt) %>% arrange(date)
  
  # Indices of start dates
  t.start_ind <- which(is.na(t.min_dat[-1, "n2o.ha.day"])) + 1
  # End dates are weighted averages of their adjoining dates
  t.min_dat[t.start_ind, 2] <- t.min_dat[t.start_ind-1, 2] +
    (t.min_dat[t.start_ind+1,2]-t.min_dat[t.start_ind-1,2])/
    as.numeric(t.min_dat[t.start_ind+1,1]-t.min_dat[t.start_ind-1,1])
  # Exception in case the first date has nothing going for it
  if(is.na(t.min_dat[1,2])){t.min_dat[1,2] <- 0}

  # Identify the years for which we have data
  t.years <- setdiff(year(t.min_dat$date), "2015")
  
  # Make a dataframe to hold the output
  t.out <- data.frame(plot=t.dat[1,"plot"],
                      year=t.years,
                      n2o.agg=0)
  # Loop over years and aggregate
  for(t.j in seq_along(t.years)){
    t.min_dat.sub <- 
      filter(t.min_dat, date >= as.Date("0-01-01")+years(t.years[t.j]),
             date <= as.Date("0-01-01")+years(t.years[t.j]+1))
    t.n <- 2:nrow(t.min_dat.sub)
    t.out[t.j, "n2o.agg"] <- 
      with(t.min_dat.sub, as.numeric(date[t.n]-date[t.n-1]) %*%
                          (n2o.ha.day[t.n]+n2o.ha.day[t.n-1])/2)
    }
  n2o.agg.dat <- rbind(n2o.agg.dat, t.out)
}
rm(list=ls(pattern="t\\."))

# Merge flux data with ancillary data
n2o.agg.dat <- n2o.dat %>% select(site:block, plot:avp.fert) %>% 
  unique %>% merge(n2o.agg.dat) %>% mutate(year=factor(year))

### Figure 1 ###
ggplot(n2o.agg.dat, aes(y=n2o.agg, x=avp, fill=site.fert)) +
  geom_boxplot(color="black") + scale_y_log10() +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Cropping system", 
       y=expression(paste("Cumulative annual ", N[2], "O emissions (",
                          "g-N ", ha^{-1}, ")", sep="")),
       fill="Site and fertilization")+
  main.theme +
  theme(legend.position="bottom")

ggsave("../figures/Figure 1.png", width=7.5, height=4.5, units="in")

### Figure S1 ###
n2o.agg.dat %>% 
  subset((trt=="G10" & fert=="unfrt") |
           (trt == "G01") |
           (trt %in% c("G05", "G06", "G07", "G08") & fert=="fert")) %>%
  droplevels %>%
  mutate(trt=factor(trt, labels=c("Corn", "Switchgrass", "Miscanthus", 
                    "Native grass mix", "Poplar", "Prairie"))) %>%
  ggplot(aes(y=n2o.agg, x=year, color=site)) +
  geom_point(size=3, position=position_dodge(width=0.4), alpha=0.5) +
  labs(x="Year", 
       y=expression(paste("Cumulative annual   ",N[2], "O emissions (",
                          "g-N ", ha^{-1}, ")", sep="")),
       color="Site")+
  scale_y_log10() + scale_color_brewer(palette="Dark2")+
  facet_wrap(~trt, nrow=3)+
  main.theme 
ggsave("../figures/Figure S1.png", width=7.5, units="in")
  