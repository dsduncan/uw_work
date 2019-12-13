# Cumulative and peak fluxes ----

# Libraries ----
library(dplyr)
library(ggplot2)
library(lme4)
library(lsmeans)
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

# 15-level allegedly dichromat-friendly palette
cf.pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d")
site.pal <- cf.pal[c(11, 2)]

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

# Data entry ----

# Fix date format
# Filter date range
# Remove 2009 ARL Miscanthus (remember the replant!)
# Create some additional factors for useful comparisons
n2o.dat <- read.delim("../data/R outputs/field measurements.txt") %>%
  mutate(date=as.Date(date)) %>% 
  filter(date > "2008-12-31", date < "2015-01-01",
         !(trt=="G06" & site=="ARL" & date<"2010-01-01"))
n2o.dat$avp.fert <- 
  factor(n2o.dat$avp.fert, 
         levels=levels(n2o.dat$avp.fert)[c(1, 11, 2, 4, 8, 6, 9,
                                           12, 3, 5, 7, 10)])
n2o.dat$avp <- 
  factor(n2o.dat$avp,
         levels=levels(n2o.dat$avp)[c(1, 7, 2, 3, 5, 4, 6)])
n2o.dat$site.fert <-
  factor(n2o.dat$site.fert,
         levels=levels(n2o.dat$site.fert)[c(1,3,2,4)])

# 3.1 cumulative annual N2O fluxes ----

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

write.table(n2o.agg.dat, "../data/R outputs/n2o agg.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)

### Figure 2 ###
# Aggregate fluxes

f2<- ggplot(n2o.agg.dat, aes(y=n2o.agg, x=site.fert, fill=site.fert)) +
  geom_boxplot(color="black") + scale_y_log10() +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Cropping system", 
       y=expression(Aggregate~annual~N[2]*O~emissions~(g~N~ha^{-1})),
       fill="Site and fertilization")+
  facet_grid(.~avp, scale="free_x", space="free_x") +
  main.theme +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.margin.x=unit(0, "lines"),
        strip.text=element_text(size=12, color="black"))

ggsave("../figures/Figure 2.tiff", f2, width=7.5, height=4.5, units="in",
       dpi=600, compression="lzw")

### Figure 3 ###
f3.lab <- as_labeller(c(G01="Continuous corn",
                         G02="G02",
                         G03="G03",
                         G04="G04",
                         G05="Switchgrass",
                         G06="Miscanthus",
                         G07="Native grass",
                         G08="Poplar",
                         G09="Old field",
                         G10="Restored prairie"))
f3 <- n2o.agg.dat %>% group_by(site.fert, year, system) %>%
  summarise(med=median(n2o.agg)) %>%
  ungroup %>%
  mutate(year=as.numeric(year)) %>%
  ggplot(aes(x=year, y=med, color=site.fert)) +
  geom_point(data=n2o.agg.dat, aes(y=n2o.agg), alpha=0.3, size=2) + 
  geom_path(size=1) +
  scale_y_log10() + scale_color_brewer(palette="Dark2")+
  facet_wrap(~system, nrow=2, labeller=f3.lab) +
  labs(x="Year",
       y=expression(paste("Aggregate annual ", N[2], "O emissions (g-N ",
                          ha^{-1}, ")", sep="")),
       fill="Site & Fertilization") +
  main.theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="bottom",
        legend.title=element_blank())

ggsave("../figures/Figure 3.tiff", f3, width=9.5, units="in",
       dpi=600, compression="lzw")

# 3.2 Peak fluxes ----

### FIGURE 4 ###
# Dates of peak fluxes for each plot-year, with site & trt data
n2o.peaks <- n2o.dat %>% group_by(plot.year) %>% 
  top_n(1, n2o.ha.day) %>% ungroup %>% data.frame %>%
  mutate(status="Peak")

# When did fertilization events happen?
fert.grp <- read.delim("../data/R inputs/fert timing groups.txt") %>%
  mutate(date=as.Date(date))

levels(fert.grp$fert.time) <- c("Early", "Late", "Poplar")
levels(n2o.peaks$fert.time) <- c("Early", "Late", "Poplar", "Unfertilized")

n2o.peaks %>% 
  mutate(year=factor(year)) %>%
  ggplot(aes(x=date, fill=fert.time)) +
  geom_blank() + 
  geom_vline(aes(color=fert.time, xintercept=as.numeric(date), linetype=fert.time), 
             data=fert.grp, size=1) + 
  geom_histogram(bins=20) + 
  facet_grid(site ~ year, scale="free_x") +
  scale_fill_brewer(palette="Paired") + scale_color_brewer(palette="Paired") + 
  coord_cartesian(ylim=c(0, 60)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(y="Number of peak events",
       x="Sampling date",
       fill="Fertilization timing",
       color="", linetype="") +
  main.theme +
  theme(legend.position="bottom")

ggsave("../figures/Figure 4.tiff", width=10, units="in",
       dpi=600, compression="lzw")

### Figure 5 ###
# Mash things back together so I can sort by status
n2o.peak.plot <- rbind(
  n2o.peaks,
  n2o.dat[!(interaction(select(n2o.dat, plot.year, date)) %in%
              interaction(select(n2o.peaks, plot.year, date))), ] %>%
    mutate(status="Nonpeak"))

f5.xlabs <- levels(n2o.peak.plot$avp.fert)
f1.lab <- as_labeller(c(n2o.ha.day="N[2]*O~flux~'('*g~N~ha^{-1}~day^{-1}*')'",
                        no3="NO[3]^'-'~'('*mu*g~N~g^{-1}*')'",
                        nh4="NH[4]^'+'~'('*mu*g~N~g^{-1}*')'",
                        soilT="Soil~temp~(degree~C)",
                        wfps="WFPS~('%')"),
                      default=label_parsed)

f5 <- n2o.peak.plot %>%
  melt(measure.vars=c("n2o.ha.day", "no3", "nh4", "soilT", "wfps")) %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=avp.fert, y=value, fill=status)) +
  geom_boxplot(position="dodge") +
  geom_hline(yintercept=0) +
  facet_grid(variable ~ site, scales="free", space="free_x", 
             labeller=f1.lab, switch="y") +
  scale_y_continuous(trans="asinh") +
  scale_x_discrete(labels=f5.xlabs) +
  scale_fill_brewer(palette="Set1") +
  labs(y="", fill="",
       x="Cropping system") +
  main.theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x=element_text(margin=margin(t=5)),
        legend.position=c(0.1, -0.18),
        legend.direction="horizontal")

f5.a <- f5 + scale_y_continuous()
f5 <- ggplotGrob(f5)
f5.a <- ggplotGrob(f5.a)

f5$grobs[c(13:14, 18:19, 8:9)] <- f5.a$grobs[c(13:14, 18:19, 8:9)]
grid::grid.draw(f5)

ggsave("../figures/Figure 5.tiff", f5, width=6.5, height=8, units="in",
       dpi=600, compression="lzw")

# Weighted averages, for export
# Data for dates on which at least one plot in the treatment peaked
n2o.peaks.day <- n2o.peaks %>% select(trt, site, date) %>%
  merge(n2o.dat) %>% 
  select(plot.year, site, avp, avp.fert, n2o.ha.day, soilT, wfps, no3, nh4) %>%
  group_by(plot.year, site, avp, avp.fert) %>%
  summarise(n2o.ha.day=mean(n2o.ha.day, na.rm=TRUE),
            no3=mean(no3, na.rm=TRUE),
            nh4=mean(nh4, na.rm=TRUE),
            wfps=mean(wfps, na.rm=TRUE),
            soilT=mean(soilT, na.rm=TRUE)) %>%
  ungroup %>% data.frame
write.table(n2o.peaks.day, "../data/R outputs/n2o.peaks.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
