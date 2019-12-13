# Understanding drivers of N2O fluxes
#
# We have a large dataset of N2O fluxes with accompanying environmental drivers.
# Let's look at N2O flux patterns, and how the drivers relate to these.
# I'm asking 4 major questions here:
# - Do *cumulative annual emissions* differ systematically by site?
# - Do distributions of daily fluxes and drivers differ among systems?
# - Do drivers impose limits on fluxes, and do these constraints differ?
# - Are peak fluxes coinciding with high levels of drivers?
# 
# Note that our VWC data are *bad*. A lot of ARL values are too high to make
# any sense. For general modeling, I guess this isn't terrible, since we mostly
# care about relative differences, but if I try to give any kind of "real" 
# interpretation, I'll probably need to toss in some filters. 
#
# 2016-05-27
# David Duncan

# Libraries ----
library(plyr)  # Need to load plyr before dplyr for masking purposes
library(dplyr)
library(ggplot2)
library(grid)
library(lme4)
library(lsmeans)
library(lubridate)
library(quantreg)
library(reshape2)
library(scales)

# Plotting ----
main.theme <- theme(panel.grid=element_blank(),
                    panel.background=element_rect(color="black", fill="white"),
                    panel.margin.x=unit(0, "lines"),
                    axis.ticks.x=element_blank(),
                    axis.text=element_text(color="black", size=8),
                    axis.title=element_text(size=12),
                    axis.title.x=element_blank(),
                    legend.key=element_blank(),
                    legend.position="bottom",
                    legend.title=element_blank(),
                    strip.background=element_blank(), axis.text.x=element_blank(),
                    strip.text.x=element_text(size=11),
                    plot.margin=unit(c(0.01, 0.02, 0.0, 0.01), "npc"))

# 15-level allegedly dichromat-friendly palette
cf.pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d")

cb.pal <- c("#db6d00","#004949", "#ff6db6","#920000","#490092", "#006ddb", "#924900",
            "#009292", "#ffb6db","#b66dff", "#b6dbff", "#ffff6d")
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
n2o.dat$avp <- factor(n2o.dat$avp, levels=levels(n2o.dat$avp)[c(1:3, 5, 4, 6:7)],
                      labels = c("Annual", "Switchgrass", "Miscanthus", "Poplar", 
                         "Native grass", "Old field", "Prairie"))

# Useful to compare fertilized/unfertilized by site
n2o.dat <- n2o.dat %>% 
  mutate(site.fert=factor(interaction(site, fert), 
                          labels=c("ARL fertilized", "KBS fertilized", 
                                   "ARL unfertilized", "KBS unfertilized")),
         avp.fert=factor(interaction(fert, avp),
                         labels=c("Annual", "Switchgrass fert", "Switchgrass unfrt",
                                  "Miscanthus fert", "Miscanthus unfrt", "Poplar",
                                  "Native grass fert", "Native grass unfrt",
                                  "Old field fert", "Old field unfrt", "Prairie fert",
                                  "Prairie unfrt")))

# Correct fertilization dates for unfertilized treatments
n2o.dat[n2o.dat$fert=="unfrt", "fert.time"] <- "unfrt"

# 3.1 Dataset overview ----
### Table 1 ###
# Looking at how many observations we have for various measurements across the 
# sample period
# N2O
n2o.dat %>% transmute(show=interaction(site, year)) %>% table
# Soil N
n2o.dat %>% filter(!is.na(no3), !is.na(nh4)) %>% 
  transmute(show=interaction(site, year)) %>% table
# Soil Temperature
n2o.dat %>% filter(!is.na(soilT)) %>% 
  transmute(show=interaction(site, year)) %>% table
# WFPS
n2o.dat %>% filter(!is.na(wfps)) %>% 
  transmute(show=interaction(site, year)) %>% table

# 3.3 Distribution of fluxes and parameters ----
# I want to show the distribution of raw values for N2O fluxes and key 
# environmental drivers, to show the extent to which the systems differ.

# First, melt down the data to make an uber table
n2o.mlt <- n2o.dat %>% 
  select(n2o.ha.day, no3, nh4, soilT, wfps, avp, system, site.fert) %>%
  melt %>% filter(!is.na(value))
# Count observations for each point
# Note: not actually used anymore
n2o.cnt <- merge(n2o.mlt %>% count(avp, site.fert, variable),
                 aggregate(value ~ site.fert + avp + variable, 
                           n2o.mlt, median))

### Figure 1 ###
# Fluxes, nitrate by system
f1.lab <- as_labeller(c(n2o.ha.day="N[2]*O~flux~'('*g~N~ha^{-1}~day^{-1}*')'",
                        no3="NO[3]^'-'~'('*mu*g~g^{-1}*')'",
                        nh4="NH[4]^'+'~'('*mu*g~g^{-1}*')'",
                        soilT="Soil~temp~(degree~C)",
                        wfps="WFPS~('%')"),
                      default=label_parsed)

n2o.mlt %>% 
  filter(variable %in% c("n2o.ha.day", "no3", "nh4")) %>% 
  droplevels %>%
  ggplot(aes(y=value, x=site.fert, fill=site.fert)) +
  geom_hline(yintercept=0) +
  geom_boxplot(outlier.size=1, width=1.28) +
  labs(y="") +
  scale_fill_brewer(palette="Dark2") + 
  scale_y_continuous(trans="asinh") +
  facet_grid(variable~avp, scale="free", space="free", switch="both", 
             labeller=labeller(.rows=f1.lab, .cols=label_value)) +
  main.theme

ggsave("../figures/Figure 1.tiff", width=7.5, height=7, units="in",
       dpi=600, compression="lzw")

### Figure S1 ###
# Annual systems aren't different from each other. See?
n2o.mlt %>% 
  filter(variable %in% c("n2o.ha.day", "no3", "nh4"), avp=="Annual") %>% 
  droplevels %>%
  ggplot(aes(x=system, y=value, fill=site.fert), color="black") +
  geom_hline(yintercept=0) +
  geom_boxplot() +
  labs(x="Cropping system", y="") +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(trans="asinh") +
  facet_grid(variable~., scale="free", switch="y", labeller=f1.lab) +
  main.theme +
  theme(axis.text.x=element_text(size=10))

ggsave("../figures/Figure S1.tiff", width=7.5, units="in",
       dpi=600, compression="lzw")

### Figure S2 ###
# Soil temp & WFPS? Who needs 'em...
n2o.mlt %>% filter(variable %in% c("soilT", "wfps")) %>% 
  droplevels %>%
  ggplot(aes(y=value, x=site.fert, fill=site.fert)) +
  geom_hline(yintercept=0) +
  geom_boxplot(outlier.size=1) +
  labs(y="") +
  scale_fill_brewer(palette="Dark2") + 
  facet_grid(variable~avp, scale="free", space="free_x", switch="both", 
             labeller=labeller(.rows=f1.lab, .cols=label_value)) +
  main.theme
ggsave("../figures/Figure S2.tiff", width=7.5, units="in",
       dpi=600, compression="lzw")

# Are NO3 and NH4 concentrations related?

### Figure S3 ###
# NO3 vs NH4 
n2o.dat %>%
  mutate(fert=factor(fert, labels=c("Fertilized", "Unfertilized"))) %>%
  filter(!is.na(nh4)) %>%
  ggplot(aes(x=no3, y=nh4, color = avp)) +
  geom_point(alpha=0.25) + stat_smooth(se=FALSE) + 
  scale_color_brewer(palette="Dark2") +
  facet_grid(fert ~ site)+
  scale_x_continuous(trans="asinh", expand=c(0,0)) +
  scale_y_continuous(trans="asinh", expand=c(0,0)) +
  labs(color="Cropping system",
       x=expression(NO[3]^'-'~'('*mu*g~g^{-1}*')'),
       y=expression(NH[4]^'+'~'('*mu*g~g^{-1}*')')) +
  main.theme +
  theme(aspect.ratio = 1,
        axis.title.x=element_text(vjust=-0.3),
        axis.text.x=element_text(),
        axis.ticks.x=element_line(),
        panel.margin.x=unit(0.3, "lines"),
        panel.margin.y=unit(0.3, "lines"))
ggsave("../figures/Figure S3.tiff", units="in", height=6.5, width=6,
       dpi=600, compression="lzw")

# 3.4 Environmental constraints on N2O fluxes ----

# QR for Nitrate + Ammonium
# Does adding ammonium data improve the nitrate model at ARL?
na.a.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3) + asinh(nh4), tau=0.95,
           data=filter(n2o.dat, site=="ARL"))
summary(na.a.rq)
n.a.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3), tau=0.95,
           data=filter(n2o.dat, site=="ARL"))
anova(na.a.rq, n.a.rq)

# Does adding ammonium data improve the nitrate model at KBS?
na.k.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3) + asinh(nh4), tau=0.95,
              data=filter(n2o.dat, site=="KBS"))
summary(na.k.rq)
n.k.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3), tau=0.95,
             data=filter(n2o.dat, site=="KBS"))
anova(na.k.rq, n.k.rq)

# Melt the data so I can loop over everything and make a nice table
qr.dat <- melt(n2o.dat, idvars=c("site", "avp.fert"), 
               measure.vars=c("no3", "nh4", "wfps", "soilT")) %>%
  filter(!is.na(value))

# Transform predictors like this
t.preds <- c("asinh(no3)", "asinh(nh4)", "wfps", "soilT")
# Table to hold things
qr.coefs <- c()

for(t.i in seq_along(t.preds)){
  # Identify variable being used 
  t.var <- t.preds[t.i]
  # Make a formula with that variable
   t.form <- as.formula(paste("asinh(n2o.ha.day) ~", t.var))
   
   for(t.j in seq_along(levels(n2o.dat$site))){
    t.site <- levels(n2o.dat$site)[t.j]
    # Subset the data
    t.qr <- filter(n2o.dat, site==t.site)
    
    # QR, model
    t.rq <- rq(t.form, tau=0.95, data=t.qr)
    # QR, null
    t.null <- rq(asinh(n2o.ha.day) ~ 1, data=t.qr)
    # Summary (for parameters)
    rq.sum <- summary(t.rq)
    
    R1 <- 1 - (t.rq$rho/t.null$rho)

    # Unconventional dataframe... Why don't we do it like this more often?
    data.frame(site=t.site, variable=t.var, int=rq.sum$coefficients[1,1],
          slope=rq.sum$coefficients[2,1], slope.p=rq.sum$coefficients[2,4],
          gof=R1) %>%
      rbind(qr.coefs) -> qr.coefs
  }
}
qr.coefs <- mutate(qr.coefs,
                   variable=revalue(variable, c("asinh(no3)"="no3",
                                        "asinh(nh4)"="nh4")))
rm(list=ls(pattern="t\\."))

### Figure 2 ###
f2.a <- qr.dat %>%
  ggplot(aes(x=value, y=n2o.ha.day)) +
  geom_hline(yintercept=0) +
  geom_point(aes(color=site), alpha=0.08, size=3) +
  geom_abline(aes(slope=slope, intercept=int), data=qr.coefs) +
  facet_grid(site~variable, scales="free_x", switch="x",
             labeller=labeller(.rows=label_value, .cols=f1.lab)) +
  scale_color_manual(values=site.pal) +
  scale_y_continuous(trans="asinh") +
  labs(y=expression(g~N[2]*O-N~flux~ha^{-1}~day^{-1})) +
  main.theme +
  theme(aspect.ratio=1,
        axis.text.x=element_text(),
        axis.ticks.x=element_line(),
        legend.position="none",
        panel.margin.y=unit(2, "lines"))

f2.b <- f2.a +
  scale_x_continuous(trans="asinh")
f2.a <- ggplotGrob(f2.a)
f2.b <- ggplotGrob(f2.b)

f2.a$grobs[c(4:7, 14:15)] <- f2.b$grobs[c(4:7, 14:15)]
grid::grid.draw(f2.a)
ggsave("../figures/Figure 2.tiff", f2.a, width=10.5, units="in",
       dpi=600, compression="lzw")
ggsave("../figures/Figure 2.pdf", f2.a, width=10.5, units="in")
,
       dpi=600, compression="lzw")

# Using AIC to look whether separate slopes/intercepts by system (crop + 
# fertilization) improves the model over a universal fit

# ARL: Nitrate
n.a.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3), tau=0.95,
             data=filter(n2o.dat, site=="ARL"))
n.a.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+asinh(no3), tau=0.95,
               data=filter(n2o.dat, site=="ARL"))
n.a.ss.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3) + avp.fert:asinh(no3), tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
n.a.sb.rq <- rq(asinh(n2o.ha.day) ~ fert*asinh(no3), tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
AIC(n.a.rq); AIC(n.a.si.rq); AIC(n.a.ss.rq); AIC(n.a.sb.rq)
anova(n.a.rq, n.a.si.rq, n.a.sb.rq)

# ARL: Ammonium
a.a.rq <- rq(asinh(n2o.ha.day) ~ asinh(nh4), tau=0.95,
             data=filter(n2o.dat, site=="ARL"))
a.a.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+asinh(nh4), tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
a.a.ss.rq <- rq(asinh(n2o.ha.day) ~ asinh(nh4) + avp.fert:asinh(nh4), tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
a.a.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*asinh(nh4), tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
AIC(a.a.rq); AIC(a.a.si.rq); AIC(a.a.ss.rq); AIC(a.a.sb.rq)
anova(a.a.rq, a.a.si.rq, a.a.sb.rq)

# ARL: WFPS
w.a.rq <- rq(asinh(n2o.ha.day) ~ wfps, tau=0.95,
             data=filter(n2o.dat, site=="ARL"))
w.a.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+wfps, tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
w.a.ss.rq <- rq(asinh(n2o.ha.day) ~ wfps + avp.fert:wfps, tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
w.a.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*wfps, tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
AIC(w.a.rq); AIC(w.a.si.rq); AIC(w.a.ss.rq); AIC(w.a.sb.rq)
anova(w.a.rq, w.a.si.rq, w.a.sb.rq)

# ARL: SoilT
t.a.rq <- rq(asinh(n2o.ha.day) ~ soilT, tau=0.95,
             data=filter(n2o.dat, site=="ARL"))
t.a.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+soilT, tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
t.a.ss.rq <- rq(asinh(n2o.ha.day) ~ soilT + avp.fert:soilT, tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
t.a.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*soilT, tau=0.95,
                data=filter(n2o.dat, site=="ARL"))
AIC(t.a.rq); AIC(t.a.si.rq); AIC(t.a.ss.rq); AIC(t.a.sb.rq)
anova(t.a.rq, t.a.si.rq, t.a.sb.rq)

# KBS: Nitrate
n.k.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3), tau=0.95,
             data=filter(n2o.dat, site=="KBS"))
n.k.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+asinh(no3), tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
n.k.ss.rq <- rq(asinh(n2o.ha.day) ~ asinh(no3) + avp.fert:asinh(no3), tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
n.k.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*asinh(no3), tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
AIC(n.k.rq); AIC(n.k.si.rq); AIC(n.k.ss.rq); AIC(n.k.sb.rq)
anova(n.k.rq, n.k.si.rq, n.k.sb.rq)

# KBS: Ammonium
a.k.rq <- rq(asinh(n2o.ha.day) ~ asinh(nh4), tau=0.95,
             data=filter(n2o.dat, site=="KBS"))
a.k.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+asinh(nh4), tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
a.k.ss.rq <- rq(asinh(n2o.ha.day) ~ asinh(nh4) + avp.fert:asinh(nh4), tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
a.k.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*asinh(nh4), tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
AIC(a.k.rq); AIC(a.k.si.rq); AIC(a.k.ss.rq); AIC(a.k.sb.rq)
anova(a.k.rq, a.k.si.rq, a.k.sb.rq)

# KBS: WFPS
w.k.rq <- rq(asinh(n2o.ha.day) ~ wfps, tau=0.95,
             data=filter(n2o.dat, site=="KBS"))
w.k.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+wfps, tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
w.k.ss.rq <- rq(asinh(n2o.ha.day) ~ wfps + avp.fert:wfps, tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
w.k.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*wfps, tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
AIC(w.k.rq); AIC(w.k.si.rq); AIC(w.k.ss.rq); AIC(w.k.sb.rq)
anova(w.k.rq, w.k.si.rq, w.k.sb.rq)

# KBS: SoilT
t.k.rq <- rq(asinh(n2o.ha.day) ~ soilT, tau=0.95,
             data=filter(n2o.dat, site=="KBS"))
t.k.si.rq <- rq(asinh(n2o.ha.day) ~ avp.fert+soilT, tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
t.k.ss.rq <- rq(asinh(n2o.ha.day) ~ soilT + avp.fert:soilT, tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
t.k.sb.rq <- rq(asinh(n2o.ha.day) ~ avp.fert*soilT, tau=0.95,
                data=filter(n2o.dat, site=="KBS"))
AIC(t.k.rq); AIC(t.k.si.rq); AIC(t.k.ss.rq); AIC(t.k.sb.rq)
anova(t.k.rq, t.k.si.rq, t.k.sb.rq)
summary(t.k.ss.rq)

# Distribution of coefficients
# Making a loop to go over all (yes all) models so I can graph it later

# We're making a generalized formula, because we're cool like that!
t.variable <- c("asinh(no3)", "asinh(nh4)", "wfps", "soilT")
qr.crop.coefs <- data.frame()

for(t.i in seq_along(t.variable)){
  # Generalized formula, allows us to cycle predictors
  # Note that we've suppressed intercepts so coefficients are treatment means
    t.form <- formula(paste("asinh(n2o.ha.day)~avp.fert-1+avp.fert:", 
                          t.variable[t.i]))
    t.form.p <- formula(paste("asinh(n2o.ha.day)~avp.fert*", t.variable[t.i]))
    for(t.j in seq_along(levels(n2o.dat$site))){
      t.dat <- filter(n2o.dat, site==levels(n2o.dat$site)[t.j])
      
      # Parameter + upper and lower confidence intervals
      t.rq <- rq(t.form, tau=0.95, data=t.dat) %>% 
        summary(se="ker") %>% coefficients %>% data.frame %>% 
        select(1:2)
      # P-values for comparison to annual systems
      t.pva <- (rq(t.form.p, tau=0.95, data=t.dat) %>%
        summary(se="ker") %>% coefficients)[,4] 

      # Dataframe to hold output (includes parameter, whether it's a slope,
      # site, avp.fert)
      # Note the clever regex to re-extract cropfert from the table
      t.rq.n <- nrow(t.rq)/2
      t.pva.n <- length(t.pva)/2
      qr.crop.coefs <- 
        data.frame(intercept=t.rq$Value[1:t.rq.n],
                   intercept.se=t.rq$Std..Error[1:t.rq.n],
                   intercept.p=t.pva[1:t.pva.n],
                   slope=t.rq$Value[(t.rq.n+1):(2*t.rq.n)],
                   slope.se=t.rq$Std..Error[(t.rq.n+1):(2*t.rq.n)],
                   slope.p=t.pva[(t.pva.n+1):(2*t.pva.n)]) %>%
        mutate(site=levels(n2o.dat$site)[t.j],
               variable=factor(t.variable[t.i], levels=t.variable),
               avp.fert=factor(gsub("avp.fert([[:alpha:]//. ]+)","\\1",
                 row.names(t.rq)[1:t.rq.n]), 
                 levels=levels(n2o.dat$avp.fert))) %>%
        rbind(qr.crop.coefs)
      print(c(t.i, t.j))
    }
}
rm(list=ls(pattern="t\\."))

f3.lab <- as_labeller(c("asinh(no3)"="'['*NO[3]^'-'*']'",
                        "asinh(nh4)"="'['*NH[4]^'+'*']'",
                        wfps="WFPS",
                        soilT="Soil~temp"),
                      default=label_parsed)

levels(qr.crop.coefs$avp.fert) <- gsub("\\.", " ", levels(qr.crop.coefs$avp.fert))

### Figure 3 ###
# Slopes
ggplot(qr.crop.coefs, aes(x=avp.fert, y=slope, fill=avp.fert)) +
  geom_crossbar(aes(ymin=slope-slope.se, ymax=slope+slope.se), 
                position="dodge") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=cf.pal[-1]) + 
  facet_grid(variable~site, scales="free", space="free_x",
             labeller=labeller(.rows=f3.lab)) +
  labs(x="Cropping system",
       y="Quantile regression slope") +
  main.theme +
  theme(legend.position="none",
        axis.title.x=element_text(),
        axis.text.x=element_text(angle=45, size=10, vjust=1, hjust=1),
        plot.margin=unit(c(0.01, 0.02, 0.02, 0.01), "npc"))
ggsave("../figures/Figure 3.tiff", width=7.5, height=8, units="in",
       dpi=600, compression="lzw")

### Figure 4 ###
# Intercepts
ggplot(qr.crop.coefs, aes(x=avp.fert, y=sinh(intercept), fill=avp.fert)) +
  geom_crossbar(aes(ymin=sinh(intercept-intercept.se), 
                    ymax=sinh(intercept+intercept.se)), position="dodge") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=cf.pal[-1]) + 
  scale_y_continuous(trans="asinh") +
  facet_grid(variable~site, scales="free", space="free_x",
             labeller=labeller(.rows=f3.lab)) +
  labs(x="Cropping system",
       y=expression(Intercept~'('*g~N[2]*O-N~ha^{-1}~
                      day^{-1}*')')) +
  main.theme +
  theme(legend.position="none",
        axis.title.x=element_text(),
        axis.text.x=element_text(angle=45, size=10, vjust=1, hjust=1),
        plot.margin=unit(c(0.01, 0.02, 0.02, 0.01), "npc"))
ggsave("../figures/Figure 4.tiff", width=7.5, height=8, units="in",
       dpi=600, compression="lzw")
