# Yet another attempt at makint rStan work for me...
# TODO:
# 1. Load data
# 2. Form stan-like data structure
# 3. Try basic stanscript
# 4. ?
# 5. Profit!

load("G:/post Stan.RData")

# Libraries ----
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(rstan)
library(shinystan)

# Functions ----
stanicize <- function(df, site.pick){
  # Given a site to subset by, create a list that will serve as Stan input data
  # This is only for the NOE thing...
  require(dplyr)
  t.sub <- df %>% filter(site==site.pick) %>% droplevels
  model.dat <- list(N=nrow(t.sub),  # Number of observations
                    P=t.sub$plot.year %>% unique %>% length, # Number of plots
                    M=model.matrix(n2o.ha.day ~ plot.year - 1, t.sub),
                    n2o=t.sub$n2o.ha.day,
                    no3=t.sub$no3,
                    nh4=t.sub$nh4,
                    wfps=sapply(t.sub$wfps, 
                                FUN=function(x){
                                  return(max(min(x/100, 0.9999), 0.0001))}),
                    soilT=t.sub$soilT)
  return(model.dat)
}

quantFind <- function(df, tsite, quantile, target){
  # Given a dataframe and a site to subset by, return the plot whose target
  # is the specified quantile
  t.df <- filter(df, site==tsite)
  t.q <- quantile(t.df[,target], quantile)
  t.names <- c()
  for(t.i in seq_along(t.q)){
    t.names[t.i] <- t.df[which.min(abs(t.df[,target] - t.q[t.i])),
                         "plot.year"] %>% as.character
  }
  names(t.names) <- quantile
  t.names <- data.frame(plot.year=t.names, quantile=factor(quantile))
  return(t.names)
}

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
# 2 colors, for ARL vs KBS
site.pal <- cf.pal[c(11, 2)]

# Data import ----
# Remember, we only want samples with *all* of the data and positive fluxes
field <- read.delim("../data/R outputs/field measurements.txt") %>%
  filter(!is.na(soilT), !is.na(wfps), !is.na(no3), n2o.ha.day > 0) %>% 
  droplevels %>%
  mutate(year=factor(year),
         date=as.Date(date),
         avp.fert=factor(avp.fert, 
                         levels=c("Annual fertilized", "Switchgrass fertilized",
                                  "Miscanthus fertilized", "Native grass fertilized",
                                  "Poplar fertilized", "Old field fertilized",
                                  "Prairie fertilized", "Switchgrass unfertilized",
                                  "Miscanthus unfertilized", "Native grass unfertilized",
                                  "Old field unfertilized", "Prairie unfertilized")))

# Further, we want at least 5 observations for a plot year (for stats!)
t.keep <- (field %>% select(plot.year) %>% 
  table %>% data.frame %>%
  filter(Freq >= 5))[,1] %>% as.character 

field <- field %>% filter(plot.year %in% t.keep) %>% droplevels %>%
  arrange(plot.year, date)

# Stanalysis ----
# Setup input data
arl.standat <- stanicize(field, "ARL")
kbs.standat <- stanicize(field, "KBS")

# Actually run thing
# NOTE: Need to run kbs.stan FIRST 
kbs.stan <- stan(file="noe.stan", data=kbs.standat, iter=1000, warmup=500, 
                 chains=1)
arl.stan <- stan(file="noe.stan", data=arl.standat, iter=200, warmup=50, 
                 chains=4, fit=kbs.stan)

# 
pairs(kbs.stan, pars=c("DKm", "DWa", "DQ10", "NKm", "NWa", "NWb", "NQ10"))
pairs(arl.stan, pars=c("DKm", "DWa", "DQ10", "NKm", "NWa", "NWb", "NQ10"))

print(arl.stan, pars=c("DKm", "DWa", "DWb", "DQ10", "NKm", "NWa", "NWb", "NQ10"))

# Evaluate site-level parameters ----
# Not actually used...
stan.sum <- rbind(summary(arl.stan,pars=c("DKm", "DWa", "DQ10", "NKm", "NWa", "NWb", 
                                   "NQ10"))[[1]] %>% data.frame %>%
                    mutate(site="ARL", varname=row.names(.)),
                  summary(kbs.stan,pars=c("DKm", "DWa", "DQ10", "NKm", "NWa", "NWb", 
                                   "NQ10"))[[1]] %>% data.frame %>%
                    mutate(site="KBS", varname=row.names(.))) %>%
  select(site, varname, mean, se_mean) %>%
  arrange(varname, site)

### Figure 6 (left) ###
# Posterior distributions of parameters
# Make a dataframe full of parameter estimates, for density plots
t.ks <- extract(kbs.stan, pars=c(unique(stan.sum$varname)))
t.ar <- extract(arl.stan, pars=c(unique(stan.sum$varname)))

pars.dist <- rbind(
  data.frame(par=rep(names(t.ar), each=length(t.ar[[1]])), site="ARL",
             values=unlist(t.ar)),
  data.frame(par=rep(names(t.ks), each=length(t.ks[[1]])), site="KBS",
             values=unlist(t.ks)))

f6.l.lab <- pars.dist %>% group_by(par) %>%
  summarise(min=floor(min(values))) %>% ungroup %>% data.frame %>%
  cbind(data.frame(par=unique(stan.sum$varname),
                   max=c(2100, 650, 2500, 6600, 650, 2500, 5500)))

f6.l <- ggplot(pars.dist) + 
  geom_histogram(aes(x=values, fill=site), alpha=0.8, position="identity", 
                 bins=50) +
  geom_text(data=f6.l.lab, aes(x=min, y=max, label=par), vjust=0.8, hjust=0,
            fontface="italic") +
  scale_fill_manual(values=site.pal) +
  facet_wrap(~par, ncol=1, scales="free") +
  labs(x="Model parameter value", y="Frequency") +
  main.theme +
  theme(strip.text=element_blank(),
        legend.position="none")

### Figure 6 (right) ###
# 
yvals <- c("Dn", "Dt", "Dw", "Nn", "Nt", "Nw")
xvals <- c("no3", "soilT", "wfps", "nh4", "soilT", "wfps")

pred.plot <- data.frame(site=factor(),
                        pred=factor(),
                        xvals=vector("numeric"),
                        yvals=vector("numeric"))

for(t.i in seq_along(yvals)){
  t.a <- extract(arl.stan, pars=yvals[t.i])
  pred.plot <- rbind(pred.plot,
    data.frame(site="ARL", pred=yvals[t.i],
               xvals=arl.standat[[xvals[t.i]]],
               yvals=exp(colMeans(t.a[[1]])))) %>%
    unique
  t.a <- extract(kbs.stan, pars=yvals[t.i])
  pred.plot <- rbind(pred.plot,
                     data.frame(site="KBS", pred=yvals[t.i],
                                xvals=kbs.standat[[xvals[t.i]]],
                                yvals=exp(colMeans(t.a[[1]])))) %>%
    unique
}

f6.r.lab <-pred.plot %>% group_by(pred) %>%
  summarise(xvals=floor(min(xvals)),
           yvals=ceiling(max(yvals))) %>%
  ungroup %>% data.frame

f6.r <- ggplot(pred.plot, aes(x=xvals, y=yvals)) +
  geom_line(aes(color=site), size=1) +
  geom_text(data=f6.r.lab, aes(label=pred), vjust=0.8, fontface="italic") + 
  scale_color_manual(values=site.pal) +
  labs(x="Environmental variable value", y="Multiplier", color="Site") +
  facet_wrap(~pred, nrow=7, drop=FALSE, scales="free") +
  main.theme + 
  guides(color=guide_legend(override.aes=list(size=4))) +
  theme(strip.text=element_blank(),
        legend.position="bottom")

f6.leg <- ggplotGrob(f6.r)$grobs[[28]]

f6 <- grid.arrange(f6.l, (f6.r+theme(legend.position="none")), f6.leg, 
                   ncol=2, 
                   layout_matrix=cbind(rep(1, 7), c(rep(2, 6), 3)))
ggsave("../figures/Figure 6.tiff", f6, width=7.5, height=10.5, units="in",
       dpi=600, compression="lzw")

# Evaluate model fit ----
# Set up a table with all of the model estimates
stan.est <- rbind(filter(field, site=="ARL"), filter(field, site=="KBS")) %>%
  select(plot.year, n2o.ha.day, date) %>% 
  cbind(rbind(summary(arl.stan, pars="est")[[1]] %>% 
                data.frame %>% mutate(site="ARL"),
              summary(kbs.stan, pars="est")[[1]] %>% 
                data.frame %>% mutate(site="KBS")))

# Create data frame to hold output of a loop
stan.rmse <- cbind(stan.est %>% select(plot.year, site) %>% unique,
                   rmse=0)

# Loop by plot-year to calculate rmse
for(t.i in seq_along(stan.rmse$plot.year)){
  t.sub <- filter(stan.est, plot.year==stan.rmse$plot.year[t.i], 
                  site==stan.rmse$site[t.i])
  stan.rmse[t.i, "rmse"] <- 
    sqrt(sum((t.sub$n2o.ha.day - t.sub$mean)^2)/
           nrow(t.sub))
}

stan.eval <-  rbind(
  merge(stan.est, quantFind(stan.rmse, "ARL", c(0.1, 0.5, 0.9), "rmse")),
  merge(stan.est, quantFind(stan.rmse, "KBS", c(0.1, 0.5, 0.9), "rmse")))

### Figure 7 ###
ggplot(stan.eval, aes(x=date)) +
  geom_point(aes(y=n2o.ha.day), color="red") + 
  geom_line(aes(y=n2o.ha.day), color="red", size=1) +
  geom_point(aes(y=mean), color="black", size=3)+ 
  labs(x="Date",
       y=expression(N[2]*O~flux~(g-N~ha^{-1}~day^{-1}))) +
  scale_y_continuous(expand=c(0.02,0.02)) +
  facet_wrap(~ site + quantile, nrow=2, scales="free", 
             labeller=label_wrap_gen(multi_line=FALSE)) +
  main.theme
ggsave("../figures/Figure 7.tiff",
       dpi=600, compression="lzw")

# Potential rates ----
stan.pot <- arrange(field, site) %>%
  select(plot.year, block, avp.fert) %>% unique %>%
  cbind(rbind(summary(arl.stan, pars="PDR")[[1]] %>% 
                data.frame %>% mutate(site="ARL", pars="PDR"),
              summary(kbs.stan, pars="PDR")[[1]] %>% 
                data.frame %>% mutate(site="KBS", pars="PDR"),
              summary(arl.stan, pars="PNR")[[1]] %>% 
                data.frame %>% mutate(site="ARL", pars="PNR"),
              summary(kbs.stan, pars="PNR")[[1]] %>% 
                data.frame %>% mutate(site="KBS", pars="PNR"))) %>%
  arrange(pars, site, avp.fert, X50.) %>%
  mutate(ranks=factor(as.numeric(row.names(.))))
levels(stan.pot$ranks) <- rep(stan.pot$ranks[1:(nrow(stan.pot)/2)],2)

### Figure 8 ###
# Posterior estimates of N2O from de/nitrification, by site and treatment
stan.pot %>%
  mutate(pars=factor(pars, labels=c("Denitrification", "Nitrification"))) %>%
  ggplot(aes(x=ranks, fill=avp.fert)) +
  geom_boxplot(stat="identity", position=position_dodge(width=1),
               aes(lower=X25., middle=X50., upper=X75., ymax=X97.5., 
                   ymin=X2.5.), lwd=0.25) +
  facet_grid(pars ~ site, scales="free", space="free_x") +
  ylab(expression(Potential~N[2]*O~flux~(g~N~ha^{-1}~day^{-1}))) +
  scale_fill_manual(values=cf.pal) +
  scale_y_continuous(expand=c(0,0)) +
  labs(fill="System and fertilization") +
  main.theme +
  guides(fill=guide_legend(override.aes=list(color=NA))) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../figures/Figure 8.tiff", width=9.5, height=5.75, units="in",
       dpi=600, compression="lzw")

# Export data ----
stan.pot %>%
  select(plot.year:avp.fert, X50., site, pars) %>%
  reshape(v.names="X50.", idvar="plot.year", timevar="pars", direction="wide") %>%
  rename(PDR=X50..PDR, PNR=X50..PNR) %>%
  write.table("../data/R outputs/N2O potentials from RStan.txt",
            sep="\t", row.names=FALSE, quote=FALSE)
