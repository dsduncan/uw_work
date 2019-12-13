# Plotting ----
main.theme <- theme(panel.grid=element_blank(),
                    panel.background=element_rect(color="black", fill="white"),
                    axis.text.y=element_text(color="black", size=12),
                    axis.title.y=element_text(size=16),
                    axis.text.x=element_text(color="black", size=16, angle=45,
                                             hjust=1, vjust=1),
                    axis.title.x=element_blank(),
                    legend.position="none",
                    strip.background=element_blank(), 
                    strip.text=element_text(size=16))

n2o.dat <- filter(n2o.dat,
                  !(site=="KBS" & avp.fert %in% c("Switchgrass unfrt", "Miscanthus unfrt", 
                                                  "Old field unfrt"))) %>%
  droplevels


frz_days <- read.delim("../data/R inputs/Freeze dates.txt") %>%
  mutate(date=as.Date(date), 
         n2o.ha.day=0)

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
#rm(list=ls(pattern="t\\."))

# Merge flux data with ancillary data
n2o.agg.dat <- n2o.dat %>% select(site:block, plot:avp.fert) %>% 
  unique %>% merge(n2o.agg.dat) %>% mutate(year=factor(year))

# Agg flux plot ----

ggplot(n2o.agg.dat, aes(y=n2o.agg, x=avp, fill=avp, alpha=fert)) +
  geom_boxplot(color="black", width=1.1) + scale_y_log10() +
  geom_vline(xintercept=which(levels(n2o.agg.dat$avp)=="Switchgrass")-0.5) +
  geom_vline(xintercept=which(levels(n2o.agg.dat$avp)=="Native grass")-0.5, 
             linetype=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_alpha_discrete(range=c(1, 0.5), guide=FALSE) +
  labs(x="Cropping system", 
       y=expression(Aggregate~annual~N[2]*O~emissions~(g~N~ha^{-1})),
       fill="Site and fertilization")+
  facet_grid(.~site, scale="free_x", space="free_x") +
  guides(fill=guide_legend(nrow=1)) +
  main.theme

ggsave("C://My Program Data/Dropbox/Articles/Dissertation/images/aggflux.png",
       height=5.6, width=9.5, units="in")

# Agg flux by year ----
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
n2o.agg.dat %>% group_by(site.fert, year, system) %>%
  filter((system %in% c("G01", "G05", "G06", "G07", "G08", "G09")), fert=="fert") %>%
  summarise(med=median(n2o.agg)) %>%
  ungroup %>%
  mutate(year=as.numeric(year)) %>%
  ggplot(aes(x=year, y=med, color=site.fert)) +
  geom_point(data=filter(n2o.agg.dat, (system %in% c("G01", "G05", "G06", "G07", "G08", "G09")), fert=="fert"),
             aes(y=n2o.agg), alpha=0.7, size=5, position=position_dodge(width=0.2)) + 
  geom_path(size=1) +
  scale_y_log10(breaks=c(100, 500, 1000, 5000)) + 
  scale_color_manual(values=site.pal)+
  facet_wrap(~system, nrow=2, labeller=f3.lab) +
  labs(x="Year",
       y=expression(paste("Aggregate annual ", N[2], "O emissions (g-N ",
                          ha^{-1}, ")", sep="")),
       fill="Site & Fertilization") +
  main.theme + 
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key=element_blank(),
        legend.text=element_text(size=14),
        strip.text=element_text(size=16))
ggsave("C://My Program Data/Dropbox/Articles/Dissertation/images/by year.png",
       height=7, width=9.5, units="in")

# Daily fluxes ----
f1.lab <- as_labeller(c(n2o.ha.day="N[2]*O~flux~'('*g~N~ha^{-1}~day^{-1}*')'",
                        no3="NO[3]^'-'~'('*mu*g~g^{-1}*')'",
                        nh4="NH[4]^'+'~'('*mu*g~g^{-1}*')'",
                        soilT="Soil~temp~(degree~C)",
                        wfps="WFPS~('%')"),
                      default=label_parsed)

n2o.mlt <- n2o.dat %>% 
  select(n2o.ha.day, no3, nh4, soilT, wfps, avp, fert, avp.fert, site) %>%
  melt %>% filter(!is.na(value))

n2o.mlt %>% 
  filter(variable %in% c("no3")) %>% 
  droplevels %>%
  ggplot(aes(y=value, x=avp, fill=avp, alpha=fert)) +
  geom_blank() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=which(levels(n2o.agg.dat$avp)=="Switchgrass")-0.5) +
  geom_vline(xintercept=which(levels(n2o.agg.dat$avp)=="Native grass")-0.5, 
             linetype=2) +
  geom_boxplot(outlier.size=1, width=1.1) +
  scale_alpha_discrete(range=c(1, 0.5), guide=FALSE) +
  scale_fill_brewer(palette="Dark2") + 
  scale_y_continuous(trans="asinh") +
  labs(y="") +
  facet_grid(variable~site, scale="free", switch="y", 
             labeller=labeller(.rows=f1.lab, .cols=label_value)) +
  guides(fill=guide_legend(nrow=1)) +
  main.theme +
  theme(strip.text=element_text(size=14, color="black"))

ggsave("G://Program Data/Dropbox/Articles/Dissertation/images/daily nit.png",
       height=5.6, width=9.5, units="in")

# Near-final figure ----
mod.pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d")
n2o.agg.dat %>%
  filter(!(system %in% c("G02", "G03", "G04", "G09")),
         !(avp.fert %in% c("Miscanthus unfrt", "Native grass unfrt")),
         (site=="ARL" & year %in% c("2010", "2011", "2012")) |
           (site=="KBS" & year == "2012"),
         !(site=="ARL" & system=="G07" & year=="2011")) %>%
  ggplot(aes(y=n2o.agg, x=avp, fill=avp, alpha=fert)) +
  geom_boxplot(color="black", width=1.1) + scale_y_log10() +
  geom_vline(xintercept=which(levels(n2o.agg.dat$avp)=="Switchgrass")-0.5) +
  geom_vline(xintercept=which(levels(n2o.agg.dat$avp)=="Native grass")-0.5, 
             linetype=2) +
  scale_fill_manual(values=mod.pal) +
  scale_alpha_discrete(range=c(1, 0.5), guide=FALSE) +
  labs(x="Cropping system", 
       y=expression(Aggregate~annual~N[2]*O~emissions~(g~N~ha^{-1})),
       fill="Site and fertilization")+
  facet_grid(.~site, scale="free_x", space="free_x") +
  guides(fill=guide_legend(nrow=1)) +
  main.theme

ggsave("C://My Program Data/Dropbox/Articles/Dissertation/images/mod fluxes.png",
       height=5.6, width=9.5, units="in")

