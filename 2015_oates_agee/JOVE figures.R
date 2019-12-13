# Libraries ----
library(plyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(lme4)
library(lsmeans)
library(grid)
library(reshape2)


cbbPalette <- c("#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#E69F00", "#999999")
#c("gray92", "gray55", "gray30")
shapes <- c(21,22,23,24,25)

main.theme <- theme(panel.grid=element_blank(),
                    panel.background=element_rect(color="black", fill="white"),
                    panel.margin.x=unit(0, "lines"),
                    axis.text.x=element_text(color="black", size=12),
                    axis.text.y=element_text(color="black", size=10),
                    axis.title=element_text(size=14),
                    legend.key=element_blank(),
                    legend.position="none",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.text=element_text(size=12),
                    strip.background=element_blank(),
                    strip.text.x=element_text(size=11),
                    plot.margin=unit(c(0.01, 0.02, 0.0, 0.01), "npc"))



# Data import ----
plfa.data <- read.delim("../Data/AARS split depth data for R.txt") %>%
  filter(Depth=="0 to 10",
         System %in% c("Corn", "Prairie")) %>%
  droplevels %>%
  mutate(trt=factor(interaction(System, Nfert),
                    levels=c("Corn.Yes", "Prairie.Yes", "Prairie.No"),
                    labels=c("Corn", "Fertilized prairie", "Unfertilized prairie")),
         fb.ratio=Fungi/Bacteria) 

# Isolate plfa concentrations, replace blanks with 0
plfa.markers <- plfa.data[,15:83]
plfa.markers[is.na(plfa.markers)] <- 0

# arcsine-square root transform PLFA markers to rein in outliers somewhat
plfa.markers <- asin(sqrt(plfa.markers))/(pi/2)
# remove lipids 20 carbons and longer (likely to be plant contaminants)
plfa.markers <- plfa.markers[,c(1:51,60:69)]

# Biomass ----
ggplot(plfa.data, aes(x=trt, y=Biomass, color=trt)) +
  geom_point(size=6, alpha=0.8) +
  labs(x="Cropping system",
       y=expression(FAME~biomass~'('*nmol~g^{-1}~soil*')'),
       color="Cropping system") +
  scale_color_brewer(palette="Dark2") +
  main.theme
ggsave("../JOVE Figs/JOVE F1.png")

plfa.data %>% select(trt, Bacteria, Fungi) %>% melt %>%
  group_by(trt, variable) %>%
  summarise(mval=mean(value)) %>%
  ggplot(aes(x=trt, y=mval, fill=trt, alpha=variable)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette="Dark2") +
  scale_alpha_discrete(range=c(1,0.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 0.3)) +
  labs(x="Cropping system",
       y=expression(FAME~biomass~'('*nmol~g^{-1}~soil*')')) +
  main.theme +
  theme(legend.position="bottom")
ggsave("../JOVE Figs/JOVE F2 alt.png")

plfa.mds <- metaMDS(plfa.markers, distance="euclidian", trymax=200) %>%
  scores %>% as.data.frame %>% cbind(trt=plfa.data$trt) 

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls <- ddply(plfa.mds, "trt", find_hull)

ggplot(plfa.mds, aes(x=NMDS1, y=NMDS2, color=trt)) +
  geom_polygon(data=hulls, aes(fill=trt), alpha=0.5) +
  geom_point(size=6) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  main.theme +
  theme(legend.position="bottom")
ggsave("../Jove Figs/Jove F3.png")
