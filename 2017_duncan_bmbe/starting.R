### Housekeeping----
# function for extracting seconds out of time in data
toSeconds <- function(x){
  if (!is.character(x)) stop("x must be a character string of the form H:M:S")
  if (length(x)<=0)return(x)
  
  unlist(
    lapply(x,
           function(i){
             i <- as.numeric(strsplit(i,':',fixed=TRUE)[[1]])
             if (length(i) == 3) 
               i[1]*3600 + i[2]*60 + i[3]
             else if (length(i) == 2) 
               i[1]*60 + i[2]
             else if (length(i) == 1) 
               i[1]
           }  
    )  
  )  
}

# libraries
library(ggplot2)
library(lme4)

# Constants
chamb_vol <- 0.389 # ml
soil_mass <- 127 # g
C_molec_mass <- 12.01 # ug per umol
m_per_day <- 1440 # i.e. final output is in day-1
gas_const <- 0.0820575
air_temp <- 294.26 #K
kPa_per_atm <- 101.325


# loop setup

t.folders <- paste0("../data/runs/",dir("../data/runs/"))
n.files <- length(list.files(t.folders, recursive=TRUE))

flux.df <- data.frame(chamber=rep("00", n.files), date=as.Date("2014-01-01"),
                      IRGA="neither", mode="none", l.r2=0, l.flux=0,
                      q.term=0, q.p = 1, q.r2=0, q.flux=0, 
                      co2.flux=0, co2.flux_day=0,
                      stringsAsFactors = FALSE)

t.pos <- 1 # Position in the output dataframe

for(t.i in 1:length(t.folders)){
  t.samples <- dir(t.folders[t.i])
  t.conv <- t.pos-1 # useful number, I'm hoping
  t.range <- t.pos:(t.pos+length(t.samples)-1) # use range instead of subset
  
  t.name <- t.folders[t.i]
  # Extract date from folder name
  flux.df[t.range,"date"] <- substr(t.name, (nchar(t.name)-9), nchar(t.name))

  for(t.j in t.range){
    t.path <- paste(t.folders[t.i], t.samples[t.j-t.conv], sep="/")
    t.data <- read.table(t.path, header=TRUE, skip=1, sep=";",
                         stringsAsFactors=FALSE)
    
    # Get time in units of minutes, then count minutes from incubation start
    t.data$clock <- toSeconds(t.data$Time.H.M.S)
    t.data$mins <- with(t.data, (1+clock-clock[1])/60)
    
    # Headspace CO2 in umol
    t.data$co2.umol <- with(t.data, CO2.ppm.*
                              ((CellPres.kPa./kPa_per_atm)*chamb_vol)
                               / (gas_const*air_temp))
    # Headspace CO2 in ug per g soil (penultimate units)
    t.data$co2.ug_g <- t.data$co2.umol*C_molec_mass/soil_mass
    
    # subset out data to remove early noise
    t.data_sub <- subset(t.data, mins>=0.25)
    
    t.l <- lm(co2.ug_g~mins, data=t.data_sub)
    flux.df[t.j,"l.r2"] <- summary(t.l)$adj.r.squared
    flux.df[t.j,"l.flux"] <- t.l$coefficients[2]
          
    t.q <- lm(co2.ug_g~mins + I(mins^2), data=t.data_sub)
    flux.df[t.j,"q.term"] <- t.q$coefficients[3]
    flux.df[t.j,"q.p"] <- summary(t.q)$coefficients[3,4]
    flux.df[t.j,"q.r2"] <- summary(t.q)$adj.r.squared
    flux.df[t.j,"q.flux"] <- t.q$coefficients[2]
    
    if(t.q$coefficients[3] < 0){
      flux.df[t.j,"mode"] <- "quad"
      flux.df[t.j,"co2.flux"] <- flux.df[t.j,"q.flux"]
    } else{
      flux.df[t.j,"mode"] <- "lin"
      flux.df[t.j,"co2.flux"] <- flux.df[t.j,"l.flux"]
    }
    
    flux.df[t.j,"co2.flux_day"] <- flux.df[t.j,"co2.flux"]*m_per_day
    
    flux.df[t.j, "chamber"] <- paste0("C",substr(t.samples[t.j-t.conv],8,9))
    flux.df[t.j,"IRGA"] <- if(mean(t.data_sub$CellTemp.c.) < 51){
                              "Left"}else{"Right"}
  }
  
  t.pos <- max(t.range)+1 # advance range for the next date
}

write.table(flux.df, "temp.txt", row.names=FALSE)

t.key <- read.table("../data/zip_chamber_data.txt", header=TRUE, sep="\t")
names(t.key)[5] <- "chamber"
t.key$plant <- factor(t.key$plant)
int.df <- merge(flux.df, t.key)

t.sub <- int.df

t.sub$group <- factor(with(t.sub,paste0(genotype,tissue)))

#lmer1 <- lmer(co2.flux_day~group + (1|rep/chamber) + (1|date), data=t.sub)
#summary(lmer1)
#t.sub$resid2 <- residuals(lmer1)^2

t.avg_resid <- by(t.sub$resid2, t.sub$chamber, mean)
exclude.df <- data.frame(chamber=names(t.avg_resid), resid=as.vector(t.avg_resid))
exclude.df <- merge(exclude.df, t.key)

write.table(exclude.df,"temp.txt", row.names=FALSE)


t.plant <- c(1:18)
t.cont <- subset(t.sub, tissue=="")
t.root <- subset(t.sub, tissue=="ROOT")
t.leaf <- subset(t.sub, tissue=="LEAF")

ggplot(t.cont, aes(x=date, y=co2.flux_day,color=genotype, linetype=tissue))+
  geom_point()+stat_smooth(se=FALSE)+
  geom_point(data=t.root)+stat_smooth(data=t.root,se=FALSE)+
  geom_point(data=t.leaf)+stat_smooth(data=t.leaf,se=FALSE)+
  ylim(0,200)

t.sub2 <- subset(t.sub, plant==t.plant)
ggplot(t.sub2, aes(x=date, y=resid2, color=chamber))+
  geom_point()+stat_smooth(aes(color=plant))
  geom_line()#+ylim(0,500)


  
t.sub4 <- subset(t.sub, chamber%in%c("C05","C06"))
ggplot(t.sub4, aes(x=date, y=co2.flux_day, color=chamber))+
  geom_point()+geom_line()#+ylim(0,500)

t.sub5 <- subset(int.df, genotype=="SOIL")
ggplot(t.sub5, aes(x=date, y=co2.flux_day, color=genotype))+
  geom_point()+stat_smooth()+ylim(0,55)
