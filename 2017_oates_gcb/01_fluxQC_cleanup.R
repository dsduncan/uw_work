# Libraries ----

source("00_useful_functions.R")
require(data.table)
require(lubridate)


# Read in AARS data ----
# Data are in a single file, also including second round of visual inspections.
# For raw data, `study` == 'aars bcse'. Data from 2009-12 were already 
# evaluated and don't need to be run again (except for microplots in 2012).

# Read data
t_aars <- fread("data/fluxQC_downloads/old_fluxQC_dump.csv")

# Standardize
t_aars[, ':=' (
  sampled_on = as.Date(sampled_on, "%Y-%m-%d"),
  site = "AARS",
  replicate = sub("R", "A", replicate)
  )]

# Filter to exclude earlier data
t_aars <- t_aars[
  study == "aars bcse" & 
    (sampled_on > '2012-12-31' | grepl("M", treatment)
     )]



# Read in KBS data ----
# KBS data had to be read in as individual files. These files have to be
# joined into one data table.
t_kbs <- batchRead("data/staging_KBS_data/")
# Modify data
t_kbs[, ':=' (
  sampled_on = as.Date(sampled_on, format="%Y-%m-%d"),
  site = "KBS",
  replicate = sub("R", "K", replicate)
)]
# Standardize column names
setnames(t_kbs, "compound", "name")




# Standardize KBS microbplot treatment names
t_kbs <- subset(t_kbs, !(treatment %in% c("L1","L2","L3","M1","M2","M3","M4")))  
t_kbs$treatment <- gsub("M1","M", t_kbs$treatment)
t_kbs$treatment <- gsub("no_cover","M", t_kbs$treatment)

# No need to re-analyze data from Oates et al. 2015, so only keep data from
# 2012+, or from microplots before then.
t_kbs_mirco <- grepl("M", t_kbs$treatment)
t_kbs <- subset(t_kbs, (sampled_on > "2011-12-31" | t_kbs_mirco))

# 3.3 Combine AARS and KBS datasets ----
# Align columns and merge data frames
t_cols <- c("site", "sampled_on", "treatment", "replicate",  
                    "lid", "avg_height_cm", "minutes", "ppm", "name")
t_aars <- t_aars[,t_cols]
t_kbs <- t_kbs[,t_cols]
conc.temp <- rbind(t_aars, t_kbs)

# Rename columns for simplicity and consistency, then sort
names(conc.temp)[c(2:7)] <- c("date", "trt", "block", "bucket", "height.cm", "d.min")
conc.temp <- conc.temp[with(conc.temp, order(site, date, trt, block, d.min, name)),]

## Unstack GHGs
# Note: There will be errors due to duplicated measurements. These have been
# reviewed and the duplicated measurements come from samples lacking data in
# the d.min field. 
conc <- reshape(conc.temp, v.names = "ppm", timevar = "name",
                idvar = c("date", "trt", "block", "d.min", "height.cm"), direction = "wide")

## Remove unusable data
# No or very low N2O concentration values
conc <- conc[(!is.na(conc$ppm.n2o)) & conc$ppm.n2o > 0.2,]
# Remove data without bucket height
conc <- conc[conc$height.cm > 0,]

# Create a unique series name to store metadata during the HMR process
# Also gives a unique identifier for counting observations per sample
conc$series <- factor(with(conc, paste(date, site, block, trt ,sep = "_")))
# Use only series with 3+ data points. It makes sense to attach this
# information to conc so I can use it later on for identifying samples that
# cannot be nonlinear.
counts.temp <- table(conc$series)
conc$count <- counts.temp[match(conc$series, names(counts.temp))]
conc <- subset(conc, count >= 3)

# It's always good to sort
conc <- conc[with(conc, order(site, date, trt, block, d.min)),]

## Cleanup
rm(list=ls(pattern="temp"))

# 3.4 Set bucket dimensions----
# Different buckets were used at different times in each site, so we need to
# calculate a separate volume for each observation.

conc$height.m <- 0
conc$surfRad.m <- 0
conc$vol.m3 <- 0

### Plastic (type Z) buckets
# Slightly but noticeably conical, growing larger near the top of the bucket.
# Cross sectional area at surface will thus be a function of insertion height.
# All constants were obtained by measuring a set of 3 "representative" buckets
# we had lying around.
kZLid <- 1.25 # Distance from lid bottom to bucekt lip, cm
kZBot <- 12.765 # Bottom radius, cm
kZHeight <- 25.45 # Total height of plastic bucket, cm
kZDelta <- 1.385 # Increase in radius at base of lid, cm
kZTop <- (kZBot+kZDelta)/100 # Radius at base of lid in m

temp.Z.height <- conc[conc$bucket == "Z", "height.cm"] - kZLid
conc[conc$bucket == "Z", "height.m"] <- temp.Z.height/100

conc[conc$bucket == "Z", "surfRad.m"] <- 
  (kZBot + (kZHeight - temp.Z.height)/kZHeight*kZDelta)/100

# Because plastic buckets are truncated cones, we calculate their volume as:
# V = 1/3 pi (r^2+R^2+r*R) h
# where r and R are radii at the surface and lid bottom
temp.Z.rad <- conc[conc$bucket == "Z", "surfRad.m"]
{conc[conc$bucket == "Z", "vol.m3"] <-
    pi/3 * (temp.Z.rad^2 + kZTop^2 + kZTop*temp.Z.rad)* #R^2+r^2+r*R
    conc[conc$bucket == "Z", "height.m"]} # h

### Metal (Type Y) buckets
# Assumed to be perfect cylinders.
# All constants derived from the manufacturing specifications
kYLid <- 0.3 # Distance from lid bottom to bucket lip, cm
kYRad <- 14.225 # Radius in cm

temp.Y.height <- conc[conc$bucket == "Y", "height.cm"]-kYLid
conc[conc$bucket == "Y", "height.m"] <- temp.Y.height/100

conc[conc$bucket == "Y", "surfRad.m"] <- kYRad/100

conc[conc$bucket == "Y", "vol.m3"] <- 
  pi*conc[conc$bucket == "Y", "surfRad.m"]^2* # radius
  conc[conc$bucket == "Y", "height.m"]        # height

### Area at soil surface (both bucket types)
conc$surfA.m2 <- pi*conc$surfRad.m^2

### Remove temporary objects
rm(list=ls(pattern="k"))
rm(list=ls(pattern = "temp"))

#### Backup conc data ###
# write.table(conc, "../data/R intermediate/conc_backup.txt", sep = "\t", row.names=FALSE, quote=FALSE)
####---###---###---###

# 3.5 Calculate fluxes----
# We calculate fluxes at this stage primarily to determine the "shape of the 
# time series data prior to the second round of visual inspections. 
# The HMR package requires a discreet input *file* with a specific format. It
# will thus be necessary to organize the data into this format and write
# the file prior to analysis.
# Fluxes with >= 4 observations wil be written as a file and sent to HMR for
# analysis. Fluxes with 3 observations will have separate linear estimates.
# The datasets will eventually be combined.

#### Restore from backups if necessary
# conc <- read.table("../data/R intermediate/conc_backup.txt", sep = "\t", header = TRUE, as.is=TRUE)
# conc$date <- as.Date(conc$date, "%Y-%m-%d")

# Organize data for export to HMR
# Note that to get fluxes for a different gas, you'd change the Concentration
# term
conc.out <- with(conc, data.frame(Series = series, V = vol.m3, A = surfA.m2,
                                  Time = d.min, Concentration = ppm.n2o))

# Break data into 2 groups: those with 4+ timepoints, those with 3 or fewer.
# First class will have HMR run on them, rest will have linear fits.
conc.out.nlin <- conc.out[conc$count >= 4,]
conc.out.lin <- conc.out[conc$count < 4,]

#### Write out table for HMR to read
write.table(conc.out.nlin, "../data/R intermediate/HMR_input.txt", sep=";", 
            quote=FALSE, row.names=FALSE)

# Note that the directory needs to be changed for HMR to work.
# ProgHMR is a wrapper for HMR that is in useful funcitons. It doesn't change
# anything substantive, but makes the process slightly more fun.
setwd("../data/R intermediate")
temp.nlr <- ProgHMR("HMR_input.txt", step = 500)
setwd("../../N2O for GCB")

#### Restore from backups if necessary
# temp.nlr <- read.table("../data/R intermediate/HMR_backup-HMR_input.txt", sep = ";", header = TRUE, as.is=TRUE)

# Data with 3 points are given linear fits and given a special warning
temp.lr <- LRintoHMR(conc.out.lin, warn = "Fit with <4 points")

# Both data types are combined and sorted
flux <- rbind(temp.nlr, temp.lr)
flux <- flux[order(flux$Series),]

rm(list=ls(pattern="temp"))
rm(list=ls(pattern="conc."))

#### BACKUP ###
# write.table(flux, "../data/R intermediate/flux_backup.txt", sep = "\t", row.names=FALSE, quote=FALSE)

# 4.0 Secondary visual inspection ----
# As per our protocol, initial flux estimation with HMR is followed by secondary
# visual inspection, again using the FluxQC pipeline. Samples are analyzed in
# groups based on the estimation method (linear, nonlinear, no flux), so we can
# qualitatively evaluate the evidence for that flux type.

#### Restore backups if necessary
# conc <- read.table("../data/R intermediate/conc_backup.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
# conc$date <- as.Date(conc$date, "%Y-%m-%d")
# flux <- read.table("../data/R intermediate/flux_backup.txt", sep = "\t", header = TRUE, as.is = TRUE)

## Prep uploads ##
# First, if a nonlinear flux's 95% confidence interval (assuming fixed K value
# for HMR v0.4.1+) includes the linear estimate, the linear estimate is used
# instead. This imposes a higher standard of evidence on the use of the 
# nonlinear estimate, as it is more error prone.
flux$HvL <- with(flux, (f0-LR.f0)/f0.se)
kCrit <- qt(p = 0.975, df = 4-2) # SE radius of 95% confidence interval with 4 data points

temp.switch <- flux$Method == "HMR" & abs(flux$HvL) < kCrit
flux$Method[temp.switch] <- "LR"

# 5.0 Write out concentrations for second round of visual inspections in FluxQC ----
# Write out nonlinear fluxes for upload to FluxQC
PrepFQC(flux, conc, method="HMR", datadir="../data/FluxQC uploads/Nonlinear visual")

# Write out linear fluxes for upload to FluxQC
PrepFQC(flux, conc, method="LR", datadir="../data/FluxQC uploads/Linear visual")

# Write out No flux data for upload to FluxQC
PrepFQC(flux, conc, method="No flux", datadir="../data/FluxQC uploads/No flux visual")

# 6.0 Second round of visual inspections ----
# At this point, you need to manually upload data into FluxQC, run the visual
# inspection to see if any points are better off being removed. Then, download
# the re-inspected concentration data.
