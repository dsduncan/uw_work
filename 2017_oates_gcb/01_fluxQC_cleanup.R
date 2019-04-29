#### Libraries ####

source("00_useful_functions.R")
require(data.table)
require(lubridate)

#### Read in AARS data ####
# Data are in a single file, also including second round of visual inspections.
# For raw data, `study` == 'aars bcse'. Data from 2009-12 were already 
# evaluated and don't need to be run again (except for microplots in 2012).

# Read data
t_aars <- fread("data/fluxQC_downloads/old_fluxQC_dump.csv")

# Filter to exclude earlier data
t_aars <- t_aars[
  study == "aars bcse" & 
    (sampled_on > '2012-12-31' | grepl("M", treatment)
    )]

# Standardize data
t_aars[, ':=' (
  sampled_on = as.Date(sampled_on, "%Y-%m-%d"),
  site = "AARS",
  replicate = sub("R", "A", replicate),
  id = NULL,
  study = NULL
  )]


#### Read in KBS data ####
# KBS data had to be read in as individual files. These files have to be
# joined into one data table.
t_kbs <- batchRead("data/staging_KBS_data/")

# Standardize data
t_kbs[, ':=' (
  sampled_on = as.Date(sampled_on, format="%Y-%m-%d"),
  site = "KBS",
  replicate = sub("R", "K", replicate),
  treatment = sub("[LM][1-4]|no_cover", "M", treatment),
  study = NULL
)]

# Rename column that identifes the trace gas
setnames(t_kbs, "compound", "name")

# Microplots not associated with a treatment
t_kbs <- t_kbs[(treatment != "M") & (!is.na(ppm))]

# Only keep main plots after 2011 or any mircoplots
t_kbs <- t_kbs[sampled_on > '2011-12-31' | grepl("M", treatment)] 

#### Combine datasets ####

# Align columns and merge data frames
t_cols <- c(
  "site", "sampled_on", "treatment", "replicate", "lid", "avg_height_cm", 
  "minutes", "ppm", "name")
setcolorder(t_aars, t_cols)
setcolorder(t_kbs, t_cols)
t_dat <- rbind(t_aars, t_kbs)

# Rename columns (names dictated by convention)
t_oldnames <- c(
  "sampled_on", "treatment", "replicate", "lid", "avg_height_cm", "minutes")
t_newnames <- c("date", "trt", "block", "bucket", "height_cm", "d_min")
setnames(t_dat, t_oldnames, t_newnames)
setorder(t_dat, site, date, trt, block, d_min, name, -ppm)

# Drop duplicated columns (both of them!)
t_dup <- duplicated(t_dat, by=c('site','date','trt','block','d_min','name'))
t_dat <- t_dat[!t_dup]

# Reshape data to put GHGs in same row
conc <- reshape(
  t_dat, v.names="ppm", timevar="name",
  idvar = c("date", "trt", "block", "d_min", "height_cm"), direction = "wide")
# Convert conc from dataframe to data.table
setDT(conc)

# Remove unusable data
conc <- conc[(!is.na(ppm.n2o)) & ppm.n2o > 0.2 & height_cm > 0]

# Create a unique series name to store metadata during the HMR process
# (Also weed out measurements with < 3 data points)
conc[, series := factor(paste(date, site, block, trt, sep="_"))]
t_cnts = conc[, list(nobs = .N), by = series]
conc <- conc[series %in% t_cnts[nobs >= 3, series]]
# Sort
setorder(conc, site, date, trt, block, d_min)

# Cleanup
rm(list=ls(pattern="t_"))

#### Calculate bucket dimensions ####
# Different buckets were used at different times in each site and inserted to
# different depths, so we need to calculate a separate volume for each 
# observation.

## Plastic (type Z) buckets ##
# Conical, larger near the top, so area at surface depends on insertion height.
# Constants based on measurement of 3 representative buckets.
t_lid_z    <- 0.0125           # Distance from lid bottom to bucket lip, m
t_bot_z    <- 0.12765          # Bottom radius, m
t_height_z <- 0.2545           # Total height of plastic bucket, m
t_delta_z  <- 0.01385          # Increase in radius at base of lid, m
t_top_z <- t_bot_z + t_delta_z # Radius at base of lid, m

conc[bucket == "Z", height_m := (height_cm/100) - t_lid_z]
conc[bucket == "Z", surfRad_m := 
       (t_bot_z + t_delta_z * (t_height_z - height_m) / t_height_z)]
# V = 1/3 * pi * (r^2 + R^2 + r*R) * h
conc[bucket == "Z", vol_m3 := 
       pi/3 * (surfRad_m^2 + t_top_z^2 + surfRad_m*t_top_z) * height_m]

## Metal (Type Y) buckets ##
# Assumed to be perfect cylinders. Constants derived from manufacturing
# specifications.
t_lid_y <- 0.003   # Distance from lid bottom to bucket lip, m
t_rad_y <- 0.14225 # Bucket radius, m

conc[bucket == "Y", height_m := (height_cm/100) - t_lid_y]
conc[bucket == "Y", surfRad_m := t_rad_y / 100]
conc[bucket == "Y", vol_m3 := pi * surfRad_m^2 * height_m]

## Area at soil surface ##
conc[, surfA_m2 := pi * surfRad_m ^2]

## Remove temporary objects ##
rm(list=ls(pattern = "t_"))

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
