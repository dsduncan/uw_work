###---###---###---###
# These are some useful functions that automate a lot of the work necessary to
# juggle data between the FluxQC visual inspection platform and the HMR package
# for estimating nonlinear fluxes.
# 
# David Duncan
# 2017-09-17
###---###---###---###


FluxAgg <- function(dataset, date.table, site = "site", plot = "plot", year = "year", trt = "trt", 
                    block = "block", flux = "n2o.ha.day", date = "date", date.format = "%Y-%m-%d"){
  # FluxAgg 
  # Takes daily trace gas fluxes and linearly interpolates them across a year. Requires
  # a specifically formatted dataset as well as a table of start and end dates for emissions
  
  # Basic formatting checks
  dataset[,date] <- as.Date(dataset[,date], format = date.format)
  dataset[,year] <- as.character(dataset[,year])
  
  date.table[,"year"] <- as.character(date.table[,"year"])
  date.table[,c("start", "end")] <- sapply(date.table[,c("start", "end")], as.character)
  date.table[,"start"] <- as.Date(date.table[,"start"])
  date.table[,"end"] <- as.Date(date.table[,"end"])
  
  # Sort data by site, year, treatment, block, date
  dataset <- dataset[order(dataset[,site], dataset[,year], dataset[,trt], 
                           dataset[,block],dataset[,date]),]
  
  # Create a dataset to hold output
  # Will consist of all of the unique plot-years
  out.df <- unique(dataset[,c(site, year, trt, block)]) 
  
  # Empty columns for holding output
  out.df <- cbind(out.df, flux.Agg = NA)
  
  # Dates at which soil first thawed and froze permanently. In the absence of
  # contradictory observations, we assume these dates demarcate the period over
  # which trace gasses are emitted from the soil.
  
  out.site.year <- paste(out.df[,site],out.df[,year],sep = "")
  date.site.year <- paste(date.table[,"site"],date.table[,"year"], sep = "")
  
  out.df$start <- date.table[match(out.site.year, date.site.year), "start"]
  out.df$end <- date.table[match(out.site.year, date.site.year), "end"]
  
  
  # Lists for filtering dataset for unique plot-years
  ds.hash <- paste(dataset[,site], dataset[,year], dataset[,trt], dataset[,block],sep = "")
  out.hash <- paste(out.df[,site], out.df[,year], out.df[,trt], out.df[,block],sep = "")
  
  # Iterate over each plot-year in out.df to generate yearly aggregates
  for (i in 1:nrow(out.df)){
    temp.data <- dataset[ds.hash == out.hash[i],]
    # Number of measurements in dataset
    temp.meas <- nrow(temp.data)
    
    # Important to ensure there's at least a couple of days of data
    if(temp.meas >= 2){
      # Determine whether "official" start and end dates bracket our observations
      out.df[i,"start"] <- min(out.df[i, "start"], temp.data[1,date])
      out.df[i,"end"] <- max(out.df[i, "end"], temp.data[temp.meas,date])
    
    
      # Days between measurements
      # When measurements were taken outside of the start and end date, 
      # the distance will be 0
      temp.days <- vector("numeric", temp.meas+1)
      temp.days[1] <- temp.data[1,"date"] - out.df[i,"start"]
      temp.days[2:temp.meas] <- temp.data[2:temp.meas, date] - 
       temp.data[1:temp.meas-1, date]
      temp.days[temp.meas+1] <- out.df[i,"end"] - temp.data[temp.meas, date]
    
      # Official fluxes, averaged over consecutive dates
      temp.flux <- vector("numeric", temp.meas+1)
      temp.flux[c(1,temp.meas+1)] <- temp.data[c(1,temp.meas),flux]/2
      temp.flux[2:temp.meas] <- (temp.data[1:temp.meas-1,flux] + 
                                 temp.data[2:temp.meas,flux])/2
    
      # Multiply days between sample times by average flux between sample times
      # to get cumulative flux between times
      out.df$flux.Agg[i] <- temp.days %*% temp.flux
    } else{
      out.df$flux.Agg[i] <- 0
    }
  }
  
  return(out.df)
}

LRintoHMR <- function(dataset, warn = 'Fit outside of HMR'){
  # LRintoHMR takes timeseries concentration data, generates a linear estimate of
  # concentration change over time (flux), and formats the data to resemble the
  # output of the HMR function, with LR.alwys = TRUE. All data, including 
  # concentrations and dimensions, need to be set as if you were about to load
  # the data into HMR. I am doing this because the process is so slow and it makes
  # no sense to run it for data we were going to treat as linear anyway.
  
  out <- data.table(Series = unique(dataset$Series), 
                   f0 = 0, f0.se = 0, f0.p = 0, 
                   Method = 'LR', Warning = warn,
                   LR.f0 = 0, LR.f0.se = 0, LR.f0.p = 0,
                   LR.Warning = '',
                   IG.f0 = 0, IG.f0.se = 0, IG.f0.p = 0)
  # HMR adjusts for effective chamber height, so we must do it here
  dataset[, Concentration := Concentration*V/A]
  # Progress bar for fun
  pb <- txtProgressBar(min = 0, max = nrow(out), style = 3)
  for (i in seq_along(out$Series)){
    # Linear flux model
    lm_flux <- dataset[Series == out[i, Series]] %>%
      lm(formula = Concentration ~ Time, data = .) %>%
      summary
    # Extract values
    out[i, c('f0', 'f0.se', 'f0.p', 'LR.Warning') := list(
      lm_flux$coefficients["Time", "Estimate"],
      lm_flux$coefficients["Time", "Std. Error"],
      lm_flux$coefficients["Time", "Pr(>|t|)"],
      case_when(
        lm_flux$coefficients["(Intercept)", "Estimate"] < 0 ~ 'Negative intercept',
        TRUE ~ 'None')
    )]

    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # For a linear regression, the linear values will, of course, be the same
  out[, c('LR.f0', 'LR.f0.se', 'LR.f0.p', 'IG.f0', 'IG.f0.se', 'IG.f0.p') := 
        list(f0, f0.se, f0.p, f0, f0.se, f0.p)]
  
  close(pb)
  return(out)  
}

PrepFQC <- function(hmr.dat, conc.dat, method, datadir=""){
  # PrepFQC takes a dataset formatted like an HMR output with a matching 
  # dataset that could be used to build an HMR input and takes a subset from 
  # one of the three estiamtion methods (specified), and formats it for upload
  # into FluxQC for subsequent visual inspection.
  
  sub.series <- as.character(hmr.dat[hmr.dat$Method == method, "Series"])
  conc.dat <- conc.dat[conc.dat$series %in% sub.series,]
  len <- nrow(conc.dat)
  kStep <- 150 # number of samples to run at once, based on FluxQC capacities
  
  ## Setup file
  # Follows the style of GLBRC 2011 inputs
  # Be *very* careful with data formatting, as this can be very finicky
  
  # Need to include date and id information to keep things unique
  temp.treat <- paste(conc.dat$trt, conc.dat$date, sep=";")
  temp.block <- paste(conc.dat$block, conc.dat$id, sep=";")
  
  setup <- data.frame(Treatment = temp.treat, Rep = temp.block)
  
  # Bucket sample order (T0-T3), assumes series are in order, which they should be 

  t1 <- c(FALSE, 
          conc.dat$series[1:(len-1)] == conc.dat$series[2:len])
  t2 <- c(FALSE, FALSE,
          conc.dat$series[1:(len-2)] == conc.dat$series[3:len])
  t3 <- c(FALSE, FALSE, FALSE,
          conc.dat$series[1:(len-3)] == conc.dat$series[4:len])
  setup$T0 <- t1+t2+t3
  
  # Microplot data should be embedded in the block or treatment name
  setup$micro <- ""
  setup$chamber <- ""

  # Vial # matches up the setup and data sheets
  setup$vial <- 1:len

  # I think FluxQC expects a simple range of bucket types
  setup$lid <- conc.dat$bucket

  # FluxQC expects 4 heights
  # Because we are clever, we shall give them effective heights, which is all
  # HMR really needs. This will prevent us from needing to calcluate them from
  # measured heights again. It is very important to remember this though, or else
  # we'll keep making the buckets shorter and shorter each time.
  setup$h1 <- conc.dat$vol.m3/conc.dat$surfA.m2
  setup$h2 <- setup$h1
  setup$h3 <- setup$h1
  setup$h4 <- setup$h1

  # Time and temperature, pretty easy
  setup$temp <- 0
  setup$hours <- 0
  setup$min <- conc.dat$d.min
  setup$seconds <- 0
  setup$time <- setup$min
  setup$atime <- setup$min

  # Clever, hide the series name in the notes for use further down
  setup$notes <- conc.dat$series
  
  ## Data file
  #The data file requires only a vial number and GHG concentrations
  # It is important that these be in PPM, based on the way FluxQC graphs data
  data <- with(conc.dat, data.frame(vial=1:len, n2o=ppm.n2o, co2=ppm.co2, ch4=ppm.ch4))
  
  ## Write files
  # We may need to swap directories back and forth
  old.dir <- getwd()
  if(nchar(datadir) > 0){
    setwd(datadir)
  }
  
  # Make file names
  # NOTE: Dates are fictitious and used only for sorting 
  # These hard-coded dates build off of the last set I uploaded.
  # This would need to be redone for every data set. Clearly not a good way
  # to do this.
  if(method == "HMR"){
    prefix <- "NLR"
    s.date <- "2000-09-10"
  } else if(method == "LR"){
    prefix <- "LR"
    s.date <- "2000-04-04"
  } else {
    prefix <- "MV"
    s.date <- "2000-01-06"
  }
  
  #Helper function to simplify some redundant code
  helpPrep <- function(setup, data, i, prefix, s.date){
    ## Helper function to create a separate file for a subset of data
    ## Necessary because FluxQC chokes on big data files
    
    # Use the date to order the files
    s.date2 <- strftime(as.character(as.Date(s.date,"%Y-%m-%d")+i-1), "%Y-%m-%d")
    setup.name <- paste(prefix, i, "setup.csv", sep = "_")
    data.name <- paste(prefix, i, "data.csv", sep = "_")
    
    # Writeout Setup file
    fileCon <- file(setup.name)
    writeLines(c("GLBRC 2011 Series 1,,,,,,,,,,,,,,,,,",
                 "data entered by: AUTO,,,,,,,,,,,,,,,,,",
                 "data checked by: hope,,,,,,,,,,,,,,,,,",
                 paste("sample date: ",s.date2,",,,,,,,,,,,,,,,,,", sep = ""),
                 ",,,,,,,Box,Box,Box,Box,,,,,,,",
                 "Treatment,Rep,T0-T3,Microplots,Chamber #,Vial #,Lid Type,Ht (cm), Ht (cm),Ht (cm),Ht (cm),Soil temp,Hours,Minutes,Seconds,Time ,Actual time,Notes"), 
               fileCon)
    close(fileCon)
    write.table(setup, setup.name, append=TRUE, quote=FALSE,sep = ",", 
                col.names=FALSE, row.names=FALSE)
    
    # Write out Data file
    fileCon <- file(data.name)
    writeLines(c("Z:\\GLBRC SERIES1 040711\\GLBRC SERIES 1 040711 2011-04-08 10-19-17\\374Fhm01.D,,,",
                 "No flux visual inspection,,,",
                 "All Values in PPM,,NA=invalid data,",
                 "Sample,,,",
                 ",N20,CO2,CH4"), fileCon)
    close(fileCon)
    write.table(data, data.name, append=TRUE, quote=FALSE, sep = ",", 
                col.names=FALSE, row.names=FALSE)
  }

  # Write out files with "kStep" series in them, hopefully of manageable size
  for(i in 1:(length(sub.series)%/%kStep)){
    ser.sub <- sub.series[(((i-1)*kStep)+1):(i*kStep)]
    set.sub <- setup[setup$notes %in% ser.sub,]
    dat.sub <- data[setup$notes %in% ser.sub,]
    helpPrep(set.sub, dat.sub, i, prefix, s.date)
  }
  
  i <- i + 1
  ser.sub <- sub.series[(((i-1)*kStep)+1):length(sub.series)]
  set.sub <- setup[setup$notes %in% ser.sub,]
  dat.sub <- data[setup$notes %in% ser.sub,]
  helpPrep(set.sub, dat.sub, i, prefix, s.date)

  setwd(old.dir)
}

ProgHMR <- function(filename, step=100, LR.always=TRUE, FollowHMR=TRUE){
  # ProgHMR
  # Function for running HMR piecewise and providing some indication of how it's progressing
  # filename: name of properly formatted text file to generate HMR data
  # step: number of samples to run at a time (for reporting progress)
  # LR.always: HMR argument, always include linear estimate
  # FollowHMR: HMR argument, use HMR recommendations, necessary to automate procedure
  
  require(data.table)
  require(HMR)
    
  # Input file needs same formatting as normal HMR
  infile <- fread(filename, sep=";", header=TRUE)
  # Sort the input file, just in case
  setorder(infile, Series, Time)
  
  # Total number of series to be run
  ser <- infile$Series %>% unique
  lser <- length(ser)
  
  # Create empty dataframe for output. This should be faster than appending
  # an existing one via rbind
  hmr.out <- data.table(Series = ser, 
                        f0 = 0, f0.se = 0, f0.p = 0, 
                        Method = "", Warning = "", 
                        LR.f0 = 0, LR.f0.se = 0, LR.f0.p = 0, 
                        LR.Warning="",
                        IG.f0 = 0, IG.f0.se = 0, IG.f0.p = 0)
  # Want to write out the table, then append the data to the file
  out.filename <- paste0("HMR_backup-", filename)
  fwrite(hmr.out[0], out.filename, sep=";")
  
  # Main loop
  # Have HMR run step series at a time, add the output to the output table, also
  # append it to a backup file being written out, just in case.
  # Structuring it as a while loop (rather than for) to handle the case where
  # the number of series is smaller than step.
  min.ser <- 1
  while (min.ser <= lser){
    # Take time to calculate how long this iteration took
    t0 <- proc.time()[3]
    
    # Range of series to calculate at once, last set may be a bit short
    max.ser <- min.ser + step - 1
    step.range <- min.ser:min(max.ser, lser)
    ser.list <- ser[step.range]
    
    # Copy HMR's output so we can make a nice table out of it.
    # Need to do this to keep appending the backup data file.
    # Note that HMR output is apparently all text (not sure why), so it is
    # necessary to turn numbers into numbers.
    # As of v0.4.1, HMR has started incorporating uncertainty about the K 
    # parameter in its estimate of error for HMR fits. In my hands, this causes
    # problems with estimating that standard error. I'm disabling that feature
    # for now (kappa.fixed).
    HMR(filename, series=ser.list, FollowHMR=FollowHMR, LR.always=LR.always, 
        kappa.fixed=TRUE)
    
    # Reading in file fixes unit conversion issues
    step_out <- fread(paste0("HMR - ", filename))
    step_out[, c('f0.lo95', 'f0.up95', 'LR.f0.lo95', 'LR.f0.up95') := NULL]

    step_out[, c('IG.f0', 'IG.f0.se', 'IG.f0.p') := list(f0, f0.se, f0.p)]
    
    
    for(i in seq_along(ser.list)) {
      if(step_out[i, Method] == 'HMR') {
        lm_flux <- infile[Series == ser.list[i]][1:3] %>%
          lm(formula = (Concentration * V / A)~ Time, data = .) %>%
          summary
        step_out[i, c('IG.f0', 'IG.f0.se', 'IG.f0.p') := list(
          lm_flux$coefficients["Time", "Estimate"],
          lm_flux$coefficients["Time", "Std. Error"],
          lm_flux$coefficients["Time", "Pr(>|t|)"]
        )]
      }
    }
    
    # Temporary file gets added into exisiting data table
    hmr.out[step.range] <- step_out

      
    # Append results to file, just in case there's a crash or something
    fwrite(step_out, out.filename, append=TRUE, sep=";")
    
    # Index advances forward
    min.ser <- max.ser + 1

    # Finally, calculate and report progress and time elapsed
    pct <- min(round((max.ser/lser)*100, digits=1), 100)
    dt <- round(proc.time()[3]-t0, digits=1)
    feedback <- paste0(pct,"% finished, ", dt, " seconds.")
    print(feedback)    
  }
  # Returning full table, but output written out differs slightly from this, so
  # it's best to use the "officially written" file instead.
  return(hmr.out)
}

batchRead <- function(path){

  # Read all of the files in the folder in the supplied path and append them
  # into a single table. Called sufficiently often to warrant a standalone
  # function.
  require(data.table)
  
  # List of all files in directory
  files <- dir(path)
  # Path to each file
  file_paths <- paste(path, files, sep = "")
  # Use fread to read in files, then group together into a single object
  out <- rbindlist(lapply(file_paths, fread))

  return(out)  
}

asinhAvg <- function(x){
  return(sinh(mean(asinh(x))))
}

plotPrep <- function(trt, block, ec_days, frFlux, df){
  m1.temp <- which(df$trt == trt & df$block == block)
  
  # Subset of data matching (unique combination of) trt and block
  sub.temp <- subset(df, trt==trt & block==block)
  
  # Slightly sloppy. Site is redundant with block but used frequently, so I 
  # still need to have a site variable.
  site.temp <- sub.temp[1,"site"]
  
  ###---###---###---###
  # No gapfill approach
  
  # Make a dataframe with all days for which we have data, observed or assumed,
  # and fill all with winter flux estimates (or no flux)
  days.temp <- sort(unique(c(sub.temp$date, 
                             subset(ec_days, site==site.temp)$date)))
  flux.temp <- data.frame(date=days.temp, n2o.0=0,
                          n2o=subset(frFlux, site==site.temp)$n2o)
  
  # Add in observed data where we have it 
  # Index of dates with field data
  ind1.temp <- match(sub.temp$date, days.temp)
  flux.temp[ind1.temp,c("n2o","n2o.0")] <- sub.temp$n2o.ha.day
  
  yr_end.temp <- data.frame(date=yr_end[!(yr_end %in% days.temp)], 
                            n2o=0, n2o.0=0)
  
  # For the first day of each year, interpolate the fluxes of the two
  # adjoining days. This simplifies calendar-year integration.
  for(j in 1:nrow(yr_end.temp)){
    ti.temp <- yr_end.temp[j,"date"]
    t0.temp <- days.temp[max(which(days.temp <= ti.temp))]
    t1.temp <- days.temp[max(which(days.temp <= ti.temp))+1]
    f0.temp <- as.vector(flux.temp[flux.temp$date==t0.temp, c("n2o","n2o.0")])
    f1.temp <- as.vector(flux.temp[flux.temp$date==t1.temp, c("n2o","n2o.0")])
    fi.temp <- f0.temp + (f1.temp-f0.temp) * 
      (as.numeric(ti.temp-t0.temp)/as.numeric(t1.temp-t0.temp))
    yr_end.temp[j,c("n2o","n2o.0")] <- fi.temp
  }
  
  flux.temp <- rbind(flux.temp, yr_end.temp)
  flux.temp <- flux.temp[order(flux.temp$date),]
  
  ###---###---###---###---###
  # Gapfill approach
  
  # Subset all days from a site for which we have any data. They already have
  # site average fluxes recorded!
  flux2.temp <- gf_flux[gf_flux$site==site.temp,c("date","n2o","n2o.0")]
  
  # Add in observed data where we have it 
  # Index of dates with field data
  ind1.temp <- match(sub.temp$date, flux2.temp$date)
  flux2.temp[ind1.temp,c("n2o","n2o.0")] <- sub.temp$n2o.ha.day
  
  yr_end2.temp <- data.frame(date=yr_end[!(yr_end %in% flux2.temp$date)], 
                             n2o=0, n2o.0=0)
  
  # For the first day of each year, interpolate the fluxes of the two
  # adjoining days. This simplifies calendar-year integration.
  for(j in 1:nrow(yr_end2.temp)){
    ti.temp <- yr_end2.temp[j,"date"]
    t0.temp <- flux2.temp[max(which(flux2.temp$date <= ti.temp)),"date"]
    t1.temp <- flux2.temp[max(which(flux2.temp$date <= ti.temp)) + 1,"date"]
    f0.temp <- as.vector(flux2.temp[flux2.temp$date==t0.temp, c("n2o","n2o.0")])
    f1.temp <- as.vector(flux2.temp[flux2.temp$date==t1.temp, c("n2o","n2o.0")])
    fi.temp <- f0.temp + (f1.temp-f0.temp) * 
      (as.numeric(ti.temp-t0.temp)/as.numeric(t1.temp-t0.temp))
    yr_end2.temp[j,c("n2o","n2o.0")] <- fi.temp
  }
  
  flux2.temp <- rbind(flux2.temp, yr_end2.temp)
  flux2.temp <- flux2.temp[order(flux2.temp$date),]
  
}