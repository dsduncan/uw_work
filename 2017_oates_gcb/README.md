**Environmental factors function as constraints on soil nitrous oxide fluxes in 
bioenergy feedstock cropping systems**  
David S. Duncan, Lawrence G. Oates, Ilya Gelfand, Neville Millar, G. Philip 
Robertson, Randall D. Jackson  
_GCB Bioenergy, 2018_ 
[https://doi.org/10.1111/gcbb.12572](https://doi.org/10.1111/gcbb.12572)

## Overview ##

## Scripts ##

### `00_useful_functions.R` ###

Contains useful functions for handling data processing.

### `01_fluxQC_cleanup.R` ###
 
We created a largely manual/visual inspection process for trace gas fluxes. This
 process employed an online tool designed to simplify the process of visualizing
 concentration time-series. Individual concentration values could be removed 
using the tool's GUI. The censored data could then be downloaded and processed 
normally. The workflow was as follows:

1. Upload raw data to FluxQC
2. Visually inspect data for obvious outliers (more detailed criteria in the 
   document "Protocol for Gas Flux Determination")
3. Download inspected data and calculate fluxes
4. Force data to be linear if:
	1. The HMR package's algorithm identifies the data as linear
 	2. There are only 3 valid observations in the dataset
	3. The HMR package identifies the data as HMR, but the slope for HMR is 
       within the 95% confidence interval for the linear fit
5. Re-upload all data to FluxQC
6. Conduct second round of visual inspection, cognizant of how the fluxes were 
   estimated the first time. In essence, this is a second opportunity to rule 
   against the use of HMR by eliminating a data point with too much apparent 
   leverage.

There are two quirks in this process:

1. There was a batch download process that pulled down *all* data that David 
   Duncan had initially uploaded. This would include raw data and second-round 
   data. Raw data uploaded by KBS, on the other hand, needed to be downloaded 
   one sampling date at a time (data are typically uploaded by sampling dates).
2. This process had been carried out for the earlier part of the dataset 
   (2009-2012) for a previous manuscript, and as such did not need to be run 
a second time.

This entire process is not easily recreated/re-run, so I'm pulling aside the
code for this for documentation purposes, rather than to attempt to re-run
the data.