# Setting up initial data
#
# I use KO counts from minimally-assembled samples, straight from JGI, with
# relativization to average copy numbers via the MUSiCC algorithm (site).
#
# David Duncan
# 2016-05-09
 
# Libraries ----
library(dplyr)

# Functions ----
columnizer <- function(dat, col){
  # Convenience function to return the col'th element in a list
  out <- sapply(dat, function(x) x[col])
  return(out)
}

# Initial screen ---- 
### NOTE ###
# These data were pulled from JGI. There were some errors in sample assignments
# that were manually corrected. These changes are documented in README.txt, in
# the same folder.
# I don't need to keep the function name data (I can always look it up later)

counts <- read.delim("../data/starting/JGI KO counts.txt") %>%
  select(-Func_name)

# Remove samples that are missing USiCCs, as these were misannotated by IMG
# NOTE: I'm hoping to correct this at a later date!
keepers <- which((filter(counts, Func_id == "K02357") %>% select(-Func_id)) > 0)
counts <- select(counts, Func_id, keepers+1)

# Remove KOs with fewer than 17 median reads across all samples. This is the
# most stringent cutoff that allowed us to retain all nitrification/
# denitrification KOs
# This screens out ~72% of KOs, hopefully removing a lot of potential noise
counts <- counts[(apply(counts[,-1], 1, median) >= 17), ]

# Old code for seeing effects of cutoffs
# "K10944" %in% counts$Func_id
# median(t(counts[counts$Func_id == "K02567", -1]))

# Now that we have a smaller, better dataset, let's write it out and run it
# through MUSiCC (http://elbo.gs.washington.edu/software_musicc.html)
# We'll do both inter and intra sample variation correction
write.table(counts, "../data/R outputs/KO preMUSiCC.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

# Cleaning up normalized KO data ----
counts.rel <- read.delim("../data/starting/KO post MUSiCC.txt")
kos <- as.character(counts.rel$Func_id)
samples <- names(counts.rel)[-1]

# Transpose and stick column names back on
ko.mat <- t(counts.rel[ ,-1])
row.names(ko.mat) <- 1:nrow(ko.mat)
colnames(ko.mat) <- kos

# Parse sample names
t.start <- regexpr("G[0-9]{2}", samples)
t.end <- regexpr("[aw1-5]_1[0-2]", samples)
t.names <- substr(samples, t.start, t.end + 3)
t.nchar <- nchar(t.names)
t.bstart <- regexpr("[AK][1-5]", t.names)  # Start position of the block

# Impose uniform naming for cropping systems, and place in a logical order
t.crops <- c("Corn", "Switchgrass", "Miscanthus", "Native.Grass", 
              "Poplar", "Old.Field", "Prairie")

# Gather all of the metadata together, then bind to the COG counts
ko.dat <- data.frame(orig = samples,  # Original name, to verify the process
                      # Only need the one letter to denote site
                      site = factor(substr(t.names, t.bstart, t.bstart), 
                                    labels=c("ARL", "KBS")),
                      year = paste0("20", 
                                    substr(t.names, t.nchar - 1, t.nchar)),
                      crop = factor(substr(t.names, 1, 3),
                                    labels=t.crops),
                      rep = substr(t.names, t.bstart, t.bstart + 1),
                      # All redos have an underscore before block
                      is.redo = grepl("[0-9]_[0-9][AK]", t.names),
                      # Microplots denoted with "w"
                      is.micro = grepl("w", t.names)) %>%
  mutate(orig = as.character(orig)) %>%
  cbind(ko.mat)

row.names(ko.dat) <- 1:nrow(ko.dat)

# Add in sequencing information
# On the suspicion that metagenome sequencing depth may impact its position, I
# am including sequence length in this set of metadata.
ko.dat <- merge(read.delim("../data/starting/sample meta.txt") %>%
                mutate(year = factor(year)),
                ko.dat[ , -1],
                all.y = TRUE)

# Write it all out!
write.table(ko.dat, "../data/R outputs/KO copy numbers.txt", 
            row.names=FALSE, quote=FALSE, sep="\t")

