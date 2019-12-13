# Chapter 4.01: Initial data setup
#
# I'm using COG counts from minimally-assembled samples, straight from JGI.
# This script takes the JGI files, sets up treatment metadata and performs
# relativization by COG length and housekeeping gene abundance.
#
# David Duncan
# 2016-04-02

# Libraries ----
library(dplyr)

# Functions ----
columnizer <- function(dat, col){
  # Convenience function to return the col'th element in a list
  out <- sapply(dat, function(x) x[col])
  return(out)
}

# Reading in data ---- 
### NOTE ###
# These data were pulled from JGI. There were some errors in sample assignments
# that were manually corrected. These changes are documented in README.txt, in
# the same folder.

counts <- read.csv("../data/starting/abundance_cog.csv")
# Make row names (future column names) COGs names
row.names(counts) <- counts$Func_id  
# Transpose so COGs are columns and samples are rows
counts <- counts %>% select(-(1:2)) %>% t() %>% as.data.frame()

# Extract and parse sample names
# Names are irregular, but they all contain the information needed to identify
# the sample. Note that the name appears twice, with the first instance being
# corrected and used here.
t.full <- row.names(counts)
t.start <- regexpr("G[0-9]{2}", t.full)
t.end <- regexpr("[aw1-5]_1[0-2]", t.full)
t.names <- substr(t.full, t.start, t.end + 3)
t.nchar <- nchar(t.names)
t.bstart <- regexpr("[AK][1-5]", t.names)  # Start position of the block

# Impose uniform naming for cropping systems, and place in a logical order
t.crops <- c("Corn", "Switchgrass", "Miscanthus", "Native.Grass", 
              "Poplar", "Old.Field", "Prairie")

# Gather all of the metadata together, then bind to the COG counts
cogs <- data.frame(orig = t.full,  # Original name, to verify the process
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
  cbind(counts)

row.names(cogs) <- 1:nrow(cogs)

write.csv(cogs, "../data/R outputs/clean_cog_counts.csv", row.names=FALSE)

# Add in sequencing information ----
# On the suspicion that metagenome sequencing depth may impact its position, I
# am including sequence length in this set of metadata.
cogs <- merge(read.delim("../data/starting/sample meta.txt") %>%
                mutate(year = factor(year)),
              cogs %>% select(-1), all = TRUE)
  
# Relativizing ----
# There's a two-step relativization process:
# 
# 1. Relativize everything to COG model length (longer sequences should be hit
# more often). 
# This gives average read depth for each COG.
# 
# 2. Relativze everything to average read depth for 37 single-copy housekeeping
# genes (based on He et al., YEAR? PUB?). 
# This gives average copies per cell.
### READ THIS ###
# This assumes all reads are unassembled. Depths for assembled reads are
# already reported by JGI as average read depth, and as such will be 
# massively undercounted. Given my minimal assembly, this isn't a huge issue
# for me (I hope).

# Read in metadata about COGs, sort out the ones removed in the 2014 update, 
# then make sure everything is sorted properly
cog.meta <- read.delim("../data/starting/COG meta.txt") %>% 
  filter(!outdated) %>%
  arrange(match(names(cogs)[-(1:10)], cog))

# There *has* to be a more elegant way to do this, I just can't figure it out
# Division is element-wise (i.e. down by columns). For some reason, I can't
# just transpose the second matrix/vector, so I need to transpose the first
# so element-wise division works, then transpose back to my original structure.
cogs.rel <- cogs %>% select(1:10) %>% cbind(
  ((cogs %>% select(-(1:10)) %>% t()) / 
     cog.meta$dna_len) %>% t())

# Mean read depth for single copy COGs
t.sngl <- cogs.rel[ , 
                   (cog.meta %>% filter(sngl_cpy))$cog %>% as.character()] %>%
  rowMeans()

# Read depth of all COGS/single-copy read depth == copies per cell!
cogs.rel[ , -(1:10)] <- cogs.rel[ , -(1:10)] / t.sngl

write.csv(cogs.rel, "../data/R outputs/cog_relativized.csv", row.names=FALSE)
