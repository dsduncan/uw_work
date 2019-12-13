# Intro ----
# Correlating microbial community measures to cumulative N2O fluxes
#
# This will be my second chapter. The focus is on correlating individual
# markers (think individual PLFAs, 16S OTUs or COGs) to N2O emissions within
# the AARS BCSE corn and prairie systems (G01 and G10). Due the the small
# number of independent observations (2 systems x 3 years x 5 reps), I'm using
# elastic net modeling to shrink/eliminate predictors.
# NOte that there was an earlier attempt to use PCA on the data first, then
# feed that into the elastic net model. There's an old script looking at that
# (01a N correlations (lasso).R).
# 
# 2016-04-28
# David Duncan

# Libraries ----
library(dplyr)
library(ggplot2)
library(glmnet)
library(lsmeans)
library(nlme)
library(permute)
library(vegan)

# Functions ----
glmnet.wrap <- function(mat, resp, nfold, alpha) {
  # Automates extraction of deviance ratios and SSE from glmnet fits
  # resp    Response variable
  # mat     Model matrix
  # nfold   Number of folds for cross-validation
  # alpha   Parameter balancing ridge (0.0) and lasso (1.0) approaches
  
  suppressWarnings(
    fit <- cv.glmnet(mat, resp, nfolds=nfold, alpha=alpha, 
                   family = "gaussian", nlambda = 1000))
  
  # Identify the deviance ratio at the closest value to the "optimal" lambda
  # Note we're using the largest value within 1 se of the true "optimum"
  dev <- fit$glmnet.fit$dev.ratio[abs(fit$glmnet$lambda - fit$lambda.1se) == 
                                min(abs(fit$glmnet.fit$lambda - fit$lambda.1se))]
  
  # I don't use SSE much, so why calculate, right?
  #pred <- predict(fit, mat)
  #sse <- sum((resp - pred) ^ 2)
  out <- list(deviance = dev)
  return(out)
}

glmnet.alphas <- function(mat, resp, nfold, alphas=seq(0.5, 1, 0.1)){
  # Run glmnets over a range of alphas, returning deviance/SSE
  # alphas    Range of alphas to calculate
  out2 = data.frame(alpha = rep(0, length(alphas)), deviance = 0)

    for(i in seq_along(alphas)){
    out <- glmnet.wrap(mat, resp, nfold, alphas[i])
    out2[i, ] <- c(alphas[i], out)
    }
  return(out2)
}

glmnet.perm <- function(mat, resp, nfold = NULL, group=NULL, nperm=200, 
                        alphas=seq(0.5, 1, 0.1)){
  # Permutation testing of glmnetfits over a range of alphas
  # Note that permutations assume year as a factor
  # group   Grouping factor used for permutations
  # nperm   Number of permutations
  
  nalpha <- length(alphas)
  # Hold output
  out3 <- data.frame(alpha = rep(0, nalpha * nperm),
                     deviance = 0)
  
  # Number of items
  if(is.null(nfold)){
    nfold = length(resp)
  }
  
  # Rules for permutation, mostly dealing with blocking factors
  if(is.null(group)){
    ctrl <- how()
  } else{
    ctrl <- how(within = Within(type = "free"), blocks = group)
  }
  
  # Create set of random draws for permutation
  draws <- shuffleSet(length(resp), nperm, ctrl)
  
  # Loop over random draws (with progress bar!)
  pb <- txtProgressBar(0, nperm, style=3)
  for(i in 1:nperm){
    t.resp <- resp[draws[i, ]]
    out3[(nalpha * (i - 1) + 1:nalpha), ] <- glmnet.alphas(mat, t.resp, nfold, alphas)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out3)
}

quant.calc <- function(fits, perms, metric = "deviance"){
  # Convenience function for calculating the quantile of a fit in a set of 
  # permutations
  # fits    Output of glmnet.alphas with real data
  # perms   Output of glmnet.perm (permuted dataset)
  # metric  Metric to comapre ("deviance" or "SSE")
  
  out <- rep(0, nrow(fits))
  for(i in 1:nrow(fits)){
    alpha <- fits$alpha[i]
    out[i] <- sum(perms[perms$alpha == alpha, metric] >= fits[i, metric]) / 
      sum(perms$alpha == alpha)
  }
  names(out) <- fits$alpha
  return(out)
}

raw.wrap <- function(df, exclude=1:6, alphas=seq(0.5, 1.0, 0.05), resp="lflux",
                     trts, dat){
  
  # Automate glmnet fitting and permutation testing in one pretty package
  # df        Input dataframe
  # exclude   Vector of columns that DO NOT contain biomarker data
  # alphas    Vector of alpha (mixing parameter) values
  # resp      Name of the response vector in the dataframe
  # trts      List of names of treatments (used for subsetting)
  # dat       Data type (e.g. PLFA, COGs, etc)
  out.dev <- cbind(expand.grid(alphas, trts), dat=dat, terms=0, deviance=0, 
                   quant=0, pctl_80=0, pctl_85=0, pctl_90 = 0, pctl_95=0)
  names(out.dev)[1:2] <- c("alpha", "trt")
  out.coef <- list()
  
  for(i in 1:length(trts)){
    # Subset data based on treatment
    dat.sub <- filter(df, trt == trts[i])

    # Filter data to only include markers present it all 3 years
    # Greatest value for each year
    t.filt1 <- aggregate(. ~ dat.sub$year, dat.sub[ , -exclude], 
                         function(x) x[rank(x, ties.method = "first") == 4])[ , -1]
    # Is that greater than 0 in all years?
    t.filt2 <- colSums(t.filt1 > 0) == nrow(t.filt1)
    dat.sub <- dat.sub[ , c(rep(TRUE, max(exclude)), t.filt2)]

    # Create model matrix (note that this should be much faster than the 
    # normal model matrix function)
    t.mat <- model.matrix(as.formula(paste(resp, "~ year")), 
                          data = dat.sub)[ , -1] %>%
      cbind(as.matrix(dat.sub[, -exclude]))

    # Response variable and an index range of output cells to fill
    t.resp <- dat.sub[ , resp]
    t.indx <- length(alphas)*(i-1)+1:length(alphas)
    
    t.coef <- list()
    
    for(j in 1:length(alphas)){
      # Get glmnet fit
      t.fit <- suppressWarnings(
        cv.glmnet(t.mat, t.resp, nfolds=nrow(t.mat), alpha=alphas[j],
                         family="gaussian", nlambda=1000))

      # Extract deviance
      t.dev <- t.fit$glmnet.fit$dev.ratio
      out.dev[t.indx[j], "deviance"] <- 
        t.dev[abs(t.fit$glmnet$lambda - t.fit$lambda.1se) == 
                min(abs(t.fit$glmnet.fit$lambda - t.fit$lambda.1se))]
      
      # Add coefficients to a separate list
      t.coef[[j]] <- coef(t.fit, lambda = "lambda.1se")
      t.coef[[j]]@Dimnames[[2]] <- paste("alpha", alphas[j])
    }
    
    # Combine the lists of coefficients
    t.coef2 <- Reduce(function(x,y) {cbind(x, y)}, t.coef)
    # Keep only those terms that actually have coefficients (most won't)
    out.coef[[i]] <- suppressWarnings(
      t.coef2[unique(t.coef2@i) + 1, ] %>% as.matrix)
    
    # Permute data within year to see how things are looking
    out.perms <- glmnet.perm(t.mat, t.resp, group=dat.sub$year, alphas=alphas,
                             nfold=15, nperm=1000)
    out.dev[t.indx, "terms"] <- 
      if(dim(out.coef[[i]])[2] == 1){1} else
      {(abs(out.coef[[i]]) > 0) %>% colSums}
    out.dev[t.indx, "quant"] <- quant.calc(out.dev[t.indx, ], out.perms)
    out.dev[t.indx, 7:10] <- 
      aggregate(deviance ~ alpha, out.perms, FUN = quantile, 
                probs = seq(0.8, 0.95, 0.05)) %>% select(-alpha) %>% as.matrix
  }
  
  out <- list(deviances = out.dev, coef = out.coef)
  return(out)
}

# Data entry ----

# n2o   Aggregate annual N2O fluxes for corn (G01) and prairie (G10)
# plfa  Abolute abundances of PLFAs (contains fungal markers, with plant and 
#       rare markers removed)
# ssu   16S rRNA SSU (v6-8) amplicon OTUs (relative abundances)
# nos   nosZ amplicon OTUs (relative abundances)
# cog   COGs from Illumina shotgun sequencing (relativized to 30 single-copy
#       housekeeping genes)

n2o <- read.delim("../data/n2o.txt") %>% 
  mutate(year = factor(year), lflux = log(flux.Agg))
plfa <- merge(n2o, 
              read.delim("../data/plfa.txt") %>% mutate(year = factor(year)))
ssu <- merge(n2o, 
             read.delim("../data/ssu.txt") %>% mutate(year = factor(year)))
nos <- merge(n2o, 
             read.delim("../data/nos.txt") %>% mutate(year = factor(year)))
cog <- merge(n2o,
             read.delim("../data/cog.txt") %>% mutate(year = factor(year)))

# Raw markers----
# Path to save workspace so I can run this from home computer
t.path <- "G:/Program Data/Dropbox/Articles/Dissertation/Duncan Ch2 community-
flux correlations/Ch2 Analysis/lasso permutations.RData"

plfa.r <- raw.wrap(plfa, trts=unique(plfa$trt), dat="plfa")
save.image(t.path)
ssu.r <- raw.wrap(ssu, trts=unique(ssu$trt), dat="ssu")
save.image(t.path)
nos.r <- raw.wrap(nos, trts=unique(nos$trt), dat="nos")
save.image(t.path)
cog.r <- raw.wrap(cog, trts=unique(cog$trt), dat="cog")
save.image(t.path)

# Figure 1 N2O by crop & year ----
# In which we compare crop X year effects on annual n2o flux
n2o.lme <- lme(lflux ~ trt * year, random = ~ 1 | block, data = plfa,
               weights = varIdent(form = ~1|trt*year))
n2o.sum <- summary(lsmeans(n2o.lme, ~year*trt))
n2o.sum$letters <- gsub(" ", "", 
                        cld(lsmeans(n2o.lme, pairwise ~ year * trt),
                            sort = FALSE, Letters = letters)$.group)
                 
# Plotting (so much for so little...)
n2o.err <- aes(ymax=exp(lsmean + SE), ymin = exp(lsmean - SE))
n2o.let <- aes(y=exp(lsmean + SE + 0.2), x = trt, group = year, label = letters)
n2o.ylab <- expression(paste("Aggregate emissions (g ", N[2], "O-N ", ha^{-1}, ")"))
n2o.xlab <- "Cropping system"
n2o.llab <- "Year"

n2o.y <- scale_y_log10(expand = c(0,0), breaks = c(100, 500, 1000, 5000, 10000))
n2o.x <- scale_x_discrete(labels = c("Corn", "Prairie"))

ggplot(n2o.sum, aes(x = trt, y = exp(lsmean), fill = year)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  geom_errorbar(n2o.err, width = 0.15, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept=0) +  # This should not have worked...
  geom_text(n2o.let, position = position_dodge(width = 0.9), fontface = "italic",
            size = 6) +
  n2o.y + n2o.x + coord_cartesian(ylim = c(100, 20000)) + 
  labs(x = n2o.xlab, y = n2o.ylab, fill = n2o.llab) +
  scale_fill_brewer(palette="PuOr") +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = "black", fill = "transparent"),
        axis.text = element_text(color = "black", size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20))
ggsave("../figures/Ch 2 Fig 1.tiff", width = 7.5, height = 6, units = "in", 
       compression="lzw", dpi=600)

# Figure 2 modeling results ----

all.r <- rbind(plfa.r[[1]], ssu.r[[1]], nos.r[[1]], cog.r[[1]]) %>% 
  mutate(dat = factor(dat, levels = c("plfa", "ssu", "nos", "cog"), 
                      ordered=TRUE))
levels(all.r$dat) <- c("Lipids", "16S rRNA", "nosZ", "COGs")
levels(all.r$trt) <- c("Corn", "Prairie")

ggplot(all.r, aes(x=alpha, y=deviance)) +
  geom_point(aes(size=terms)) + geom_line() +
  geom_area(aes(y=pctl_80), stat = "identity", alpha = 0.12) +
  geom_area(aes(y=pctl_85), stat = "identity", alpha = 0.12) +
  geom_area(aes(y=pctl_90), stat = "identity", alpha = 0.12) +
  geom_area(aes(y=pctl_95), stat = "identity", alpha = 0.12) +
  facet_grid(dat ~ trt) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "alpha value", y = "Deviance ratio", size = "Number of terms") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color="black", fill="transparent"),
        axis.text = element_text(color = "black", size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill="transparent"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.position = "bottom",
        legend.key=element_blank())
ggsave("../figures/Ch 2 Fig 2.tiff", width = 7, units = "in",height = 9, 
       compression="lzw", dpi=600)
          
# Figure S1 norB vs N2O ----

nor.xlab <- expression("Estimated COG3256 copies" ~ cell^{-1})

ggplot(cog, aes(x = COG3256, y = flux.Agg, color = trt, shape = year)) +
  geom_point(size = 4) + 
  geom_smooth(aes(color = trt, group = trt), method = 'lm', se = FALSE) +
  scale_y_log10() +
  scale_color_brewer(palette="Set1") +
  labs(y = n2o.ylab, x = nor.xlab, color = "System", shape = "Year") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color="black", fill="transparent"),
        axis.text = element_text(color = "black", size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill="transparent"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.key=element_blank())

ggsave("../figures/Ch 2 Fig S1.tiff", width = 7.5, height = 6, units = "in",
       compression="lzw", dpi=600)

