
load("glmnet.RData")

# Libraries ----
library(dplyr)
library(ggplot2)
library(glmnet)
library(permute)
library(reshape2)

# Plotting ----
main.theme <- theme(panel.grid=element_blank(),
                    panel.background=element_rect(color="black", fill="white"),
                    axis.text=element_text(color="black", size=8),
                    axis.title=element_text(size=16),
                    legend.key=element_blank(),
                    legend.position="bottom",
                    strip.background=element_blank(),
                    strip.text=element_text(size=14))

# 15-level allegedly dichromat-friendly palette
cf.pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d")
# 2 colors, for ARL vs KBS
site.pal <- cf.pal[c(11, 2)]

# Functions ----
enet.r <- function(mat, resp, nfold, alpha, ltype="lambda.1se", 
                   rtype="lambda"){
  # Get a lambda value from a cv.glmnet call
  suppressWarnings(
    fit <- cv.glmnet(mat, resp, nfolds=nfold, alpha=alpha, family="gaussian"))
  lambda <- unlist(fit[ltype], use.names=FALSE)
  if(grepl("lambda", rtype)){
    return(lambda)
  } else if(grepl("dev", rtype)){
    rind <- abs(fit$glmnet$lambda - lambda) == 
      min(abs(fit$glmnet.fit$lambda - lambda))
    out <- fit$glmnet.fit$dev.ratio[rind]
    return(out)
  }
}
  
enet.pdev <- function(mat, resp, nfold, alphas, group, nperm=100, ...){
  # Get a bunch of deviance ratios out
  require(permute)
  
  out <- data.frame(alpha=rep(alphas, nperm),
                    pdevs=0)
  draws <- shuffleSet(length(resp), nperm, 
                      how(within = Within(type = "free"), blocks = group))
  pb <- txtProgressBar(0, nrow(out), style=3)
  for(i in 1:nperm){
    t.resp <- resp[draws[i,]]
    for(j in seq_along(alphas)){
      out[length(alphas)*(i-1) + j, "pdevs"] <- 
        enet.r(mat, t.resp, nfold, alphas[j], rtype="dev")
      setTxtProgressBar(pb, length(alphas)*(i-1) + j)
    }
  }
  close(pb)
  return(out)
}

enet.real <- function(mat, resp, nfold, alphas, ltype="lambda.1se", niter=20){
  # General wrapper for glmnet procedure
  
  # Hold a big table o'deviances
  out.dev <- data.frame(alpha=alphas, terms=0, devs=0)
  # Placeholder for coefficients
  t.coef <- list()
  
  # Loop over levels of alpha first
  pb <- txtProgressBar(0, length(alphas)*niter, style=3)
  for(i in seq_along(alphas)){
    
    # We'll run multiple versions of the real cross-validated glmnet and
    # average their lambda values to get a more robust estimate
    lambdas <- vector("numeric", niter)
    for(j in seq_along(lambdas)){
      lambdas[j] <- enet.r(mat, resp, nfold, alphas[i], ltype)
      setTxtProgressBar(pb, niter*(i-1)+j)
    }
    # Once we have a good lambda, we can directly fit a glmnet
    fit <- glmnet(mat, resp, family="gaussian", alpha=alphas[i])
    
    # Contribute to deviance table
    out.dev[i, "devs"] <- fit$dev.ratio[which.min(abs(fit$lambda - mean(lambdas)))]
    
    # Contribute to coefficients table
    t.coef[[i]] <- coef(fit, s=mean(lambdas))
    # Rename the table in the list, to track alphas
    t.coef[[i]]@Dimnames[[2]] <- paste("alpha", alphas[i])
  }
  close(pb)
  t2.coef <- Reduce(function(x,y) {cbind(x, y)}, t.coef)
  out.coef <- suppressWarnings(
    t2.coef[unique(t2.coef@i)+1, ] %>% as.matrix)
  
  out.dev[, "terms"] <- (abs(out.coef) > 0) %>% colSums

  out <- list(out.coef, out.dev)
  return(out)
}

# Data entry ----
all.dat <- read.delim("../data/R outputs/Pot and microb.txt") %>%
  mutate(year=factor(year))

cog.met <- read.delim("../data/R inputs/COG meta.txt", 
                      stringsAsFactors=FALSE) %>%
  filter(!outdated) %>% select(cog:name)

funct.met <- read.delim("../data/R inputs/COG functions.txt") %>%
  filter(funct!="Y") %>% droplevels

# 3.6 Elastic net models! ----
# Run it all together

# Define terms for the run
sites <- levels(all.dat$site)
pots <- c("agg.flux", "peak.flux", "PDR", "PNR")
alphas <- seq(0, 1, 0.1)
nfold <- 14  # leave 2 out for KBS, leave ~5 out for ARL
path <- "glmnet.RData"

# Containers for output
out.coef <- vector("list", length(sites)*length(pots))
out.dev <- c()

# Huge loop, would normally be a function, but I need to be able to save output
# as I go along
for(i in seq_along(sites)){
  # Subset out by site first, to make the model matrix
  # Note that the column range is hard-coded. This is not ideal...
  t.site <- sites[i]
  t.sub <- filter(all.dat, site==t.site)
  t.mat <- model.matrix(as.formula(paste("~ year + trt")), 
                      data=t.sub) %>%
    cbind(as.matrix(t.sub %>% select(-site:-is.redo)))
    
  for(j in seq_along(pots)){
    t.meas <- pots[j]
    t.resp <- log(t.sub[, t.meas])
      
    # Permutations
    print(paste("Permuting", sites[i], pots[j]))
    t.perm <- enet.pdev(t.mat, t.resp, nfold, alphas, nperm=200,
                        group=interaction(t.sub$year, t.sub$trt)) %>%
      aggregate(pdevs~alpha, data=., FUN=quantile, 
                probs=seq(0.8, 0.975, 0.025))
    t.perm <- cbind(site=t.site, measure=t.meas, alpha=t.perm$alpha,
                    data.frame(t.perm$pdevs))
    save.image(path)
    
    # Real data
    print(paste("Fitting", sites[i], pots[j]))
    t.real <- enet.real(t.mat, t.resp, nfold, alphas)
      
    out.coef[[length(pots)*(i-1)+j]] <- 
      data.frame(t.real[[1]], site=t.site, measure=t.meas)
      
    out.dev <- rbind(out.dev, cbind(t.real[[2]][,-1], t.perm))
    save.image(path)  
  }
}

# Figure making ----

### Figure 10 ###
out.dev %>%
  mutate(terms=sapply(terms, function(x){min(x, 132)})) %>%
  filter(measure=="agg.flux") %>%
  ggplot(aes(x=alpha, y=devs)) +
  geom_area(aes(y=X95., fill=site), stat="identity", alpha=0.3) +
  geom_point(aes(size=terms, color=site)) + 
  geom_line(aes(color=site)) +
  scale_color_manual(values=site.pal, guide=FALSE) +
  scale_fill_manual(values=site.pal, guide=FALSE) +
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(y="Proportion of deviance explained",
       size="Number of terms",
       x="") +
  facet_grid(measure ~ site, scales="fixed",
             labeller=as_labeller(c(ARL="ARL~'('*3~years*','~105~samples*')'",
                                    KBS="ARL~'('*1~year*','~27~samples*')'", 
                                    PDR="PDR", PNR="PNR", agg.flux="",
                                   peak.flux="Peak fluxes"),
                                  default=label_parsed)) +
  main.theme +
  theme(legend.position="right",
        axis.text=element_text(size=12),
        axis.text.x=element_text(color="white"),
        axis.ticks.x=element_blank())
ggsave("../../images/enet mod.png", width=9.5, height=5.6, units="in")

ggsave("../figures/Figure 10.tiff", height=9.5, units="in",
       dpi=600, compression="lzw")

### Figure 11 ###
# Terms by function category

out.counts <- c()

for(t.i in 1:length(out.coef)){
  t.dat <- out.coef[[t.i]] %>% mutate(cog=row.names(.)) %>%
    merge(x=cog.met, y=.) %>% droplevels
  t.name <- t.dat %>% select(alpha.0:alpha.1) %>% names
  
  for(t.j in seq_along(t.name)){
    t.dat2 <- cbind(
      select(t.dat, cog:name, site:measure),
      select_(t.dat, t.name[t.j])) %>% 
      filter(abs(.[6])>0)
    
    t.out <- data.frame(funct.met,
                        site=t.dat2$site[1],
                        measure=t.dat2$measure[1],
                        alpha=t.name[t.j]) %>%
      merge(aggregate(funct.met$funct, by=list(funct.met$funct),
                      FUN=function(x){
                        sum(grepl(x, t.dat2$funct))}),
            by.x="funct", by.y="Group.1", sort=FALSE) %>%
      rename(freq=x)  %>%
      mutate(pfreq=freq/sum(freq))
    out.counts <- rbind(out.counts, t.out)
  }
}

out.counts %>% 
  mutate(alpha=(as.numeric(alpha)-1)/10,
         funct=factor(funct, levels=funct.met$funct[25:1])) %>%
  filter(!is.na(site), measure=="peak.flux") %>%
  ggplot(aes(y=funct, x=alpha, size=pfreq, color=site)) +
  geom_point(alpha=0.8) +
  scale_size_continuous(range=c(0, 9)) +
  scale_color_manual(values=site.pal, guide=FALSE) +
  labs(y="Function category",
       size="Proportion of retained COGs") +
  facet_grid(.~site) +
  main.theme +
  theme(axis.text.y=element_text(size=11))
ggsave("../figures/Figure 11.tiff", 
       dpi=600, compression="lzw")
 
# Denitrification COGs ----
den.cog <- read.delim("../../Duncan Ch3 microbial community/data/starting/denit COGs.txt")

out.counts <- c()

rbind(
  out.coef[[1]] %>% mutate(cog=row.names(.)),
  out.coef[[5]] %>% mutate(cog=row.names(.))) %>%
  merge(x=den.cog, y=.) %>% droplevels %>%
  melt(id.vars=c("cog", "step", "site", "measure")) %>%
  mutate(variable=factor(variable,
                         labels=c(paste0("0.", 0:9), "1.0"))) %>%
  ggplot(aes(x=variable, y=value, color=cog, group=cog)) +
  geom_point(size=3) + geom_line() + geom_hline(yintercept=0) +
  facet_grid(.~site) +
  labs(x="Alpha", y="Coefficient value", color="COG") +
  main.theme + 
  theme(legend.position="right")

ggsave("../../images/enet denit cogs.png", width=9.5, height=5.6, units="in")

for(t.i in 1:length(out.coef)){
  t.dat <- out.coef[[t.i]] %>% mutate(cog=row.names(.)) %>%
    merge(x=den.cog, y=.) %>% droplevels
  t.name <- t.dat %>% select(alpha.0:alpha.1) %>% names
  
  for(t.j in seq_along(t.name)){
    t.dat2 <- cbind(
      select(t.dat, funct, site:measure),
      select_(t.dat, t.name[t.j])) %>% 
      filter(abs(.[4])>0)
    
    t.out <- data.frame(funct.met,
                        site=t.dat2$site[1],
                        measure=t.dat2$measure[1],
                        alpha=t.name[t.j]) %>%
      merge(aggregate(funct.met$funct, by=list(funct.met$funct),
                      FUN=function(x){
                        sum(grepl(x, t.dat2$funct))}),
            by.x="funct", by.y="Group.1", sort=FALSE) %>%
      rename(freq=x)  %>%
      mutate(pfreq=freq/sum(freq))
    out.counts <- rbind(out.counts, t.out)
  }
}



