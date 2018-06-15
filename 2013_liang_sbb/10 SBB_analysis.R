
#### Libraries ####
library(vegan)

#### Data import ####
env_dat <- read.delim("data/chaot.tsv")


setwd("C:/Users/David/Dropbox/Projects/Articles/Chao short comm")
chaot <- read.table("chaot.tsv", sep = "\t", header = TRUE)
chaor <- read.table("chaor.tsv", sep = "\t", header = TRUE)
chao.resp <- chaot[,4:6] # response (amino sugars)
chao.plfa <- chaot[,15:20] # functional-group plfas
chao.dist <- vegdist(chao.resp, method = "euclidian") # amino sugar distance matrix
chao.dist2 <- vegdist(chao.plfa, method = "euclidian") # plfa distance matrix

# RDA: Amino Sugars
rda.env <- rda(chao.resp ~ TN + Clay + Age + pH + Burning + CropSp, scale = TRUE, data = chaot)
# Note TC and Sand removed, compare to below:
# rda.env1 <- rda(chao.resp ~ TN + TC + Clay + Sand + Age + pH + Burning + CropSp, scale = TRUE, data = chaot)

RsquareAdj(rda.env)
# RsquareAdj(rda.env1)
summary(rda.env)
anova(rda.env, by = "margin", step = 1000)
# anova(rda.env1, by = "margin", step = 1000)
anova(rda.env, by = "axis")
plot(rda.env)

#RDA: PLFAs
rda.plfa <- rda(chao.plfa ~ TN + Clay + Burning + CropSp + Age + pH, scale = TRUE, data = chaot)
RsquareAdj(rda.plfa)
summary(rda.plfa)
anova(rda.plfa, by = "margin", step = 1000)
anova(rda.plfa, by = "axis")
plot(rda.plfa)

#RDA: Everything together
rda.all <- rda(chao.resp ~ TN + Clay + Age + Burning + CropSp + pH + GMp + GMm + SF + AMF + Act, scale = TRUE, data = chaot)
RsquareAdj(rda.all)
anova(rda.all, by = "margin", step = 1000)
anova(rda.all, by = "terms", permutation = 999)
anova(rda.all, by = "axis")
plot(rda.all)
summary(rda.all)

# GLMs
# GluN
glu.glm.0 <- glm(GluN ~ as.factor(Burning) + as.factor(CropSp) + Age + TC + TN + pH + Sand + Clay + GMp + GMm + SF + AMF + Act + Protozol, data = chaot)
# glu.glm.1 <- glm(GluN ~ 1, data = chaor)
glu.glm.1 <- glu.glm.2
glu.glm.2 <- glm(GluN ~ TN + Clay-1, data = chaot)
anova(glu.glm.1, glu.glm.2, test = "F")
add1(glu.glm.2, glu.glm.0, k = log(29))

RsquareAdj(glu.glm.2)
summary(glu.glm.2)

# GalN
gal.glm.0 <- glm(GalN ~ as.factor(Burning) + as.factor(CropSp) + Age + TC + TN + pH + Sand + Clay + GMp + GMm + SF + AMF + Act + Protozol, data = chaot)
#gal.glm.1 <- glm(GalN ~ 1, data = chaor)
gal.glm.1 <- gal.glm.2
gal.glm.2 <- glm(GalN ~ TC + Age-1, data = chaot)
anova(gal.glm.1, gal.glm.2, test = "F")
add1(gal.glm.2, gal.glm.0, k = log(29))

summary(gal.glm.2)
RsquareAdj(gal.glm.2)

#MurA
mur.glm.0 <- glm(MurA ~ as.factor(Burning) + as.factor(CropSp) + Age + TC + TN + pH + Sand + Clay + GMp + GMm + SF + AMF + Act + Protozol, data = chaor)
#mur.glm.1 <- glm(MurA ~ 1, data = chaot)
mur.glm.1 <- mur.glm.2
mur.glm.2 <- glm(MurA ~ Clay + TC -1, data = chaot)
anova(mur.glm.1, mur.glm.2, test = "F")
add1(mur.glm.2, mur.glm.0, k = log(29))

summary(mur.glm.2)
RsquareAdj(mur.glm.2)

# Total
tot.glm.0 <- glm(I(GluN + GalN + MurA) ~ as.factor(Burning) + as.factor(CropSp) + Age + TC + TN + pH + Sand + Clay + GMp + GMm + SF + AMF + Act + Protozol, data = chaot)
# tot.glm.1 <- glm(I(GluN + GalN + MurA) ~ 1, data = chaot)
tot.glm.1 <- tot.glm.2
tot.glm.2 <- glm(I(GluN + GalN + MurA) ~ TC , data = chaot)
anova(tot.glm.1, tot.glm.2, test = "F")
add1(tot.glm.2, mur.glm.0, k = log(29))

summary(tot.glm.2)
RsquareAdj(tot.glm.2)

# Correlation
cor.test(chaor$GlucN, chaor$Total.N, method = "pearson")
cor.test(chaor$TC, chaor$TN, method = "pearson")
cor.test(chaor$Sand, chaor$Clay, method = "pearson")
cor.test(chaor$GluMur, chaor$FB, method = "pearson")
cor.test(chaor$GluGal, chaor$FB, method = "pearson")
