# Modeling fluxes from Zip-Lignin incubations
#
# Giving this SoilR package a try...


# Libraries
library(FME)
library(lattice)
library(MASS)
library(SoilR)

# Functions
eCO2func=function(pars){
  mod=TwopFeedbackModel(
    t=days,
    ks=pars[1:2],
    a21=pars[3]*pars[1],
    a12=pars[4]*pars[2],
    C0=Ctotal*c(pars[5],1-pars[5]),
    In=0,
    pass=TRUE
  )
  AccR=getAccumulatedRelease(mod)
  return(data.frame(time=days,eCO2=rowSums(AccR)))
}

eCO2cost=function(pars){
  modelOutput=eCO2func(pars)
  return(modCost(model=modelOutput, obs=BorealCO2, err="eCO2sd"))
}


# Data import

#
BorealCO2=subset(eCO2, subset=Sample=="AK_T25", select=-Sample)
names(BorealCO2)<-c("time","eCO2","eCO2sd")

plot(BorealCO2[,1:2], xlab="Days", ylab="Evolved CO2 (mgC g-1 soil)")
arrows(BorealCO2[,1],BorealCO2[,2]-BorealCO2[,3],BorealCO2[,1], 
       BorealCO2[,2]+BorealCO2[,3],code=3,angle=90,length=0.1)

days=seq(0,42)
Ctotal=7.7
inipars=c(k1=0.5,k2=0.05,alpha21=0.5,alpha12=0.1,gamma=0.5)
eCO2fit=modFit(f=eCO2cost,p=inipars,method="Marq",
               upper=c(Inf,Inf,1,1,1),lower=c(0,0,0,0,0))
fitmod=eCO2func(eCO2fit$par)
pairs(eCO2mcmc)

eCO2mcmc <- modMCMC(eCO2cost, eCO2fit$par, niter=100)
