# Load necessary packages - smcDefaults() is a list of relevant packages
library(devtools)
library(cnpParticles)
library(popbio)
library(scales)
library(pbapply)
library(Rage)
library(popdemo)
library(ggplot2)
library(gganimate)
library(plotly)
library(tidyverse)
library(dplyr)
library(viridis)
library(seewave)
library(npreg)
library(scales)
library(gtools)

# Primary System Properties -------------------------
# System-wide properties
independentRuns <- 25
maxIteration <- 35

particle_Number <- 1000
discretisation_Par <- 25
domain_Filter <- 0.05
dilation_parameter <- rep(0.1, 50)


# General System Properties -------------------------
# IPM Boundaries

t_0 <- IPMmesh(IPMboundary(minValue = 2, maxValue = 8, jointBuffer = 0.2),
               mesh.resoution = 100)
t_1 <- IPMmesh(IPMboundary(minValue = 2, maxValue = 8, jointBuffer = 0.2),
               mesh.resoution = 100)


# Parameters -------------------------
No_1 <- solveParameter("fertInt")
No_2 <- solveParameter("fertSlope")


# Priors -------------------------
# Setting up a prior series

PriorsList <- list()

model_Parameters_1 <- modelParams(No_1, No_2) 
model_Parameters_1[[1]]$prior <- function(x){return(rnorm(x, 1.332, 0.05))}
model_Parameters_1[[2]]$prior <- function(x){return(runif(x, -5, 5))}

for(i in 1:independentRuns){
  PriorsList[[i]] <- model_Parameters_1
}

# Pre-IPM  -------------------------
# Empty dataframe
IPM_complete <- clearParamFrame()

# Assign population parameters - here pulled in as defaults
IPM_complete <- seriesA(IPM_complete)

# Extract the true values
true_IPMparams <- trueParameter(param.name = names(model_Parameters_1),
                                param.dataframe = IPM_complete)

# Dataframe w/ parameters to solve coded as NA
IPM_paramaters <- rmParameter(param.name = names(model_Parameters_1),
                              param.dataframe = IPM_complete)

# Map parameters into a list for the IPM model to use
IPM_pre <- IPMparameters(param.dataframe = IPM_paramaters)


# Reference parameterisation -------------------------
# Generate the "true" (known/a priori)  Integral Population Model
IPM_observed <- generateIPM(IPM.parameters = IPM_pre,
                            ... = true_IPMparams$`true values`)

# Generate the "true"  st/age distribution
SSD_observed <- stableStage(A = IPM_observed)


# Evaluate Param Performance -------------------------
iNestedRegister <- list()

for(k in 1:independentRuns){
  iRegister <- list()
  
  iRegister[[1]]  <- initiateSMC(model.parameters = PriorsList[[k]],
                                 particle.number = particle_Number)
  
  for(i in 1:maxIteration){
    try(iRegister[[i]][,3] <- pbmapply(IPMdistance,
                                       1:particle_Number,
                                       MoreArgs = list(IPM.pre = IPM_pre,
                                                       observed.stage.dist = SSD_observed,
                                                       evaluate.params = iRegister[[i]])))
    
    try(prob_Surface <- multiSurface(perf.dataframe = iRegister[[i]],
                                     div = discretisation_Par))
    
    try(new_Domain <- subDomain(surface.object = prob_Surface,
                                cut.threshold = domain_Filter))
    
    try(new_Sample <- sirSampling(new_Domain,
                                  dilation_parameter[[i]]))
    
    try(iRegister[[i+1]] <- new_Sample)
    try(print(paste("local", i, "|", "iteration", k)))
  }
  iNestedRegister[[k]]<- iRegister
}

save(iNestedRegister, file="confirm.RData")

