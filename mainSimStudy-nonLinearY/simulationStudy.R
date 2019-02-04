#######################################
#                                     #
#   simulationStudy.R                 #
#   Contains code to run a complete   #
#   simulation.                       #
#                                     #
#######################################

# - The commented-out code below is code that may be helpful 
# - if you are running this simulation on a computing cluster.

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
	print("No arguments supplied.")
	# specify defaults
	# n         <- 1500
	scenario  <- "C"
	numsim    <- 100
	seedVal   <- 1
	procedure <- "FS-bt" 
	numGrp    <- 10
} else {
	scenario    <- as.numeric(args[[1]])
	scenario    <- LETTERS[scenario]
	numsim      <- as.numeric(args[[2]])
	seedVal     <- as.numeric(args[[3]])
	procedure   <- as.character(args[[4]])
	numGrp      <- as.numeric(args[[5]])
	n           <- 1500
}
inputArgs <- paste0("Arguments: n = 1500, scenario = ", scenario, 
                    ", numsim = ", numsim, ", study is ", procedure, 
                    ", and seed is ", seedVal, ".")
print(inputArgs)

# - The `lib.loc` argument of library() may be helpful
# - if you are running this simulation on a computing cluster.

library(MASS)                              
library(BayesTree, lib.loc="~/R/library")  # for BART
library(twang, lib.loc="~/R/library")      # for GBM
library(cluster, lib.loc="~/R/library")    # for FS (the PAM clustering procedure)

source("helperFcns.R")                     # helper functions

set.seed(seedVal)                          # set the seed 

# - Setting the correct estimation procedure.
# - The value of the numbers assigned to the `study` variable
# -  mean nothing. Later in the script, this variable is used
# -  to identify which of the three TEH estimation procedures 
# -  to apply. 
if(procedure == "GBM")           { 
	study <- 2 
} else if(procedure == "BART")   {
	study <- 3 
} else if(procedure == "FS-bt")  { 
	study <- 8
	source("CIFT-Bootstrap.R")
} else stop("Please specify a valid procedure.")


# - # - # - # - # - # - # - # - # - # - # - # - #


# - Load the appropriate dataset. 
dataFile <- paste0("dataFiles", scenario, "/ds", scenario, numsim, ".Rdata")
load(dataFile)

# - These are variables that are used by the TEH estimation script to 
# -  which variable of the dataset is the treatment variable (`col.trt`), 
# -  which is the outcome (`col.y`), which variables are binary (`cols.cat`),
# -  which variables are continuous (`cols.interval`),  
# -  and which variables in the dataset are the covariates 
# -  (`cols.covars`). 
col.trt <- 1
col.y <- 2
cols.cat <- 7:9 + 2
cols.interval <- 1:6 + 2 
cols.covars <- c(cols.cat, cols.interval)

# - This is where the estimation procedure is applied to the data.
estim <- estimation(ds, procedure=study, numGrp=numGrp)

# - This is where we save the resumts of the estimation procedure.
fileName <- paste0("./res/simNonLinear-study-", procedure, "-n", n, "-numsim", numsim,
	               "-scen", scenario, "-seed", seedVal, "-numGrp", numGrp,
	               ".RData")
	               
save.image(fileName)