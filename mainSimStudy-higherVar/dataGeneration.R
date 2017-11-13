########################################
#                                      #
#   dataGeneration.R                   #
#   Contains code to generate *numsim* #
#   datasets of size *n*.              #
#                                      #
########################################

set.seed(1)

n      <- 1500
numsim <- 100

library(MASS)

source("helperFcns.R")

# - # - # - # - # - # - # - # - # - # - # - # - #

alpha1  <-  0.1 # weak.T
alpha2  <- -0.1 # weak.T
alpha3  <-  1.1 # strong.T
alpha4  <- -1.1 # strong.T
alpha7  <-  0.4 # moderate.T
alpha8  <- -0.1 # weak.T
alpha9  <-  1.1 # strong.T
alpha10 <- -4   # strongg.T

beta0  <- -3.85
beta1  <-  0.5  # weak.Y
beta2  <- -2    # strong.Y
beta3  <- -0.5  # weak.Y
beta4  <-  2    # strong.Y
beta5  <-  1    # moderate.Y
beta6  <- -1    # moderate.Y
beta7  <-  0
beta8  <- -2    # strong.Y
beta9  <-  1    # moderate.Y
beta10 <-  4    # strongg.Y
beta11 <- -4    # strongg.Y

gamma1 <- 5

varOutcome <- 100      # outcome variance


##  DATA GENERATION + ESTIMATION     
## 
##  In these simulations, I put all available covariates in the model.

for(scenario in LETTERS[1:4]) {
	
	##  PARAMETER VALUES         
	##  here i set fixed values for the simulation.
	
	if(scenario == "A") { 
		# confounding + no effect modification
		
		alpha <- numeric(10)
		names(alpha) <- c("constant", paste0("X",1:6), paste0("Z",1:3))
		alpha["X1"] <- alpha1
		alpha["X2"] <- alpha2
		alpha["X3"] <- alpha3
		alpha["X4"] <- alpha4
		alpha["X5"] <- alpha7
		
		beta  <- numeric(14)
		names(beta) <- c("constant", "T", paste0("X",1:6),  
		                 paste0("Z",1:3), "TZ1", "TZ2", "TZ3")
		beta["constant"] <- beta0
		beta["T"]    <- gamma1 
		beta["X1"]   <- beta1
		beta["X2"]   <- beta2
		beta["X3"]   <- beta3
		beta["X4"]   <- beta4
		beta["X6"]   <- beta5
		
	} else if(scenario == "B") { 
		# effect modification + no confounding
		
		alpha <- numeric(10)
		names(alpha) <- c("constant", paste0("X",1:6), paste0("Z",1:3))
		alpha["X5"] <- alpha7
	
		beta  <- numeric(14)
		names(beta) <- c("constant", "T", paste0("X",1:6),  
		                 paste0("Z",1:3), "TZ1", "TZ2", "TZ3")
		beta["constant"] <- beta0
		beta["T"]    <- gamma1 
		beta["X6"]   <- beta5
		beta["Z1"]   <- beta6
		beta["Z2"]   <- beta7
		beta["Z3"]   <- beta8
		beta["TZ1"]  <- beta9
		beta["TZ2"]  <- beta10
		beta["TZ3"]  <- beta11 
	} else if(scenario == "C") {
		# effect modification + confounding (not by EMs)
		
		alpha <- numeric(10)
		names(alpha) <- c("constant", paste0("X",1:6), paste0("Z",1:3))
		alpha["X1"] <- alpha1
		alpha["X2"] <- alpha2
		alpha["X3"] <- alpha3
		alpha["X4"] <- alpha4
		alpha["X5"] <- alpha7
		
		beta  <- numeric(14)
		names(beta) <- c("constant", "T", paste0("X",1:6), 
		                 paste0("Z",1:3), "TZ1", "TZ2", "TZ3")
		beta["constant"] <- beta0
		beta["T"]    <- gamma1 
		beta["X1"]   <- beta1
		beta["X2"]   <- beta2
		beta["X3"]   <- beta3
		beta["X4"]   <- beta4
		beta["X6"]   <- beta5
		beta["Z1"]   <- beta6
		beta["Z2"]   <- beta7
		beta["Z3"]   <- beta8
		beta["TZ1"]  <- beta9
		beta["TZ2"]  <- beta10
		beta["TZ3"]  <- beta11
	} else if(scenario == "D") {
		# effect modification + confounding (EMs assoc w trt)
		alpha <- numeric(10)
		names(alpha) <- c("constant", paste0("X",1:6), paste0("Z",1:3))
		alpha["X1"] <- alpha1
		alpha["X2"] <- alpha2
		alpha["X3"] <- alpha3
		alpha["X4"] <- alpha4
		alpha["X5"] <- alpha7
		alpha["Z1"] <- alpha8
		alpha["Z2"] <- alpha9
		alpha["Z3"] <- alpha10
		
		beta  <- numeric(14)
		names(beta) <- c("constant", "T", paste0("X",1:6),  
		                 paste0("Z",1:3), "TZ1", "TZ2", "TZ3")
		beta["constant"] <- beta0
		beta["T"]    <- gamma1 
		beta["X1"]   <- beta1
		beta["X2"]   <- beta2
		beta["X3"]   <- beta3
		beta["X4"]   <- beta4
		beta["X6"]   <- beta5
		beta["Z1"]   <- beta6
		beta["Z2"]   <- beta7
		beta["Z3"]   <- beta8
		beta["TZ1"]  <- beta9
		beta["TZ2"]  <- beta10
		beta["TZ3"]  <- beta11
	} 
	
	ds.temp <- lapply(rep(n, numsim), datagen, alpha=alpha, beta=beta)
	for(i in 1:numsim) { 
		ds <- ds.temp[[i]]
		save("ds", file=paste0("./dataFiles", scenario, "/ds",scenario, i,".Rdata")) 
	}
}