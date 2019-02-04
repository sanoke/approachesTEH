#######################################
#                                     #
#   resultProcessing.R                #
#   Contains helper functions needed  #
#   for simulation to run.            #
#                                     #
#######################################

# start from a clean slate
rm(list=ls())

# the number of subgroups that the data were partitioned in to
numGrp <- 10

# strings that make it easier to process the result files
prefix1 <- "res/simNonLinear-study-"
prefix2 <- "-n1500-numsim100-scen"
prefix3 <- "-seed1-numGrp"
prefix4 <- ".RData"

# - setting some data parameters that also ease the 
#   processing of the result files

# the subgroup-specific treatment effect for each of the 
#   8 treatment groups
tmtEff <- c(1,2,5,5,6,6,9,10)

# the values of (E1, E2, E3) that yield the treatment effect 
#   values above
tmtEff.lbl <- c("(0,0,1)", "(1,0,1)", "(0,1,1)", "(0,0,0)",
                "(1,1,1)", "(1,0,0)", "(0,1,0)", "(1,1,0)")
                
# subgroup labels (to use in producing figures, etc)                
tmtEff.full <- paste0(tmtEff, ":", tmtEff.lbl)
tmtEff.full.noEff <- paste0(5, ":", tmtEff.lbl)
tmtEff.lvls <- c(paste0("TE(",1:10,")"), "undef") # lump undefined TE grps together
tmtEff.lvls.FS <- c(paste0("TE(",1:15,")"))

# - result processing loop: values to loop through
scenario <- c("A","B","C","D", "D*")
method <- c("GBM", "BART", "FS-bt")
numsim <- 100

# initializing the lists that will hold the summary values
res1 <- list() # divide by numsim
res2 <- list() # divide by membership count
res3 <- list() # divide some by numsim, some by membership count

# expected cell proportions
expCell0 <- matrix(1e-20, nrow=length(tmtEff), ncol=length(tmtEff.lvls)-1)
expCell  <- matrix(0, nrow=length(tmtEff), ncol=length(tmtEff.lvls)-1)
expCell0[1, 1:2]  <- expCell[1, 1:2] <- c(150,38) / 188
expCell0[2, 2:3]  <- expCell[2, 2:3] <- c(112,76) / 188
expCell0[3, 3:6]  <- expCell[3, 3:6] <- c(37,76,73,2) / 188
expCell0[4, 3:6]  <- expCell[4, 3:6] <- c(37,76,73,2) / 188
expCell0[5, 5:8]  <- expCell[5, 5:8] <- c(2,73,76,37) / 188
expCell0[6, 5:8]  <- expCell[6, 5:8] <- c(2,73,76,37) / 188
expCell0[7, 8:9]  <- expCell[7, 8:9] <- c(76,112) / 188
expCell0[8, 9:10] <- expCell[8, 9:10] <- c(38,150) / 188

sensitivity <- list()

# for debugging (e.g., loop over just one iteration)
# m <- "GBM"
# scen <- "D"
# j <- 1
 

# # # calculate the percentage of 'expected' matches
for(m in method) {
	
	for(scen in c("A","B","C","D", "D*")) {
			
		    # initializing method-specific summary objects
			rowPerc <- matrix(0, nrow=length(tmtEff), ncol=length(tmtEff.lvls))
			simCount <- numeric(length(tmtEff.lvls))
			
			print(paste(c(m,scen)))
			
			if(scen == "D*") {
				setwd("~/GitHub/approachesTEH/mainSimStudy-nonLinearY")					
				scen <- "D"
				nonLinear <- TRUE
				prefix1 <- "res/simNonLinear-study-"
			} else {
				setwd("~/GitHub/approachesTEH/mainSimStudy")
				prefix1 <- "res/simStudy-"
				nonLinear <- FALSE
			}			
				        
			for(j in 1:numsim) { # calculate the rowPerc matrix and sum across sim iterations

				filename <- paste0(prefix1, m, "-n1500-numsim", j, "-scen", 
								   scen, prefix3, numGrp, prefix4)
				load(filename)
				
				iterCount <- numeric(length(tmtEff.lvls))				            
				            
				# count the number of undefined groups
				numUndef <- sum( is.na(estim$TE) )
				
				# assign people to their new grouping (by order stats)
				#   (observations within a subgroup w/ an undefined treatment
				#    group are reassigned; only subgroups with a defined treatment
				#    group have an integer assignment.)
				newGrp <- paste0("TE(",estim$PSgrp - numUndef,")")	
				
				# turn 'newGrp' into a factor variable, to ease averaging across simulations
				newGrp <- factor(newGrp, levels=tmtEff.lvls)
				newGrp[ is.na(newGrp) ] <- "undef"
				
				# create a new 'trueGrp' variable that shows what group each observation belongs to
				ds$trueGrp2 <- paste0(ds$trueGrp, ":(",ds$Z1, ",", ds$Z2, ",", ds$Z3, ")")	
				
				# if there's no effect modification, there's only one treatment level 
				if(scen == "A") {
					ds$trueGrp2 <- factor(ds$trueGrp2, levels=tmtEff.full.noEff)
				} else {
					ds$trueGrp2 <- factor(ds$trueGrp2, levels=tmtEff.full)
				}	
				
				# construct tabulation w/ row percentages
				rowPerc <- rowPerc + round(table(ds$trueGrp2, newGrp) / 
				                     apply(table(ds$trueGrp2, newGrp),1,sum),3)
				                     
				# count which groups had observations in this simulation iteration
				iterCount[as.numeric(unique(newGrp))] <- 1 
				simCount <- simCount + iterCount  
				
				# print(c("FS ", mean(estim[[j]]$n0 + estim[[j]]$n1)))
				
			} # - END, simulation iterations
			
			if(nonLinear) scen <- "D*"
					
			# group assignment percentage, averaged over all simulations	
			res1[[m]][[scen]] <- round(rowPerc / numsim, 3)
			res1[[m]][[scen]][ , which(simCount==0) ] <- NA  # fix division by 0
			colnames(res1[[m]][[scen]]) <- paste0(tmtEff.lvls, ":", simCount)
			#print(round(res1[[m]][[scen]]*100,1))
			
			# group assignment percentage, averaged only over the simulations where 
			#   there was membership in the group
			res2[[m]][[scen]] <- round(rowPerc %*% diag(1/simCount), 3)
			res2[[m]][[scen]][ , which(simCount==0) ] <- NA  # fix division by 0
			colnames(res2[[m]][[scen]]) <- paste0(tmtEff.lvls, ":", simCount)
			#print(round(res2[[m]][[scen]]*100,1))
			
			# group assignment percentage, averaged only over the simulations where 
			#   there was membership in the group. undefined groups are averaged
			#   over 100. but because of the way FS is coded, there are no 
			#   undefined groups (balance is forced)
			res3[[m]][[scen]] <- round(rowPerc %*% diag(1/simCount), 3)
			res3[[m]][[scen]][ , which(simCount==0) ] <- NA  # fix division by 0
			colnames(res3[[m]][[scen]]) <- paste0(tmtEff.lvls, ":", simCount)
			#print(round(res3[[m]][[scen]]*100,1))
		
	} # - END, scenarios
} # - END, methods


# colors that were used to generate the images
# the function below produces a 3-column matrix of RGB values
# (if you give it a matrix, it gives the values by column)
RGBcolors <- col2rgb(rainbow(30,end=4/6)[1:22])
colnames(RGBcolors) <- paste0("cut", round(seq(0,80,length.out=22)))
colorMapping <- colorRamp(rainbow(30,end=4/6)[1:22])

distanceFunc <- function(x,y) {
	cellDist <- apply( (x-y)^2 , 1, sum ) 
	return(mean(sqrt(cellDist)))
}

distance0   <- colorMapping( matrix(0.10, nrow=8, ncol=10) ) # homogeneity
distance100 <- colorMapping( expCell )                       # expected cell distribution

for(m in c("BART", "GBM", "FS-bt")) {
	index1 <- which(c("BART", "GBM", "FS-bt")==m)
	for(scen in c("A","B","C","D","D*")) {
		index2 <- which(c("A","B","C","D","D*") == scen)
		print(paste(m,scen))
		colorGrid <- colorMapping(res2[[m]][[scen]][,1:10])
		
		if(scen == "A") {
			# distance from 'homogeneity' 
			print( 1 - distanceFunc( colorGrid, distance0 )   / distanceFunc( distance0, distance100 ))
			# print( 1 - distanceFunc( colorGrid, distance0 ) )
		} else {
			# distance from the pattern we expect (this value used in manuscript)
			print( 1 - distanceFunc( colorGrid, distance100 ) / distanceFunc( distance0, distance100 ))
			# print( 1 - distanceFunc( colorGrid, distance100 ) )
		}
		
	}
}