#######################################
#                                     #
#   resultProcessing-forest.R         #
#   Contains helper functions needed  #
#   for simulation to run.            #
#                                     #
#######################################

library(ggplot2)

# start from a clean slate
rm(list=ls())

# function needed for the image to display properly
rotate <- function(x) t(apply(x, 2, rev))

# the number of subgroups that the data were partitioned in to
numGrp <- 10

# strings that make it easier to process the result files
prefix1 <- "res/simLargeVar-study-"
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
# scenario <- c("A","B","C","D")
scenario <- c("D")
scenarioDesc <- c("confounding + no EM",
                  "EM + no confounding",
                  "EM + confonding",
                  "EM + confounding by EMs")
method <- c("GBM", "BART", "FS-bt")
numsim <- 100

# initializing the lists that will hold the summary values
res1 <- list() # divide by numsim
res2 <- list() # divide by membership count
res3 <- list() # divide some by numsim, some by membership count

# for debugging (e.g., loop over just one iteration)
# m <- "GBM"
# scen <- "D"
# j <- 1


# by estimated treatment effect
for(scen in scenario) {

	for(m in method) {
		
		estimates <- NULL
		print(paste(c(m,scen)))
			
		for(j in 1:numsim) { # calculate the rowPerc matrix and sum across sim iterations
			
			filename <- paste0(prefix1, m, "-n1500-numsim", j, "-scen", 
								   scen, prefix3, numGrp, prefix4)
			load(filename)
			
			# count the number of undefined groups 
			numUndef <- sum(is.na(estim$TE))

			# assign people to their new grouping (by order stats)
			#   (observations within a subgroup w/ an undefined treatment
			#    group are reassigned; only subgroups with a defined treatment
			#    group have an integer assignment.)
			#   (recall -- smallest group numbers are undefined, so we are going to shift
			#    the numbering so Grp 1 is the smallest estimated TE)
			newGrp <- paste0("TE(",estim$PSgrp - numUndef,")")	
			
			# turn 'newGrp' into a factor variable, to ease averaging across simulations
			newGrp <- factor(newGrp, levels=tmtEff.lvls)
			
			# those with undefined TE, set as undefined
			if(numUndef > 0) {
				for(i in 1:numUndef) {
					# code to lump the undefined groups together
					newGrp[estim$PSgrp == i] <- "undef"
					# code for separate undefined groups
					#newGrp[estim$PSgrp == i] <- paste0("undef", i)
				}				
			}
			
			# create a new 'trueGrp' variable that shows what group each observation belongs to
			ds$trueGrp2 <- paste0(ds$trueGrp, 
			                      ":(",ds$Z1, ",", ds$Z2, ",", ds$Z3, ")")	
			# if there's no EM, there's only one treatment level 
			if(scen == "A") {
				ds$trueGrp2 <- factor(ds$trueGrp2)
			} else {
				ds$trueGrp2 <- factor(ds$trueGrp2, levels=tmtEff.full)
			}			
			
			# save TE estimates and the estimated grouping
			estimates.temp <- cbind( sort(estim$TE) , 
			                         levels(newGrp)[1:sum(!is.na(estim$TE))] )  			                         
			estimates <- rbind( estimates, estimates.temp ) 
		}

		estimates <- data.frame(TE = as.numeric(estimates[,1]), grp = estimates[,2] )
		estimates$grp <- factor(estimates$grp, levels=c(paste0("TE(",1:10,")")))
		
		
		res1[[m]][[scen]] <- estimates
	
		# count the number of sims in each group
		counts <- summary(estimates$grp)
		counts <- paste0("(", round(counts/numsim*100), "%)")
		
		index1 <- which(LETTERS[1:4]==scen) 
	
		ggplot(estimates, aes(x=grp, y=TE, fill=grp)) + 
		    geom_hline(yintercept=c(1,2,5,6,9,10), linetype="dotted") +
	        stat_boxplot(geom ='errorbar', stat_params = list(width = 0.5), lwd=0.001) +
	        geom_boxplot(lwd=0.3, outlier.size=1, fatten=1) +
	        guides(fill=FALSE) + coord_flip() + 
	        stat_summary(fun.y=mean, geom="point", shape=20, size=1, color="white") +
	        theme_classic() +
	        theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), 
	              axis.title.y = element_blank(),
	              axis.line.x = element_line(lineend="round")) +
	        # scale_fill_hue(l=40) +      # darker
	        scale_x_discrete(limits = rev(levels(estimates$grp))) +
	        ylab("Treatment Effect") +
	        ylim(0,13) +
	        ggtitle(paste0(m, " (", scen, ": ", scenarioDesc[index1],") n=1500")) +
	        #annotate("text", x=numGrp:1, y=12.5, label=counts, colour="blue", size=2) +
	        geom_hline(yintercept=c(1,2,5,6,9,10), linetype="dotted")
		ggsave(paste0("forestPlots/estTrtEff-n", n, "-numsim100-",m,"-",scen, "star.png"), 
			    width=6, height=3, units="in", dpi=800)
	}
}


# by the true treatment effect
for(scen in scenario) {

	for(m in method) {
		
		print(paste(c(m,scen)))
		
		estimates <- NULL

		# load the estimation results for this method	
		for(j in 1:numsim) { # calculate the rowPerc matrix and sum across sim iterations
			
			filename <- paste0(prefix1, m, "-n1500-numsim",j,"-scen", scen, prefix3, numGrp, prefix4)
			load(filename)
			
			ds$newGrp <- paste0("TE(",estim$PSgrp - numUndef,")")	
			
			# turn 'newGrp' into a factor variable, to ease averaging across simulations
			ds$newGrp <- factor(ds$newGrp, levels=tmtEff.lvls)
			
			# those with undefined TE, set as undefined
			if(numUndef > 0) {
				for(i in 1:numUndef) {
					# code to lump the undefined groups together
					ds$newGrp[estim$PSgrp == i] <- "undef"
					# code for separate undefined groups
					#newGrp[estim$PSgrp == i] <- paste0("undef", i)
				}				
			}
			
			# create a new 'trueGrp' variable that shows what group each observation belongs to
			ds$trueGrp2 <- paste0(ds$trueGrp, 
			                    ":(",ds$Z1, ",", ds$Z2, ",", ds$Z3, ")")	
			
			# if there's no EM, there's only one treatment level 
			if(scen == "A") {
				ds$trueGrp2 <- factor(ds$trueGrp2)
			} else {
				ds$trueGrp2 <- factor(ds$trueGrp2, levels=tmtEff.full)
			}		
			
			# list everyone's estimated treatment effect
			estTrt <- sort(estim$TE, na.last=TRUE)
			estGrp.counts <- summary(ds$newGrp)
			ds$estTE <- NA
			for(i in tmtEff.lvls[1:10]) {
				index2 <- which(tmtEff.lvls == i)
				if(estGrp.counts[index2] == 0) next
				ds$estTE[ ds$newGrp == i ] <- estTrt[index2]
			}
			
			# now calculate TE estimates based on the true grouping
			lvl <- levels(ds$trueGrp2)
			estimates.temp <- rep(NA, length(lvl))
			for(i in lvl) {
				index3 <- which(lvl == i)
				estimates.temp[index3] <- mean( ds$estTE[ ds$trueGrp2  == i ] , na.rm=TRUE )
			}
			
			# save TE estimates and the true grouping
			estimates.temp <- cbind( tmtEff.full , estimates.temp )  			                         
			estimates <- rbind( estimates, estimates.temp ) 
		}

		estimates <- data.frame(TE = as.numeric(estimates[,2]), grp = estimates[,1] )
		estimates$grp <- factor(estimates$grp, levels=tmtEff.full)
		
		
		res1[[m]][[scen]] <- estimates
	
		# count the number of sims in each group
		counts <- summary(estimates$grp)
		counts <- paste0("(", round(counts/numsim*100), "%)")
		
		index1 <- which(LETTERS[1:4]==scen)
	
		ggplot(estimates, aes(x=grp, y=TE, fill=grp)) + 
		    geom_hline(yintercept=c(1,2,5,6,9,10), linetype="dotted") +
	        stat_boxplot(geom ='errorbar', stat_params = list(width = 0.5), lwd=0.001) +
	        geom_boxplot(lwd=0.3, outlier.size=1, fatten=1) +
	        guides(fill=FALSE) + coord_flip() + 
	        stat_summary(fun.y=mean, geom="point", shape=20, size=1, color="white") +
	        theme_classic() +
	        theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), 
	              axis.title.y = element_blank(),
	              axis.line.x = element_line(lineend="round")) +
	        # scale_fill_hue(l=40) +      # darker
	        scale_x_discrete(limits = rev(levels(estimates$grp))) +
	        ylab("Treatment Effect") +
	        ylim(0,13) +
	        ggtitle(paste0(m, " (", scen, ": ", scenarioDesc[index1],") n=1500")) 
	        #annotate("text", x=8:1, y=12.5, label=counts, colour="blue", size=2) +
		ggsave(paste0("forestPlots/trueTrtEff-n", n, "-numsim100-",m,"-",scen, "star.png"), 
			    width=6, height=3, units="in", dpi=800)
	}
}