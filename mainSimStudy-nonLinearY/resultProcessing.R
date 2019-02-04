#######################################
#                                     #
#   resultProcessing.R                #
#   Contains code for processing the  #
#   simulation study results.         #
#                                     #
#######################################

# start from a clean slate
rm(list=ls())

# function needed for the image to display properly
rotate <- function(x) t(apply(x, 2, rev))

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
scenario <- c("A","B","C","D")
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


# # # calculate the percentage of 'expected' matches
for(m in method) {
	
	for(scen in c("A","B","C","D")) {
			
		    # initializing method-specific summary objects
			rowPerc <- matrix(0, nrow=length(tmtEff), ncol=length(tmtEff.lvls))
			simCount <- numeric(length(tmtEff.lvls))
			
			# 10 November 2017
			# No longer including simulation A in the figure; 
			# including nonlinear version of D
			if(scen == "A") {
				setwd("~/GitHub/approachesTEH/mainSimStudy-nonLinearY")					
				scen <- "D"
				nonLinear <- TRUE
				prefix1 <- "res/simNonLinear-study-"
			} else {
				setwd("~/GitHub/approachesTEH/mainSimStudy")
				prefix1 <- "res/simStudy-"
				nonLinear <- FALSE
			}
			
			print(paste(c(m,scen)))
				        
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
				#   (recall -- smallest group numbers are undefined, so we are going to shift
			    #    the numbering so Grp 1 is the smallest estimated TE)
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

#test  <- res[["FS"]][["A"]]
#image(rotate(test), col=rainbow(30,end=4/6)[1:20])

# put all the tables into one large matrix
matchVals1 <- matchVals2 <- matchVals3 <- matrix(NA, 
                                                 nrow=length(tmtEff)*length(method), 
                                                 ncol=4*length(tmtEff.lvls))
for(m in c("BART", "GBM", "FS-bt")) {
	index1 <- which(c("BART", "GBM", "FS-bt")==m)
	for(scen in c("B","C","D","D*")) {
		index2 <- which(c("B","C","D","D*") == scen)
		matchVals1[(8*(index1-1)+1):(8*index1) , (11*(index2-1)+1):(11*index2)] <- res1[[m]][[scen]] 
		matchVals2[(8*(index1-1)+1):(8*index1) , (11*(index2-1)+1):(11*index2)] <- res2[[m]][[scen]] 
		matchVals3[(8*(index1-1)+1):(8*index1) , (11*(index2-1)+1):(11*index2)] <- res3[[m]][[scen]] 	
	}
}

matchVals1 <- round(matchVals1*100,1)
matchVals2 <- round(matchVals2*100,1)
matchVals3 <- round(matchVals3*100,1)
rownames(matchVals1) <- rownames(matchVals2) <- rownames(matchVals3) <- rep(method, each=8)
colnames(matchVals1) <- colnames(matchVals2) <- colnames(matchVals3) <- rep(paste("Scen",c("B","C","D","D*")), each=11)

numCol <- 22 # (helps in determining the color scale)

filename <- paste0("~/GitHub/approachesTEH/mainSimStudy-nonLinearY/mainSimStudy-annot.png")
png(filename, width=21, height=9, res=800, units="in")
	# op <- par(xaxt="n", yaxt="n", mar=c(5, 4, 4, 4) + 0.1)
	op <- par(xaxt="n", yaxt="n", mar=c(5, 4, 4, 4) + 0.1)
	image(rotate(matchVals2), col=rainbow(30,end=4/6)[1:numCol])
	abline(h=c(0.98,2.02)/3, v=c(0.98,2,3.02)/4, lwd=3)
	#mtext("Scenario A", side=3, line=1.2, font=2, at=0.12, cex=1.7)
	#mtext("confounding + no EM", side=3, line=0.3, at=0.12, cex=1.2)
	mtext("Scenario B", side=3, line=1.2, font=2, at=0.12, cex=1.7)
	mtext("EM + no confounding", side=3, line=0.3, at=0.12, cex=1.2)
	mtext("Scenario C", side=3, line=1.2, font=2, at=0.38, cex=1.7)
	mtext("EM + confounding", side=3, line=0.3, at=0.38, cex=1.2)
	mtext("Scenario D", side=3, line=1.2, font=2, at=0.64, cex=1.7)
	mtext("EM + confounding by EMs", side=3, line=0.3, at=0.64, cex=1.2)
	mtext("Scenario D*", side=3, line=1.2, font=2, at=0.89, cex=1.7)
	mtext("nonlinear confounding", side=3, line=0.3, at=0.89, cex=1.2)
	mtext(c("BART", "GBM", "FS-bt"),
	      side=2,
	      line=0.1,
	      at=rev(c(0.15, 0.50, 0.85)),
	      cex=1.7)	
	mtext(rep(1:11,4), side=1, at=seq(0,1,length=44), line=0.2)
	yVals <- seq(1,0,length=nrow(matchVals2))
	for(j in 1:nrow(matchVals2)) {
		text(seq(0,1,length=44), yVals[j], round(matchVals2[j,]), cex=0.8)
	}
	par(op)
	op <- par(xpd=TRUE)
	points(1:numCol*0.0194 + 0.29,
	       rep(-0.11,numCol), 
	       pch=15, col=rainbow(30,end=4/6,alpha=0.7)[1:numCol], cex=5)
	#text(0.335, -0.14, "smallest", font=2)
	#text(0.675, -0.14, "largest", font=2)
	text(seq(0.0194+0.29, numCol*0.0194+0.29, length=6) , 
	     -0.11, round(seq(min(matchVals2,na.rm=TRUE),
	                max(matchVals2,na.rm=TRUE),
	                length=6)), font=2)
	#text(1.03, yVals, c(paste0(5,":",tmtEff.lbl),rep(tmtEff.full, 3)), cex=1, font=2)		
	text(1.035, yVals, rep(tmtEff.full, 3), cex=1, font=2)		
dev.off()


filename2 <- paste0("~/GitHub/approachesTEH/mainSimStudy-nonLinearY/mainSimStudy.png")
png(filename2, width=21, height=9, res=800, units="in")
	# op <- par(xaxt="n", yaxt="n", mar=c(5, 4, 4, 4) + 0.1)
	op <- par(xaxt="n", yaxt="n", mar=c(5, 4, 4, 4) + 0.1)
	image(rotate(matchVals2), col=rainbow(30,end=4/6)[1:numCol])
	abline(h=c(0.98,2.02)/3, v=c(0.98,2,3.02)/4, lwd=3)
	#mtext("Scenario A", side=3, line=1.2, font=2, at=0.12, cex=1.7)
	#mtext("confounding + no EM", side=3, line=0.3, at=0.12, cex=1.2)
	mtext("Scenario B", side=3, line=1.2, font=2, at=0.12, cex=1.7)
	mtext("EM + no confounding", side=3, line=0.3, at=0.12, cex=1.2)
	mtext("Scenario C", side=3, line=1.2, font=2, at=0.38, cex=1.7)
	mtext("EM + confounding", side=3, line=0.3, at=0.38, cex=1.2)
	mtext("Scenario D", side=3, line=1.2, font=2, at=0.64, cex=1.7)
	mtext("EM + confounding by EMs", side=3, line=0.3, at=0.64, cex=1.2)
	mtext("Scenario D*", side=3, line=1.2, font=2, at=0.89, cex=1.7)
	mtext("nonlinear confounding", side=3, line=0.3, at=0.89, cex=1.2)
	mtext(c("BART", "GBM", "FS-bt"),
	      side=2,
	      line=0.1,
	      at=rev(c(0.15, 0.50, 0.85)),
	      cex=1.7)	
	yVals <- seq(1,0,length=nrow(matchVals2))
	par(op)
	op <- par(xpd=TRUE)
	points(1:numCol*0.0194 + 0.29,
	       rep(-0.11,numCol), 
	       pch=15, col=rainbow(30,end=4/6,alpha=0.7)[1:numCol], cex=5)
	#text(0.335, -0.14, "smallest", font=2)
	#text(0.675, -0.14, "largest", font=2)
	text(seq(0.0194+0.29, numCol*0.0194+0.29, length=6) , 
	     -0.11, round(seq(min(matchVals2,na.rm=TRUE),
	                max(matchVals2,na.rm=TRUE),
	                length=6)), font=2)	
dev.off()
