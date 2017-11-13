#######################################
#                                     #
#   helperFcns.R                      #
#   Contains helper functions needed  #
#   for simulation to run.            #
#                                     #
#######################################

expit <- function(x) { exp(x) / (1 + exp(x)) }
logit <- function(x) { log( x / (1 - x) ) }


#   this function generates one data set, using the parameter
#     values specified in the argument
datagen <- function(n,          # sample size
                    alpha,      # treatment mean coefficients (vector of length 10)
                    beta        # outcome mean coefficients (vector of length 14)
                    ) {
	
	baseCovars <- mvrnorm(n, numeric(6), diag(6))
		
	Z1 <- rbinom(n, 1, 0.5)
	Z2 <- rbinom(n, 1, 0.5)
	Z3 <- rbinom(n, 1, 0.5)
	
	allCovars <- cbind(baseCovars, Z1, Z2, Z3)
	colnames(allCovars) <- c( paste0("X",1:6), paste0("Z", 1:3) )

	pTreat <- expit( cbind(1, allCovars) %*% alpha )
	trt <- rbinom(n, 1, pTreat)
	
	# intercept, treatment, interaction terms
	Y0mean <- (beta["constant"] + 
	           beta["X1"] * allCovars[,"X1"] + 
	           beta["X2"] * allCovars[,"X2"] + 
	           beta["X3"] * allCovars[,"X3"] + 
	           beta["X4"] * allCovars[,"X4"] +
	           beta["X6"] * allCovars[,"X6"] + 
	           beta["Z1"] * allCovars[,"Z1"] + 
	           beta["Z2"] * allCovars[,"Z2"] + 
	           beta["Z3"] * allCovars[,"Z3"] )^2 +
	           beta["T"]  * 0 + 
	           beta["TZ1"] * 0 * allCovars[,"Z1"] + 
	           beta["TZ2"] * 0 * allCovars[,"Z2"] + 
	           beta["TZ3"] * 0 * allCovars[,"Z3"] 
    Y0     <- rnorm(n, Y0mean, sqrt(varOutcome))
	
	Y1mean <- (beta["constant"] + 
	           beta["X1"] * allCovars[,"X1"] + 
	           beta["X2"] * allCovars[,"X2"] + 
	           beta["X3"] * allCovars[,"X3"] + 
	           beta["X4"] * allCovars[,"X4"] +
	           beta["X6"] * allCovars[,"X6"] + 
	           beta["Z1"] * allCovars[,"Z1"] + 
	           beta["Z2"] * allCovars[,"Z2"] + 
	           beta["Z3"] * allCovars[,"Z3"] )^2 +
	           beta["T"]  * 1 + 
	           beta["TZ1"] * 1 * allCovars[,"Z1"] + 
	           beta["TZ2"] * 1 * allCovars[,"Z2"] + 
	           beta["TZ3"] * 1 * allCovars[,"Z3"] 
	Y1     <- rnorm(n, Y1mean, sqrt(varOutcome))
	
	Y <- ifelse(trt, Y1, Y0) # observed outcome
	
	trueGrp <- round(Y1mean - Y0mean) # treatment effect group
	
	ds <- data.frame(trt, Y, allCovars, Y0, Y1, trueGrp)
	
	return( ds )
}


estimation <- function(ds,		  # data (in data frame)
					   procedure, # estimation procedure for calc TE
					   numGrp     # number of subgroups to partition the data into
					   ) {		   	
					   	
	if(procedure == 1) { # 1 = initialization
		TE <- var <- n0 <- n1 <- numeric(numGrp)
		PSgrp <- rep(1:numGrp, each=nrow(ds)/numGrp)
		return(list(TE=TE, var=var, n0=n0, n1=n1, PSgrp=PSgrp))
	}
	if(procedure == 2) { # 2 = GBM
		gbm1 <- ps(trt ~ .,                   # predicts treatment from all other vars
	               estimand = "ATE",          # estimate average treatment effect
	               verbose = FALSE,
	               data = ds[,-match(c("Y", "Y0", "Y1", "trueGrp"),
	                              names(ds), nomatch=0)],
	                                          # the dataset dropping the outcome
	               n.trees = 20000,           # runs for 20,000 iterations
	               shrinkage = 0.0005,        # sets the shrinkage parameter (see
	                                          # Appendix B)
	               stop.method = "ks.mean",
	               interaction.depth = 2,     # as recommended in McCaffery et al.
	               bag.fraction = 0.5)        # as recommended in McCaffery et al.
	            
		gbmPS <- round(as.vector(gbm1$ps)$ks.mean.ATE,4)
	
		# divide the propensity score into `numGrp`-iles
		cutpts <- sort(unique(c(0,quantile(gbmPS, seq(1/numGrp,1,1/numGrp)))))
		PSgrp <- as.numeric(cut(gbmPS, cutpts, labels=factor(1:(length(cutpts)-1))))
		
		# calculate treatment effects in each subgrpup
		TE.gbm <- var.gbm <- n0.gbm <- n1.gbm <- rep(NA, length(unique(PSgrp)))
		for(j in 1:length(unique(PSgrp))) {
			Y1data <- ds$Y[ PSgrp == unique(PSgrp)[j] & ds$trt == 1]
			Y0data <- ds$Y[ PSgrp == unique(PSgrp)[j] & ds$trt == 0]
			n0.gbm[j] <- sum(PSgrp == unique(PSgrp)[j] & ds$trt == 0)
			n1.gbm[j] <- sum(PSgrp == unique(PSgrp)[j] & ds$trt == 1)	
			Y1 <- mean( Y1data )
			Y0 <- mean( Y0data )
			TE.gbm[j] <- round(Y1 - Y0, 4)
			var.gbm[j] <- var(Y1data)/length(Y1data) + var(Y0data)/length(Y0data)
		}
		
		# assign group numbers by trt effect size 
		# (smaller effect -> smaller group number, to facilitate coherent plotting)
		# also construct the dissimilarity matrix (deprecated)
		lvls <- unique(PSgrp)[order(TE.gbm, decreasing=TRUE)]
		PSgrp2 <- numeric(nrow(ds))
		dissim <- matrix(0, nrow=n, ncol=numGrp)
		for(j in 1:length(lvls)) {
			PSgrp2[PSgrp == lvls[j]] <- numGrp
			dissim[ PSgrp == lvls[j] , numGrp ] <- 1
			numGrp <- numGrp - 1
		}
		
		return(list(TE=TE.gbm, var=var.gbm, n0=n0.gbm, n1=n1.gbm, PSgrp=PSgrp2, D=dissim, mmt=gbmPS))	   
	} else if(procedure == 3) { # 3 = BART
		# supplementary material at
			# http://amstat.tandfonline.com/doi/suppl/10.1198/jcgs.2010.08162#tabModule
		
			#   estimation of counterfactuals, w vars defined as follows:
			#   1) y.train: the outcome variable for the whole data set.
			#   2) x.train: matrix of the treatment variable (first column)
			#               and predictor variables for the whole data set.
			#   3) x.test:  same as x.train but with everyone's treatment assignment
			#               reversed.
		
		# training data (all the data) 
		xt <- as.matrix(ds[,-match(c("Y", "Y0", "Y1", "trueGrp"), 
		                names(ds), nomatch=0)])
		# data for 'prediction' (all the data, but w/ trt assnmt flipped)
		trt.rev <- ds$trt + 1
		trt.rev[trt.rev == 2] <- 0
		## the treatment var is in the first column
		xp <- cbind(trt.rev, as.matrix(ds[,-match(c("trt", "Y", "Y0", "Y1", "trueGrp"), 
		            names(ds), nomatch=0)]))
		            
		# calling bart()
		bart.tot <- bart(x.train=xt, y.train=ds$Y, x.test=xp, verbose=FALSE)
		
		# estimating ITE
		# bart.tot$yhat.train has a row for every draw from the posterior,
		#  and a column for every person in the sample.
		# bart.tot$yhat.test has a row for every draw from the posterior predictive,
		#  and a column for every person in the test data (so nrow(xp) columns)
		# (the default is 1000 posterior draws)
		diffs <- bart.tot$yhat.train - bart.tot$yhat.test
		
		# the mean of each posterior draw (so the average is taken across people,
		#   so should have k posterior draws of the ATE)
		# (rounding because it makes some data subsetting easier, and because
		#   I only care about TEs up to a certain precision anyway.)
		mndiffs <- round(apply(diffs, 2, mean),3) 
		
		# have to reverse the order here because TE is Y(1) - Y(0). for the control
		#   group, their Y(1) is contained in the test.
		mndiffs[!as.logical(ds$trt)] <- -mndiffs[!as.logical(ds$trt)]
			     
		# next, we want to assign each person to a stratum based on their treatment
		#   effect value. 
		# (20 groups are hard-coded)
		cutpts <- sort(unique(c( min(mndiffs)-0.1 ,
		                quantile(mndiffs, seq(1/numGrp,1,1/numGrp)))))
		BARTgrp <- as.numeric(cut(mndiffs, cutpts, labels=factor(1:(length(cutpts)-1))))
		
		TE.BART <- var.BART <- n0.BART <- n1.BART <- rep(NA, length(unique(BARTgrp)))
		for(j in 1:length(TE.BART)) {
			n0.BART[j] <- sum(BARTgrp == unique(BARTgrp)[j] & ds$trt == 0)
			n1.BART[j] <- sum(BARTgrp == unique(BARTgrp)[j] & ds$trt == 1)	
			TE.BART[j] <- mean( mndiffs[ BARTgrp == unique(BARTgrp)[j] ] )
			var.BART[j] <- var( mndiffs[ BARTgrp == unique(BARTgrp)[j] ] )
		}	
		
		lvls <- unique(BARTgrp)[order(TE.BART, decreasing=TRUE)]
		BARTgrp2 <- numeric(nrow(ds))
		dissim <- matrix(0, nrow=n, ncol=numGrp)
		for(j in 1:length(lvls)) {
			BARTgrp2[BARTgrp == lvls[j]] <- numGrp
			dissim[ BARTgrp == lvls[j] , numGrp ] <- 1
			numGrp <- numGrp - 1
		}
		
		return(list(TE=TE.BART, var=var.BART, n0=n0.BART, n1=n1.BART, PSgrp=BARTgrp2, D=dissim, mmt=mndiffs))	
			            		            
	} else if(procedure == 8) { # 8 = FS, bootstrapped
		
		# make a copy of the dataset; going to make some changes 
		# (would have been better to code the dataset like this in the first place,
		#    but oh well.)
		dat <- ds
		
		# FS code needs each observation to have an ID number
		dat$id <- 1:nrow(ds)
		
		# rename outcome to match CIT coding
		colnames(dat)[match("Y", names(dat), nomatch=0)] <- "y"
		
		col.group <- ncol(dat) + 1
		
		# number of folds for the cross-validation, 
		# or the number of bootstrap samples.
		V <- 100
				
		# proportion of a sample for learning (rest is test)
		# (same fraction used in the paper)
		ratio <- 2/3
		
		# minumum terminal node size
		N0 <- 20
		
		# minimum number of observations in each treatment group
		n0 <- 5
		
		n <- nrow(dat); id <- 1:n
		
		# column indices of variables to split on (depends on structure of dataset)
		split.var <- cols.covars
		
		# column indices of categorical covars (depends on structure of dataset)
		ctg <- cols.cat
		
		method <- "bootstrap"                   # "CV" or "bootstrap"
		bsize <- NULL							# size of "best tree"
		D <- matrix(0, n, n)					# dissimilarity matrix
		concord <- rep(0, (V-1))                # ???
		node.assign <- as.list(1:V)             # ???
		
		calcBootstrap <- function() {
			out <- tryCatch({
				message(paste0("Attempting to grow tree #", v , " on a bootstrap sample."))
				
				tre0 <- grow.CIFT(data=L1, test=L2, min.ndsz=N0, N0=n0, split.var=cols.covars, 
		                          ctg=cols.cat, max.depth=8, mtry=length(cols.covars))
		   		prn <- prune.size.testsample(tre0)
				bsize <- c(bsize, prn$size[4])
		    	btre.BIC <- prn$btree[[4]]
		        senddown.1(data=dat, tre=btre.BIC, char.var=cols.cat)
			},
			error = function(cond) {
				# message(paste0("There was an error growing tree #", v, "."))
				# message("Original message:")
				# message(cond)
				error <- error + 1
				return(NA)
			},
			warning = function(cond) {
				# message(paste0("There was a warning growing tree #", v, "."))
				# message("Original message:")
				# message(cond)
				warning <- warning + 1
				return(NA)
			},
			finally = {})
			
			return(out)
		}
		
		error <- warning <- 0
		counter <- 0
		for (v in 1:V) {
		    print(paste0("v = ", v)) 
		    if (method == "CV") {
		        id.L1 <- sample(id, size=floor(n*ratio), replace=FALSE) 
		        L1 <- dat[sort(id.L1), ]
		        L2 <- dat[!is.element(id, id.L1), ]
		    } else if (method == "bootstrap") {
		        id.L1 <- sample(id, size=n, replace=TRUE)  
		        L1 <- dat[id.L1, ]
		        L2 <- dat
		    } else stop("Not a Legitimate value for {method}!")
		    
		    ##
		    node.assign[[v]] <- calcBootstrap()
		    if( !is.na(node.assign[[v]]) ) { counter <- counter + 1 }
		}
		
		# constructing the D matrix for each bootstrap tree
		D <- as.list(1:V)
		counter2 <- nodeIndex <- 1
		while(nodeIndex <= V) {
			if( !is.na(node.assign[[nodeIndex]]) ) {
				D.v <- as.matrix(dist(node.assign[[nodeIndex]], method = "euclidean", 
			                      diag = FALSE, upper = FALSE)) 
			    D[[counter2]] <- (D.v != 0) 
			    counter2 <- counter2 + 1
			}
			nodeIndex <- nodeIndex + 1
			print(paste("nodeIndex is", nodeIndex, "and counter is", counter2, "."))
		}
		
		
		d <- as.list(1:V)
		D0 <- matrix(0, n, n)
		for (v in 1:V) {
			print(v)
		    # summing the dissimilarity matrices together
		    D0 <- D[[v]] + D0
		}
		
		
		# now that we have the D matrix, we cluster!
		pamRes  <- pam(round(D0/numsim,4), numGrp, diss=TRUE, cluster.only=TRUE, 
		               keep.diss=FALSE, keep.data=FALSE)

		TE.pam <- var.pam <- n0.pam <- n1.pam <- rep(NA, numGrp)

		for(j in 1:numGrp) {
			Y1data <- ds$Y[ pamRes == j & ds$trt == 1]
			Y0data <- ds$Y[ pamRes == j & ds$trt == 0]
			n0.pam[j] <- sum(pamRes == j & ds$trt == 0)
			n1.pam[j] <- sum(pamRes == j & ds$trt == 1)	
			Y1 <- mean( Y1data )
			Y0 <- mean( Y0data )
			TE.pam[j] <- round(Y1 - Y0, 4)
			var.pam[j] <- var(Y1data)/length(Y1data) + var(Y0data)/length(Y0data)
		}
		
		# assign group numbers by trt effect size 
		# (smaller effect -> smaller group number, to facilitate coherent plotting)
		lvls <- seq(1,numGrp,1)[order(TE.pam, decreasing=TRUE)]
		PSgrp2 <- numeric(nrow(ds))
		numGrp0 <- numGrp
		for(j in 1:length(lvls)) {
			PSgrp2[pamRes == lvls[j]] <- numGrp0
			numGrp0 <- numGrp0 - 1
		}
		
		return(list(TE=TE.pam, var=var.pam, n0=n0.pam, 
		            n1=n1.pam, PSgrp=PSgrp2, D=round(D0/numsim,4)))

	} else {
		stop("Please specify a valid procedure.")
	}
}