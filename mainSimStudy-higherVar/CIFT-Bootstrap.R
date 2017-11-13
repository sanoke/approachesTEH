#  - - - - - - - - - - - - - - - - - - 
#  FACILITATIONG SCORE CODE 
#  as adapted from code emailed at
#  request by X. Su, used in the 
#  analysis of Su et al (2012)
#  - - - - - - - - - - - - - - - - - - 

# ==========================================================================================  #
#  TREE-RELATED FUNCTIONS FOR SUBGROUP ANALYSIS - THE MAIN PART
# ==========================================================================================  #


# ==================================================================
#  loglik() 
#           node : data frame that contains observations in the node
#           N0   : the minimum number of observations in each tmt grp
#  RETURNS THE LOG-LIKELIHOOD SCORE FOR EACH NODE (a scalar)
#  (see equation (8) in paper for derivation)
# ==================================================================

loglik <- function(node, N0=5) {
	# if() statement checks each element of the first argument to see
	#     if it's in colnames() -- returns a two-dimensional vector
	#     with TRUE/FALSE 
    if (sum(is.element(c("y", "trt"), colnames(node)))!=2) 
        stop("Please make sure both y and trt columns in data!!!") 
    y <- node$y; trt <- node$trt
    loglik <- NULL
    n <- length(y); n1 <- sum(trt); n0 <- n-n1 
    if (min(n1, n0) >= N0){
        sse <-  ((anova(lm(y~trt)))$"Sum Sq")[2]
        loglik <- -n/2*log(n*sse) + n1*log(n1) + n0*log(n0)
    }  
    return(loglik)  
}



# ==================================================================
#  power.set() 
#           x : a categorical varialbe.
#  RETURNS THE POWER SET FOR A CATEGORICAL VARIABLE
#  (as an object where each element is a set from the larger
#   power set)
# ==================================================================

power.set <- function(x) {
	if(length(x) == 0) {
		return(vector(mode(x), 0))	
    }
	x <- sort(unique(x)); n <- length(x); K <- NULL
    for(m in x) K <- rbind(cbind(K, F), cbind(K, T))
    # creation of the complete power set
    out <- apply(K, 1, function(x, s) s[x], s = x)
    # drop the null set, and the set with everything
    out <- out[-c(1, length(out))]
    l <- length(out); i <- 1
    while (i <= l) {
    	if (length(out[[i]]) >= ceiling(n/2+.5)) {out[[i]] <- NULL; i <- i-1}
        i <- i+1 
        l <- length(out) 
	}
	# an object where each element is of the 
	#    power set {levels of x}
    out
}



# ==============================================================
# partition.CIFT() 
#           dat : data frame to bisect
#           test : test data to apply split to; LRT test statistics 
#                  are calculated.
#           name : nchar(name) is the maximum tree depth 
# 	 			   (i.e., number of splits)
#           split.var : variables to consider for splitting on
#           ctg : vector indicating column indices of categorical variables
#           min.ndsz: minimum node size
#           N0 : minumum number of observations in each treatment group
#           max.depth : 
#           mtry : number of variables (out of split.var) to consider
#                  for splitting (take a random sample of split.var)
#  COMPLETES ONE SINGLE SPLIT USING LRT
#  (produces the bisected data)
# ==============================================================

partition.CIFT <- function(dat, test, name="1", 
             split.var, ctg=NULL, 
             min.ndsz=20, N0=5, max.depth=10, mtry=NULL)
{   
	# remember the call to this function
    call <- match.call(); out <- match.call(expand = F)
    # add arguments to the call above
    out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
    n <- nrow(dat); n.test <- nrow(test); var <- varname <- NA; cut <- NA; 
    max.score <- -1e20; 
    vnames <- colnames(dat)
    trt <- dat$trt; y <- dat$y; score.test <- NA;    
    # depth is the number of splits
    depth <- nchar(name) 
    # calculate the loglik of ...?
    LOGLIK0 <- loglik(dat, N0=N0)
    # calculate the loglik of ...?
    LOGLIK.0 <- loglik(test, N0=N0)
    if(depth < max.depth && min(n, n.test) >= min.ndsz & 
            !is.null(LOGLIK0) & !is.null(LOGLIK.0)) {
        # m.try is the number of variables (out of split.var) that we want to consider
        #    for potential splitting; to consider the full list, would set equal to
        #    length(split.var). 
        m.try <- ifelse(is.null(mtry), length(split.var), mtry) 
        # this for() loop looks at all possible splits of all split.var and chooses the
        #    split that yields the largest test statistic.
        # looks at the selected split vars in some random order.
        for(i in sample(split.var, size=m.try, replace=F)) {
            vname <- vnames[i]
            x <- dat[,i]; temp <- sort(unique(x));  
            if(length(temp) > 1) { 
            	# if the var is categorical, calculate the power set;
            	#    otherwise, consider all observed values as a place to cut.
                if (is.element(i,ctg)) zcut <- power.set(temp)
                else zcut <- temp[-length(temp)]
                # for all possible cut points...
                for(j in zcut) {
                    score <- NA
                    # if i denotes the column of a categorical variable...
                    if (is.element(i,ctg)) {
                    	# grp selects out the potential group to be bisected to the left
                    	grp <- sign(is.element(x, j))
                    	cut1 <- paste(j, collapse=" ")}      
                    else  {
                    	grp <- sign(x <= j)
                    	cut1 <- as.character(j)
                    }  
                    node.L <- dat[grp==1, c("y", "trt")]
                    node.R <- dat[grp!=1, c("y", "trt")]
                    loglik.L <- loglik(node.L, N0=N0)
                    loglik.R <- loglik(node.R, N0=N0)                                     
                    # if the LL from each node are valid numbers, then set
                    #   the score equal to the sum.
                    if ((is.null(loglik.L) + is.null(loglik.R)) == 0)  score <- loglik.L + loglik.R
					# if the score is not null and it's bigger than the max, then set
					#   all variables to this new cut.
                    if (!is.na(score) && score >= max.score) {
                    	max.score <- score; var <- i; cut <- cut1; best.cut<-j; varname <- vname
                    } 
    			}
    		}
    	}
    }
    # if we were able to select a split var above, then calculate the test statistic
    #    based on this split, for the test data
    if (!(is.na(var))) {
        if (is.element(var,ctg)) grp.test <- sign(is.element(test[,var], best.cut))
        else  grp.test <- sign(test[,var] <= best.cut) 
        nod.L <- test[grp.test==1, c("y", "trt")]
        nod.R <- test[grp.test!=1, c("y", "trt")]
        loglik.L <- loglik(nod.L, N0=N0)
        loglik.R <- loglik(nod.R, N0=N0)
        if ((is.null(loglik.L) + is.null(loglik.R)) == 0)  score.test <- loglik.L + loglik.R  
    }
    # if we were able to calculate a test statistic for the test data, 
    #    bisect the training and test data 
    if (!is.na(score.test)) {
        out$name.l <- paste(name, 1, sep="")       
        out$name.r <- paste(name, 2, sep="")        
        out$left.test <- test[grp.test==1,  ]
        out$right.test <- test[grp.test==0,  ]
        if (is.element(var,ctg)) {
        	out$left  <- dat[is.element(dat[,var], best.cut),]      
            out$right <- dat[!is.element(dat[,var], best.cut), ]
        } else {
            out$left  <- dat[dat[,var] <= best.cut, ]
            out$right <- dat[dat[,var] > best.cut, ]
        }
        LRT <- 2*(max.score - LOGLIK0)
        LRT.test <- 2*(score.test - LOGLIK.0)
    } else { 
    	var <- varname <- NA; cut <- NA;  LRT <- LRT.test <- NA 
    }
    out$info <- data.frame(node=name, size = n, var = var, cut = cut, 
            score=ifelse(max.score==-1e10, NA, LRT), score.test=LRT.test, 
            size.test=n.test, varname=varname)
    out 
}




# ==============================================================
# grow.CIFT() 
#           datA : data frame to bisect
#           test : test data to apply split to; LRT test statistics 
#                  are calculated.
#           split.var : variables to consider for splitting on
#           ctg : vector indicating column indices of categorical variables
#           min.ndsz: minimum node size
#           N0 : minumum number of observations in each treatment group
#           max.depth : 
#           mtry : number of variables (out of split.var) to consider
#                  for splitting (take a random sample of split.var)
#  CONSTRUCTS TREE
# ==============================================================

grow.CIFT <- function(data, test, min.ndsz, N0=10, split.var, ctg=NULL, 
	max.depth=10, mtry=length(split.var))
{
    out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
    # coerces data and test (both data sets) into lists
    list.nd <- list(data); list.test <- list(test)
    name <- 1
    while (length(list.nd)!=0) {    
      for (i in 1:length(list.nd)) {
        if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1) { 
        split <- partition.CIFT(list.nd[[i]], list.test[[i]], name[i], min.ndsz=min.ndsz, 
                      N0=N0, split.var=split.var, ctg=ctg, max.depth=max.depth, mtry=mtry)
        print(split$info)
        out <- rbind(out, split$info)
        # out <- rbind(out, c(as.numeric(paste(rep(depth.tre,depth.tre), collapse="")), n.split, split$info))
        if (!is.null(split$left) && !is.null(split$left.test)){
            temp.list <- c(temp.list, list(split$left, split$right))
              temp.test <- c(temp.test, list(split$left.test, split$right.test))
            temp.name <- c(temp.name, split$name.l, split$name.r)
        }}}
        list.nd <- temp.list; list.test <- temp.test; name <- temp.name
        temp.list <- temp.test <- temp.name <- NULL
    }   
    out$node <- as.character(out$node)
    out <- out[order(out$node), ]
    out
}





# ====================================================================
# Pruning and Size Selection Based on LeBlanc and Crowley (JASA, 1992)
# ====================================================================

# THE PRUNING ALGORITHM GOOD FOR THE BOOTSTRAP TREE SIZE SELECTION METHOD
# -----------------------------------------------------------------------
prune.size <- function(tre)
{
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre[,4])); subtree <- 1            
     while (n.tmnl > 1 ) {
            # if (n.tmnl==5) {btre <- tre; print(btre)}
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
            G.test <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T)
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, G.test=G.test))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c(3:6, 8)] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, G.test=0))    
      result <- as.data.frame(result)
      result
}




de <- function(x, tree)
{
    if(length(x) != 1) stop("The length of x in function de must be 1.")    
    y <- tree$node;  de <- NA
    if(sum(match(x, y), na.rm = T) != 0) {
        temp <- 1:length(y)
        start <- match(x, y) + 1    
        end <- length(y)
        if(start <= length(y) & nchar(y[start]) > nchar(x)) {
            temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
            if(!is.na(temp1)) end <- temp1
            de <- y[start:end]
    }}
    de
}




# -----------------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE TREE SIZE SELECTION VIA TEST SAMPLE METHOD
# -----------------------------------------------------------------------------

prune.size.testsample <- function(tre)
{
     call <- match.call(); out <- match.call(expand = F)
     out$result <- out$size  <- out$... <- NULL
     ntest <- as.numeric(tre$size.test[1])
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre[,4])); subtree <- 1            
     a <- cbind(Ga.2=2, Ga.3=3, Ga.4=4, Ga.BIC=log(ntest))*4         ###################################### 4 times for CIFT
     max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4) 
     while (n.tmnl > 1 ) {
            # print(tre)
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
            Ga <- G - a*l 
            for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tre}}                        
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, Ga))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c(3:6, 8)] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            if (n.tmnl ==1) {for (k in 1:4){if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tre}}}
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, Ga=cbind(Ga.2=0, Ga.3=0, Ga.4=0, Ga.BIC=0)))     
      result <- as.data.frame(result)
      out$result <- result; out$size <- size; out$btree <- btree
      out 
}






# ==========================================================================================  #
#  FUNCTIONS RELATED TO THE BOOTSTRAP FUNCTION
# ==========================================================================================  #


# OPTION sample.id IS TO CONTROL WHETHER WE SAMPLE ID OR OBSERVATIONS IN BOOTSTRAPPING FOR LONGITUDINAL DATA
# OPTION LeBlanc IS TO APPLY THE WHOLE SAMPLE OR THE OUT-OF-BAD SAMPLE IN THE BOOTSTRAP PROCEDURE
# OPTION min.boot.tree.size IS TO MAKE SURE A NON-NULL TREE IS GROWN AT EACH BOOTSTRAP 11/1/2007

boottrap.grow.prune <- function(B=30, data, N0=20, n0=5, split.var, ctg=1000, 
	max.depth=6, mtry=NULL, sample.id=T, LeBlanc=F, min.boot.tree.size=1)  
{
	# call contains an object with the boottrap.grow.prune() function call.
	# out is the same as call, but omitting arguments matching "..."
	# (so I think that call and out will be the same here)
    call <- match.call(); out <- match.call(expand = F)
    # some sort of initialization
    out$boot.tree <- out$boot.prune <- out$... <- NULL
    # a string that contains the current date and time
    time.start <- date()
    tree0 <- grow.CIFT(data=data, test=data, min.ndsz=N0, N0=n0, split.var=split.var, 
	ctg=ctg, max.depth=max.depth, mtry=mtry)  
    print(tree0);  
    prune0 <- prune.size(tree0); 
    boot.tree <- list(tree0); boot.prune <- list(prune0) 
    b <- 1
    while (b <= B) {
        print(paste("###################### b = ", b, " ###########################", sep=""))
        if (sample.id) {
            # SAMPLING ID
            ID <- sort(unique(data$id))
            samp <- sample(x=ID, size=length(ID), replace=T)
            # samp.oob <- ID[!is.element(ID, samp)]
            dat.oob <- data[!is.element(data$id,unique(samp)), ];
            dat <- NULL
            for (i in samp) dat <- rbind(dat, data[data$id==i,])
        }
        else{
            # SAMPLING OBSERVATION
            samp <- sample(1:nrow(data), size=nrow(data), replace=T) 
            dat <- data[samp, ];     
            dat.oob <- data[-unique(samp),]
        }
        n.oob <- nrow(dat.oob); 
	  # print(n.oob)   ##############################      
        if (LeBlanc) {tre <- grow.CIFT(data=dat, test=data, min.ndsz=N0, N0=n0, split.var=split.var, 
		ctg=ctg, max.depth=max.depth, mtry=mtry)}  
        else {tre <- grow.INT(dat, dat.oob, min.ndsz=N0, N0=n0, split.var=split.var, 
		ctg=ctg, max.depth=max.depth, mtry=mtry)}
        print(tre)        
        if (nrow(tre)> min.boot.tree.size) {
            boot.tree <- c(boot.tree, list(tre)); 
            prune <- prune.size(tre); # print(prune)
            boot.prune <- c(boot.prune, list(prune));
            b <- b+1
        }
    }
    time.end <- date(); 
    print(paste("The Start and End time for ", B, "bootstrap runs is:"))
    print(rbind(time.start, time.end))
    out$boot.tree <- boot.tree
    out$boot.prune <- boot.prune
    out
}   



bootstrap.size <- function(boot.prune, n)
{   
    #  COMPUTE THE ALPHA PRIME'S
    prune0 <- boot.prune[[1]] 
    n.subtree <- nrow(prune0)
    alpha <- as.numeric(as.vector(prune0$alpha));
    temp <- c(alpha[1], alpha[-length(alpha)])
    alpha.prime <- sqrt(alpha*temp)  
    # cbind(alpha,  alpha.prime=prune0$alpha.prime)
    b <- length(boot.prune)
    G <- as.numeric(as.vector(prune0$G)); 
    size.tmnl <- as.numeric(as.vector(prune0$size.tmnl)); 
    subtree <- as.numeric(as.vector(prune0$subtree)); 
    # tree.penalty <- log(nrow(teeth))
    G.a <- matrix(0, n.subtree, 5)
    penalty <- c(0, 2:4, log(n))*4    ########################## four times with CIFT
    for (i in 1:n.subtree) {
        a <- alpha.prime[i]
        bias <- 0
        for (j in 2:b){
            prune.bs <- boot.prune[[j]]
            alpha.bs <- as.numeric(as.vector(prune.bs$alpha)); 
            g <- as.numeric(as.vector(prune.bs$G)); 
            g.test <- as.numeric(as.vector(prune.bs$G.test)); 
            indx <- 1
            if (sum(alpha.bs <=a)>0) {          
                temp1 <- which.max(which(alpha.bs<=a))
                indx <- ifelse(is.null(temp1), 1, temp1)
            }
            temp2 <- (g-g.test)[indx]
            bias <- bias + temp2 
            # print(cbind(i, a, j, bias, indx, temp2))
        }
        G.honest <- G[i] - bias/(b-1) 
        G.a[i,] <- G.honest - penalty*(size.tmnl[i]-1)
    }
    out <- data.frame(cbind(size.tmnl, G.a))
    colnames(out) <- c("tmnl.size", "G", "G.2", "G.3", "G.4", "G.log(n)")
    out
}



# =================================
# SEVERAL OTHER RELATED FUNCTIONS
# =================================

obtain.btree <- function(tre, bsize=6)
{
     n.tmnl <- sum(is.na(tre[,4])); indicator <- T            
     while (n.tmnl > 1 && indicator ==T) {
            if (n.tmnl==bsize) {btre <- tre; print(btre); indicator=F}
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c(3:6, 8)] <- NA
            n.tmnl <- sum(is.na(tre$cut))          
      }
      #btre
      tre
}



# THIS send.down FUNCTION RUNS A TREE STRUCTURE DOWN A DATA SET
send.down <- function(data, tre, char.var=1000)
{
    call <- match.call(); out <- match.call(expand = F)
    out$tree <- out$data <- out$... <- NULL
    dat <- cbind(data, node=1); tre <- cbind(tre, n.test=NA)
    cut.point <- as.vector(tre$cut); 
    split.var <- as.numeric(as.vector(tre$var)); 
    for (i in 1:nrow(tre)){
        in.node <- (dat$node)==(tre$node[i]);
        tre$n.test[i] <- sum(in.node)                       
        if (!is.na(split.var[i])){
            # print(cbind(i, var=tre$var[i], cut=tre$cut[i]))
            var.split <- dat[,split.var[i]]; 
            cut <- cut.point[i]
            if (!is.element(split.var[i], char.var)) { 
                cut1 <- as.numeric(cut)    
                l.nd <- dat$node[in.node & var.split <= cut1] 
                r.nd <- dat$node[in.node & var.split > cut1]
                dat$node[in.node & var.split <= cut1] <- paste(l.nd, 1, sep="")
                dat$node[in.node & var.split >  cut1] <- paste(r.nd, 2, sep="")  
            }
            else {
                var.split <- as.character(var.split)
                cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
                l.nd <- dat$node[in.node & is.element(var.split, cut1)] 
                r.nd <- dat$node[in.node & !is.element(var.split, cut1)]                  
                dat$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")  
                dat$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="")}                   
    }}
    # print(tre)
    out$data <- dat
    out$tree <- tre
    out 
}





# THIS senddown FUNCTION IS WRITTEN FOR THE RANDOM FORESTS PART
 
senddown <- function(data, tre, char.var=1000, n0=5)
{
    call <- match.call(); out <- match.call(expand = F)
    out$tree <- out$data <- out$... <- NULL
    dat <- cbind(data, node=1); tre <- cbind(tre, n2=NA, chi=NA)
    cut.point <- as.vector(tre$cut); 
    split.var <- as.numeric(as.vector(tre$var)); 
    for (i in 1:nrow(tre)){
        in.node <- (dat$node)==(tre$node[i]);
        y <- dat$y[in.node]; trt <- dat$trt[in.node]; id <- dat$id[in.node]; 
        n.0 <- length(y)
        dat.0 <- data.frame(y=y, trt=trt, id=id)
        tre$n2[i] <- sum(in.node); z <- NA                      
        if (n.0 > n0 && !is.na(split.var[i]) ){
            # print(cbind(i, var=tre$var[i], cut=tre$cut[i]))
            var.split <- dat[,split.var[i]]; 
            cut <- cut.point[i]
            if (!is.element(split.var[i], char.var)) { 
                cut1 <- as.numeric(cut)    
                l.nd <- dat$node[in.node & var.split <= cut1] 
                r.nd <- dat$node[in.node & var.split > cut1]
                z <- sign(var.split[in.node] <= cut1)
                dat$node[in.node & var.split <= cut1] <- paste(l.nd, 1, sep="")
                dat$node[in.node & var.split >  cut1] <- paste(r.nd, 2, sep="")  
            }
            else {
                cut1 <- unlist(strsplit(as.character(cut), split=" "))  ############################
                l.nd <- dat$node[in.node & is.element(var.split, cut1)] 
                r.nd <- dat$node[in.node & !is.element(var.split, cut1)]   
                z <- sign(is.element(var.split[in.node], cut1))  
                dat$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")  
                dat$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="")}                   
        }
        # print(dim(dat.0)); print(z); print(length(z))   
        if (!is.na(z)) tre$chi[i] <- ttest(dat.0, z, n0=n0)
    }
    # print(tre)
    out$data <- dat
    out$tree <- tre
    out 
}


# THIS senddown.1 FUNCTION WAS WRITTEN FOR AGGREGATED GROUPING METHOD - 2/23/2011
senddown.1 <- function(data, tre, char.var=NULL)
{
    call <- match.call(); out <- match.call(expand = F)
    out$tree <- out$data <- out$... <- NULL
    dat <- cbind(data, node=1); tre <- cbind(tre, n.test=NA)
    cut.point <- as.vector(tre$cut); 
    split.var <- as.numeric(as.vector(tre$var)); 
    for (i in 1:nrow(tre)){
        in.node <- (dat$node)==(tre$node[i]);
        tre$n.test[i] <- sum(in.node)                       
        if (!is.na(split.var[i])){
            # print(cbind(i, var=tre$var[i], cut=tre$cut[i]))
            var.split <- dat[,split.var[i]]; 
            cut <- cut.point[i]
            if (!is.element(split.var[i], char.var)) {  # for continuous var
                cut1 <- as.numeric(cut)    
                l.nd <- dat$node[in.node & var.split <= cut1] 
                r.nd <- dat$node[in.node & var.split > cut1]
                dat$node[in.node & var.split <= cut1] <- paste(l.nd, 1, sep="")
                dat$node[in.node & var.split >  cut1] <- paste(r.nd, 2, sep="")  
            } else {   # for categorical vars
                var.split <- as.character(var.split)
                cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
                l.nd <- dat$node[in.node & is.element(var.split, cut1)] 
                r.nd <- dat$node[in.node & !is.element(var.split, cut1)]                  
                dat$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")  
                dat$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="")}   
                print(dat$node[1])                
    }}
    # print(tre)
    # if (is.null(data$id)) {id <- 1:nrow(data)} else {id <- data$id}
    # out <- data.frame(id =id, node = dat$node)
    # return(out) 
    return(dat$node)
}






# =============================
# LaTeX Plot of Tree Structure 
# =============================

latex.plot <- function(file="tree-code.tex", tree, group=rep("I", n.node), cat.var=c(3) )   {
  n.node <- nrow(tre)
  sink(file=file)
#  cat("\n \\begin{Tree} \n")
  for (i in n.node:1) {
    node <- tre[i, ] 
    if (is.na(node$var))  {
      cat(paste("\\node{\\external\\type{circle}\\cntr{", node$size, "}\\lft{\\textbf{", group[i], "}}} \n", sep="")) 
    } else {
      if (!is.element(node$var, cat.var)) { 
        cat(paste("\\node{\\type{frame}\\cntr{$\\texttt{V", node$var, "} \\leq ", node$cut, "$}} \n", sep=""))
      } else {
        # cut <- node$cut
        cut1 <- unlist(strsplit(as.character(node$cut), split=" "))
        cat(paste("\\node{\\type{frame}\\cntr{\\texttt{V", node$var, "} $ \\in \\{", paste(cut1, collapse=","), "\\} $}} \n", sep="")) 
    }}
  }
  sink()  
} 

# latex.plot(file="tree-code.tex", tree=tree0, group=rep("I", n.node), cat.var=c(3))

  

# ================================================================
# FUNCTION concordance COMPUTES THE PROPORTION OF CONCORDANT PAIRS
# ================================================================

concordance <- function(y1, y2)
{
   s <- 0;  n <- 0
   for(i in seq(along=y1))
   {
      for(j in seq(along=y1))
      {
     		if(i != j)
     		{
        		if(y1[i] > y1[j])
        		{
           			s <- s + (y2[i] > y2[j]) + 0.5 * (y2[i] == y2[j])
           			n <- n + 1
        		}
     		}
     }
  }
  s / n
}



# =============================
# MODIFY THE Kendall FUNCTION
# =============================

Kendall.S <- function (x, y) 
{
    tau <- 0
    ptau <- 0
    sltau <- 0
    score <- 0
    varscore <- 0
    denom <- 0
    iws <- numeric(length(x))
    ifault <- 0
    outF <- .Fortran("tauk2", as.single(x), as.single(y), as.integer(length(x)), 
        as.single(tau), as.single(ptau), as.single(sltau), as.single(score), 
        as.single(varscore), as.single(denom), as.integer(iws), 
        as.integer(ifault), PACKAGE = "Kendall")
    S <- outF[[7]]
    S
}


# =================================================================================
# FUNCTION extract() EXTRACTS FREQUENCY AND MEAN INFO FROM DATA SET WITH NODE GROUPS. 
# Node groups are indicated by a categorical variable that's part of the dataset.
# =================================================================================

extract <- function(data, y.col=10, trt.col=1, group.col, effect.est=F){
    options(digits=7)
    dat <- data
    if (!is.data.frame(dat)) dat <- as.data.frame(dat)
    colnames(dat)[y.col] <- "y"
    colnames(dat)[trt.col] <- "trt"
    colnames(dat)[group.col] <- "group"
    groups <- sort(unique(dat$group))
    out0 <- NULL
    for (t in groups ){
        dat.t <- dat[dat$group==t, ]
        nt <- nrow(dat.t)
        n1 <- sum(dat.t$trt==1); 
        n0 <- nt - n1
        mean1 <- mean(dat.t$y[dat.t$trt==1])
        sd1 <- sd(dat.t$y[dat.t$trt==1])
        mean0 <- mean(dat.t$y[dat.t$trt==0])
        sd0 <- sd(dat.t$y[dat.t$trt==0])
        if (effect.est==T) {
            # fit <- lm(y~trt + ps  + educ, data=dat.t) 
		fit <- lm(y~trt, data=dat.t) 
            b <- coef(fit)[2]
            se <- sqrt((summary(fit)$cov.unscaled)[2,2]) * summary(fit)$sigma 
            t0 <- b/se
            df0 <- fit$df.residual
            p.value <- 2*pt(abs(t0), df=df0, lower.tail =F) 
            print(paste("===================== Node ", t, " ====================", sep="")) 
            print(cbind(b, se, t0, df0, p.value))
        }
        out0 <- rbind(out0, cbind(t=t, n1=n1, mean1=mean1, sd1=sd1, 
                n0=n0, mean0=mean0, sd0=sd0))  
    }
    out0 <- as.data.frame(out0)
    colnames(out0) <- c("group", "n1", "mean1", "sd1", "n0", "mean0", "sd0")
    return(out0)
}


# =========================================== END ===============================================  #
