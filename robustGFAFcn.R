

##########       robustGFA        ##########
#####  Author: Henry Yeh, Ph.D.   #####
#####  The following functions are developed to combine multiple   #####
#####  replicates of GFA and obtain robust components and loadings #####


### 1. robustLoadings() ###
# The GFA package includes a function robustComponents() that aims to match components 
# across replicates based on "recontructed values" or projections (X times t(W)) and 
# adjusts signs to achieve positive correlations in projections across replicates.  
# Consequently, the signs of loadings are not guarantteed to be consistent across replicates.

# Our focus is not on prediction or reconstruction but on loadings.  
# Below is a function, revised from robustComponent(), to determine robust components by  

# * matching components across replicates based on loadings rather than reconstruction
# * adjusting signs of loadings so that signs of loadings of the same components are as 
#   consistent as possible across replicates.

robustLoadings <- function(models, corThr=0.95, matchThr=0.9){
  n.rep <- length(models)
  corList <- list()
  rob <- list(Krobust=0, indices=matrix(NA, n.rep, 0), cor=matrix(NA, n.rep, 0))
  
  maxK <- 0
  ps.mean.W.rep <- list()
  compStorage <- vector("list", length=n.rep) #Store the components that can still be used
  for (r in 1:n.rep){
    res <- models[[r]]
    maxK <- max(maxK, res$K)
    ps.mean.W.rep[[r]] <- list()
    tmp <- apply(res$posterior$W, 2:3, mean)
    for (k in 1:res$K){
      ps.mean.W.rep[[r]][[k]] <- tmp[,k]
    }
    compStorage[[r]] <- which(sapply(ps.mean.W.rep[[r]], function(x) {sum(abs(x)) > 0}))
  }
  
  for (r1 in 1:n.rep){
    # matching <- vector(length = maxK)
    sim <- array(NA, dim = c(maxK, n.rep)) # similarity
    matchInd <- w.sign <- matrix(0, nrow = maxK, ncol = n.rep)
    
    for (k1 in compStorage[[r1]]){
      for (r2 in which(sapply(compStorage, length)>0)){
        d <- sapply(ps.mean.W.rep[[r2]][compStorage[[r2]]], 
                    function(x) cor(x, ps.mean.W.rep[[r1]][[k1]]))
        w.sign[k1,r2] <- sign(d)[which.max(abs(d))]
        sim[k1, r2] <- round(max(abs(d))*w.sign[k1,r2], 3)
        matchInd[k1, r2] <- compStorage[[r2]][which.max(abs(d))]*w.sign[k1,r2]
      }
      
      if (mean(abs(sim[k1,])>corThr, na.rm=T) >= matchThr){ # Robust component found!
        avg.ps.mean.W <- rep(0, models[[1]]$D)
        for (r2 in which(matchInd[k1,] != 0)){
          tmp <- sapply(ps.mean.W.rep[[r2]][[abs(matchInd[k1, r2])]], '*', w.sign[k1, r2])
          avg.ps.mean.W <- avg.ps.mean.W + tmp
          compStorage[[r2]] <- compStorage[[r2]][!compStorage[[r2]]==matchInd[k1,r2]]
        }
        avg.ps.mean.W <- avg.ps.mean.W/sum(matchInd[k1,]!=0)
        
        rob$Krobust <- rob$Krobust + 1
        rob$effect <- array(c(rob$effect, avg.ps.mean.W),
                            dim=c(length(avg.ps.mean.W), rob$Krobust),
                            dimnames=list(colnames(models[[1]]$posterior$W), 
                                          paste0("K",1:rob$Krobust)))
        rob$indices <- cbind(rob$indices, matchInd[k1,])
        rob$cor <- cbind(rob$cor, sim[k1,])
      }
    }
  }
  rob$effect <- round(rob$effect, 3)
  rob$cor <- round(rob$cor, 3)
  rownames(rob$indices) <- rownames(rob$cor) <- paste0("rep",1:n.rep)
  colnames(rob$indices) <- colnames(rob$cor) <- paste0("K",1:rob$Krobust)
  return(rob)
}

### 2. rob.var.exp() ###
# the following function 
# models: a list of replicates of GFA
# rob:    an object obtained from robustLoadings() function
rob.var.exp <- function(models, rob){
  n.rep <- length(models)
  K.rob <- rob$Krobust
  ps.mean.W.rep <- array(NA, dim=c(models[[1]]$D, K.rob, n.rep))
  for (r in 1:n.rep){
    indices <- rob$indices[r,]
    tmp.w <- apply(models[[r]]$posterior$W[,, abs(indices)], 2:3, mean)
    ps.mean.W.rep[,,r] <- sweep(tmp.w, MARGIN=2, sign(indices), '*')
  }
  ve.rep <- apply(ps.mean.W.rep^2, 2:3, mean)*100
  order.r1 <- order(-ve.rep[,1])
  for (r in 1:n.rep){
    ve.rep[,r] <- ve.rep[order.r1,r]
    rob$indices[r,] <- rob$indices[r,][order.r1]
    rob$cor[r,] <- rob$cor[r,][order.r1]
    # rob$indices[r,] <- rob$indices[r,][order(-ve.rep[,r])]
    # rob$cor[r,] <- rob$cor[r,][order(-ve.rep[,r])]
  }
  # ve.rep <- ve.rep[order(-rowMeans(ve.rep)),]
  colnames(ve.rep) <- paste0('rep.', 1:n.rep)
  rownames(ve.rep) <- paste0('Component ', 1:K.rob)
  ve <- data.frame(matrix(NA, K.rob, 3))
  names(ve) <- c('Component', 'Mean', 'SE')
  ve$Component <- 1:K.rob
  ve$Mean <- apply(ve.rep, 1, mean)
  ve$SE <- apply(ve.rep, 1, sd)/sqrt(n.rep)
  # ve$Comp.r <- rob$Krobust+1-order(ve$Mean)
  
  print(round(ve.rep, 3))
  
  # ve <- ve[order(ve$Comp.r),]
  p <- ggplot(ve, 
              aes(x=Component, y=Mean, ymin=Mean-SE, ymax=Mean+SE)) + 
        geom_pointrange() + 
        xlab('Robust components') + ylab('Percent variance explained') + 
        theme_bw() # use a white background 
  print(p + ggtitle(paste0('The ', length(indices), ' robust components explain ',
                           round(sum(ve$Mean),1), '+/-', round(sd(colSums(ve.rep)), 1),
                           '% variance of original variables')))
  print(ve)
  return(list(rob=rob,ve.rep=ve.rep))
}

### 3. a function to create a heat map of loadings
w.plot <- function(w, D, K, gr1, conf.level, replicate){
  mar <- c(6,4,4,6)
  par(mar=mar)
  cols <- colorRampPalette(c("orange","red","white","blue","cyan"))(19)
  if(any(is.na(w))) cols <- colorRampPalette(c("orange","red","#DDDDDD","blue","cyan"))(19)
  M <- max(abs(w),na.rm=T)
  breaks <- seq(-M,M,length=20)
  
  title <- c("Matrix W^T","Components","Features")
  if (!is.null(replicate)){
    title[1] <- paste0('Replicate ', replicate,': ', title[1])
  } else {
    if (is.null(conf.level)){
      title[1] <- paste0(title[1], ' (all components & ', round(sum(w!=0)), ' loadings)')
    } else if (!is.null(conf.level)){
      title[1] <- paste0(title[1], ' (', sum(w!=0),' non-zero loadings at ', conf.level*100, '% confidence)')
    }
  }
  image(1:D, 1:K, w[,K:1], col=cols, breaks=breaks, axes=F, main=title[1],
        xlab="",ylab="")
  title(xlab=title[3],line=mar[1]-1)
  title(ylab=title[2],line=mar[2]-1)
  box()
  par(las=2)
  axis(1, 1:D, rownames(w), cex.axis=D^(-1/5))
  axis(2, K:1, colnames(w), cex.axis=K^(-1/5))
  
  #Grouping
  par(xpd=T)
  mu <- gr1[-1]/2+gr1[-length(gr1)]/2
  N <- K
  for(i in 1:length(mu)) {
    if(i!=length(mu)) lines(rep(gr1[i+1]+1/2,2), c(.5, N*1.03+.5), lwd=2)
    text(mu[i],N*1.03+.5,names(gr1)[i])
  }
  #Colorbar
  n <- length(cols)
  cba <- D + 1/2 + D/60; cbw <- D/40
  for(i in 1:n){
    polygon(c(0,cbw,cbw,0)+cba, c(0,0,N/n,N/n)+N*(i-1)/n+1/2,
            col=cols[i], border=NA)
  }
  #Colorbar: axis
  lines(rep(cba+cbw,2),c(0,N)+1/2)
  m <- 10^floor(log10(M)); m <- floor(M/m)*m
  for(l in c(-m,0,m)) {
    ly <- N*(l/M/2+.5)+1/2
    lines(cba+cbw-c(cbw,-cbw)/5, rep(ly,2))
    text(cba+cbw*2.5,ly,l)
  }
  par(xpd=F)
}

### 4. a function taht excutes previous functions and makes a forest plot for 
###    each robust component
w_rob <- function(models, corThr=0.9, matchThr=0.9, var.exp=TRUE, conf.level=0.95, 
                  heatmap.rep=TRUE, sparse=FALSE, forest=TRUE, block.names){
  variables <- rownames(models[[1]]$W)
  D <- length(variables)
  n.rep <- length(models)
  rob <- robustLoadings(models, corThr, matchThr)
  
  # Percent variance explained across replicates
  if (var.exp) ve <- rob.var.exp(models, rob)
  
  # create a dataframe to store credible intervals of loadings across replicates
  w.ci <- data.frame(matrix(NA, n.rep*rob$Krobust*D, 6))
  names(w.ci) <- c('Replicate', 'Component', 'variable', 'Lower', 'Median', 'Upper')
  w.ci$Replicate <- rep(1:n.rep, each=length(variables)*rob$Krobust)
  w.ci$Component <- rep(rep(1:rob$Krobust, each=length(variables)), n.rep)
  w.ci$variable <- as.factor(rep(variables, n.rep*rob$Krobust))
  
  # compute credible intervals for loadings in each replicate
  for (r in 1:n.rep){
    w.ps <- round(models[[r]]$posterior$W[,,abs(ve$rob$indices[r,])], 3)
    w.ps <- sweep(w.ps, MARGIN=3, sign(ve$rob$indices[r,]), '*')
    p.lo <- (1-conf.level)/2
    p.up <- 1 - p.lo
    tmp <- round(apply(w.ps, 2:3, function(x) quantile(x, c(p.lo, 0.5, p.up))), 3)
    for (k in 1:rob$Krobust){
      w.ci[w.ci$Replicate==r & w.ci$Component==k, 4:6] <- t(tmp[,,k])
    }
  }
  w.ci$contain.0 <- (w.ci$Lower * w.ci$Upper < 0 | w.ci$Lower==0 & w.ci$Upper==0)*1
  
  # compute medians of credible intervals across replicates
  w.ci.med <- data.frame(matrix(NA, D*rob$Krobust, 5))
  names(w.ci.med) <- c('Component', 'variable', 'Lower', 'Median', 'Upper')
  w.ci.med$Component <- rep(1:rob$Krobust, each=D)
  w.ci.med$variable <- rep(variables, rob$Krobust)
  tmp <- aggregate(Lower ~ Component + variable, median, data=w.ci)
  tmp <- merge(tmp, aggregate(Median ~ Component + variable, median, data=w.ci), by=c('Component','variable'))
  tmp <- merge(tmp, aggregate(Upper ~ Component + variable, median, data=w.ci), by=c('Component','variable'))
  for (k in 1:rob$Krobust){
    for (d in 1:D){
      w.ci.med[w.ci.med$Component==k & w.ci.med$variable==variables[d], 3:5] <-
        tmp[tmp$Component==k & tmp$variable==variables[d], 3:5]
    }
  }
  w.ci.med$contain.0 <- (w.ci.med$Lowe * w.ci.med$Upper<0 | 
                           w.ci.med$Lower==0 & w.ci.med$Upper==0)
  w.ci.med$all.0 <- (w.ci.med$Lower==0 & w.ci.med$Upper==0)*1
  
  # heat maps
  gr <- models[[1]]$groups; M <- length(gr)
  if (is.null(block.names)) { 
    names(gr) <- paste("Source",1:M) 
  } else { names(gr) <- block.names }
  gr1 <- c(0,cumsum(sapply(gr,length))); names(gr1) <- c(names(gr),"NA")
  if (heatmap.rep){
    for (r in 1:n.rep){
      w.tmp <- matrix(w.ci$Median[w.ci$Replicate==r]*(1-w.ci$contain.0[w.ci$Replicate==r]), 
                      D, rob$Krobust)
      colnames(w.tmp) <- 1:rob$Krobust
      rownames(w.tmp) <- variables
      w.plot(w.tmp, D=nrow(w.tmp), K=ncol(w.tmp), gr1, conf.level, r) 
    }
  }
  w.med <- matrix(w.ci.med$Median*(1-w.ci.med$contain.0), D, rob$Krobust)
  colnames(w.med) <- 1:rob$Krobust
  rownames(w.med) <- variables
  # print('Robust heat map')
  w.plot(w.med, D=nrow(w.med), K=ncol(w.med), gr1, conf.level, replicate=NULL) 
  
  if (forest){
    for (k in 1:rob$Krobust){
      # Exclude variables whose credible intervals of loadings contain zero
      w.ci.med.k <- w.ci.med[w.ci.med$Component==k,]
      w.ci.med.k$variable <- factor(w.ci.med.k$variable, as.character(variables))
      if (sparse) { w.ci.med.k <- w.ci.med.k[!w.ci.med.k$contain.0,] }
      rob.w <- ggplot(data=w.ci.med.k, 
                      aes(x=variable, y=Median, ymin=Lower, ymax=Upper)) +
        geom_pointrange() + 
        geom_hline(yintercept=0, lty=2) + # add a dotted line at 0 as reference
        coord_flip() + # flip coordinates for horizontal lines
        xlab('Variable') + ylab('Median and 95% Credible Interval') + 
        theme_bw() # use a white background
      print(rob.w + ggtitle(paste0('Robust component ', k)))
    }
  }
  
  # compute median of posterior mean of factor scores across replicates
  x.rep <- array(NA, dim=c(dim(models[[1]]$X)[1], ve$rob$Krobust, n.rep))
  for (r in 1:n.rep){
    x.rep[,,r] <- round(models[[r]]$X[, abs(ve$rob$indices[r,])], 3)
    x.rep[,,r] <- sweep(x.rep[,,r], MARGIN=2, sign(ve$rob$indices[r,]), '*')
  }
  x.rob <- apply(x.rep, 1:2, median)
  return(list(rob=ve$rob, w.ci.med=w.ci.med, x.rob=x.rob))
}


### 0. save GFA results in fewer digits to save space
#gfa_digits <- function(gfa.res, digits){
#  gfa.res$W <- round(gfa.res$W, digits)
#  gfa.res$X <- round(gfa.res$X, digits)
#  gfa.res$Z <- round(gfa.res$Z, digits)
#  gfa.res$posterior$W <- round(gfa.res$posterior$W, digits)
#  gfa.res$posterior$X <- round(gfa.res$posterior$X, digits)
#  return(gfa.res)
#}