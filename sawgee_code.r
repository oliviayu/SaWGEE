# library(expm)

getLeftTerm = function(A, Xmat, cc.cum, k, count_bycluster){
  cluster_mu <- A[(cc.cum[k]+1):cc.cum[k+1]]
  Design = Xmat[(cc.cum[k]+1):cc.cum[k+1], ]
  
  # Compute the row sums
  row_sum_result = (Design%*%t(Design)%*%cluster_mu)*cluster_mu
  diag_cluster.M <- apply(Design**2, 1, sum)*cluster_mu**2
  rm(Design); rm(cluster_mu)
  cluster.M_sum <- sum(row_sum_result)
  
  ### 
  First <- diag_cluster.M 
  Second  <-  -2*(row_sum_result - First)/(count_bycluster[k]-1)
  Third <-  (cluster.M_sum-2*row_sum_result+diag_cluster.M -
               sum(diag_cluster.M) + diag_cluster.M)/(
                 (count_bycluster[k]-1)*(count_bycluster[k]-2))
  
  if(any(First+Second+Third < 0)){
    stop("Negative results in sampling probabilities...")
  }
  # Left.term <- c(Left.term, sqrt(First+Second+Third))
  # rm(First); rm(Second); rm(Third)
  sqrt(First+Second+Third)
}


getBlockDiag <- function(len, xvec=NULL){
  K <- length(len)
  
  if(is.null(xvec)){
    xvec <- rep.int(1, sum(len^2))
  }
  
  row.vec <- col.vec <- vector("numeric", sum(len^2))
  add.vec <- cumsum(len) - len
  if(K == 1){
    index <- c(0, sum(len^2))
  }else{
    index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  }
  
  for(i in 1:K){
    row.vec[(index[i] + 1):(index[i+1])] <- rep.int( (1:len[i]) + add.vec[i], len[i])
    col.vec[(index[i] + 1):(index[i+1])] <- rep( (1:len[i]) + add.vec[i], each=len[i])
  }	
  BlockDiag <- sparseMatrix(i = row.vec, j = col.vec, x = xvec)
  
  if(!is.null(xvec)){
    testsymm <- abs(sum(skewpart(BlockDiag)))
    if(testsymm != 0) {
      warning("Correlation matrix is not computed to be exactly symmetric. Taking only the symmetric part.")
    }
  }
  return(list(BDiag = symmpart(BlockDiag), row.vec =row.vec, col.vec=col.vec))
}


## Inverse of (1-\alpha)I+\alpha*1*1^T

getAlphaInvEX_modified <- function(alpha.new, diag.vec, BlockDiag){
  return(as(BlockDiag %*% Diagonal(x = (-alpha.new/((1-alpha.new)*(1+(diag.vec-1)*alpha.new))))
            + Diagonal( x = ((1+(diag.vec-2)*alpha.new)/((1-alpha.new)*(1+(diag.vec-1)*alpha.new))
                             + alpha.new/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)))), "symmetricMatrix"))
}


updateAlphaEX_modified <- function(YY, mu, VarFun, phi, id, len, StdErr,
                                   Resid, p, BlockDiag, included,
                                   includedlen, sqrtW, useP,sampling_W){
  W <- sqrtW^2
  # Resid <- StdErr %*% included %*% sqrtW %*% inv_sqrtmZ %*% Diagonal(x = YY - mu)
  
  denom <- phi*(sum(includedlen*(includedlen-1))/2 - useP * p)
  
  ## Modified (add sampling_W: (z_{ij}z_{ij'})/(\pi_{ij}\pi_{ij'}*))
  
  BlockDiag <- sampling_W*(StdErr  %*%Diagonal(x = YY - mu) %*%  included %*% (BlockDiag) %*% W %*% Diagonal(x = YY - mu)  %*% StdErr)
  
  
  ## Original method
  alpha <- sum(triu(BlockDiag, k=1))
  # numpos <- sum(triu(BlockDiag, k=1) != 0)
  
  alpha.new <- alpha/denom
  
  return(alpha.new)
}


getfam <- function(family){
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))  
  }
  
  if(is.function(family)){
    family <- family()
    return(family)
  }else if(inherits(family, "family")){
    return(family)
  }else if(is.list(family)){
    if(length(match(names(family), c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv"))) == 4){
      famname <- "custom"
      LinkFun <- family$LinkFun
      InvLink <- family$InvLink
      VarFun <- family$VarFun
      InvLinkDeriv <- family$InvLinkDeriv
    }else{
      famname <- "custom"
      LinkFun <- family[[1]]
      VarFun <- family[[2]]
      InvLink <- family[[3]]
      InvLinkDeriv <- family[[4]]
    }
    
    FunList <- list("family"= famname, "LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv) 
    return(FunList)
  }else{
    stop("problem with family argument: should be string, family object, or list of functions")
  }
}


updateBeta = function(YY, XX, beta, off, InvLinkDeriv, InvLink,
                      VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, W, included){
  beta.new <- beta
  conv=F
  for(i in 1:10){
    eta <- as.vector(XX%*%beta.new) + off
    
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    mu <- InvLink(eta)
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    hess <- crossprod( StdErr %*% dInvLinkdEta %*%XX, R.alpha.inv %*% W %*% StdErr %*%dInvLinkdEta %*% XX)
    esteq <- crossprod( StdErr %*%dInvLinkdEta %*%XX , R.alpha.inv %*% W %*% StdErr %*% as.matrix(YY - mu))
    
    #hess <- crossprod( StdErr %*% dInvLinkdEta %*%XX, included %*% R.alpha.inv  %*% W %*% StdErr %*%dInvLinkdEta %*% XX)
    #esteq <- crossprod( StdErr %*%dInvLinkdEta %*%XX , included %*% R.alpha.inv %*% W %*% StdErr %*% as.matrix(YY - mu))
    
    
    update <- solve(hess, esteq)
    
    
    beta.new <- beta.new + as.vector(update)
  }
  return(list(beta = beta.new, hess = hess))
}


updatePhi_modified <- function(YY, mu, VarFun, p, StdErr, included, includedlen, 
                               sqrtW, useP,sampling_W,original_includedlen){
  nn <- sum(includedlen)
  
  ## Original
  resid <- diag(StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu))
  
  ## Modified
  # resid <- diag(StdErr %*% included %*% sqrtW %*% inv_Z %*% Diagonal(x = YY - mu))
  
  # phi <- (1/(sum(included)- useP * p))*crossprod(resid,resid)
  
  phi <- (1/(sum(original_includedlen)- useP * p))*crossprod(resid, diag(sampling_W)*resid)
  
  ## Modified (Need Tingting's comment), check the denominator
  
  #phi <- (1/(sum(original_includedlen)- useP * p))* resid %*% inv_Z %*%matrix(resid,ncol=1)
  
  return(as.numeric(phi))
}


getSandwich_modified = function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv,
                                InvLink, VarFun, hessMat, StdErr, dInvLinkdEta,
                                BlockDiag, W, included){
  
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
  scoreDiag <- Diagonal(x= YY - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  
  numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% XX,  R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*%  StdErr %*% dInvLinkdEta %*% XX))
  #numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% XX, included %*% R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*% included %*% StdErr %*% dInvLinkdEta %*% XX))
  
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  
  return(list(sandvar = sandvar, numsand = numsand))
}

getSandwich = function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv,
                       InvLink, VarFun, hessMat, StdErr, dInvLinkdEta,
                       BlockDiag, W, included){
  
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
  scoreDiag <- Diagonal(x= YY - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  
  numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% XX,  R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*%  StdErr %*% dInvLinkdEta %*% XX))
  #numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% XX, included %*% R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*% included %*% StdErr %*% dInvLinkdEta %*% XX))
  
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  
  return(list(sandvar = sandvar, numsand = numsand))
}


predict.geem <- function(object, newdata = NULL,...){
  coefs <- object$beta
  if(is.null(newdata)){
    return(as.vector(object$X %*% object$beta))
  }else{
    if(dim(newdata)[2] != length(coefs)){warning("New observations must have the same number of rows as coefficients in the model")}
    return(as.vector(newdata %*% object$beta))
  }
}

# library(expm)


geem.modified <- function(formula, original_id, id, waves=NULL, data = parent.frame(),
                          family = gaussian, corstr = "independence", Mv = 1,
                          weights = NULL, corr.mat = NULL, 
                          init.beta=NULL,
                          init.alpha=NULL, init.phi = 1, scale.fix = FALSE, nodummy = FALSE,
                          sandwich = FALSE, useP = TRUE, maxit = 20, tol = 0.00001, DF){
  
  # dat=simdatPerm
  # familiy=gaussian
  # id=simdatPerm$idvar
  # formula=yvar~tvar
  call <- match.call()
  
  famret <- getfam(family)
  
  if(inherits(famret, "family")){
    LinkFun <- famret$linkfun
    InvLink <- famret$linkinv
    VarFun <- famret$variance
    InvLinkDeriv <- famret$mu.eta
  }else{
    LinkFun <- famret$LinkFun
    VarFun <- famret$VarFun
    InvLink <- famret$InvLink
    InvLinkDeriv <- famret$InvLinkDeriv
  }
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix=TRUE, then init.phi must be supplied")
  }
  
  useP <- as.numeric(useP)
  
  ### First, get all the relevant elements from the arguments
  dat <- model.frame(formula, data, na.action=na.pass)
  nn <- dim(dat)[1]
  cor.match=3
  if(typeof(data) == "environment"){
    id <- id
    weights <- weights
    if(is.null(call$weights)) weights <- rep(1, nn)
    waves <- waves
  }  else{
    if(length(call$id) == 1){
      subj.col <- which(colnames(data) == call$id)
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      }else{
        id <- eval(call$id, envir=parent.frame())
      }
    }else if(is.null(call$id)){
      id <- 1:nn
    }
    
    if(length(call$weights) == 1){
      weights.col <- which(colnames(data) == call$weights)
      if(length(weights.col) > 0){
        weights <- data[,weights.col]
      }else{
        weights <- eval(call$weights, envir=parent.frame())
      }
    }else if(is.null(call$weights)){
      weights <- rep.int(1,nn)
    }
    
    if(length(call$waves) == 1){
      waves.col <- which(colnames(data) == call$waves)
      if(length(waves.col) > 0){
        waves <- data[,waves.col]
      }else{
        waves <- eval(call$waves, envir=parent.frame())
      }
    }else if(is.null(call$waves)){
      waves <- NULL
    }
  }
  dat$id <- id
  dat$weights <- weights
  dat$waves <- waves
  
  
  if(!is.numeric(dat$waves) & !is.null(dat$waves)) stop("waves must be either an integer vector or NULL")
  
  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  # Note that we need to assign weight 0 to rows with NAs
  # in order to preserve the correlation structure
  na.inds <- NULL
  
  if(any(is.na(dat))){
    na.inds <- which(is.na(dat), arr.ind=T)
  }
  
  #SORT THE DATA ACCORDING TO WAVES
  if(!is.null(waves)){
    dat <- dat[order(id, waves),]
  }else{
    dat <- dat[order(id),]
  }
  
  
  # Figure out the correlation structure
  if(!is.null(dat$waves)){
    wavespl <- split(dat$waves, dat$id)
    idspl <- split(dat$id, dat$id)
    
    maxwave <- rep(0, length(wavespl))
    incomp <- rep(0, length(wavespl))
    
    for(i in 1:length(wavespl)){
      maxwave[i] <- max(wavespl[[i]]) - min(wavespl[[i]]) + 1
      if(maxwave[i] != length(wavespl[[i]])){
        incomp[i] <- 1
      }
    }
    
    #If there are gaps and correlation isn't independent or exchangeable
    #then we'll add some dummy rows
    if( !is.element(cor.match,c(1,3)) & (sum(incomp) > 0) & !nodummy){
      dat <- dummyrows(formula, dat, incomp, maxwave, wavespl, idspl)
      id <- dat$id
      waves <- dat$waves
      weights <- dat$weights
    }
  }
  
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0
    for(i in unique(na.inds)[,2]){
      if(is.factor(dat[,i])){
        dat[na.inds[,1], i] <- levels(dat[,i])[1]
      }else{
        dat[na.inds[,1], i] <- median(dat[,i], na.rm=T)
      }
    }
  }
  includedvec <- weights>0
  
  inclsplit <- split(includedvec, id)
  
  dropid <- NULL
  allobs <- T
  if(any(!includedvec)){
    allobs <- F
    for(i in 1:length(unique(id))){
      if(all(!inclsplit[[i]])){
        dropid <- c(dropid, unique(id)[i])
      }
    }
  }
  
  dropind <- c()
  
  if(is.element(cor.match, c(1,3))){
    dropind <- which(weights==0)
  }else if(length(dropid)>0){
    dropind <- which(is.element(id, dropid))
  }
  if(length(dropind) > 0){
    dat <- dat[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- id[-dropind]
  }
  nn <- dim(dat)[1]
  K <- length(unique(id))
  # cat("K=", K, "\n")
  
  modterms <- terms(formula)
  
  X <- model.matrix(formula,dat)
  Y <- model.response(dat)
  offset <- model.offset(dat)
  
  p <- dim(X)[2]
  
  
  
  ### if no offset is given, then set to zero
  if(is.null(offset)){
    off <- rep(0, nn)
  }else{
    off <- offset
  }
  
  # Is there an intercept column?
  interceptcol <- apply(X==1, 2, all)
  
  ## Basic check to see if link and variance functions make any kind of sense
  linkOfMean <- LinkFun(mean(Y[includedvec])) - mean(off)
  
  if( any(is.infinite(linkOfMean) | is.nan(linkOfMean)) ){
    stop("Infinite or NaN in the link of the mean of responses.  Make sure link function makes sense for these data.")
  }
  if( any(is.infinite( VarFun(mean(Y))) | is.nan( VarFun(mean(Y)))) ){
    stop("Infinite or NaN in the variance of the mean of responses.  Make sure variance function makes sense for these data.")
  }
  
  if(is.null(init.beta)){
    if(any(interceptcol)){
      #if there is an intercept and no initial beta, then use link of mean of response
      init.beta <- rep(0, dim(X)[2])
      init.beta[which(interceptcol)] <- linkOfMean
    }else{
      stop("Must supply an initial beta if not using an intercept.")
    }
  }
  
  
  # Number of included observations for each cluster
  
  #### Added 
  # original_K <- length(unique(original_id))
  # 
  # original_includedlen <- rep(0, original_K)
  # original_len <- rep(0,original_K)
  # original_uniqueid <- unique(original_id)
  # original_includedvec <- rep(TRUE,length(original_id))
  # 
  # original_tmpwgt <- as.numeric(original_includedvec)
  # original_idspl <-ifelse(original_tmpwgt==0, NA, original_id)
  # oroignal_includedlen <- as.numeric(summary(split(original_Y, original_idspl, drop=T))[,1])
  original_includedlen <-table(original_id)
  ###############
  
  includedlen <- rep(0, K)
  len <- rep(0,K)
  uniqueid <- unique(id)
  
  tmpwgt <- as.numeric(includedvec)
  idspl <-ifelse(tmpwgt==0, NA, id)
  includedlen <- as.numeric(summary(split(Y, idspl, drop=T))[,1])
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x=(as.numeric(weights>0)))
  
  # Get vector of cluster sizes... remember this len variable
  #len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  alpha.new <- 0.2
  
  #if no initial overdispersion parameter, start at 1
  if(is.null(init.phi)){
    phi <- 1
  }else{
    phi <- init.phi
  }
  
  beta <- init.beta
  
  #Set up matrix storage
  StdErr <- Diagonal(nn)
  dInvLinkdEta <- Diagonal(nn)
  Resid <- Diagonal(nn)
  #  if( (max(len)==1) & cor.match != 1 ){
  #    warning("Largest cluster size is 1. Changing working correlation to independence.")
  #    cor.match <- 1
  #    corstr <- "independence"
  #  }
  
  # Initialize for each correlation structure
  
  tmp <- getBlockDiag(len)
  BlockDiag <- tmp$BDiag
  rm(tmp)
  
  #Create a vector of length number of observations with associated cluster size for each observation
  n.vec <- vector("numeric", nn)
  index <- c(cumsum(len) - len, nn)
  for(i in 1:K){
    n.vec[(index[i]+1) : index[i+1]] <-  rep(includedlen[i], len[i])
  }
  
  stop <- F
  converged <- F
  count <- 0
  beta.old <- beta
  unstable <- F
  phi.old <- phi
  
  
  # Main fisher scoring loop
  while(!stop){
    count <- count+1
    
    eta <- as.vector(X %*% beta) + off
    
    mu <- InvLink(eta)
    
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    if(!scale.fix){
      phi <- updatePhi_modified(Y, mu, VarFun, p, StdErr, included, includedlen, sqrtW, 
                                useP, corr.mat, original_includedlen)
    }
    if (is.nan(phi) | phi>1e+4){
      phi<-init.phi
    }
    phi.new <- phi
    # print(phi)
    
    #EXCHANGEABLE
    alpha.new <- updateAlphaEX_modified(Y, mu, VarFun, phi, id, len, StdErr,
                                        Resid, p, BlockDiag, included,
                                        original_includedlen, sqrtW, useP, corr.mat)
    
    #R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
    
    ## modified version 
    # R.alpha.inv2 <- (getAlphaInvEX_modified(alpha.new, n.vec, BlockDiag)/phi)* corr.mat
    n.vec2 = rep(original_includedlen, len)
    R.alpha.inv <- (getAlphaInvEX_modified(alpha.new, n.vec2, BlockDiag)/phi)* corr.mat
    
    beta.list <- updateBeta(Y, X, beta, off, InvLinkDeriv, InvLink, VarFun, 
                            R.alpha.inv, StdErr, dInvLinkdEta, tol, W, included)
    beta <- beta.list$beta
    
    phi.old <- phi
    if( max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol ){converged <- T; stop <- T}
    if(count >= maxit){stop <- T}
    beta.old <- beta
  }
  
  # biggest <- which.max(len)[1]
  # index <- sum(len[1:biggest])-len[biggest]
  # 
  # if(K == 1){
  #   biggest.R.alpha.inv <- R.alpha.inv
  #   if(cor.match == 6) {
  #     biggest.R.alpha <- corr.mat*phi
  #   }else{
  #     biggest.R.alpha <- solve(R.alpha.inv)
  #   }
  # }else{
  #   biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
  #   if(cor.match == 6){
  #     biggest.R.alpha <- corr.mat[1:len[biggest] , 1:len[biggest]]*phi
  #   }else{
  #     biggest.R.alpha <- solve(biggest.R.alpha.inv)
  #   }
  # }
  
  eta <- as.vector(X %*% beta) + off
  mu <- InvLink(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
  
  if(sandwich){
    sandvar.list <- getSandwich(Y, X, eta, id, 
                                R.alpha.inv,
                                phi, InvLinkDeriv, InvLink, VarFun, 
                                beta.list$hess, StdErr, dInvLinkdEta, BlockDiag, W, included)
    
  } else {
    sandvar.list$sandvar = NULL
  }
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
  scoreDiag <- Diagonal(x= Y - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  Mmat <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% X,  
                                R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*%  StdErr %*% dInvLinkdEta %*% X))
  
  if(!converged){warning("Did not converge")}
  if(unstable){warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")}
  
  
  # Create object of class geem with information about the fit
  dat <- model.frame(formula, data, na.action=na.pass)
  X <- model.matrix(formula, dat)
  
  if(is.character(alpha.new)){alpha.new <- 0}
  results <- list()
  results$beta <- as.vector(beta)
  results$phi <- phi
  results$alpha <- alpha.new
  if(cor.match == 6){
    results$alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat,1)!=0)])
  }
  results$coefnames <- colnames(X)
  results$niter <- count
  results$converged <- converged
  results$naiv.var <- solve(beta.list$hess)  ## call model-based
  results$var <- sandvar.list$sandvar
  results$Hmat = beta.list$hess
  results$Mmat = Mmat
  results$call <- call
  results$corr <- "Exchange"
  results$clusz <- len
  results$FunList <- famret
  results$X <- X
  results$offset <- off
  results$eta <- eta
  results$dropped <- dropid
  results$weights <- weights
  results$terms <- modterms
  results$y <- Y
  # results$biggest.R.alpha <- biggest.R.alpha/phi
  results$formula <- formula
  results$original=original_includedlen
  class(results) <- "geem"
  return(results)
}



#' Adjust sampling probabilities such that max(p)<=1 and sum(p)=M.
#' 
#' @param Prob.N crude Prob.N calculated from proposed formula
#' @param M expected number of samples been selected
#' @param adjustP adjust probabilities 
#' @param group a vector specifying the clusters
#' @param maxC maximum number of samples allowed within a cluster 
adjust_p = function(Prob.N, M, adjustP=T, group, maxC=NULL){
  
  if(adjustP==T){
    Prob.OS.sampling = Prob.N*M
    
    if(!is.null(maxC)){
      df = data.frame(p=Prob.OS.sampling, group=group) 
      ez = group_by(df, group) %>% summarise(expect = sum(p)) %>%
        mutate(exceed = as.numeric(round(expect) > maxC))
      # print(summary(ez))
      while(any(ez$exceed > 0)){
        df2=left_join(df, ez, by="group")
        df2$newp = df2$p
        df2$newp[df2$exceed==1] = (df2$p/df2$expect*maxC)[df2$exceed==1]
        left_exp = M - maxC*sum(ez$exceed)
        org_exp = sum(df2$p[df2$exceed==0])
        df2$newp[df2$exceed==0] = (df2$p/org_exp*left_exp)[df2$exceed==0]
        
        df = data.frame(p=df2$newp, group=group) 
        ez = group_by(df, group) %>% summarise(expect = sum(p)) %>%
          mutate(exceed = as.numeric(round(expect) > maxC))
        # print(summary(ez))
      }
      Prob.OS.sampling = df$p
    }
    
    while(max(Prob.OS.sampling) > 1){
      onep = which(Prob.OS.sampling >= 1)
      normp = which(Prob.OS.sampling < 1)
      Prob.OS.sampling[onep] = 1
      Prob.OS.sampling[normp] =  (M-length(onep))*
        Prob.OS.sampling[normp]/sum(Prob.OS.sampling[normp])
    }
    
  } else {
    Prob.OS.sampling <- pmin(Prob.N*M, 1)
  }
  Prob.OS.sampling
}

#' Calculate sampling probabilities given a strategy
#' 
get_sampling_prob = function(fm, DF, mu, M, strategy, bycluster,
                             adjustP=T, cc_m=1, track=0, maxC, method=1, rho_ini=0){
  
  count_bycluster <- table(DF$idvar)
  
  if(strategy=='os'){
    n_id = length(count_bycluster)
    cc.cum <- c(0, cumsum(count_bycluster))
    
    A <- sqrt(((mu)*(1-mu)))
    Xmat = model.matrix(fm, data=DF)
    
    if(method == 2){
      # calculate Hmat= sum(DVD')
      Lambda_mat <- lapply(1:n_id, function(j){
        kk = unique(DF$idvar)[j]
        D_mat = crossprod(Xmat[DF$idvar==kk, ], diag(A[DF$idvar==kk]))
        mi=count_bycluster[j]
        if(rho_ini != 0){
          R_inv = matrix(-rho_ini/((1-rho_ini)^2+mi*rho_ini*(1-rho_ini)), mi,mi)
          diag(R_inv) = 1/(1-rho_ini) - rho_ini/((1-rho_ini)^2+mi*rho_ini*(1-rho_ini))
        } else {
          R_inv = diag(1, mi,mi)
        }
        D_mat%*%R_inv%*%t(D_mat)/n_id
      }) %>% Reduce("+",.)
      Xmat = Xmat%*%solve(Lambda_mat)
    }
    
    ymat = model.response(model.frame(fm, data=DF))
    Res.Err <- abs(ymat-mu)/A
    rm(mu)
    
    Left.term <- c()
    for(k in 1:n_id){
      Left.term = c(Left.term, getLeftTerm(A, Xmat, cc.cum, k, count_bycluster))
      # cat("k=", k, "...\n")
    }
    rm(A)
    
    Prob <- Res.Err*Left.term
    Prob.N <- Prob/sum(Prob)
    rm(Prob); rm(Left.term); rm(Res.Err)
    
    if(bycluster==F){
      Prob.OS.sampling = adjust_p(Prob.N, M, adjustP, group=DF$idvar, maxC)
    } else {
      group_prob = data.frame(prob=Prob.N, idvar = DF$idvar,
                              Mk=rep(M, count_bycluster)) %>%
        group_by(idvar) %>%
        mutate(groupsum=sum(prob), stdP = prob/groupsum) %>%
        mutate(p = adjust_p(stdP, Mk, adjustP))
      Prob.OS.sampling = group_prob$p
      rm(group_prob)
    }
    rm(Prob.N)
    
  } else if (strategy=='random'){
    
    if(bycluster==F){
      Prob.OS.sampling = rep(M/nrow(DF), nrow(DF))
    } else {
      Prob.OS.sampling = rep(M/count_bycluster, count_bycluster)
    }
    
  } else if (strategy == 'cc'){
    
    ymat = model.response(model.frame(fm, data=DF))
    Prob.OS.sampling = rep(NA, nrow(DF))
    ncase=sum(ymat)
    Prob.OS.sampling[ymat==1] = pmin(M/(cc_m+1)/ncase, 1)
    Prob.OS.sampling[ymat==0] = (M-sum(Prob.OS.sampling[ymat==1]))/(nrow(DF)-ncase)
    # Prob.OS.sampling[ymat==0] = M/(cc_m+1)*cc_m/(nrow(DF)-ncase)
    
  } else if (strategy=='ccc'){
    
    ymat = model.response(model.frame(fm, data=DF))
    ncase = aggregate(ymat, list(DF$idvar), sum) %>%
      rename(idvar=Group.1, ncase=x) %>%
      mutate(M=M, size = count_bycluster)
    new_DF = left_join(DF, ncase, by="idvar")
    py1 = pmin(new_DF$M/(cc_m+1)/new_DF$ncase, 1)
    # py0 =  new_DF$M/(cc_m+1)*cc_m/(new_DF$size-new_DF$ncase)
    new_DF$py1=py1
    new_DF = filter(new_DF, ymat==1) %>% group_by(idvar) %>% summarise(sumpy1=sum(py1)) %>%
      left_join(new_DF, ., by='idvar')
    py0 =  (new_DF$M-new_DF$sumpy1)/(new_DF$size-new_DF$ncase)
    
    Prob.OS.sampling = py1*ymat + py0*(1-ymat)
    
    while(max(Prob.OS.sampling) > 1){
      onep = which(Prob.OS.sampling >= 1)
      normp = which(Prob.OS.sampling < 1)
      Prob.OS.sampling[onep] = 1
      Prob.OS.sampling[normp] =  (sum(Mccc)-length(onep))*
        Prob.OS.sampling[normp]/sum(Prob.OS.sampling[normp])
    }
    
  }
  
  # cat("range(p)=", range(Prob.OS.sampling), '\n')
  # cat("sum(p)=", sum(Prob.OS.sampling), ", M=", sum(M), '\n')
  if(max(Prob.OS.sampling) > 1 | min(Prob.OS.sampling) < 0) stop("range(Prob.OS.sampling) is outside [0,1]...")
  if(round(sum(Prob.OS.sampling)) < round(sum(M))) stop("sum(Prob.OS.sampling) < sum(M)")
  
  Prob.OS.sampling
}

#' Fit Weighted GEE for sub-sampled data
#' 
# mu <- predict.glm(m.pilot, DF, type="response")
# ini.beta = coef(m.pilot)
# fm = yvar~x1+x2
subWGEE <- function(mu, DF, M, T_n, ini.beta, fm, id="idvar", 
                    bycluster=F, adjustP=T, sandwich=T, 
                    strategy=c("os", "random", "cc", "ccc"), hybrid=F,
                    maxC = NULL,
                    cc_m=1, track = 0, silent=T, parallel=F, mc.cores=1,
                    method, rho_ini) {
  
  cat("Always scale the continuous variables before running model!!!\n")
  if(id!="idvar"){
    DF$idvar = DF[[id]]
  }
  if(!class(DF$idvar) %in% c('numeric', 'integer')){
    stop('Please convert the idvar to numeric/integer.')
  }
  if(strategy=="ccc"){
    bycluster=T
  }
  if(bycluster == T & length(M) == 1){
    stop("Length of M must equal to number of clusters.")
  }
  if(!is.null(maxC)){
    if(sum(M) > maxC*length(unique(DF$idvar))) stop("M is too large!")
  }
  
  if(hybrid==F){
    Prob.OS.sampling = get_sampling_prob(fm, DF, mu, M, strategy, bycluster,
                                         adjustP, cc_m, track, maxC, method, rho_ini)
  } else {
    if(strategy!='os') stop('"hybrid" option is for OS only...')
    ymat = model.response(model.frame(fm, data=DF))
    Prob.OS.sampling = rep(NA, nrow(DF))
    Prob.OS.sampling[ymat==1] <- 1
    control_idx = which(ymat==0)
    Prob.OS.sampling[control_idx] = get_sampling_prob(fm, DF[control_idx, ], 
                                                      mu[control_idx], 
                                                      M - sum(ymat), 
                                                      strategy, bycluster,
                                                      adjustP, cc_m, track, maxC, method, rho_ini)
  }
  print("calculating prob is done...")
  print(summary(Prob.OS.sampling))
  if(parallel==F){  
    
    output = c()
    sim = 1; nfail = 0
    Hmat = Mmat = matrix(0, length(ini.beta), length(ini.beta))
    TT=0
    
    while(sim <= T_n & nfail < T_n*10){
      
      sample_subject <- rbinom(nrow(DF), 1, Prob.OS.sampling)
      sample_idx <- which(sample_subject==1) %>% sort()
      rm(sample_subject)
      # sample_idx <- sample(nrow(DF), sum(M), prob=Prob.OS.sampling, replace = F) %>% sort()
      if(track>0) print(length(sample_idx))
      
      DF.sub <- DF[sample_idx, ]
      Prob.OS <- Prob.OS.sampling[sample_idx]
      sub.count_bycluster <- table(DF.sub$idvar)
      if(track>0){
        cat("cluster size of subdata:", range(as.numeric(sub.count_bycluster)), "\n")
      }
      if(track >1){
        group_by(DF.sub, idvar) %>% summarise(case=sum(yvar), control=n()-sum(yvar),
                                              rate=mean(yvar), n=n()) %>%
          summary() %>% print()
      }
      # 
      sub.cc.cum <-c(0, cumsum(sub.count_bycluster))
      # Z = lapply(1:length(sub.count_bycluster), function(k){
      #   clus.idx <- (sub.cc.cum[k]+1):sub.cc.cum[k+1]
      #   w_k = 1/Prob.OS[clus.idx]
      #   outer(w_k, w_k)
      # }) %>% do.call(bdiag, .) %>% as.matrix()
      # diag(Z) <- 1/Prob.OS
      
      Zvec = lapply(1:length(sub.count_bycluster), function(k){
        clus.idx <- (sub.cc.cum[k]+1):sub.cc.cum[k+1]
        w_k = 1/Prob.OS[clus.idx]
        res= outer(w_k, w_k)
        diag(res) = w_k
        as.numeric(res)
      }) %>% do.call(c, .)
      Z2 = getBlockDiag(sub.count_bycluster, xvec=Zvec)

      rm(Prob.OS)      
      OS.modified <- try(geem.modified(fm, original_id=DF$idvar, id=idvar, 
                                       data = DF.sub, family = binomial(link = "logit"), 
                                       corr.mat = Z2$BDiag,
                                       init.beta = ini.beta, sandwich = sandwich,
                                       DF=DF), 
                         silent=silent)
      # OS.modified$var
      
      if(class(OS.modified) != 'try-error'){
        print(OS.modified$beta)
        
        if(OS.modified$converged){
          res = c(nrow(DF.sub), OS.modified$beta, OS.modified$phi, OS.modified$alpha)
          
          if(sandwich==T){
            res = c(res, as.numeric(as.matrix(OS.modified$var)))
          }
            Hmat = Hmat + OS.modified$Hmat
            Mmat = Mmat + OS.modified$Mmat
            TT = TT+1

          output=rbind(output, res)
          
          cat(paste("Resampling", sim, "is done for", strategy, " strategy ...\n"))
          sim = sim+1
        } else {
          cat(paste("Resampling", sim, "didn't converge for", strategy, " strategy ...\n"))
          nfail = nfail+1
        }
      } else {
        cat(paste("Resampling", sim, "is failed for", strategy, " strategy ...\n"))
        nfail = nfail+1
      }
    }
    
  } else {
    
    output = parallel::mclapply(1:T_n, function(sim){
      
      success = F
      nfail=0
      res = c()
      while(success == F & nfail <= 10){
        
        sample_subject <- rbinom(nrow(DF), 1, Prob.OS.sampling)
        sample_idx <- which(sample_subject==1)
        if(track>0) print(length(sample_idx))
        rm(sample_subject)
        
        DF.sub <- DF[sample_idx, ]
        Prob.OS <- Prob.OS.sampling[sample_idx]
        sub.count_bycluster <- table(DF.sub$idvar)
        if(track>0){
          cat("cluster size of subdata:", range(as.numeric(sub.count_bycluster)), "\n")
        }
        if(track >1){
          group_by(DF.sub, idvar) %>% summarise(case=sum(yvar), control=n()-sum(yvar),
                                                rate=mean(yvar), n=n()) %>%
            summary() %>% print()
        }
        
        sub.cc.cum <-c(0, cumsum(sub.count_bycluster))
        Z = lapply(1:length(sub.count_bycluster), function(k){
          clus.idx <- (sub.cc.cum[k]+1):sub.cc.cum[k+1]
          w_k = 1/Prob.OS[clus.idx]
          outer(w_k, w_k)
        }) %>% do.call(bdiag, .) %>% as.matrix()
        diag(Z) <- 1/Prob.OS
        rm(Prob.OS)
        
        OS.modified <- try(geem.modified(fm, original_id=DF$idvar, id=idvar, 
                                         data = DF.sub, family = binomial, corr.mat = Z,
                                         init.beta = ini.beta, sandwich = sandwich), 
                           silent=silent)
        
        res = c()
        if(class(OS.modified) != 'try-error'){
          print(OS.modified$beta)
          
          if(OS.modified$converged){
            res = c(nrow(DF.sub), OS.modified$beta, OS.modified$phi, OS.modified$alpha)
            
            if(sandwich==T){
              res = c(res, as.numeric(as.matrix(OS.modified$var)))
            }
            success = TRUE
            cat(paste("Resampling", sim, "is done for", strategy, " strategy ...\n"))
          } else {
            cat(paste("Resampling", sim, "didn't converge for", strategy, " strategy ...\n"))
            nfail=nfail+1
          }
        } else {
          cat(paste("Resampling", sim, "is failed for", strategy, " strategy ...\n"))
          nfail=nfail+1
        }
      }
      res
    }, mc.cores = mc.cores) %>% do.call(rbind,.)
  }
  
  Hmatbar = Hmat/TT
  Mmatbar = Mmat/TT
  # Sigma_bar = solve(Hmatbar)%*%Mmatbar%*%solve(Hmatbar)
  sandvar <- t(solve(Hmatbar, Mmatbar))
  Sigma_bar <- t(solve(t(Hmatbar), sandvar))

  list(output=output, Sigma_bar=Sigma_bar, Hmatbar=Hmatbar, Mmatbar=Mmatbar, nfail=nfail)
}

## Summarize simulation results
summary_sim = function(output_glm, true_par){
  p = length(true_par)
  mu = apply(output_glm[,2:(p+1)], 2, mean)
  bias =  abs((mu - true_par)/true_par)
  sd = apply(output_glm[,2:(p+1)], 2, sd)
  msd = apply(output_glm[, (p+2):(2*p+1)], 2, mean)
  
  rmse = sqrt(apply(apply(output_glm[,2:(p+1)], 1, 
                          function(x){(x-true_par)^2}), 1, mean))/abs(true_par)
  cp = sapply(2:(p+1), function(k){
    mean(output_glm[,k] - 1.96*output_glm[,k+p] <= true_par[k-1] &
           output_glm[,k] + 1.96*output_glm[,k+p] >= true_par[k-1])
  })
  list(mu=mu, Rbias=bias, sd=sd, msd=msd, Rrmse=rmse, cp=cp) %>%
    do.call(cbind,.)
}

summary_sim_subWGEE = function(output, os_sd, true_par, rmNA=T, trim=0.025){
  p = length(true_par)
  xx  = as.data.frame(output)
  xx_sd = as.data.frame(os_sd) %>% group_by(b) %>% 
    summarise_all(mean, trim=trim)
  
  xxx_sd = lapply(unique(xx$b), function(b){
    T_n = sum(xx$b==b)
    xxx=scale(xx[xx$b==b, 2:(p+1)], center=T, scale = F)
    ggg=c()
    for(i in 1:(nrow(xxx))){
      ggg=rbind(ggg,  c(diag((xxx[i,])%*%t(xxx[i,]))))
    }
    cov_est = matrix(as.numeric(xx_sd[xx_sd$b==b, -1]), ncol=p, nrow=p) -
      apply(ggg, 2, function(x){mean(x,trim=trim)})
    
    sqrt(diag(cov_est))
  }) %>% do.call(rbind,.)
  
  names(xx) = c("b", paste0("V", 2:ncol(xx)))
  xxx = group_by(xx, b) %>% 
    summarise(mu1=mean(V2, trim=trim), mu2=mean(V3, trim=trim), mu3=mean(V4, trim=trim)) %>%
    cbind(xxx_sd)
  if(rmNA==T){
    xxx=xxx[complete.cases(xxx),]
    cat("NA removed.", nrow(xxx), "repetitions left.\n")
  }
  summary_sim(xxx, true_par)
}

summary_subWGEE = function(output_os, true_par, rmNA=T, trim=0){
  p = length(true_par)
  T_n = nrow(output_os)  
  if(is.null(T_n)) T_n=1
  
  if(T_n>1){
    if(rmNA==T){
      output_os=output_os[complete.cases(output_os),]
      cat("After removing NAs", nrow(output_os), "iterations left.\n")
    }
    beta_hat = apply(output_os[,1:p], 2, mean)
    
    if(trim ==0){
      cov_est = matrix(apply(output_os[, -c(1:(p+2))], 2, mean), 
                       ncol=p, nrow=p) - 
        cov(output_os[, 1:(p)])*(T_n-1)/T_n
      
    } else if (trim>0){
      
      cv1 = matrix(apply(output_os[, -c(1:(p+2))], 2, function(x){mean(x, trim=trim)}), 
                   ncol=p, nrow=p)
      xxx = scale(output_os[,1:p], center=T, scale = F)
      ggg=c()
      for(i in 1:(nrow(xxx))){
        ggg=rbind(ggg,  c(diag((xxx[i,])%*%t(xxx[i,]))))
      }
      cov_est = cv1 - apply(ggg, 2, function(x){mean(x,trim=trim)})
    }
  } else {
    beta_hat = as.numeric(output_os[1:p])
    cov_est = matrix(output_os[-c(1:(p+2))], 
                     ncol=p, nrow=p) 
  }
  
  sd_est = sqrt(as.numeric(diag(cov_est)))
  res = cbind(beta_hat, sd_est, beta_hat-1.96* sd_est,beta_hat+1.96* sd_est)
  colnames(res) = c('mu', 'sd', 'ci_l', 'ci_u')
  res
}
