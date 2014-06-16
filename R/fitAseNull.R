#' @title fit QuASAR for a single sample
#'
#' @description
#' \code(fitAseNull) for a single sample conduct genotyping and estimate sample error
#'
#' @param ref reference read count
#' @param alt alternate read count
#' @param eps starting value
#' @param log.gmat log of genotype priors
#' @param max.it maximum number of iterations for algorithm
#' @param tol tolerance to assess convergance
#' @param fixGprior logical to fix genotype priors across all steps of the algorithm
#' @return list of genotypes, log genotypes, error estimate, log-likelihood, and the sum of
#' log likelihoods 
fitAseNull <- function(ref,alt,eps=0.1,log.gmat,max.it=100,tol=1E-16,fixGprior=TRUE){
  L <- length(ref);
  stopifnot(L == length(alt));
  ## Parameter initialization
  ##log.gmat <- matrix(log(1/3),L,3)
  if(missing(log.gmat)){
    log.gmat <- matrix(log(c(0.25,0.5,0.25)),1,3)
    colnames(log.gmat) <- c("g0","g1","g2");
  }
  it.num <- 0
  logit.eps <- qlogis(0.1);
  converged <- FALSE;
  logliksum <- 0;
  while(it.num <= max.it){
  ############ E-step #############
    log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps) + log.gmat[,"g0"],
                    g1 = (ref + alt) * log(0.5) + log.gmat[,"g1"],
                    g2   = ref * log(eps) + alt * log(1 - eps) + log.gmat[,"g2"]
                    )
    ##log.gt[log.gt < (-200)] <- (-200)
    ##browser()
    ## This normalizes marginal posterior probabilities to add 1
	log.gt.max <- pmax(log.gt[,"g0"],log.gt[,"g1"],log.gt[,"g2"])
	log.gt <- apply(log.gt,2,function(col){col-log.gt.max})	
	loglik <- log(exp(log.gt) %*% rep(1,3))+log.gt.max;	
	loglik2 <- log(exp(log.gt) %*% rep(1,3));		
	
    new.logliksum <- sum(loglik);
    if(abs(new.logliksum-logliksum)<tol)
      converged <- TRUE;
    log.gt <- apply(log.gt,2,function(col){col-loglik2})
    gt <- exp(log.gt);
    if((it.num == max.it) | (converged==TRUE))
      break;
    ###### M-STEP  #####
    ## w6=w7 and w5=w4 , w2=w1, 
    ## epsilon
    converged <- TRUE;
    new.logit.eps <- log(sum( gt[,"g2"] * ref + gt[,"g0"] * alt)) - log(sum( gt[,"g0"] * ref + gt[,"g2"] * alt))
   ## browser()
    if((logit.eps-new.logit.eps)>tol)
      converged <- FALSE;
    logit.eps <- new.logit.eps;
    eps <- plogis(logit.eps)
    ## ## ## opt 2
    if(!fixGprior){
      log.gmat[,"g0"] <- log.gt[,"g0"]
      log.gmat[,"g1"] <- log.gt[,"g1"]
      log.gmat[,"g2"] <- log.gt[,"g2"]
    }
    ## ##
    it.num <- it.num +1;
    cat("#it:",it.num,"eps=",eps,"Post:",colMeans(gt),"loglik",logliksum,"DeltaLogLik",abs(new.logliksum-logliksum),"\n");
    stopifnot(!is.na(eps))
    logliksum <- new.logliksum;
  }
  colnames(log.gt) <- c("g0","g1","g2");
  ##log.gt[log.gt < (-200)] <- (-200)
  log.gt.max <- pmax(log.gt[,"g0"],log.gt[,"g1"],log.gt[,"g2"])
  log.gt <- apply(log.gt,2,function(col){col-log.gt.max})	
  ## This normalizes marginal posterior probabilities to add 1
  loglik2 <- log(exp(log.gt) %*% rep(1,3));
  log.gt <- apply(log.gt,2,function(col){col-loglik2})
  gt <- exp(log.gt);
    
  invisible(list(gt=gt,log.gt=log.gt,eps=eps,loglik=loglik,logliksum=logliksum))
}
