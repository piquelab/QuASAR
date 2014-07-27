#' @title fit QuASAR for a multiple samples
#'
#' @description
#' joint geontype across samples and estimate sample error
#'
#' @param ref matrix of reference read count
#' @param alt matrix alternate read count
#' @param eps vector of starting values
#' @param log.gmat log of genotype priors
#' @param max.it maximum number of iterations for algorithm
#' @param tol tolerance to assess convergance
#' @param fixGprior logical to fix genotype priors across all steps of the algorithm
#' @param verbose logical turn on reporting during algorithm iterations
#' @return list of genotypes, log genotypes, vector oferror estimates, log-likelihood, and the sum of
#' log likelihoods 
#' @export
fitAseNullMulti <- function(ref,alt,eps=rep(0.1,ncol(ref)),log.gmat,max.it=100,tol=1E-16,fixGprior=TRUE,verbose=TRUE){
	L <- nrow(ref);
	S <- ncol(ref);
	##browser()
	stopifnot(L == nrow(alt));
	## Parameter initialization
	##log.gmat <- matrix(log(1/3),L,3)
	if(missing(log.gmat)){
		log.gmat <- matrix(log(c(0.25,0.5,0.25)),1,3)
		colnames(log.gmat) <- c("g0","g1","g2");
        }

        it.num <- 0
	logit.eps <- qlogis(eps);
	converged <- FALSE;
	logliksum <- 0;
	rs <- rowSums((ref + alt)) * log(0.5);
	while(it.num <= max.it){
		############ E-step #############
		log.gt <- cbind(ref %*% log(1 - eps) + alt %*% log(eps) + log.gmat[,"g0"],
				rs  + log.gmat[,"g1"],
				ref %*% log(eps) + alt %*% log(1 - eps) + log.gmat[,"g2"]
		)
		colnames(log.gt) <- c("g0","g1","g2");
		log.gt.max <- pmax(log.gt[,"g0"],log.gt[,"g1"],log.gt[,"g2"])
		log.gt <- apply(log.gt,2,function(col){col-log.gt.max})	
		## This normalizes marginal posterior probabilities to add 1
		loglik <- log(exp(log.gt) %*% rep(1,3))+log.gt.max;	
		loglik2 <- log(exp(log.gt) %*% rep(1,3));		
		##browser()		
		## This normalizes marginal posterior probabilities to add 1
		##loglik <- log(exp(log.gt) %*% rep(1,3));
		#norm.genotypes<-t(sapply(1:L,function(ii){
		#	exp(log.gt[ii,])*(exp(loglik[ii])^-1)
		#
		#	}))
		#cat(norm.genotypes[1,],'\n')
		new.logliksum <- sum(loglik);
		#browser()
		if(abs(new.logliksum-logliksum)<tol)
			converged <- TRUE;
		## Normalize such that rows of exp(log.gt) add to 1 
		log.gt <- apply(log.gt,2,function(col){col-loglik2})
		gt <- exp(log.gt);
		if((it.num == max.it) | (converged==TRUE))
			break;
		###### M-STEP  #####
		## w6=w7 and w5=w4 , w2=w1, 
		## epsilon
		converged <- TRUE;
		num <- t(gt[,"g2"]) %*% ref + t(gt[,"g0"]) %*% alt + 0.01;  ## 0.01 added pseudo-counts to avoid 0,0 reads
		den <- t(gt[,"g0"]) %*% ref + t(gt[,"g2"]) %*% alt + 0.01;
		new.logit.eps <- as.vector(log(num) - log(den));
		##browser()
		if(max(abs(logit.eps-new.logit.eps))>tol)
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
		if(verbose){
			cat("#it:",it.num,"eps=",eps,"Post:",colMeans(gt),"loglik",logliksum,"DeltaLogLik",abs(new.logliksum-logliksum),"\n");
		}
		stopifnot(!is.na(eps))
		logliksum <- new.logliksum;
	}
	log.gt <- cbind(ref %*% log(1 - eps) + alt %*% log(eps) + log.gmat[,"g0"],
			rs  + log.gmat[,"g1"],
			ref %*% log(eps) + alt %*% log(1 - eps) + log.gmat[,"g2"]
	)
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
