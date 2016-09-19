##' @title UnionExtractFields
##' @description
##' From a list of files (fileList), find the union of all loci and return a list with the data prepared to run QuASAR.
##' @param fileList List of files *.quasar.in.gz
##' @param combine Collapses all samples into a single sample if true. Default is false
##' @return returns a R list with the following elements:
##' Ref: Matrix with number of reads matching the reference allele.
##' Alt: Matrix with number of reads matching the alternate allele.
##' Err: Matrix with number of reads matchine neither the ref. or alt. allele.
##' anno: Object with the annotation for the SNPs. 
##' @author Chris Harvey
#' @export
UnionExtractFields <- function(fileList, combine=FALSE){			
	#browser()
	tmpFile <- scan(pipe("mktemp -t"),character(0))
	system(paste("less ", paste(fileList, collapse=" "), " | grep -v -w '^chr\\|^NA' | cut -f 1-7 | sortBed -i stdin | uniq | gzip > ", tmpFile))
	anno <- read.table(gzfile(tmpFile),sep="\t",as.is=T)	
	aux <- sapply(fileList,function(fn){			
				#fn <- fileList[1]
				cat("Processing:",fn,"\n")
				command=paste("intersectBed -a ",tmpFile," -b ",fn," -wao | cut -f 1-3,15-17 ",sep="")
				aa <- read.table(pipe(command),sep="\t",as.is=T,na.strings=".")
				aa[is.na(aa)] <- 0
				stopifnot(identical(aa[,1:3],anno[,1:3]))
				aa[,-(1:3)]
			})
	colnames(anno) = c("chr","pos0","pos","ref","alt","rsID", "af")
	Ref <- as.matrix(do.call(cbind,aux[1,]))
	Alt <- as.matrix(do.call(cbind,aux[2,]))
	Err <- as.matrix(do.call(cbind,aux[3,]))
	return.list<-list(ref=Ref,alt=Alt,err=Err,anno=anno);
	if(combine==TRUE){
		allRef<-apply(Ref, MARGIN=1, sum)
		allAlt<-apply(Alt, MARGIN=1, sum)
		allErr<-apply(Err, MARGIN=1, sum)
		return.list$all<-as.matrix(cbind(allRef, allAlt, allErr))

	}
	return(return.list)
}
