#' @title UnionExtractFields
#'
#' @description
#' from a list of clean bed files, find the intersection of all loci and return a list with reference, alternate, and error counts for each loci   
#' across all samples
#' @param fileList list of clean bed files as described in the QuASAR vingette
#' @return  
#'
#' @export
UnionExtractFields <- function(fileList, combine=FALSE){			
	#browser()
	tmpFile <- scan(pipe("mktemp -t"),character(0))
	system(paste("less ", paste(fileList, collapse=" "), " | grep -v -w '^chr\\|^NA' | cut -f 1-3,6-7 | sortBed -i stdin | uniq | gzip > ", tmpFile))
	##TODO a better way to do this is to have the directory and the file list as separate arguments.
	#sNames <- strsplit(gsub("/nfs/hpnfs/groups/piquelab/charvey/ASE/prod2/", "", gsub(".pileup.clean.bed.gz", "", fileList)), "/")
	#sNames <- gsub("_NEB_Rep1.q10.p10", "", sapply(sNames, `[`, 1))	

	anno <- read.table(gzfile(tmpFile),sep="\t",as.is=T)	

	aux <- sapply(fileList,function(fn){			
				#fn <- fileList[1]
				cat("Processing:",fn,"\n")
				command=paste("intersectBed -a ",tmpFile," -b ",fn," -wao | cut -f 1-3,13-15 ",sep="")
				aa <- read.table(pipe(command),sep="\t",as.is=T,na.strings=".")
				aa[is.na(aa)] <- 0
				stopifnot(identical(aa[,1:3],anno[,1:3]))
				aa[,-(1:3)]
			})

	colnames(anno) = c("chr","pos0","pos", "rsID", "af")

	Ref <- as.matrix(do.call(cbind,aux[1,]))
	#colnames(Ref) <- sNames
	Alt <- as.matrix(do.call(cbind,aux[2,]))
	#colnames(Alt) <- sNames
	Err <- as.matrix(do.call(cbind,aux[3,]))
	#colnames(Err) <- sNames

	return.list<-list(ref=Ref,alt=Alt,err=Err,anno=anno);

	if(combine==TRUE){
		allRef<-apply(Ref, MARGIN=1, sum)
		allAlt<-apply(Alt, MARGIN=1, sum)
		allErr<-apply(Err, MARGIN=1, sum)

		return.list$all<-as.matrix(cbind(allRef, allAlt, allErr))

	}
	return(return.list)
}
