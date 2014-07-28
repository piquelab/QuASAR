###
library(devtools)

document("~/QuASAR",clean=T)


##install_github('QuASAR', 'piquelab')
install_git('https://github.com/piquelab/QuASAR')

##install.package
library('QuASAR')
install.packages("~/QuASAR_1.0.tar.gz")


###
library(QuASAR)



#source(paste(LPG, '/charvey/source/ASE/fitAseModels.v4.R', sep=''))
## qqplot functions
##source(paste(LPG, '/gmb/AI/results/qqman.r', sep=''))
##library('ggplot2')

##require(qvalue)
#qvalue(p = runif(100))

#x11(display="localhost:10.0" ,type="Xlib")


folderName="~/QuASAR/sampleinput/"
fileNames=paste0(folderName,dir(folderName,"Et.*gz"))


##################################################################    
## prepare the relevant samples
ase.dat <- UnionExtractFields(fileNames, combine=TRUE)

## Renaming nice but not necessary
## for(ii in seq_along(1:3)){colnames(ase.dat[[ii]])<-paste("HT", barcodes,sep='')}

ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=5)
str(ase.dat.gt)

################################################################## 
## QuASAR Model Fitting
## ase.joint ~ object for joint genotyoping across all samples
ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))

## Saving the output genotype probabilities
out_dat <- data.frame(ase.dat.gt$annotations[, -5], map=ase.joint$gt)
write.table(out_dat, file='genotypes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")






## Mandatory objects for inference
ase.dat.final <- list(ref=ase.dat.gt$ref, alt=ase.dat.gt$alt, gmat=ase.dat.gt$gmat, annotations=ase.dat.gt$annotations)
ase.joint.eps <- ase.joint$eps

n.eps <- length(ase.joint.eps)

sample.names <- colnames(ase.dat.final$ref)


test <- aseInference(gts=ase.joint$gt, priors=ase.dat.final$gmat, ref.mat=ase.dat.final$ref, alt.mat=ase.dat.final$alt, eps.vect=ase.joint.eps, min.cov=10, sample.names=sample.names)

###############################################
## save the data list as an R object
min.cov <- 10
dat.name <- paste(output.folder, "/",plate, "_",cell.line, "_cov", min.cov, '_inference.RData', sep='')
#save(inference.data, file=dat.name)

##########################################
## 1.) QQ-plots for all treatments
##########################################
for(ii in seq_along(1:length(inference.data))){
	treatment <- names(inference.data)[ii]
	pvals <- inference.data[[ii]]$dat$pval3
	pi0 <- round(inference.data[[ii]]$meta.dat[2], 2)
	hets <- inference.data[[ii]]$meta.dat[1]
	qv.2 <- inference.data[[ii]]$meta.dat[5]
	png.file <- paste('./plots/QQ/', plate, '_', cell.line, '_', treatment, "_cov", min.cov, '_', ii, '_QQ', '.png', sep='')
	png(file=png.file)
		qq(pvals)
		title <- paste(plate, " Cell ", cell.line, ": ", treatment, ' | Pi0=', pi0, ' | #hets=', hets, ' | #qv.2=', qv.2, sep='')
		title(main=title)
	dev.off()

}

##########################################
## 2.) Expression table across all treatments
##########################################
asetable <- t(sapply(seq_along(1:length(inference.data)), FUN=function(ii){
			sapply(c(.01, .05, .1, .2), FUN=function(jj){sum(inference.data[[ii]]$dat$qvals.qv3<jj)})
	}))
rownames(asetable) <- names(inference.data)
colnames(asetable) <- c('Q<.01', 'Q<.05', 'Q<.1', 'Q<.2')

outfile <- paste('./output/', plate, "_",cell.line, '_Qhits.txt', sep='')
write.table(asetable, file=outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t") 

###### THE END ######
###### THE END ######
###### THE END ######
