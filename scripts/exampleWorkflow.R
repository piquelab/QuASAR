##################################################################
## QuASAR Example Workflow
##################################################################

##################################################################    
## 1.) Download example pre-processed data
##################################################################


urlData="http://genome.grid.wayne.edu/quasar/sampleinput/"
fileNames <- paste0("t",c(2,4,6,12,18,24),"hr_Huvec_Rep1.quasar.in.gz")
sapply(fileNames,function (ii) download.file(paste0(urlData,ii),ii))


##################################################################    
## 2.) Load data into an R object
##################################################################

ase.dat <- UnionExtractFields(fileNames, combine=TRUE)
ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=5)
str(ase.dat.gt)
sample.names <- colnames(ase.dat.gt$ref)

################################################################## 
## 2.) QuASAR Model Fitting
## ase.joint ~ object for joint genotyoping across all samples
##################################################################

ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))
str(ase.joint)

## The ase.joint$gt object contains the posterior probability for each genotype: g0) is homozygote reference, g1) heterozygote, g2) homozygote alternate.
head(ase.joint$gt)

## The base-calling error parameters are stored in ase.joint$eps

## Saving the output genotype probabilities
out_dat <- data.frame(ase.dat.gt$annotations[, -5], map=ase.joint$gt)
write.table(out_dat, file='genotypes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")

##################################################################
## 3.) ASE inference
##################################################################
ourInferenceData <- aseInference(gts=ase.joint$gt, eps.vect=ase.joint$eps, priors=ase.dat.gt$gmat, ref.mat=ase.dat.gt$ref, alt.mat=ase.dat.gt$alt, min.cov=10, sample.names=sample.names, annos=ase.dat.gt$annotations)
str(ourInferenceData)
