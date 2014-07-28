##################################################################
## QuASAR Example Workflow
##################################################################
folderName="~/QuASAR/sampleinput/"
fileNames=paste0(folderName,dir(folderName,"Et.*gz"))

##################################################################    
## 1.) Prepare RNA-Seq samples
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
## Saving the output genotype probabilities
out_dat <- data.frame(ase.dat.gt$annotations[, -5], map=ase.joint$gt)
write.table(out_dat, file='genotypes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")

##################################################################
## 3.) ASE inference
##################################################################
ourInferenceData <- aseInference(gts=ase.joint$gt, eps.vect=ase.joint$eps, priors=ase.dat.gt$gmat, ref.mat=ase.dat.gt$ref, alt.mat=ase.dat.gt$alt, min.cov=10, sample.names=sample.names, annos=ase.dat.gt$annotations)
str(ourInferenceData)
