##################################################################  
## prepare the ASE data for genotyping
################################################################## 
PrepForGenotyping <- function(ase.dat, min.coverage, dampen.priors=TRUE){
	ref.all <- ase.dat[[1]]
	alt.all <- ase.dat[[2]]
	phi.all <- ase.dat[[4]]$af 
	n.samples <-dim(ref.all)[2]

	##################################################################  
	## reads.floor ~ minumum coverage across all samples
	## collapsed.indicator ~ loci with minimum coverage
	## dat.collapssed.final ~ dat.collapsed with >reads.floor minimum coverage
	## ref.all.final ~ sample wise reference counts for loci with sufficient covergae
	## alt.all.final ~ sample wise alternate counts for loci with sufficient covergae
	## annotations.included.samples ~ annotations for loci with sufficient coverage
	## gmat ~ genotype priors from 1k genomes for loci with combined coverage >15
	## gmat.collapsed ~ genotype priors from 1k genomes for all loci
	reads.floor <- min.coverage
	collapsed.indicator <- (rowSums(ref.all + alt.all) > reads.floor)
	ref.all.final <- ref.all[collapsed.indicator, ]		
	alt.all.final <- alt.all[collapsed.indicator, ]
	phi.all.final <- phi.all[collapsed.indicator]
	annotations.included.samples <- ase.dat$anno[collapsed.indicator, ]
	gmat <- cbind(g0=(1-phi.all.final)^2, g1=2*phi.all.final*(1-phi.all.final), g2=phi.all.final^2)

	if(dampen.priors){
		gmat <- (gmat+0.0001)/1.0003
	}

	list(ref=ref.all.final, alt=alt.all.final, gmat=gmat, annotations=annotations.included.samples)
}
