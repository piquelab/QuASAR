# QuASAR: Quantitative Allele Specific Analysis of Reads
QuASAR ([Harvey et al, 2015]) is an R package, that implements a statistical method for: i) genotyping from next-generation sequencing reads, and ii) conducting inference on allelic imbalance at heterozygous sites. The sequencing data can be RNA-seq, DNase-seq, ATAC-seq or any other type of high-throughput sequencing data. The input data to QuASAR is a processed pileup file (as detailed later). Here, we do not cover in depth important pre-processing steps such as choice of the aligner, read filtering and duplicate removal.

We also want to emphasize that the current software is still in development, we would kindly appreciate any comments and bug reports.
<!---
Prior to analsyis, RNA-Seq data must undergo alignment with a modern aligner, quality filtering, duplicate removal, and the creation of pileups. There are many tools and tutorials available for preprocessing Next Generation Sequencing data, but we will only describe the tools we used and expect the user to have basic familiarity with standard bioinformatics command-line tools. Our goal with this tutorial is to cover the following:

1. Installing QuASAR
2. Preprocessing 
   * Alignment, filtering, and removing duplicates. (Description of, not a tutorial how)
   * Pileups and clean pileups
3. QuASAR analyis pipeline
   * Genotyping single or multiple samples
   * Inference on ASE
   * Sample workflow

**Quick-start**: Users comfortable processing RNA-Seq data to the level of pileups should skip to the second step of preprocessing. 
-->

## 1. Installation

To install from within an R session:

```R
require(devtools)
install_github('piquelab/QuASAR')
library('QuASAR')
```

However, this method is occasionally problematic. Alternatively, you can clone/fork this repository and then build the package:
```C
git clone https://github.com/piquelab/QuASAR.git
R CMD build QuASAR
```
then in R,
```R
install.packages('QuASAR_x.y.tar.gz')
library(QuASAR)
```

## 2. Preprocessing
### Alignment & filtering
Raw reads can be aligned to the reference genome using your favorite aligner. Because allele-specific analysis is extremely sensitive to read biases and mapping errors, we strongly recommend adding steps to remove PCR duplicates and to remove reads aligning to areas with known mappability issues (e.g., [Degner et al, 2009]).


### Pileups & cleaned pileups
Note: These steps require [samtools] and [bedtools].

Using the samtools mpileup command, create a pileup file from aligned reads. Provide a fasta-formatted reference genome (hg19.fa) and a bed file of positions you wish to pileup on (e.g., 1K genomes SNP positions [1KG snp file]):

```C
samtools mpileup -f hg19.fa -l snps.af.bed input.bam | gzip > input.pileup.gz
```

Next, convert the pileup file into bed format and use intersectBed to include the allele frequencies from a bed file. The bed file with allele frequencies should be seven columns: 1-3) coordinate, 4) SNP ID, 5) reference allele, 6) alternate allele, 7) allele frequency. This can be the same [1KG snp file] used in the pileup stage. The awk filter step (below) removes positions not covered by a read, positions covered by indels, and reference skips:

```C
less input.pileup.gz | awk -v OFS='\t' '{ if ($4>0 && $5 !~ /[^\^][<>]/ && $5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $5 !~ /-[0-9]+[ACGTNacgtn]+/ && $5 !~ /[^\^]\*/) print $1,$2-1,$2,$3,$4,$5,$6}' | sortBed -i stdin | intersectBed -a stdin -b snps.af.bed -wo | cut -f 1-7,11-14 | gzip > input.pileup.bed.gz
```

Finally, get the read counts at each position, and, if desired, perform any additional filtering. The result will be the input file for QuASAR. An example processing script is provided here: [scripts/convertPileupToQuasar.R].

```C
R --vanilla --args input.pileup.bed.gz < convertPileupToQuasar.R
```

Here is an example of how the QuASAR infput file should look:

```C
zless input.quasar.in.gz | head -5
chr1	879910	879911	G	A	rs143853699	0.02	21	0	0
chr1	892379	892380	G	A	rs150615968	0.0041	22	0	0
chr1	893384	893385	G	A	rs140972868	0.01	6	0	0
chr1	894101	894102	A	T	rs188691615	0.01	6	0	0
chr1	894430	894431	G	A	rs201791495	9e-04	9	0	0
```

The fields are as follows: 
1. Chromosome 
2. Start position 
3. End position 
4. Reference allele 
5. Alternate allele 
6. SNP ID 
7. SNP allele frequency 
8. Number of reads mapping to the reference allele 
9. Number of reads mapping to the alternate allele 
10. Number of reads not mapping to either allele

## 3. Running QuASAR

### Prepare the input samples 
For a test run we provide a small sample dataset containing 6 samples from the same individual. 
The following commands will download the data to the current folder:

```R
urlData="http://genome.grid.wayne.edu/quasar/sampleinput/"
fileNames <- paste0("t",c(2,4,6,12,18,24),"hr_Huvec_Rep1.quasar.in.gz")
sapply(fileNames,function (ii) download.file(paste0(urlData,ii),ii))
```

To run the sample data, or any data, we provide a few helper functions to merge samples across the union of all annotated sites (`UnionExtractFields`), and to filter sites with insufficient coverage across all samples (`PrepForGenotyping`). Note: these functions utilize calls to [bedtools].

```R
ase.dat <- UnionExtractFields(fileNames, combine=TRUE)
ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=5)
sample.names <- colnames(ase.dat.gt$ref)
```

### Genotyping an individual from multiple samples
Genotyping an individual using `fitAseNullMulti` requires a matrix of reference counts and a matrix of alternate counts where where the columns are ordered by sample. The final argument is a matrix of priors for the minor allele frquency, for which we use the 1K genomes MAFs assumed to be at Hardy-Weinberg equilibrium.  
```R
ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))
```
This function returns a list with the following members:
```R
names(ase.joint)
[1] "gt"        "log.gt"    "eps"       "loglik"    "logliksum"
```
where the posterior probability of the genotypes, `gt`, across all samples are accessed as follows:
```C
head(ase.joint$gt)
               g0           g1           g2
[1,] 2.870026e-98 1.000000e+00 2.939460e-70
[2,] 1.465195e-27 7.773259e-04 9.992227e-01
[3,] 3.732811e-61 4.308038e-07 9.999996e-01
[4,] 9.992226e-01 7.774208e-04 1.714236e-27
[5,] 9.435425e-87 9.726281e-10 1.000000e+00
[6,] 9.999863e-01 1.372351e-05 6.274482e-46
```

g0=homozygous reference, g1=heterozygous, & g2=homozygous alternate. To save the output genotype probabilities together with the SNP annotation, we do:                                                                                                
```R
out_dat <- data.frame(ase.dat.gt$annotations[, -5], map=ase.joint$gt)
write.table(out_dat, file='genotypes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
```

Estimates of error parameters `eps` for each sample are:

```C
ase.joint$eps
[1] 0.0008748778 0.0007617141 0.0008152132 0.0007819780 0.0008956686
[6] 0.0007597717
```


### Inference on ASE
Using `aseInference` to conduct inference on ASE for an individual requires the posterior probabilities of each genotypes from the previous step `"gt"`, estimates of sequencing error for each sample `"eps"`, the same priors used in the previous step, reference counts, alternate counts, minimum coverage, sample names, and variant annotations. 
```R
ourInferenceData <- aseInference(gts=ase.joint$gt, eps.vect=ase.joint$eps, priors=ase.dat.gt$gmat, ref.mat=ase.dat.gt$ref, alt.mat=ase.dat.gt$alt, min.cov=10, sample.names=sample.names, annos=ase.dat.gt$annotations)
```
This function returns a list where each element corresponds to an input sample:
```R
names(ourInferenceData[[1]])
[1] "dat"        "n.hets"     "dispersion"
```
where `dat` contains estimates of allelic imbalance `betas`, standard errors `betas.se`, & pvalues from an LRT for ASE detailed in [Harvey et al, 2014]. Note that the number of rows (SNPs) in each sample corresponds to the the number of heterozygous SNPs passing a minimum coverage filter. 
```R
 head(ourInferenceData[[1]]$dat)
 annotations.rsID annotations.chr annotations.pos0       betas  betas.se    pval3 
1        rs2272757            chr1           881626  0.15175892 0.6005410 0.80049721
2        rs2465128            chr1           981930  0.17948875 0.6445723 0.78065789
3        rs9442391            chr1           984301 -0.15175892 0.6005410 0.80049721
4       rs12142199            chr1          1249186 -0.43478406 0.4845478 0.36955958
5           rs7290            chr1          1477243 -0.99328368 0.5969363 0.09611857
6           rs7533            chr1          1479332 -0.09853221 0.3981711 0.80455070
```
The final members of the list are the number of heterozygotes and the esimtate of dispersion for each sample.
```R
head(ourInferenceData[[1]]$n.hets)
[1] 2856
head(ourInferenceData[[1]]$dispersion)
[1] 64.07152
```

The code for this sample workflow is located here:
[scripts/exampleWorkflow.R]

<!-- links -->
[Harvey et al, 2015]:http://bioinformatics.oxfordjournals.org/content/31/8/1235
[Degner et al, 2009]:http://www.ncbi.nlm.nih.gov/pubmed/19808877
[samtools]:http://samtools.sourceforge.net/
[bedtools]:https://github.com/arq5x/bedtools2
[scripts/convertPileupToQuasar.R]:scripts/convertPileupToQuasar.R
[scripts/exampleWorkflow.R]:scripts/exampleWorkflow.R
[1KG snp file]:http://genome.grid.wayne.edu/centisnps/1kgSnps.html
