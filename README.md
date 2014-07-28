# QuASAR: Quantitative allele specific analysis of reads
QuASAR is an R package, that implements a statistical method for: i) joint genotyping across sequencing datasets of the same individual, ii) identifying heterozygous loci, and iii) conducting inference on allelic imbalance. 
The sequencing data can be RNA-seq, DNase-seq, ATAC-seq or any other type of high-throroughput sequencing data. 
The input data to QuASAR is already a clean up pileup as it will be detailed later. 
Here, we do not cover in important pre-processing steps such as choice of the aligner, read filtering and duplicate removal. 

We also want to emphisize that the current software is still in development, we would kindly appreciate any comments and bug reports. 

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

```R
require(devtools)
install_github('QuASAR', 'piquelab')
library('QuASAR')
```

Installing R packages from GitHub within an R session has sometimes problems. Alternatively, you can clone/fork this repository then and then build the package:

```C
git clone git@github.com:piquelab/QuASAR.git
R CMD build QuASAR
```

## 2. Preprocessing
### Alignment & filtering


### Pileups & cleaned pileups
Using the samtools mpileup command, create a pileup file from aligned reads. Provide a fasta-formatted reference genome (hg19.fa) and a bed file of positions you wish to pileup on (e.g., 1KG SNP positions):

```C
samtools mpileup -f hg19.fa -l snps.af.bed input.bam | gzip > input.pileup.gz
```

Next, convert the pileup file into bed format and use intersectBed to include the allele frequencies from a bed file. The awk filter step removes positions not covered, positions covered by indels, and reference skips:

```C
less input.pileup.gz | awk -v OFS='\t' '{ if ($4>0 && $5 !~ /[^\^][<>]/ && $5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $5 !~ /-[0-9]+[ACGTNacgtn]+/ && $5 !~ /[^\^]\*/) print $1,$2-1,$2,$3,$4,$5,$6}' | sortBed -i stdin | intersectBed -a stdin -b snps.af.bed -wo | cut -f 1-7,11-14 | gzip > input.pileup.bed.gz
```

Finally, run convertPileupToQuasar.R to generate a file of read counts at each SNP position:

```C
R --vanilla --args input.pileup.bed.gz < convertPileupToQuasar.R
```

The final file should look something like this:

```C
zless input.pileup.clean.bed.gz | head -5
chr1	879910	879911	G	A	rs143853699	0.02	21	0	0
chr1	892379	892380	G	A	rs150615968	0.0041	22	0	0
chr1	893384	893385	G	A	rs140972868	0.01	6	0	0
chr1	894101	894102	A	T	rs188691615	0.01	6	0	0
chr1	894430	894431	G	A	rs201791495	9e-04	9	0	0
```

The final fields are as follows:
1. Chromosome
2. Start position
3. End position
4. Reference allele
5. Alternate allele
6. SNP ID
7. SNP allele frequency
8. Number of reference reads
9. Number of alternate reads
10. Number of reads not matching reference or alternate

## 3. Running QuASAR
### Genotyping a single or multiple samples
```R
ase.joint <- fitAseNull(finalref, finalalt, log.gmat=log(ase.dat.gt$gmat))
```
```R
ase.joint <- fitAseNullMulti(finalref, finalalt, log.gmat=log(ase.dat.gt$gmat))
```


### Inference on ASE
### Sample workflow



