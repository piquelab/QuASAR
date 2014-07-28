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
<!---
Installing R packages from GitHub within an R session has sometimes problems. Alternatively, you can clone/fork this repository then and then build the package:
-->

```C
git clone git@github.com:piquelab/QuASAR.git
R CMD build QuASAR
```

## 2. Preprocessing
### Alignment & filtering


### Pileups & cleaned pileups
```C
zless SomeSample_1.pileup.clean.bed.gz | head -5
chr1	879910	879911	G	A	rs143853699	0.02	21	0	0
chr1	892379	892380	G	A	rs150615968	0.0041	22	0	0
chr1	893384	893385	G	A	rs140972868	0.01	6	0	0
chr1	894101	894102	A	T	rs188691615	0.01	6	0	0
chr1	894430	894431	G	A	rs201791495	9e-04	9	0	0
```

## 3. Genotyping with QuASAR
### Genotyping a single or multiple samples
```R
ase.joint <- fitAseNull(finalref, finalalt, log.gmat=log(ase.dat.gt$gmat))
```
```R
ase.joint <- fitAseNullMulti(finalref, finalalt, log.gmat=log(ase.dat.gt$gmat))
```


### Inference on ASE
### Sample workflow



