###################################################################
## Generate read counts covering reference and alternate alleles ##
## from a bed-formatted pilup file 															 ##
###################################################################

# Get the numerical base call quality score form the reported ascii value
qual <- function(char) { strtoi(charToRaw(char),16L)-33 };

## Default values
mincov <- 4				# Minimum coverage
maxcov <- 200000 	# Max coverage

## Get inputs from command line
cargs <- commandArgs(trail=TRUE);
if(length(cargs) >= 1)
	pileupFile <- cargs[1];
if(length(cargs) >= 2)
	mincov <- cargs[2];
if(length(cargs) >= 3)
	maxcov <- cargs[3];

# Get the input data
command <- paste("less ", pileupFile,
		"| awk ' $5 >=", mincov, " && $5 <=", maxcov, "'")
pileup <- read.table(file=pipe(command), header=F, quote="", comment.char="",
	as.is=T, sep="\t") 
names(pileup) <- c("chr", "pos-1", "pos", "ref", "num.reads", "read.alleles",
	"read.quality", "rsID", "TKG.Ref", "alt", "af")

# See if the ref allels match, then discard uncesessary columns
indMatch <- (toupper(pileup$ref) == pileup$TKG.Ref)
pileup <- pileup[indMatch, c("chr", "pos-1", "pos", "ref", "alt", "rsID",
	"num.reads", "read.alleles", "read.quality", "af")]
stopifnot(mean(indMatch)>0.8) ## Stop if too many errors
rm(indMatch)

# Duplicates can arise for (at least) 3 reasons: tri+ alleleic SNPs,
# indels (should already be filtered), and incongruencies between 
# genome assembly reference allele and 1KG reference allele
d1 <- duplicated(paste(pileup$chr, pileup$pos, sep=":"))
d2 <- duplicated(paste(pileup$chr, pileup$pos, sep=":"), fromLast=T)
pileup <- pileup[!(d1 | d2), ]
rm(d1,d2)

## Filter the alleles to remove those at the beginning 
## and end of a mapped read, then remove records with no reads
pileup$read.alleles.filt <- mapply(gsub, '[a-zA-Z., ]\\$', '$', 
																		pileup$read.alleles)
pileup$read.alleles.filt <- mapply(gsub, '\\^[[:punct:][:alnum:]][a-zA-Z., ]',
																		'^', pileup$read.alleles.filt)
pileup <- pileup[nchar(pileup$read.alleles.filt)>0, ]

## Examine the base quality scores. Output summaries, but 
## don't filter anything
qual <- sapply(1:nrow(pileup), function(ii){qual(pileup$read.quality[ii])})
qual <- unlist(qual)
qtr <- quantile(qual, seq(0, 1, 0.1))
qual.table <- table(qual);
qtr
qual.table

## Clean up the read alleles and count the matches
pileup$read.alleles.clean <- mapply(gsub, '[\\.\\, ]', pileup$ref, 
																		pileup$read.alleles.filt)
pileup$read.alleles.clean <- toupper(pileup$read.alleles.clean)
pileup$ref <- toupper(pileup$ref)
pileup$ref.matches <- as.integer(nchar(mapply(gsub, paste('[^', 
	as.character(pileup$ref), ']', sep=""), '', pileup$read.alleles.clean)))
pileup$alt.matches <- as.integer(nchar(mapply(gsub, paste('[^', 
	as.character(pileup$alt), ']', sep=""), '', pileup$read.alleles.clean)))

## Log the number of reads not matching either ref or alt allele
pileup$errors <- as.integer(nchar(mapply(gsub, paste('[^ACGT]', sep=""), 
	'', pileup$read.alleles.clean))) - (pileup$alt.matches + pileup$ref.matches)

## Reorder by position so chr names appear right
pileup <- pileup[order(pileup$chr, pileup$pos), c("chr", "pos-1", "pos", 
	"ref", "alt", "rsID", "af", "ref.matches", "alt.matches", "errors")];

## Output the clean pileup file
oName <- gsub(".*/", "", gsub(".pileup.bed.gz", "", pileupFile));
outFile <- paste(oName, ".quasar.in.gz", sep="");
write.table(pileup, gzfile(outFile), quote=F, col.names=F, 
	row.names=F, sep="\t")
