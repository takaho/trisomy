# trisomy
R script to examine chromosomal abnormality using RNA-seq data

This script uses output file of cuffdiff. Expression changes were grouped for each chromosome and significance would be calculated against total mean.
Box plot should be generated in PDF format and saved.

## Run the script as independent program
If you run the script from command line, following command let use obtain boxplot graph for each chromosome.

Rscript trisomy.R [filename] <[pdf filename]>

First argument must be the output file of cuffdiff version 2, usually named 'genes_exp.diff.' Second argument is optional. If you use non-default filename in boxplot PDF, input the name here.

## R interactive mode
The function script can be used interactively in R command line.

'''
source(“path/to/trisomy.R”)
drawExpressionChangesByChromosome("test.diff", outputFilename="sample.pdf", horizontal=FALSE, excludesXY=FALSE, ylimit=6)
'''

Argument|Value
------|-----
outputFilename|PDF filename (default:boxplot.pdf)
horizontal|Graph direction (default:TRUE)
excludesXY|Ignore or include genes on sex chromosomes (default:TRUE)
ylimit|Y axis limit (default:4)

