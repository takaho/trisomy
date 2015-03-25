#!/usr/bin/env/Rscript

if (FALSE) {
"
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"

"
Box plot drawing program written in R.
This program can be run from command line or from interactive mode.

If you run this program from command line, use following format.
Rscript trisomy.R [input filename] [output filename (default:boxplot.pdf)]

Or graphs can be drawn in interactive mode as follows.

> drawExpressionChangesByChromosome(inputFilename, outputFilename)

Written by Takaho A. End <takaho.endo@riken.jp>
"
}

getChromosomeIndex <- function(chromosome) { 
"Converter function of chromosome name into number to order chromosome correctly"
    if (setequal(chromosome, "chrX")) {
        index <- 100;
    } else if (setequal(chromosome, "chrY")) {
        index <- 101;
    } else if (setequal(chromosome, "chrM")) {
        index <- 200;
    } else if (setequal(chromosome, "chrW")) {
        index <- 102;
    } else if (setequal(chromosome, "chrZ")) {
        index <- 103;
    } else {
        m <- regexpr("chr[0-9]+", chromosome);
        if (m) {
            as.numeric(substr(chromosome, 4, attr(m, "match.length") + 1));
            index <- as.numeric(substr(chromosome, 4, attr(m, "match.length") + 1));
        } else {
            index <- 100000;
        }
    }
    index;
}

getChromosomeIndices <- function(chromosomes) {
"Give a vector of chromosome indices"    
    indices <- c();
    for (i in 1:length(chromosomes)) {
        indices <- c(indices, getChromosomeIndex(chromosomes[[i]]));
    }
    indices;
}

drawExpressionChangesByChromosome <- function(inputFilename, outputFilename="boxplot.pdf", minimumRPKM=1e-2, excludesXY=TRUE, horizontal=TRUE, ylimit=4, ttest=1.0) {
    "
Read cuffdiff output file (gene_exp.diff) and draw a box plot.
Arguments:
  inputFilename  : Cuffdiff output file. 
  outputFilename : Filename of graph PDF. (default : boxplot.pdf)
  minimumRPKM    : Ignore genes having very little expression. (default : 0.01)
  excludesXY     : Excludes sex chromosome or not. (default : TRUE)
  horizontal     : Graph direction (default: TRUE)
  ylimit         : Y-axis range (default: 4)
  ttest          : Threshold of significance calculated using Student t-test (default: 1.0 = no test)
"

    data<-read.delim(inputFilename, header=TRUE, sep="\t", colClasses=c(rep("character", 7), rep("numeric", 6), "character")); # read tab-separated text

    # Column number of Cuffdiff format. Revise numbers if the format is altered.
    geneColumn <- 3;
    locusColumn <- 4;
    stateColumn <- 7;
    valueColumn1 <- 8;
    valueColumn2 <- 9;
    fcColumn <- 10;
    expressionData <- list();
    minimumRPKM <- 1e-2;
    chromosomes <- c();
    sumlfc <- 0.0;
    num <- 0;
    for (i in 1:nrow(data)) {
        row <- data[i,];
        locus <- row[[locusColumn]];
        m <- regexpr("chr[0-9MWXYZ]{1,2}:", locus);
        if (m[[1]] > 0) {
            pos <- m[[1]];
            chromosome <- substr(locus, pos, pos + attr(m, "match.length") - 2);
            if (excludesXY & getChromosomeIndex(chromosome) >= 100) {
                next;
            }
            state <- row[[stateColumn]];
            value1 <- row[[valueColumn1]];
            value2 <- row[[valueColumn2]];
            if (setequal(state, "OK") & value1 > minimumRPKM & value2 > minimumRPKM) {
                log2fc <- row[[fcColumn]];
                sumlfc <- sumlfc + log2fc;
                num <- num + 1;
                if (is.null(expressionData[[chromosome]])) {
                    chromosomes <- c(chromosomes, chromosome);
                    expressionData[[chromosome]] <- c(log2fc);
                } else {
                    expressionData[[chromosome]] <- c(expressionData[[chromosome]], log2fc);
                }
            }
        }
    }

    # sort chromosomes
    if (horizontal) { 
        chromosomes <- chromosomes[order(-getChromosomeIndices(chromosomes))];
        labelDirection <- 1;
    } else {
        chromosomes <- chromosomes[order(getChromosomeIndices(chromosomes))];
        labelDirection <- 2;
    }

    dl <- list();
    avg <- sumlfc / num; # normalization coefficient with linear correction
    for (i in 1:length(chromosomes)) {
        chrom <- chromosomes[[i]];
        dl[[chrom]] <- expressionData[[chrom]] - avg;
    }

    # calculate p-values
    pvalues <- list();
    if (ttest >= 1.0 | ttest <= 0.0) {
        colors <- "lightgray";
    } else {
        pthr <- ttest / length(chromosomes); # Bonferroni correction
        colors <- c();
        for (i in 1:length(chromosomes)) {
            chrom <- chromosomes[[i]];
            values <- expressionData[[chrom]];
            color <- "lightgray";
            if (length(values) >= 5) { # ignore chromosomes having few values
                pvalue <- t.test(values, mu=avg);
                if (pvalue[["p.value"]] < pthr) {
                    color <- "red";
                }
                pvalues[[chrom]] <- c(length(values), mean(values), pvalue[["p.value"]]);
            } else {
                pvalues[[chrom]] <- c(length(values), 1.0);
            }
            colors <- c(colors, color);
        }
    }

    # output statistics
    chromosomes <- chromosomes[order(getChromosomeIndices(chromosomes))];
    if (ttest < 1.0 & ttest > 0.0) {
        cat("#Chr\tNum\tMean\tpVal\n");
        for (i in 1:length(chromosomes)) {
            chrom <- chromosomes[[i]];
            results <- pvalues[[chrom]];
            cat(sprintf("%s\t%d\t%.3f\t%.3e\n", chrom, results[1], results[2] - avg, results[3]));
        }
    }

    # draw box plot
    pdf(outputFilename);
    boxplot(dl, ylab="Chromosome", xlab="log2(FC)", horizontal=horizontal, las=labelDirection, col=colors, topupper=inputFilename, ylim=c(-abs(ylimit),abs(ylimit)), cex.axis=0.75, cex.lab=0.75);
    title(inputFilename);
}

    
if (!interactive()) {
    "Interface for running mode from command line.
Rscript trisomy.R [input filename] [output filename (default:boxplot.pdf)]"
    options(echo=FALSE);
    args<-commandArgs(trailingOnly=TRUE);
    if (length(args) == 0) {
        stop("no cuffdiff output file given");
    }

    inputFilename <- args[[1]];
    if (length(args) >= 2) {
        outputFilename <- args[[2]];
    } else {
        outputFilename <- "boxplot.pdf";
    }

    drawExpressionChangesByChromosome(inputFilename, outputFilename=outputFilename, horizontal=TRUE, excludesXY=TRUE, ttest=0.01);
}
