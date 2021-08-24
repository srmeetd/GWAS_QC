#!/bin/bash


Before = read.table("Missing_SNPs.dir/Before_filtering_missigness.lmiss", header = TRUE)

tiff("Missing_SNPs.dir/Before_filtering_histogram.tiff", units="in", width=5, height=5, res=600)

hist (Before$F_MISS, col =  "Blue", main = "Before filtering missigness", xlab = "Fraction of missigness",
      ylab = "Number of SNPs")
dev.off()

tiff("Missing_SNPs.dir/After_filtering_histogram.tiff", units="in", width=5, height=5, res=600)
AFter = read.table("Missing_SNPs.dir/After_filtering_missigness.lmiss", header = TRUE)
hist (AFter$F_MISS, col =  "Blue", main = "After filtering missigness", xlab = "Fraction of missigness",
      ylab = "Number of SNPs")
dev.off()

tiff("Missing_SNPs.dir/After_filtering_barplot.tiff", units="in", width=8, height=5, res=600)

bar = aggregate( F_MISS ~ CHR, AFter, mean )
barplot(height=bar$F_MISS, names=bar$CHR, col="#69b3a2",xlab = "Chromosome Names", 
        ylab  = "Average missingness fraction per chrom")
dev.off()