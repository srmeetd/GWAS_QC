#!/bin/bash


args <- commandArgs(trailingOnly = TRUE)
GT_before = data.table::fread ("args[1]",header = TRUE)
tiff("Samples_missing_GT_rate.dir//BEFORE_filtering.tiff", units="in", width=5, height=5, res=600)

hist (GT_before$F_MISS, col =  "Blue", main = "Before filtering missigness", xlab = "Missigness_fraction",
      ylab = "Number of samples")
dev.off()

GT_After = data.table::fread ("Samples_missing_GT_rate.dir/After_filter.imiss",header = TRUE)

tiff("Samples_missing_GT_rate.dir//After_filtering.tiff", units="in", width=5, height=5, res=600)

hist (GT_After$F_MISS, col =  "Blue", main = "After filtering missigness", xlab = "Missigness_fraction",
      ylab = "Number of samples")
dev.off()



