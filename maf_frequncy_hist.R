#!/bin/bash

MAF_orginial = data.table::fread ("MAF_filter.dir/Before_genotype_filter.frq",header = TRUE)
tiff("Missing_SNPs.dir/Before_GEnotype_filtering_histogram.tiff", units="in", width=5, height=5, res=600)

hist (MAF_orginial$MAF, col =  "Blue", main = "Before filtering missigness", xlab = "MAF",
      ylab = "Number of SNPs")
dev.off()

MAF_AFTER = data.table::fread ("MAF_filter.dir/After_genotype_filter.frq",header = TRUE)

tiff("Missing_SNPs.dir/AFTER_GEnotype_filtering_histogram.tiff", units="in", width=5, height=5, res=600)

hist (MAF_AFTER$MAF, col =  "Blue", main = "After filtering missigness", xlab = "MAF",
      ylab = "Number of SNPs")
dev.off()


MAF = data.table::fread ("MAF_filter.dir/After_MAF_filter.frq",header = TRUE)

tiff("Missing_SNPs.dir/AFTER_MAF_filtering_histogram.tiff", units="in", width=5, height=5, res=600)

hist (MAF$MAF, col =  "Blue", main = "After filtering MAF", xlab = "MAF",
      ylab = "Number of SNPs")
dev.off()