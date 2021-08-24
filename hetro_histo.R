#!/bin/bash

args <- commandArgs(trailingOnly = TRUE)

hetro = read.table(args[1], h = T)

tiff("Hetrozygosity.dir/Before_filtering.tiff", units="in", width=5, height=5, res=600)

hist (hetro$F, col =  "Blue", main = "After filtering missigness", xlab = "Missigness_fraction",
      ylab = "Number of samples")
dev.off()


missing = read.table(args[2], h= T)

a = merge ( hetro,missing, by= "IID")

tiff("Hetrozygosity.dir/misisgness_and_hetro.tiff", units="in", width=5, height=5, res=600)
plot (a$F_MISS,a$F, xlab = "Missigness fraction", ylab = "Hetrozygosity")
dev.off()
